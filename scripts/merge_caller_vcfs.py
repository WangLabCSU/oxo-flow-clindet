# Ported verbatim from upstream zyllifeworld/clindet @ 582a9131
# workflow/common/scripts/merge_caller_vcfs.py (attribution in NOTICE.md).
# Driven under oxo-flow via scripts/smk.py.
from collections import defaultdict
from typing import List, Optional, Tuple
import pysam
from pathlib import Path

vcf_inputs = list(snakemake.input.vcfs)
caller_names = [ Path(vcf).stem for vcf in vcf_inputs]
caller_to_vcf = dict(zip(caller_names, vcf_inputs))


def _classify_caller(caller_name: str) -> str:
    """Map caller name to a strategy tag for AD/DP extraction."""
    c = caller_name.lower()
    if c.startswith("muse"):
        return "muse"
    if c.startswith("mutect"):
        return "mutect"
    if c.startswith("varscan"):
        return "varscan"
    if c.startswith("cgppindel"):
        return "cgppindel"
    if c.startswith("caveman"):
        return "caveman"
    if "strelka" in c:
        return "strelka"
    if "somaticsniper" in c:
        return "somaticsniper"
    return "generic"


caller_strategies = {caller: _classify_caller(caller) for caller in caller_names}

merged_vcf = snakemake.output.merged_vcf

tumor_sample_name = None
normal_sample_name = None

unifed_tumor_sample_name = None
unifed_normal_sample_name = None

template = pysam.VariantFile(vcf_inputs[0])
header = template.header.copy()
for vcf_path in vcf_inputs:
    vcf = pysam.VariantFile(vcf_path)
    for filt_id in vcf.header.filters:
        if filt_id not in header.filters:
            filt_rec = vcf.header.filters[filt_id]
            description = getattr(filt_rec, "description", "")
            header.filters.add(filt_id, None, None, description)


def ensure_header():
    info_defs = {
        "CALLERS": (".", "String", "Variant callers supporting this variant"),
        "NCALLERS": (1, "Integer", "Number of supporting callers"),
        "CALLER_FILTERS": (".", "String", "Per-caller FILTER values, format caller:filter"),
        "CALLER_QUALS": (".", "String", "Per-caller QUAL values, format caller:qual"),
        "CONSENSUS_PASS": (1, "String", "Whether this variant passed consensus rule"),
        "TUMOR_SAMPLE": (1, "String", "Tumor sample name used in merged record"),
        "NORMAL_SAMPLE": (1, "String", "Normal sample name used in merged record"),
    }
    for key, (num, typ, desc) in info_defs.items():
        if key not in header.info:
            header.info.add(key, num, typ, desc)

    # 最终输出的标准 AD_ALT/DP
    if "AD" not in header.formats:
        header.formats.add(
            "AD", "R", "Integer",
            "Final allele depths using component-wise maximum across callers"
        )
    if "AD_ALT" not in header.formats:
        header.formats.add("AD_ALT", number=1, type="Integer", description="variant allele depth only")
    if "DP" not in header.formats:
        header.formats.add(
            "DP", 1, "Integer",
            "Final depth using maximum across callers"
        )

    if "TDP" not in header.info:
        header.info.add("TDP", 1, "Integer", "Tumor depth (max across callers)")

    if "TVAF" not in header.info:
        header.info.add("TVAF", 1, "Float", "Tumor variant allele frequency")

    if "NDP" not in header.info:
        header.info.add("NDP", 1, "Integer", "Normal depth (max across callers)")

    if "NVAF" not in header.info:
        header.info.add("NVAF", 1, "Float", "Normal variant allele frequency")

ensure_header()

variant_master = {}
variant_callers = defaultdict(set)
variant_filters = defaultdict(dict)
variant_quals = defaultdict(dict)
variant_samples = defaultdict(dict)        # key -> caller -> {"tumor": name, "normal": name}
variant_sampledata = defaultdict(dict)     # key -> caller -> {"tumor": {...}, "normal": {...}}


def variant_key(rec):
    alts = tuple(rec.alts) if rec.alts else ()
    return (rec.chrom, rec.pos, rec.ref, alts)


_CHROM_ORDER = {str(i): i for i in range(1, 23)}
_CHROM_ORDER.update({"X": 23, "Y": 24, "M": 25, "MT": 25})


def chrom_sort_key(chrom):
    c = chrom.replace("chr", "")
    return _CHROM_ORDER.get(c, 1000), chrom


def get_filter_string(rec):
    filt = list(rec.filter.keys())
    return "." if not filt else ";".join(filt)


def get_qual_string(rec):
    return "." if rec.qual is None else str(rec.qual)


def infer_tumor_normal_names(vcf):
    sample_names = list(vcf.header.samples)
    tumor = tumor_sample_name if tumor_sample_name in sample_names else None
    normal = normal_sample_name if normal_sample_name in sample_names else None

    if tumor is None:
        for s in sample_names:
            if "tumor" in s.lower() or "tumour" in s.lower():
                tumor = s
                break
    if normal is None:
        for s in sample_names:
            if "normal" in s.lower() or "control" in s.lower():
                normal = s
                break

    # 已识别出一个样本时，另一个自动取剩余样本
    if tumor is not None and normal is None:
        for s in sample_names:
            if s != tumor:
                normal = s
                break
    if normal is not None and tumor is None:
        for s in sample_names:
            if s != normal:
                tumor = s
                break

    # 最后兜底：按位置推断（大多数 somatic caller 是 normal 在前，tumor 在后）
    if tumor is None and len(sample_names) >= 2:
        tumor = sample_names[1]
    elif tumor is None and len(sample_names) >= 1:
        tumor = sample_names[0]
    if normal is None and len(sample_names) >= 2:
        normal = sample_names[0]

    return tumor, normal


def _as_int_list(value):
    if value is None:
        return None
    if isinstance(value, str):
        if value in (".", ""):
            return None
        parts = value.split(",")
        out = []
        for p in parts:
            if p in (".", ""):
                return None
            out.append(int(float(p)))
        return out
    if isinstance(value, (tuple, list)):
        out = []
        for x in value:
            if x is None or x == ".":
                return None
            out.append(int(x))
        return out
    return [int(value)]


def _get_sample_field(sample, key):
    return sample.get(key)


def extract_ad_dp_from_sample(rec, sample_name, caller):
    if sample_name is None or sample_name not in rec.samples:
        return None, None

    sample = rec.samples[sample_name]
    ref = rec.ref
    alts = list(rec.alts or [])
    n_alt = len(alts)
    strategy = caller_strategies.get(caller, "generic")

    # 1) generic AD + DP (all callers)
    ad = _as_int_list(_get_sample_field(sample, "AD"))
    dp = _get_sample_field(sample, "DP")
    dp = None if dp in (None, ".") else int(dp)
    if ad is not None and len(ad) >= 1 + n_alt:
        if strategy in ("muse", "mutect"):
            return ad[1], dp if dp is not None else int(sum(ad))
        return tuple(ad[:1 + n_alt]), dp if dp is not None else int(sum(ad))

    # 2) VarScan2: RD + AD
    if strategy == "varscan":
        rd = _get_sample_field(sample, "RD")
        ad_alt = _get_sample_field(sample, "AD")
        if rd not in (None, ".") and ad_alt not in (None, "."):
            try:
                ad2 = (int(rd), int(ad_alt))
                return ad2, int(sum(ad2))
            except (ValueError, TypeError):
                pass

    # 3) SomaticSniper: DP4
    if strategy == "somaticsniper":
        dp4 = _as_int_list(_get_sample_field(sample, "DP4"))
        if dp4 is not None and len(dp4) >= 4:
            ref_count = dp4[0] + dp4[1]
            alt_count = dp4[2] + dp4[3]
            return (ref_count, alt_count), int(ref_count + alt_count)

    # 4) Strelka2
    if strategy == "strelka":
        tar = _as_int_list(_get_sample_field(sample, "TAR"))
        tir = _as_int_list(_get_sample_field(sample, "TIR"))
        if tar is not None and tir is not None and n_alt == 1:
            return tir[0], int(tir[0] + tar[0])
        if n_alt == 1:
            ref_tag = f"{ref}U"
            alt_tag = f"{alts[0]}U"
            refu = _as_int_list(_get_sample_field(sample, ref_tag))
            altu = _as_int_list(_get_sample_field(sample, alt_tag))
            if refu is not None and altu is not None:
                return altu[0], dp

    # 5) CgpPindel
    if strategy == "cgppindel":
        pp = _as_int_list(_get_sample_field(sample, "PP"))
        np = _as_int_list(_get_sample_field(sample, "NP"))
        pr = _as_int_list(_get_sample_field(sample, "PR"))
        nr = _as_int_list(_get_sample_field(sample, "NR"))
        if pp is not None and np is not None and pr is not None and nr is not None and n_alt == 1:
            return int(pp[0] + np[0]), int(pr[0] + nr[0])

    # 6) CaVEMan
    if strategy == "caveman":
        faz = _as_int_list(_get_sample_field(sample, "FAZ"))
        fcz = _as_int_list(_get_sample_field(sample, "FCZ"))
        fgz = _as_int_list(_get_sample_field(sample, "FGZ"))
        ftz = _as_int_list(_get_sample_field(sample, "FTZ"))
        raz = _as_int_list(_get_sample_field(sample, "RAZ"))
        rcz = _as_int_list(_get_sample_field(sample, "RCZ"))
        rgz = _as_int_list(_get_sample_field(sample, "RGZ"))
        rtz = _as_int_list(_get_sample_field(sample, "RTZ"))
        if all(x is not None for x in (faz, fcz, fgz, ftz, raz, rcz, rgz, rtz)):
            base_depths = {
                "A": faz[0] + raz[0],
                "C": fcz[0] + rcz[0],
                "G": fgz[0] + rgz[0],
                "T": ftz[0] + rtz[0],
            }
            ad2 = base_depths.get(alts[0])
            dp2 = int(sum(base_depths.get(a, 0) for a in ("A", "C", "G", "T")))
            return ad2, dp2

    if dp is not None:
        return None, dp

    return None, None


def pick_max_vaf_caller(
    ad_list: List[Optional[int]],
    dp_list: List[Optional[int]],
) -> Tuple[Optional[int], Optional[int], Optional[int]]:
    """
    从多个 caller 的 AD/DP 中选出 VAF(=AD/DP) 最大的 caller。

    参数
    ----------
    ad_list : List[Optional[int]]
        各 caller 的肿瘤 ALT AD
    dp_list : List[Optional[int]]
        各 caller 的肿瘤 DP

    返回
    ----------
    (best_index, best_ad, best_dp)
        - best_index: VAF 最大的 caller 下标
        - best_ad:    该 caller 的 AD
        - best_dp:    该 caller 的 DP

    规则
    ----------
    - ad_list 和 dp_list 必须等长
    - 跳过以下无效值：
        - AD is None
        - DP is None
        - DP <= 0
    - 如果所有值都无效，返回 (None, None, None)
    - 若 VAF 相同，保留最先出现的 index
    """
    if len(ad_list) != len(dp_list):
        raise ValueError("ad_list and dp_list must have the same length")

    best_index = None
    best_ad = None
    best_dp = None
    best_vaf = -1.0

    for i, (ad, dp) in enumerate(zip(ad_list, dp_list)):
        # 只跳过 None，不跳过 dp=0
        if ad is None or dp is None:
            continue

        # 归一化 ad：tuple (ref_AD, alt_AD, ...) → alt_AD
        if isinstance(ad, tuple):
            ad = ad[1] if len(ad) > 1 else ad[0]

        # 特殊情况：AD=0 且 DP=0
        if ad == 0 and dp == 0:
            # 如果当前还没有有效结果，记录它
            if best_index is None:
                best_index = i
                best_ad = 0
                best_dp = 0
            continue

        # 避免除零
        if dp == 0:
            continue

        vaf = ad / dp

        if vaf > best_vaf:
            best_vaf = vaf
            best_index = i
            best_ad = ad
            best_dp = dp

    return best_index, best_ad, best_dp



def pick_gt(sample_payloads):
    for payload in sample_payloads:
        gt = payload.get("GT")
        if gt is not None:
            return gt
    return (None, None)


def create_output_record(out_vcf, src):
    new_rec = out_vcf.new_record(
        contig=src["chrom"],
        start=src["start"],
        stop=src["stop"],
        id=src["id"],
        alleles=(src["ref"],) + tuple(src["alts"] or ()),
        qual=src["qual"],
    )

    for k, v in src["info"].items():
        if k in out_vcf.header.info and k not in {
            "CALLER", "CALLERS", "NCALLERS",
            "CALLER_FILTERS", "CALLER_QUALS",
            "CONSENSUS_PASS", "TUMOR_SAMPLE", "NORMAL_SAMPLE"
        }:
            try:
                new_rec.info[k] = v
            except Exception:
                pass

    return new_rec


# 读取所有 caller VCF
for caller, path in caller_to_vcf.items():
    vcf = pysam.VariantFile(path)
    tumor_name, normal_name = infer_tumor_normal_names(vcf)

    if unifed_tumor_sample_name is None:
        unifed_tumor_sample_name = tumor_name
        unifed_normal_sample_name = normal_name

    for rec in vcf:
        key = variant_key(rec)
        variant_callers[key].add(caller)
        variant_filters[key][caller] = get_filter_string(rec)
        variant_quals[key][caller] = get_qual_string(rec)
        variant_samples[key][caller] = {"tumor": unifed_tumor_sample_name, "normal": unifed_normal_sample_name}

        tumor_ad, tumor_dp = extract_ad_dp_from_sample(rec, tumor_name, caller)
        normal_ad, normal_dp = extract_ad_dp_from_sample(rec, normal_name, caller)

        tumor_gt = rec.samples[tumor_name].get("GT") if tumor_name is not None and tumor_name in rec.samples else None
        normal_gt = rec.samples[normal_name].get("GT") if normal_name is not None and normal_name in rec.samples else None

        variant_sampledata[key][caller] = {
            "tumor": {"AD": tumor_ad, "DP": tumor_dp, "GT": tumor_gt},
            "normal": {"AD": normal_ad, "DP": normal_dp, "GT": normal_gt},
        }

        if key not in variant_master:
            variant_master[key] = {
                "chrom": rec.chrom,
                "start": rec.start,
                "stop": rec.stop,
                "id": rec.id,
                "ref": rec.ref,
                "alts": rec.alts,
                "qual": rec.qual,
                "info": {k: rec.info[k] for k in rec.info.keys()},
            }


sorted_keys = sorted(
    variant_master.keys(),
    key=lambda x: (chrom_sort_key(x[0]), x[1], x[2], x[3])
)

with pysam.VariantFile(merged_vcf, "w", header=header) as out_all:
    for key in sorted_keys:
        src_rec = variant_master[key]
        callers = sorted(variant_callers[key])
        ncallers = len(callers)

        rec = create_output_record(out_all, src_rec)
        rec.info["CALLERS"] = callers
        rec.info["NCALLERS"] = ncallers

        master_caller = callers[0]
        tumor_name = variant_samples[key][master_caller]["tumor"] or "."
        normal_name = variant_samples[key][master_caller]["normal"] or "."

        tumor_ads, tumor_dps, tumor_payloads = [], [], []
        normal_ads, normal_dps, normal_payloads = [], [], []

        for caller in callers:
            payload = variant_sampledata[key][caller]
            tumor_payloads.append(payload["tumor"])
            normal_payloads.append(payload["normal"])
            tumor_ads.append(payload["tumor"]["AD"])
            tumor_dps.append(payload["tumor"]["DP"])
            normal_ads.append(payload["normal"]["AD"])
            normal_dps.append(payload["normal"]["DP"])

        # 最终输出：所有 caller 中的最大 DP；AD 用逐分量最大值
        best_index, final_tumor_ad,final_tumor_dp = pick_max_vaf_caller(ad_list = tumor_ads,dp_list = tumor_dps)
        ## 计算肿瘤t_vaf

        final_normal_ad = normal_ads[best_index] if best_index is not None and normal_ads and best_index < len(normal_ads) else None
        final_normal_dp = normal_dps[best_index] if best_index is not None and normal_dps and best_index < len(normal_dps) else None

        # 归一化 AD：tuple (ref_AD, alt_AD) → alt_AD
        if isinstance(final_tumor_ad, tuple):
            final_tumor_ad = final_tumor_ad[1] if len(final_tumor_ad) > 1 else final_tumor_ad[0]
        if isinstance(final_normal_ad, tuple):
            final_normal_ad = final_normal_ad[1] if len(final_normal_ad) > 1 else final_normal_ad[0]

        tvaf = final_tumor_ad / final_tumor_dp if final_tumor_dp else 0
        ## 计算正常对照t_vaf
        nvaf = None
        if final_normal_ad is not None and final_normal_dp not in (None, 0):
            nvaf = final_normal_ad / final_normal_dp if final_normal_dp else 0

        sample_names = list(rec.samples.keys())

        if tumor_name in rec.samples.keys():
            out_sample = rec.samples[tumor_name]
        elif len(sample_names) >= 1:
            out_sample = rec.samples[sample_names[0]]
        else:
            out_sample = None

        if out_sample is not None:
            out_sample["GT"] = pick_gt(tumor_payloads)
            if final_tumor_ad is not None and final_tumor_dp is not None and (len(rec.alts) <= 1):
                out_sample["AD_ALT"] = final_tumor_ad
                out_sample["AD"] = (int(final_tumor_dp - final_tumor_ad), final_tumor_ad)
            if final_tumor_dp is not None:
                out_sample["DP"] = final_tumor_dp

        if normal_name and normal_name in rec.samples.keys():
            out_sample = rec.samples[normal_name]
        elif len(sample_names) >= 2:
            out_sample = rec.samples[sample_names[0]]
        else:
            out_sample = None

        if out_sample is not None:
            out_sample["GT"] = pick_gt(normal_payloads)
            if final_normal_ad is not None and final_normal_dp is not None and (len(rec.alts) <= 1):
                out_sample["AD_ALT"] = final_normal_ad
                out_sample["AD"] = (int(final_normal_dp - final_normal_ad), final_normal_ad)
            if final_normal_dp is not None:
                out_sample["DP"] = final_normal_dp

        if final_tumor_dp is not None:
            rec.info["TDP"] = int(final_tumor_dp)
        if tvaf is not None:
            rec.info["TVAF"] = round(tvaf, 6)
        if final_normal_dp is not None:
            rec.info["NDP"] = int(final_normal_dp)
        if nvaf is not None:
            rec.info["NVAF"] = round(nvaf, 6)

        if len(rec.alleles) < 3:
            out_all.write(rec)
