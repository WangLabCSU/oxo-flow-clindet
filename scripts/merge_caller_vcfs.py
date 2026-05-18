from collections import defaultdict
from typing import List, Optional, Tuple
import pysam
from pathlib import Path

# vcf_inputs = [
#     'test/hg38/results/vcf_norm/paired/CGGA_1288/cgppindel.vcf',
#     'test/hg38/results/vcf_norm/paired/CGGA_1288/caveman.vcf',
#     'test/hg38/results/vcf_norm/paired/CGGA_1288/muse.vcf',
#     'test/hg38/results/vcf_norm/paired/CGGA_1288/Mutect2.vcf',
#     'test/hg38/results/vcf_norm/paired/CGGA_1288/strelkasomaticmanta.vcf']

# vcf_inputs = [
#     'test/b37/results/vcf_norm/paired/CGGA_1288/cgppindel.vcf',
#     'test/b37/results/vcf_norm/paired/CGGA_1288/caveman.vcf',
#     'test/b37/results/vcf_norm/paired/CGGA_1288/muse.vcf',
#     'test/b37/results/vcf_norm/paired/CGGA_1288/Mutect2.vcf',
#     'test/b37/results/vcf_norm/paired/CGGA_1288/strelkasomaticmanta.vcf']

vcf_inputs = [
    'test/b37/results/vcf_norm/paired/CGGA_1287/cgppindel.vcf',
    'test/b37/results/vcf_norm/paired/CGGA_1287/caveman.vcf',
    'test/b37/results/vcf_norm/paired/CGGA_1287/muse.vcf',
    'test/b37/results/vcf_norm/paired/CGGA_1287/Mutect2.vcf',
    'test/b37/results/vcf_norm/paired/CGGA_1287/strelkasomaticmanta.vcf']
    
# vcf_inputs = list(snakemake.input.vcfs)
caller_names = [ Path(vcf).stem for vcf in vcf_inputs]
caller_to_vcf = dict(zip(caller_names, vcf_inputs))

# merged_vcf = snakemake.output.merged_vcf
# consensus_vcf = snakemake.output.consensus_vcf
summary_tsv = None

# tumor_sample_name = 'CGGA_1288_T'
# normal_sample_name = 'CGGA_1288_NC'

# unifed_tumor_sample_name = 'CGGA_1288_T'
# unifed_normal_sample_name = 'CGGA_1288_NC'

tumor_sample_name = 'CGGA_1287_T'
normal_sample_name = 'CGGA_1287_NC'

unifed_tumor_sample_name = 'CGGA_1287_T'
unifed_normal_sample_name = 'CGGA_1287_NC'


template = pysam.VariantFile(vcf_inputs[0])
header = template.header.copy()
for vcf_path in vcf_inputs:
    vcf = pysam.VariantFile(vcf_path)
    for filt_id in vcf.header.filters:
        if filt_id not in header.filters:
            filt_rec = vcf.header.filters[filt_id]
            description = getattr(filt_rec, "description", "")
            header.filters.add(filt_id, None, None, description)


def sanitize_tag(name: str) -> str:
    return "".join(ch for ch in name.upper() if ch.isalnum())


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

    # 每个 caller 的 AD / DP
    # in this version, will not out put per-caller info
    # for caller in caller_names:
    #     tag = sanitize_tag(caller)
    #     ad_tag = f"AD_{tag}"
    #     dp_tag = f"DP_{tag}"
    #     if ad_tag not in header.formats:
    #         header.formats.add(
    #             ad_tag, "R", "Integer",
    #             f"Allele depths reported by {caller}"
    #         )
    #     if dp_tag not in header.formats:
    #         header.formats.add(
    #             dp_tag, 1, "Integer",
    #             f"Read depth reported by {caller}"
    #         )
 
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


def chrom_sort_key(chrom):
    c = chrom.replace("chr", "")
    order = {str(i): i for i in range(1, 23)}
    order.update({"X": 23, "Y": 24, "M": 25, "MT": 25})
    return order.get(c, 1000), chrom


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

    if tumor is None and len(sample_names) >= 1:
        tumor = sample_names[0]
    if normal is None and len(sample_names) >= 2:
        normal = sample_names[1]

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
    try:
        return sample.get(key)
    except Exception:
        return None


def extract_ad_dp_from_sample(rec, sample_name, caller):
    """
    Return:
        ad_tuple: (ref_count, alt_count1, ...)
        dp_int
    """
    if sample_name is None or sample_name not in rec.samples:
        return None, None

    sample = rec.samples[sample_name]
    ref = rec.ref
    alts = list(rec.alts or [])
    n_alt = len(alts)
    caller_l = caller.lower()

    # 1) 通用 AD + DP
    ad = _as_int_list(_get_sample_field(sample, "AD"))
    dp = _get_sample_field(sample, "DP")
    dp = None if dp in (None, ".") else int(dp)
    if ad is not None and len(ad) >= 1 + n_alt:
        # 1.1) Muse AD + DP
        if caller_l.startswith("muse"):
            return ad[1], dp if dp is not None else int(sum(ad))
        elif caller_l.startswith("mutect"):
            return ad[1], dp if dp is not None else int(sum(ad))
        else:
            return tuple(ad[:1 + n_alt]), dp if dp is not None else int(sum(ad))

    # 2) VarScan2: RD + AD
    if caller_l.startswith("varscan"):
        rd = _get_sample_field(sample, "RD")
        ad_alt = _get_sample_field(sample, "AD")
        try:
            if rd not in (None, ".") and ad_alt not in (None, "."):
                ad2 = (int(rd), int(ad_alt))
                return ad2, int(sum(ad2))
        except Exception:
            pass

    # 3) SomaticSniper: DP4 = refF,refR,altF,altR
    dp4 = _get_sample_field(sample, "DP4")
    dp4 = _as_int_list(dp4)
    if dp4 is not None and len(dp4) >= 4:
        ref_count = dp4[0] + dp4[1]
        alt_count = dp4[2] + dp4[3]
        ad2 = (ref_count, alt_count)
        return ad2, int(sum(ad2))

    # 4) Strelka2 indel: TAR/TIR
    tar = _as_int_list(_get_sample_field(sample, "TAR"))
    tir = _as_int_list(_get_sample_field(sample, "TIR"))
    if tar is not None and tir is not None and n_alt == 1:
        ad2 = tir[0]
        dp2 =  int(sum((tir[0] ,tar[0])))
        return ad2, dp2

    # 5) Strelka2 SNV: refU / altU
    if n_alt == 1:
        ref_tag = f"{ref}U"
        alt_tag = f"{alts[0]}U"
        refu = _as_int_list(_get_sample_field(sample, ref_tag))
        altu = _as_int_list(_get_sample_field(sample, alt_tag))
        if refu is not None and altu is not None:
            ad2 = altu[0] ## use only tier 1
            return ad2, dp

    # 6) CgpPindel indel: AD = PP + NP,DP = PR + NR
    if caller_l.startswith("cgppindel"):
        pp = _as_int_list(_get_sample_field(sample, "PP"))
        np = _as_int_list(_get_sample_field(sample, "NP"))
        pr = _as_int_list(_get_sample_field(sample, "PR"))
        nr = _as_int_list(_get_sample_field(sample, "NR"))
        if pp is not None and np is not None and pr is not None and nr is not None and n_alt == 1:
            ad2 = int(sum((pp[0], np[0])))
            dp2 = int(sum((pr[0], nr[0])))
            return ad2, dp2

    # 7) caveman indel: AD = PP + NP,DP = PR + NR
    if caller_l.startswith("caveman"):
        faz = _as_int_list(_get_sample_field(sample, "FAZ"))
        fcz = _as_int_list(_get_sample_field(sample, "FCZ"))
        fgz = _as_int_list(_get_sample_field(sample, "FGZ"))
        ftz = _as_int_list(_get_sample_field(sample, "FTZ"))
        raz = _as_int_list(_get_sample_field(sample, "RAZ"))
        rcz = _as_int_list(_get_sample_field(sample, "RCZ"))
        rgz = _as_int_list(_get_sample_field(sample, "RGZ"))
        rtz = _as_int_list(_get_sample_field(sample, "RTZ"))
        alleles = ['A','T','G','C']
        if all(x is not None for x in (faz, fcz, fgz, ftz, raz, rcz, rgz, rtz)):
            base_depths = {
                "A": faz[0] + raz[0],
                "C": fcz[0] + rcz[0],
                "G": fgz[0] + rgz[0],
                "T": ftz[0] + rtz[0],
            }
            ad2 = base_depths.get(alts[0])
            dp2 = int(sum([base_depths.get(allele, "") for allele in alleles]))
            return ad2, dp2
    # 6) 只有 DP
    if dp is not None:
        return None, dp

    return None, None


def pad_ad(ad_tuple, target_len):
    if ad_tuple is None:
        return None
    if len(ad_tuple) >= target_len:
        return tuple(ad_tuple[:target_len])
    return tuple(list(ad_tuple) + [0] * (target_len - len(ad_tuple)))


from typing import List, Optional, Tuple

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


def create_output_record(out_vcf, src_rec):
    filt = list(src_rec.filter.keys())
    if not filt:
        filt = None

    new_rec = out_vcf.new_record(
        contig=src_rec.chrom,
        start=src_rec.start,
        stop=src_rec.stop,
        id=src_rec.id,
        alleles=(src_rec.ref,) + tuple(src_rec.alts or ()),
        qual=src_rec.qual,
    )

    # for f in src_rec.filter.keys():
    #     new_rec.filter.add(f)

    # 尽量保留原 INFO（前提是 header 里存在）
    for k in src_rec.info.keys():
        if k in out_vcf.header.info and k not in {
            "CALLER", "CALLERS", "NCALLERS",
            "CALLER_FILTERS", "CALLER_QUALS",
            "CONSENSUS_PASS", "TUMOR_SAMPLE", "NORMAL_SAMPLE"
        }:
            try:
                new_rec.info[k] = src_rec.info[k]
            except Exception:
                pass

    return new_rec


# 读取所有 caller VCF
for caller, path in caller_to_vcf.items():
    vcf = pysam.VariantFile(path)
    tumor_name, normal_name = infer_tumor_normal_names(vcf)

    for rec in vcf:
        key = variant_key(rec)
        variant_callers[key].add(caller)
        variant_filters[key][caller] = get_filter_string(rec)
        variant_quals[key][caller] = get_qual_string(rec)
        variant_samples[key][caller] = {"tumor": unifed_tumor_sample_name, "normal": unifed_normal_sample_name}

        tumor_ad, tumor_dp = extract_ad_dp_from_sample(rec, tumor_name, caller)
        normal_ad, normal_dp = extract_ad_dp_from_sample(rec, normal_name, caller)

        tumor_gt = rec.samples[tumor_name].get("GT") if tumor_name in rec.samples else None
        normal_gt = rec.samples[normal_name].get("GT") if normal_name in rec.samples else None

        variant_sampledata[key][caller] = {
            "tumor": {"AD": tumor_ad, "DP": tumor_dp, "GT": tumor_gt},
            "normal": {"AD": normal_ad, "DP": normal_dp, "GT": normal_gt},
        }

        if key not in variant_master:
            variant_master[key] = rec


sorted_keys = sorted(
    variant_master.keys(),
    key=lambda x: (chrom_sort_key(x[0]), x[1], x[2], x[3])
)

merged_vcf = '/public/ClinicalExam/lj_sih/projects/project_clindet/workflow/common/scripts/test_merge.vcf'
with pysam.VariantFile(merged_vcf, "w", header=header) as out_all: #, \
    #  pysam.VariantFile(consensus_vcf, "w", header=header) as out_cons:
    for key in sorted_keys:
        src_rec = variant_master[key]
        callers = sorted(variant_callers[key])
        ncallers = len(callers)

        filter_items = [f"{c}:{variant_filters[key].get(c, '.')}" for c in callers]
        qual_items = [f"{c}:{variant_quals[key].get(c, '.')}" for c in callers]

        rec = create_output_record(out_all, src_rec)
        rec.info["CALLERS"] = callers
        rec.info["NCALLERS"] = ncallers

        master_caller = callers[0]
        tumor_name = variant_samples[key][master_caller]["tumor"] or "."
        normal_name = variant_samples[key][master_caller]["normal"] or "."

        allele_count = 1 + len(rec.alts or [])
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
        t_vaf = final_tumor_ad / final_tumor_dp if final_tumor_dp else 0
        ## 计算肿瘤t_vaf

        final_normal_ad = normal_ads[best_index] if normal_ads and best_index < len(normal_ads) else None
        final_normal_dp = normal_dps[best_index] if normal_dps and best_index < len(normal_dps) else None
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
            if final_tumor_ad is not None and (len(rec.alts) <= 1):
                out_sample["AD_ALT"] = final_tumor_ad
                out_sample["AD"] = (int(final_tumor_dp - final_tumor_ad), final_tumor_ad)
            if final_tumor_dp is not None:
                out_sample["DP"] = final_tumor_dp

            # for caller in caller_names:
            #     tag = sanitize_tag(caller)
            #     ad_tag = f"AD_{tag}"
            #     dp_tag = f"DP_{tag}"
            #     payload = variant_sampledata[key].get(caller, {}).get("tumor", {})
            #     ad_val = pad_ad(payload.get("AD"), allele_count)
            #     dp_val = payload.get("DP")
            #     if ad_val is not None:
            #         out_sample[ad_tag] = ad_val
            #     if dp_val is not None:
            #         out_sample[dp_tag] = int(dp_val)

        if normal_name in rec.samples.keys():
            out_sample = rec.samples[normal_name]
        elif len(sample_names) >= 2:
            out_sample = rec.samples[sample_names[1]]
        else:
            out_sample = None

        if out_sample is not None:
            out_sample["GT"] = pick_gt(normal_payloads)
            if final_normal_ad is not None and (len(rec.alts) <= 1):
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

        if len(rec.alleles) >= 3:
            next
        else:
            out_all.write(rec)