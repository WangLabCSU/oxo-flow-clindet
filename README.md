# oxo-flow-clindet: Comprehensive Clinical Pipeline

Ported from [zyllifeworld/clindet](https://github.com/zyllifeworld/clindet). This pipeline leverages `oxo-flow` for high-performance, clinical-grade bioinformatics analysis.

## Features

- **Multi-Omics**: Supports DNA (WES/WGS) and RNA-seq.
- **Comprehensive Variant Detection**:
  - **SNV/Indel**: Mutect2, Strelka2, VarDict, Lofreq, MuSE, VarScan2.
  - **CNV**: Amber, Cobalt, Purple, ASCAT, FACETS.
  - **SV**: Manta, Delly, GRIDSS.
  - **RNA**: STAR-Fusion, Arriba, TRUST4.
- **Clinical Reporting**: Automated RMarkdown-based report generation.
- **Metadata-Driven**: Dynamic samplesheets with custom wildcard support.
- **Modular Design**: Sub-workflows can be independently included.

## Project Structure

```text
oxo-flow-clindet/
├── clindet.oxoflow      # Main workflow entry point
├── samplesheet.csv      # Sample & Metadata configuration
├── rules/               # Modular analysis rules
│   ├── common/          # QC, Alignment
│   ├── somatic/         # SNV/CNV callers & filtering
│   ├── sv/              # Structural variant detection
│   ├── unpaired/        # Tumor-only analysis
│   ├── rna/             # Transcriptome quantification & fusion
│   ├── annotation/      # VEP & functional impact
│   └── report/          # Clinical PDF/HTML generation
├── scripts/             # Native R/Python processing scripts
└── envs/                # Conda environment definitions
```

## Setup

1. **Prerequisites**: Install `conda` and build `oxo-flow`.
2. **Resource Data**: Download reference genomes and databases (hg38/b37).
3. **Samplesheet**: Prepare `samplesheet.csv` with the following minimum columns:
   - `pair_id`: Unique identifier for the experiment-control pair.
   - `tumor`: Tumor sample ID.
   - `normal`: Normal sample ID.
   - `Tumor_R1/R2`: Paths to tumor FASTQs.
   - `Normal_R1/R2`: Paths to normal FASTQs.
   - *Optional*: `Target_BED`, `Capture_Kit`, `Project`.

## Execution

### Full Pipeline (Dry-run)
```bash
oxo-flow dry-run clindet.oxoflow
```

### Run Specific Target (e.g., RNA quantification)
```bash
oxo-flow run clindet.oxoflow --target rsem_calculate_expression
```

### Resource Configuration
Modify the `[config]` section in `clindet.oxoflow` to point to your local reference paths:
```toml
[config]
reference = "/path/to/hg38.fa"
genome_build = "hg38"
vep_cache = "/path/to/vep_cache"
```

## Contributing
Optimizations for `oxo-flow-clindet` often involve feeding back improvements to the core `oxo-flow` engine.
