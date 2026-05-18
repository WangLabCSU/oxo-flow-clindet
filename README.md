# oxo-flow-clindet-oxoflow

A clinical variant detection pipeline for cancer genomics, ported from [zyllifeworld/oxo-flow-clindet](https://github.com/zyllifeworld/oxo-flow-clindet) to [oxo-flow](https://github.com/traitome/oxo-flow).

## Overview

This pipeline performs:
1. QC & Trimming (FastP)
2. Alignment (BWA-MEM)
3. Post-alignment processing (GATK MarkDuplicates)
4. Somatic Variant Calling (GATK Mutect2)
5. Annotation (Placeholder for VEP)

## Usage

1. Build `oxo-flow` from source.
2. Prepare your `samplesheet.csv`.
3. Run the pipeline:
   ```bash
   oxo-flow run oxo-flow-clindet.oxoflow
   ```

## Features

- **Metadata-driven**: Uses a flexible samplesheet where any column (like `Target_BED`, `Project`, etc.) is available as a wildcard in rules.
- **Reproducible**: Built on `oxo-flow` with provenance tracking and container/conda support.
