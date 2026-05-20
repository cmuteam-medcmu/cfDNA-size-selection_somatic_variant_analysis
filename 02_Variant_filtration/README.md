**TLOD Sweep Analysis**

**Overview**

This repository contains scripts used for post-filtering analysis of somatic variants identified by Mutect2-based workflows in both tumor-only and tumor-normal modes.

The `sweep_TLOD.py` script evaluates the effect of different TLOD thresholds on variant overlap between tumor-only and tumor-normal variant call sets.

The scripts provided in this repository are minimally modified versions of the original execution scripts used in this study. Sample identifiers and file paths have been anonymized while preserving the original computational workflow.


**Directory Structure**

```text
project/
├── scripts/
│   ├── sweep_TLOD.py
│   ├── 1_RUN_sweep_TLOD.sh
|   ├── sweep_GERMQ_per_sample.py
|   ├── 2_RUN_sweep_GERMQ.sh
|   ├── sweep_GERMQ_per_sample.py
|   ├── 3_RUN_filter_TLOD_GERMQ.sh
|   └── filter_TLOD_GERMQ.py
├── examples/
│   ├── vcf_TNM_file_list_example.txt
│   └── vcf_TO_file_list_example.txt
└── README.md
```

**Requirements**

    - Python 3.12 or later
    - Linux-based environment
    - Plain text input lists containing VCF file paths


**Input File List Format**

The `sweep_TLOD.py` script requires plain text files containing paths to VCF files, with one file path per line.

Example:

```text
/path/to/Sample01_TNM_somatic_PASS.vcf
/path/to/Sample02_TNM_somatic_PASS.vcf
/path/to/Sample03_TNM_somatic_PASS.vcf
```

Only files ending with `_somatic_PASS.vcf` are processed.

Variants are matched based on exact `CHROM:POS:REF:ALT` identity.


**Usage**

The `sweep_TLOD.py` script is executed through the accompanying bash launcher script.

Example:

```bash
TNM_VCF_LIST=/path/to/vcf_TNM_file_list.txt \
TO_VCF_LIST=/path/to/vcf_TO_file_list.txt \
OUTPUT_DIR=/path/to/output_directory \
bash 1_RUN_sweep_TLOD.sh
```
After completion of TLOD threshold evaluation, the GERMQ threshold analysis performed using the same tumor-only and tumor-normal VCF input lists.

Example:

```bash
TNM_VCF_LIST=/path/to/vcf_TNM_file_list.txt \
TO_VCF_LIST=/path/to/vcf_TO_file_list.txt \
OUTPUT_DIR=/path/to/output_directory \
bash 2_RUN_sweep_GERMQ.sh
```
The sweep_GERMQ_per_sample.py script evaluates GERMQ thresholds after applying a fixed TLOD cutoff (TLOD >= 4) to tumor-only variants.

The final post-filtering step retains only tumor-only variants passing both optimaized TLOD and GERMQ thresholds.

Example:

```bash
TO_VCF_LIST=/path/to/vcf_TO_file_list.txt \
OUTPUT_DIR=/path/to/output_directory \
bash 3_RUN_filter_TLOD_GERMQ.sh
```

This final filtering step uses only tumor-only VCF files and generates filtered VCF outputs containing variants that satisfy:

```text
    TLOD >= 4
    GERMQ >= 25
```

**Output**

The workflow generates the following output files:

**TLOD threshold analysis**

```text
TLOD_sweep_summary.tsv
```

Contains:
    - TLOD cutoff
    - Number of tumor-only variants passing the threshold
    - Number of overlapping variants with tumor-normal calls
    - Sensitivity
    - Precision

**GERMQ threshold analysis**

```text
GERMQ_sweep_per_sample.tsv
```

Contains:
    - Sample ID
    - GERMQ cutoff
    - Number of filtered tumor-only variants
    - Number of overlap with tumor-normal variants
    - Sensitivity
    - Precision

**Final filtered VCF files**

```text
*_TLOD4.vcf
*_TLOD4_GERMQ25.vcf
```

Where:
    - *_TLOD4.vcf contains variants passing TLOD >= 4
    - *TLOD4_GERMQ25.vcf contains variants passing both TLOD >= 4 and GERMQ >= 25

**Notes**

    - Only PASS variants are analyzed.
    - Multi-allelic variants are supported.
    - When duplicated variants are detected across multiple VCF files, the highest TLOD value is retained.
    - GERMQ threshold analysis is performed after applying a fixed TLOD filter.
    - Final variant filtration is applied only to tumor-only variant call sets.
    - The workflow preserves the original analysis logic used in the study.


**Reproducibility**

The scripts in this repository were used for downstream comparison of somatic variant calls generated from tumor-only and tumor-normal pipelines.

All file paths and sample identifiers shown in this repository are anonymized placeholders for reproducibility purposes.