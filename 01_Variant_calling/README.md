**Variant Calling**

**Overview**

This module contains scripts used for somatic variant calling using NVIDIA Parabricks Mutect2 workflows.

Two calling modes are included:

    - Tumor-only mode
    - Tumor-normal mode

Both workflows use:
    - Panel of Normals (1000 Genomes PON)
    - gnomAD germline resource
    - GRCh38 reference genome


**Included Scripts**

| Script | Description |
|---|---|
| 1_Mutectcal_Parabrick_TO.sh | Tumor-only variant calling |
| 2_Filtermutectcall_TO.sh | Tumor-only filtering |
| 1_Mutectcal_Parabrick_TNM.sh | Tumor-normal variant calling |
| 2_Filtermutectcall_TNM.sh | Tumor-normal filtering |


**Workflow**
```text
BAM
    ↓
Mutect2
    ↓
PostPON
    ↓
FilterMutectCalls
    ↓
PASS somatic VCF
```
