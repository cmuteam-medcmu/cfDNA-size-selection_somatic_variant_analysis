**Variant Classification**

**Overview**

This module classifies filtered somatic variants into four evidence-based categories:

    - Class A
    - Class B1
    - Class B2
    - Class C

Classification is based on:
    - COSMIC matching
    - Loss-of-function annotation
    - Ensemble pathogenicity prediction tools

Input VCF files must contain filtered somatic variants that have already been annotated using the Ensembl Variant Effect Predictor (VEP) web-based annotation platform.

The annotation configuration used in this study is provided in:

```text
VariantEffectPredictor-Homo_sapiens-Ensembl-genome-browser_115.pdf
```

Annotated VCF files should be downloaded from the 
Ensembl VEP web interface in VCF format and supplied as input for this classification workflow.

**Classification Logic**

Class A
    COSMIC matched variants

Class B1
    Strong pathogenic evidence

Class B2
    Moderate pathogenic evidence

Class C
    Weak or uncertain evidence