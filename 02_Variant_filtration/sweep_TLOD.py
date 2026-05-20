#!/usr/bin/env python3

# =========================
# Example usage:
# =========================
# TNM_VCF_LIST=/path/to/vcf_TNM_file_list.txt \
# TO_VCF_LIST=/path/to/vcf_TO_file_list.txt \
# OUTPUT_DIR=/path/to/output_directory \
# python sweep_TLOD.py

import os


def parse_info(info_str):
    info = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            info[k] = v
        else:
            info[field] = True
    return info


def load_variants_from_list(vcf_list_file, require_pass=True):
    """
    Load variants from a list of VCF files

    Return:
        dict:
            key = CHR:POS:REF:ALT
            value = dict(INFO)
    """
    variants = {}

    with open(vcf_list_file) as lf:
        for vcf_path in lf:
            vcf_path = vcf_path.strip()

            if not vcf_path:
                continue

            # Use only *_somatic_PASS.vcf
            if not vcf_path.endswith("_somatic_PASS.vcf"):
                continue

            if not os.path.exists(vcf_path):
                print(f"[WARN] Missing file: {vcf_path}")
                continue

            with open(vcf_path) as f:
                for line in f:
                    if line.startswith("#"):
                        continue

                    cols = line.rstrip().split("\t")

                    chrom, pos, _, ref, alt, _, flt, info_str = cols[:8]

                    if require_pass and flt != "PASS":
                        continue

                    info = parse_info(info_str)

                    # Support multi-allelic ALT
                    alts = alt.split(",")

                    for a in alts:
                        key = f"{chrom}:{pos}:{ref}:{a}"

                        # Keep highest TLOD if duplicated
                        if key in variants:
                            if "TLOD" in info and "TLOD" in variants[key]:
                                variants[key]["TLOD"] = str(
                                    max(
                                        float(info["TLOD"]),
                                        float(variants[key]["TLOD"])
                                    )
                                )
                        else:
                            variants[key] = info.copy()

    return variants


# -----------------------------
# Input
# -----------------------------
tnm_vcf_list = os.environ["TNM_VCF_LIST"]
to_vcf_list = os.environ["TO_VCF_LIST"]

output_dir = os.environ["OUTPUT_DIR"]

os.makedirs(output_dir, exist_ok=True)

# -----------------------------
# Load variants
# -----------------------------
tnm_vars = load_variants_from_list(tnm_vcf_list)
to_vars = load_variants_from_list(to_vcf_list)

tnm_set = set(tnm_vars.keys())

print(f"TNM variants: {len(tnm_set)}")
print(f"TO variants : {len(to_vars)}")

# -----------------------------
# Sweep TLOD thresholds
# -----------------------------
thresholds = [3, 3.5, 4, 5, 6]

results = []

for cutoff in thresholds:

    to_filtered = {
        k for k, info in to_vars.items()
        if "TLOD" in info and float(info["TLOD"]) >= cutoff
    }

    overlap = to_filtered & tnm_set

    sensitivity = len(overlap) / len(tnm_set) if tnm_set else 0
    precision = len(overlap) / len(to_filtered) if to_filtered else 0

    results.append([
        cutoff,
        len(to_filtered),
        len(overlap),
        round(sensitivity, 4),
        round(precision, 4)
    ])

# -----------------------------
# Write output
# -----------------------------
out_file = os.path.join(output_dir, "TLOD_sweep_summary.tsv")

with open(out_file, "w") as out:

    out.write(
        "TLOD_cutoff\tTO_variants\tOverlap_with_TNM\tSensitivity\tPrecision\n"
    )

    for r in results:
        out.write("\t".join(map(str, r)) + "\n")

print(f"\n[Done] Result written to: {out_file}")