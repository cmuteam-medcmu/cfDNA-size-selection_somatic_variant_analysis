#!/usr/bin/env python3

# =========================
# Example usage:
# =========================
# VCF_LIST=/path/to/vcf_TO_file_list.txt \
# OUTPUT_DIR=/path/to/output_directory \
# python filter_TLOD_GERMQ.py

import os

# -----------------------------
# Paths
# -----------------------------
vcf_list_file = os.environ["VCF_LIST"]

output_path = os.environ["OUTPUT_DIR"]

os.makedirs(output_path, exist_ok=True)

# -----------------------------
# Helper: parse INFO
# -----------------------------
def parse_info(info_str):

    info = {}

    for field in info_str.split(";"):

        if "=" in field:
            k, v = field.split("=", 1)
            info[k] = v

        else:
            info[field] = True

    return info


# -----------------------------
# Load VCF list
# -----------------------------
with open(vcf_list_file) as f:

    vcf_files = [
        line.strip()
        for line in f
        if line.strip().endswith("_somatic_PASS.vcf")
    ]

print(f"Found {len(vcf_files)} somatic PASS VCFs")

# =========================================================
# STEP 1: TLOD >= 4
# =========================================================
step1_outputs = []

for vcf in vcf_files:

    basename = os.path.basename(vcf).replace(
        ".vcf",
        "_TLOD4.vcf"
    )

    out_vcf = os.path.join(output_path, basename)

    with open(vcf) as fin, open(out_vcf, "w") as fout:

        for line in fin:

            if line.startswith("#"):
                fout.write(line)
                continue

            cols = line.rstrip().split("\t")

            info = parse_info(cols[7])

            if "TLOD" not in info:
                continue

            if float(info["TLOD"]) >= 4:
                fout.write(line)

    step1_outputs.append(out_vcf)

print(
    f"STEP 1 done: "
    f"{len(step1_outputs)} VCFs written "
    f"(TLOD >= 4)"
)

# -----------------------------
# RUN STEP 2 ONLY
# -----------------------------
step1_outputs = [
    os.path.join(output_path, f)
    for f in os.listdir(output_path)
    if f.endswith("_somatic_PASS_TLOD4.vcf")
]

print(f"Found {len(step1_outputs)} STEP1 VCFs")

# =========================================================
# STEP 2: GERMQ >= 25
# =========================================================
for vcf in step1_outputs:

    basename = os.path.basename(vcf).replace(
        "_TLOD4.vcf",
        "_TLOD4_GERMQ25.vcf"
    )

    out_vcf = os.path.join(output_path, basename)

    with open(vcf) as fin, open(out_vcf, "w") as fout:

        for line in fin:

            if line.startswith("#"):
                fout.write(line)
                continue

            cols = line.rstrip().split("\t")

            info = parse_info(cols[7])

            if "GERMQ" not in info:
                continue

            if float(info["GERMQ"]) >= 25:
                fout.write(line)

    print(f"STEP 2 done: {os.path.basename(out_vcf)}")