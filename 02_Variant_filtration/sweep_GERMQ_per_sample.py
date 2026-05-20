#!/usr/bin/env python3

# =========================
# Example usage:
# =========================
# TNM_VCF_LIST=/path/to/vcf_TNM_file_list.txt \
# TO_VCF_LIST=/path/to/vcf_TO_file_list.txt \
# OUTPUT_DIR=/path/to/output_directory \
# python sweep_GERMQ_per_sample.py

import os
import re
from collections import defaultdict

# -----------------------------
# Config
# -----------------------------
tnm_list_file = os.environ["TNM_VCF_LIST"]
to_list_file = os.environ["TO_VCF_LIST"]

output_path = os.environ["OUTPUT_DIR"]

os.makedirs(output_path, exist_ok=True)

TLOD_CUTOFF = 4
GERMQ_CUTOFFS = [10, 20, 25, 30, 40]

# -----------------------------
# Helpers
# -----------------------------
def parse_info(info_str):

    info = {}

    for field in info_str.split(";"):

        if "=" in field:
            k, v = field.split("=", 1)
            info[k] = v

    return info


def load_variants(vcf_path):
    """
    Return dict:
      key = CHR:POS:REF:ALT
      value = INFO dict
    """

    variants = {}

    with open(vcf_path) as f:

        for line in f:

            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")

            chrom, pos, _, ref, alt, _, flt, info_str = cols[:8]

            if flt != "PASS":
                continue

            info = parse_info(info_str)

            # split multi-allelic
            alts = alt.split(",")

            for a in alts:
                key = f"{chrom}:{pos}:{ref}:{a}"
                variants[key] = info

    return variants


def extract_sample_id(path):
    """
    Extract sample ID from filename
    """

    fname = os.path.basename(path)

    m = re.match(r"(Sample\d+|Tumor\d+|Case\d+)", fname)

    return m.group(1) if m else None


def read_vcf_list(list_file):

    files = []

    with open(list_file) as f:

        for line in f:

            path = line.strip()

            if path.endswith("_somatic_PASS.vcf"):
                files.append(path)

    return files


# -----------------------------
# Load file lists
# -----------------------------
tnm_files = read_vcf_list(tnm_list_file)
to_files = read_vcf_list(to_list_file)

# group by sample
tnm_by_sample = defaultdict(list)
to_by_sample = defaultdict(list)

for f in tnm_files:

    sid = extract_sample_id(f)

    if sid:
        tnm_by_sample[sid].append(f)

for f in to_files:

    sid = extract_sample_id(f)

    if sid:
        to_by_sample[sid].append(f)

samples = sorted(set(tnm_by_sample) & set(to_by_sample))

# -----------------------------
# Main analysis
# -----------------------------
results = []

for sample in samples:

    # load TNM variants (union across files)
    tnm_vars = {}

    for vcf in tnm_by_sample[sample]:
        tnm_vars.update(load_variants(vcf))

    tnm_set = set(tnm_vars.keys())

    # load TO variants (union across files)
    to_vars = {}

    for vcf in to_by_sample[sample]:
        to_vars.update(load_variants(vcf))

    # apply fixed TLOD filter first
    to_tlod_filtered = {
        k: info for k, info in to_vars.items()
        if "TLOD" in info and float(info["TLOD"]) >= TLOD_CUTOFF
    }

    for gcut in GERMQ_CUTOFFS:

        to_filtered = {
            k for k, info in to_tlod_filtered.items()
            if "GERMQ" in info and float(info["GERMQ"]) >= gcut
        }

        overlap = to_filtered & tnm_set

        sensitivity = len(overlap) / len(tnm_set) if tnm_set else 0

        precision = len(overlap) / len(to_filtered) if to_filtered else 0

        results.append({
            "Sample": sample,
            "GERMQ_cutoff": gcut,
            "TO_variants": len(to_filtered),
            "Overlap_with_TNM": len(overlap),
            "Sensitivity": round(sensitivity, 4),
            "Precision": round(precision, 4)
        })

# -----------------------------
# Write output
# -----------------------------
out_tsv = os.path.join(output_path, "GERMQ_sweep_per_sample.tsv")

with open(out_tsv, "w") as out:

    out.write(
        "Sample\tGERMQ_cutoff\tTO_variants\tOverlap_with_TNM\tSensitivity\tPrecision\n"
    )

    for r in results:

        out.write(
            f"{r['Sample']}\t"
            f"{r['GERMQ_cutoff']}\t"
            f"{r['TO_variants']}\t"
            f"{r['Overlap_with_TNM']}\t"
            f"{r['Sensitivity']}\t"
            f"{r['Precision']}\n"
        )

print(f"[DONE] GERMQ sweep written to: {out_tsv}")