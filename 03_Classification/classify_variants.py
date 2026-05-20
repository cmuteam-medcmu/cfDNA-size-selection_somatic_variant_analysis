#!/usr/bin/env python3

import os
import glob

# =========================================================
# Example usage
# =========================================================
#
# Example:
#
# export VCF_DIR="/path/to/annotated_vcfs"
# export OUT_ROOT="/path/to/output_dir"
# export COSMIC_TSV="/path/to/CancerMutationCensus_AllData_v99_GRCh38.tsv"
#
# python classify_variants.py
#
# Input VCFs should match:
# *_TLOD4-GERMQ25_annotated.vcf (annotated with VEP)
#
# =========================================================
# CONFIG
# =========================================================

VCF_DIR = os.environ.get(
    "VCF_DIR",
    "/path/to/annotated_vcfs"
)

OUT_ROOT = os.environ.get(
    "OUT_ROOT",
    "/path/to/output_dir"
)

COSMIC_TSV = os.environ.get(
    "COSMIC_TSV",
    "/path/to/CancerMutationCensus_AllData_v99_GRCh38.tsv"
)

VCF_PATTERN = "*_TLOD4-GERMQ25_annotated.vcf"

# =========================================================
# COSMIC utilities
# =========================================================

def norm_chr(c):
    return c.replace("chr", "").strip().upper()


def load_cosmic(path):

    cosmic = {}

    print("[+] Loading COSMIC...")

    with open(path) as f:

        header = [
            h.strip().replace("\ufeff", "")
            for h in f.readline().split("\t")
        ]

        idx = {k: i for i, k in enumerate(header)}

        F_GENE = "GENE_NAME"
        F_POS = "Mutation genome position GRCh38"
        F_REF = "GENOMIC_WT_ALLELE_SEQ"
        F_ALT = "GENOMIC_MUT_ALLELE_SEQ"
        F_TIER = "MUTATION_SIGNIFICANCE_TIER"

        for line in f:

            s = line.rstrip("\n").split("\t")

            gene = s[idx[F_GENE]].strip()
            pos_raw = s[idx[F_POS]].strip()
            ref = s[idx[F_REF]].strip().upper()
            alt = s[idx[F_ALT]].strip().upper()
            tier = s[idx[F_TIER]].strip()

            if ":" not in pos_raw:
                continue

            chrom_part, pos_part = pos_raw.split(":", 1)

            chrom = norm_chr(chrom_part)
            pos = pos_part.split("-")[0]

            if not pos.isdigit():
                continue

            if not (gene and chrom and pos and ref and alt):
                continue

            key = (gene, chrom, pos, ref, alt)

            cosmic[key] = tier

    print(f"[+] COSMIC loaded: {len(cosmic):,} records")

    return cosmic


def cosmic_match(vcf_chrom, vcf_pos, vcf_ref, vcf_alt, gene, cosmic_db):

    chrom = norm_chr(vcf_chrom)
    pos = str(vcf_pos)
    ref = vcf_ref.upper()
    alt = vcf_alt.upper()

    key = (gene, chrom, pos, ref, alt)

    if key in cosmic_db:
        return True, cosmic_db[key]

    return False, None


# =========================================================
# Helpers
# =========================================================

def is_snv(ref, alt):

    return len(ref) == 1 and len(alt) == 1


def is_indel(ref, alt):

    return not is_snv(ref, alt)


def is_lof(consequence):

    lof_terms = [
        "stop_gained",
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "start_lost"
    ]

    if not consequence:
        return False

    terms = consequence.split("&")

    return any(t in lof_terms for t in terms)


# =========================================================
# Tool voting
# =========================================================

def snv_tools(csq):

    results = []

    try:
        results.append(float(csq.get("REVEL_score", 0)) >= 0.75)

    except:
        results.append(False)

    results.append(csq.get("BayesDel_addAF_pred") == "D")
    results.append(csq.get("MetaLR_pred") == "D")
    results.append(csq.get("MetaRNN_pred") == "D")
    results.append(csq.get("MetaSVM_pred") == "D")
    results.append(csq.get("am_class") == "likely_pathogenic")

    return results


def indel_tools(csq):

    results = []

    try:
        results.append(float(csq.get("CADD_PHRED", 0)) >= 20)

    except:
        results.append(False)

    results.append(csq.get("PROVEAN_pred") == "D")
    results.append(csq.get("fathmm-XF_coding_pred") == "D")

    return results


# =========================================================
# Main
# =========================================================

def main():

    # -----------------------------------------------------
    # Validate inputs
    # -----------------------------------------------------

    if not os.path.isdir(VCF_DIR):
        raise FileNotFoundError(
            f"VCF_DIR not found: {VCF_DIR}"
        )

    if not os.path.exists(COSMIC_TSV):
        raise FileNotFoundError(
            f"COSMIC_TSV not found: {COSMIC_TSV}"
        )

    os.makedirs(OUT_ROOT, exist_ok=True)

    # -----------------------------------------------------
    # Load COSMIC
    # -----------------------------------------------------

    COSMIC_DB = load_cosmic(COSMIC_TSV)

    # -----------------------------------------------------
    # Load VCF files
    # -----------------------------------------------------

    vcfs = glob.glob(
        os.path.join(VCF_DIR, VCF_PATTERN)
    )

    print(f"[+] Found {len(vcfs)} VCF files")

    # -----------------------------------------------------
    # Process each sample
    # -----------------------------------------------------

    for vcf in vcfs:

        fname = os.path.basename(vcf)

        sample = fname.replace(
            "_TLOD4-GERMQ25_annotated.vcf",
            ""
        )

        out_dir = os.path.join(OUT_ROOT, sample)

        os.makedirs(out_dir, exist_ok=True)

        A, B1, B2, C = [], [], [], []

        # -------------------------------------------------
        # Read VCF headers
        # -------------------------------------------------

        with open(vcf) as f:

            headers = []
            csq_format = []

            for line in f:

                if line.startswith("##INFO=<ID=CSQ"):

                    fmt = (
                        line.split("Format:")[1]
                        .strip()
                        .strip('">')
                    )

                    csq_format = fmt.split("|")

                if line.startswith("#"):
                    headers.append(line)

                else:
                    break

        headers.append(
            '##INFO=<ID=COSMIC_MATCH,Number=1,Type=String,Description="Matched with COSMIC database">\n'
        )

        headers.append(
            '##INFO=<ID=COSMIC_TIER,Number=1,Type=String,Description="COSMIC mutation significance tier">\n'
        )

        # -------------------------------------------------
        # Process variants
        # -------------------------------------------------

        with open(vcf) as f:

            for line in f:

                if line.startswith("#"):
                    continue

                s = line.rstrip("\n").split("\t")

                chrom, pos, vid, ref, alt, qual, flt, info = s[:8]

                info_dict = {}

                for item in info.split(";"):

                    if "=" in item:

                        k, v = item.split("=", 1)

                        info_dict[k] = v

                # -----------------------------------------
                # Missing CSQ
                # -----------------------------------------

                if "CSQ" not in info_dict:

                    C.append(line)

                    continue

                csq_entries = info_dict["CSQ"].split(",")

                selected = None

                # -----------------------------------------
                # Select canonical transcript
                # -----------------------------------------

                for entry in csq_entries:

                    vals = entry.split("|")

                    csq = dict(zip(csq_format, vals))

                    if csq.get("CANONICAL") == "YES":

                        selected = csq

                        break

                if not selected:

                    selected = dict(
                        zip(
                            csq_format,
                            csq_entries[0].split("|")
                        )
                    )

                gene = selected.get("SYMBOL", "")

                # -----------------------------------------
                # CLASS A: COSMIC match
                # -----------------------------------------

                is_cosmic, tier = cosmic_match(
                    chrom,
                    pos,
                    ref,
                    alt,
                    gene,
                    COSMIC_DB
                )

                if is_cosmic:

                    new_info = (
                        info
                        + ";COSMIC_MATCH=YES"
                        + ";COSMIC_TIER="
                        + tier
                    )

                    s[7] = new_info

                    A.append("\t".join(s) + "\n")

                    continue

                # -----------------------------------------
                # LoF override
                # -----------------------------------------

                consequence = selected.get(
                    "Consequence",
                    ""
                )

                if is_lof(consequence):

                    B1.append(line)

                    continue

                # -----------------------------------------
                # SNV classification
                # -----------------------------------------

                if is_snv(ref, alt):

                    tools = snv_tools(selected)

                    pos_count = sum(tools)

                    if pos_count >= 3:

                        B1.append(line)

                    elif pos_count >= 1:

                        B2.append(line)

                    else:

                        C.append(line)

                # -----------------------------------------
                # INDEL classification
                # -----------------------------------------

                elif is_indel(ref, alt):

                    tools = indel_tools(selected)

                    pos_count = sum(tools)

                    if pos_count == 3:

                        B1.append(line)

                    elif pos_count >= 1:

                        B2.append(line)

                    else:

                        C.append(line)

                else:

                    C.append(line)

        # -------------------------------------------------
        # Write outputs
        # -------------------------------------------------

        def write(path, records):

            with open(path, "w") as o:

                for h in headers:
                    o.write(h)

                for r in records:
                    o.write(r)

        write(
            os.path.join(out_dir, f"{sample}_A.vcf"),
            A
        )

        write(
            os.path.join(out_dir, f"{sample}_B1.vcf"),
            B1
        )

        write(
            os.path.join(out_dir, f"{sample}_B2.vcf"),
            B2
        )

        write(
            os.path.join(out_dir, f"{sample}_C.vcf"),
            C
        )

        print(
            f"[OK] {sample} | "
            f"A={len(A)} "
            f"B1={len(B1)} "
            f"B2={len(B2)} "
            f"C={len(C)}"
        )

    print("\n[✓] All samples processed successfully.")


if __name__ == "__main__":

    main()