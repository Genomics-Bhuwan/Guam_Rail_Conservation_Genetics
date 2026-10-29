


```bash

set -euo pipefail

# =========================================================
# Guam Rail: Filter annotated variants by consequence
# =========================================================

WORKING_DIR="/home/bistbs/Guam_rail_SNPeff"
SNPSIFT_JAR="${WORKING_DIR}/snpEff/SnpSift.jar"
JAVA_EXEC="${WORKING_DIR}/snpEff/jdk-21.0.2/bin/java"

INPUT_VCF="${WORKING_DIR}/Guam_rail_annotated/Guam_rail_final_filtered_annotated.vcf"

OUTPUT_DIR="${WORKING_DIR}/Variants_by_Consequences"

# =========================================================
# 1. Check required files
# =========================================================

echo "Checking input files..."

if [[ ! -s "$INPUT_VCF" ]]; then
    echo "ERROR: Input VCF not found or empty:"
    echo "$INPUT_VCF"
    exit 1
fi

if [[ ! -s "$SNPSIFT_JAR" ]]; then
    echo "ERROR: SnpSift.jar not found:"
    echo "$SNPSIFT_JAR"
    exit 1
fi

if [[ ! -x "$JAVA_EXEC" ]]; then
    echo "ERROR: Java executable not found:"
    echo "$JAVA_EXEC"
    exit 1
fi

# =========================================================
# 2. Create separate folders for each consequence category
# =========================================================

mkdir -p \
    "$OUTPUT_DIR/Missense" \
    "$OUTPUT_DIR/Synonymous" \
    "$OUTPUT_DIR/LoF_and_Splicing" \
    "$OUTPUT_DIR/Intergenic"

echo "Output folders created."

# =========================================================
# 3. Extract missense variants
# =========================================================

echo "Extracting missense variants..."

"$JAVA_EXEC" -jar "$SNPSIFT_JAR" filter \
    "ANN[*].EFFECT has 'missense_variant'" \
    "$INPUT_VCF" \
    > "$OUTPUT_DIR/Missense/Guam_rail_missense.vcf"

# =========================================================
# 4. Extract synonymous variants
# =========================================================

echo "Extracting synonymous variants..."

"$JAVA_EXEC" -jar "$SNPSIFT_JAR" filter \
    "ANN[*].EFFECT has 'synonymous_variant'" \
    "$INPUT_VCF" \
    > "$OUTPUT_DIR/Synonymous/Guam_rail_synonymous.vcf"

# =========================================================
# 5. Extract LoF and splice-related variants
# =========================================================

echo "Extracting LoF and splice-related variants..."

"$JAVA_EXEC" -jar "$SNPSIFT_JAR" filter \
    "(ANN[*].EFFECT has 'transcript_ablation') | (ANN[*].EFFECT has 'splice_donor_variant') | (ANN[*].EFFECT has 'splice_acceptor_variant') | (ANN[*].EFFECT has 'stop_gained') | (ANN[*].EFFECT has 'frameshift_variant') | (ANN[*].EFFECT has 'inframe_insertion') | (ANN[*].EFFECT has 'inframe_deletion') | (ANN[*].EFFECT has 'splice_region_variant')" \
    "$INPUT_VCF" \
    > "$OUTPUT_DIR/LoF_and_Splicing/Guam_rail_lof_splicing.vcf"

# =========================================================
# 6. Extract intergenic variants
# =========================================================

echo "Extracting intergenic variants..."

"$JAVA_EXEC" -jar "$SNPSIFT_JAR" filter \
    "ANN[*].EFFECT has 'intergenic_region'" \
    "$INPUT_VCF" \
    > "$OUTPUT_DIR/Intergenic/Guam_rail_intergenic.vcf"

# =========================================================
# 7. Count and verify output variants
# =========================================================

count_variants() {
    awk '!/^#/ {n++} END {print n+0}' "$1"
}

echo ""
echo "========================================================="
echo "GUAM RAIL CONSEQUENCE FILTERING COMPLETE"
echo "========================================================="

printf "Missense variants:       "
count_variants "$OUTPUT_DIR/Missense/Guam_rail_missense.vcf"

printf "Synonymous variants:     "
count_variants "$OUTPUT_DIR/Synonymous/Guam_rail_synonymous.vcf"

printf "LoF/splicing variants:   "
count_variants "$OUTPUT_DIR/LoF_and_Splicing/Guam_rail_lof_splicing.vcf"

printf "Intergenic variants:     "
count_variants "$OUTPUT_DIR/Intergenic/Guam_rail_intergenic.vcf"

echo ""
echo "Output directory: $OUTPUT_DIR"
echo ""
find "$OUTPUT_DIR" -maxdepth 2 -type f -name '*.vcf' -print
echo ""
echo "Done."
```
