```bash

set -euo pipefail

# ==========================================
# Guam Rail SnpEff Annotation
# ==========================================

# Project paths
PROJECT="$HOME/Guam_rail_SNPeff"
SNPEFF="$PROJECT/snpEff"

# Activate Java 21
export JAVA_HOME="$SNPEFF/jdk-21.0.2"
export PATH="$JAVA_HOME/bin:$PATH"

# Genome information
GENOME_ID="bHypOws1"
GENOME_NAME="Guam_rail"

GENOME_FASTA="$PROJECT/bHypOws1_hifiasm.bp.p_ctg.fasta"
GENOME_GFF="$PROJECT/Guam_rail.withExons.gff.gz"

# Input and output files
INPUT_VCF="$PROJECT/Guam_rail_final_filtered.biallelic.vcf"
OUTPUT_VCF="$PROJECT/Guam_rail_final_filtered_annotated.vcf.gz"
REPORT_HTML="$PROJECT/Guam_rail_snpeff_report.html"

# ==========================================
# Step 1: Check Java and required tools
# ==========================================

echo "Checking Java version..."
java -version

for tool in bgzip tabix bcftools; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool is not available."
        echo "Load the required module or activate the appropriate environment."
        exit 1
    fi
done

# Check required input files
for file in \
    "$SNPEFF/snpEff.jar" \
    "$SNPEFF/snpEff.config" \
    "$GENOME_FASTA" \
    "$GENOME_GFF" \
    "$INPUT_VCF"; do
    if [[ ! -s "$file" ]]; then
        echo "ERROR: Missing or empty file: $file"
        exit 1
    fi
done

# ==========================================
# Step 2: Prepare database files
# ==========================================

echo "Preparing Guam rail database files..."

mkdir -p "$SNPEFF/data/$GENOME_ID"

cp "$GENOME_FASTA" "$SNPEFF/data/$GENOME_ID/sequences.fa"

gunzip -c "$GENOME_GFF" \
    > "$SNPEFF/data/$GENOME_ID/genes.gff"

# Register genome ID only if absent
if ! grep -q "^${GENOME_ID}\.genome[[:space:]]*:" \
    "$SNPEFF/snpEff.config"; then
    echo "${GENOME_ID}.genome : ${GENOME_NAME}" \
        >> "$SNPEFF/snpEff.config"
fi

# ==========================================
# Step 3: Build the SnpEff database
# ==========================================

cd "$SNPEFF"

echo "Building SnpEff database for Guam rail..."

java -Xmx70g -jar snpEff.jar build \
    -gff3 \
    -v \
    "$GENOME_ID"

# Stop if the expected database was not created
if [[ ! -s "$SNPEFF/data/$GENOME_ID/snpEffectPredictor.bin" ]]; then
    echo "ERROR: SnpEff database file was not created."
    echo "Check the build output before proceeding to annotation."
    exit 1
fi

echo "Database build completed successfully."

# ==========================================
# Step 4: Annotate the variants
# ==========================================

echo "Annotating Guam rail variants..."

java -Xmx32g -jar snpEff.jar ann \
    -v \
    -stats "$REPORT_HTML" \
    "$GENOME_ID" \
    "$INPUT_VCF" \
    | bgzip -c > "$OUTPUT_VCF"

# ==========================================
# Step 5: Index the annotated VCF
# ==========================================

echo "Indexing annotated VCF..."

tabix -f -p vcf "$OUTPUT_VCF"

# ==========================================
# Step 6: Verify the output
# ==========================================

echo "Checking annotated VCF..."

echo "Checking ANN header..."
bcftools view -h "$OUTPUT_VCF" | grep 'ID=ANN' || true

echo "Number of variant records:"
bcftools view -H "$OUTPUT_VCF" | wc -l

echo "First five annotated records:"
bcftools view -H "$OUTPUT_VCF" | head -5 || true

echo ""
echo "Generated files:"
ls -lh "$OUTPUT_VCF" "$OUTPUT_VCF.tbi" "$REPORT_HTML"

echo ""
echo "Guam rail annotation workflow finished."

```
