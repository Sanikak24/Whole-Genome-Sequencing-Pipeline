#!/usr/bin/env bash
# WGS Variant Calling Pipeline (SLURM + GATK Best Practices + VQSR)
# Author: Sanika Kamath
# Context: WGS/WES-style workflow on HPC
# Workflow:
#   1) FastQC (raw)
#   2) Cutadapt trimming
#   3) FastQC (trimmed)
#   4) BWA-MEM alignment (GRCh38)
#   5) Sort/Index BAM (SAMtools)
#   6) Mark duplicates (GATK MarkDuplicatesSpark)
#   7) BQSR (BaseRecalibrator + ApplyBQSR)
#   8) Alignment metrics (flagstat + mean depth)
#   9) Variant calling (HaplotypeCaller gVCF + GenotypeGVCFs)
#  10) VQSR (INDEL then SNP)
#  11) Variant calling metrics
#
# Notes:
# - This repo does NOT include FASTQs or reference resources.
# - Update the CONFIG section for your environment.

set -euo pipefail

# SLURM DIRECTIVES (edit for your cluster)
#SBATCH --job-name=wgs_gatk_vqsr
#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your_email@domain.com
# Optional (cluster-specific):
##SBATCH -p compute
##SBATCH -A your_allocation
##SBATCH -M your_cluster

# CONFIG (edit these paths)

# Input FASTQs (paired-end)
FASTQ_R1="/path/to/input/sample_R1.fastq.gz"
FASTQ_R2="/path/to/input/sample_R2.fastq.gz"

# Sample identifiers (used for read group + filenames)
SAMPLE_ID="P5"
LIB_ID="P5"
PLATFORM="ILLUMINA"

# Reference & known-sites resources (GRCh38)
REF_FASTA="/path/to/ref/Homo_sapiens_assembly38.fasta"

DBSNP_VCF="/path/to/known_sites/dbsnp_146.hg38.vcf.gz"
MILLS_INDELS_VCF="/path/to/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# VQSR resources (GRCh38)
HAPMAP_VCF="/path/to/resources/hapmap_3.3.hg38.vcf.gz"
OMNI_VCF="/path/to/resources/1000G_omni2.5.hg38.vcf.gz"
PHASE1_SNPS_VCF="/path/to/resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
AXIOM_POLY_VCF="/path/to/resources/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"

# Output locations
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
OUTDIR="${PROJECT_DIR}/results"
LOGDIR="${PROJECT_DIR}/logs"
WORKDIR="${SLURM_SCRATCH:-${PROJECT_DIR}/scratch}"

THREADS="${SLURM_CPUS_PER_TASK:-24}"

mkdir -p "${OUTDIR}" "${LOGDIR}" "${WORKDIR}"

# Helper: logging
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# Copy outputs back on exit if using scratch
cleanup() {
  log "Copying outputs back to ${OUTDIR} (if any generated in WORKDIR)..."
  rsync -av --ignore-existing "${WORKDIR}/" "${OUTDIR}/" || true
}
trap cleanup EXIT

# Load modules (edit for your environment) OR replace with conda/container
# Example module loads (cluster-specific; comment out if not applicable)
module purge || true
module load fastqc/0.11.9 || true
module load cutadapt/2.10 || true
module load gcc/8.2.0 bwa/0.7.17 samtools/1.21 || true
module load gatk/4.5.0.0 || true

#
# Record tool versions (useful for reproducibility)
# ==============================================================================
log "Recording tool versions..."
{
  echo "fastqc:   $(fastqc --version 2>/dev/null || echo 'NA')"
  echo "cutadapt: $(cutadapt --version 2>/dev/null || echo 'NA')"
  echo "bwa:      $(bwa 2>&1 | head -1 || echo 'NA')"
  echo "samtools: $(samtools --version | head -1 || echo 'NA')"
  echo "gatk:     $(gatk --version 2>/dev/null || echo 'NA')"
} > "${LOGDIR}/tool_versions.txt"

# Step 0: Sanity checks
log "Validating inputs..."
for f in "${FASTQ_R1}" "${FASTQ_R2}" "${REF_FASTA}" "${DBSNP_VCF}" "${MILLS_INDELS_VCF}"; do
  [[ -e "$f" ]] || { echo "ERROR: Missing file: $f" >&2; exit 1; }
done

log "Using WORKDIR=${WORKDIR}"
cd "${WORKDIR}"


# Step 1: Quick FASTQ checks (optional)
log "Quick FASTQ checks (first/last few lines + read count sanity)..."
zcat "${FASTQ_R1}" | head -8 > "${LOGDIR}/${SAMPLE_ID}_R1_head.txt"
zcat "${FASTQ_R2}" | head -8 > "${LOGDIR}/${SAMPLE_ID}_R2_head.txt"
zcat "${FASTQ_R1}" | tail -8 > "${LOGDIR}/${SAMPLE_ID}_R1_tail.txt"
zcat "${FASTQ_R2}" | tail -8 > "${LOGDIR}/${SAMPLE_ID}_R2_tail.txt"

# FASTQ line counts (should be divisible by 4)
R1_LINES=$(zcat "${FASTQ_R1}" | wc -l)
R2_LINES=$(zcat "${FASTQ_R2}" | wc -l)
echo "${R1_LINES}" > "${LOGDIR}/${SAMPLE_ID}_R1_wc_l.txt"
echo "${R2_LINES}" > "${LOGDIR}/${SAMPLE_ID}_R2_wc_l.txt"

# Step 2: FastQC on raw reads
log "Running FastQC on raw reads..."
fastqc "${FASTQ_R1}" -t "${THREADS}" --outdir "${OUTDIR}"
fastqc "${FASTQ_R2}" -t "${THREADS}" --outdir "${OUTDIR}"

# Step 3: Trimming (Cutadapt)
log "Running Cutadapt trimming..."
TRIM_R1="${WORKDIR}/${SAMPLE_ID}_R1.trim.fastq.gz"
TRIM_R2="${WORKDIR}/${SAMPLE_ID}_R2.trim.fastq.gz"

cutadapt -j "${THREADS}" -m 10 -q 20 \
  -a AGATCGGAAGAG -A AGATCGGAAGAG \
  -o "${TRIM_R1}" -p "${TRIM_R2}" \
  "${FASTQ_R1}" "${FASTQ_R2}" \
  > "${LOGDIR}/${SAMPLE_ID}_cutadapt.log"

# Step 4: FastQC on trimmed reads
log "Running FastQC on trimmed reads..."
fastqc "${TRIM_R1}" -t "${THREADS}" --outdir "${OUTDIR}"
fastqc "${TRIM_R2}" -t "${THREADS}" --outdir "${OUTDIR}"

# Step 5: Alignment (BWA-MEM) + sort/index (SAMtools)
log "Aligning reads with BWA-MEM (GRCh38)..."
RG="@RG\tID:${SAMPLE_ID}\tLB:${LIB_ID}\tSM:${SAMPLE_ID}\tPL:${PLATFORM}"

RAW_BAM="${WORKDIR}/${SAMPLE_ID}.sorted.bam"

bwa mem -t "${THREADS}" -R "${RG}" "${REF_FASTA}" "${TRIM_R1}" "${TRIM_R2}" \
  | samtools view -bh - \
  | samtools sort -@ "${THREADS}" -o "${RAW_BAM}"

samtools index "${RAW_BAM}"

# Step 6: Mark duplicates (GATK MarkDuplicatesSpark)
log "Marking duplicates with GATK MarkDuplicatesSpark..."
DEDUP_BAM="${WORKDIR}/${SAMPLE_ID}.dedup.bam"

gatk MarkDuplicatesSpark \
  -I "${RAW_BAM}" \
  -O "${DEDUP_BAM}" \
  --spark-master local["${THREADS}"] \
  > "${LOGDIR}/${SAMPLE_ID}_markdupsspark.log" 2>&1

samtools index "${DEDUP_BAM}"

# Step 7: Base Quality Score Recalibration (BQSR)
log "Running BaseRecalibrator..."
BQSR_TABLE="${WORKDIR}/${SAMPLE_ID}.bqsr.table"

gatk BaseRecalibrator \
  -I "${DEDUP_BAM}" \
  -R "${REF_FASTA}" \
  --known-sites "${MILLS_INDELS_VCF}" \
  --known-sites "${DBSNP_VCF}" \
  --java-options "-XX:ParallelGCThreads=${THREADS}" \
  -O "${BQSR_TABLE}" \
  > "${LOGDIR}/${SAMPLE_ID}_baserecalibrator.log" 2>&1

log "Applying BQSR..."
CLEAN_BAM="${WORKDIR}/${SAMPLE_ID}.dedup.bqsr.bam"

gatk ApplyBQSR \
  -I "${DEDUP_BAM}" \
  -R "${REF_FASTA}" \
  --bqsr-recal-file "${BQSR_TABLE}" \
  --java-options "-XX:ParallelGCThreads=${THREADS}" \
  -O "${CLEAN_BAM}" \
  > "${LOGDIR}/${SAMPLE_ID}_applybqsr.log" 2>&1

samtools index "${CLEAN_BAM}"

# Step 8: Alignment statistics (flagstat + mean depth)
log "Computing alignment statistics..."
samtools flagstat "${CLEAN_BAM}" > "${OUTDIR}/${SAMPLE_ID}_flagstat.txt"

log "Computing mean depth (simple average across covered positions)..."
samtools depth "${CLEAN_BAM}" \
  | awk '{ total += $3; n++ } END { if (n>0) print total/n; else print 0 }' \
  > "${OUTDIR}/${SAMPLE_ID}_mean_depth.txt"

# Step 9: Variant calling (HaplotypeCaller gVCF + GenotypeGVCFs)
log "Calling variants with HaplotypeCaller (gVCF mode)..."
GVCF="${WORKDIR}/${SAMPLE_ID}.g.vcf.gz"

gatk HaplotypeCaller \
  -R "${REF_FASTA}" \
  -I "${CLEAN_BAM}" \
  -O "${GVCF}" \
  -ERC GVCF \
  --native-pair-hmm-threads "${THREADS}" \
  > "${LOGDIR}/${SAMPLE_ID}_haplotypecaller.log" 2>&1

log "Genotyping gVCF with GenotypeGVCFs..."
RAW_VCF="${WORKDIR}/${SAMPLE_ID}.genotyped.vcf.gz"

gatk GenotypeGVCFs \
  -R "${REF_FASTA}" \
  -V "${GVCF}" \
  --java-options "-XX:ParallelGCThreads=${THREADS}" \
  -O "${RAW_VCF}" \
  > "${LOGDIR}/${SAMPLE_ID}_genotypegvcfs.log" 2>&1

# Step 10: VQSR (INDEL then SNP)
log "Running VQSR for INDELs..."
INDEL_RECAL="${WORKDIR}/${SAMPLE_ID}.indels.recal"
INDEL_TRANCHES="${WORKDIR}/${SAMPLE_ID}.indels.tranches"

gatk VariantRecalibrator \
  -R "${REF_FASTA}" \
  -V "${RAW_VCF}" \
  --java-options "-XX:ParallelGCThreads=${THREADS}" \
  -mode INDEL \
  -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
  --resource:axiomPoly,known=false,training=true,truth=false,prior=10 "${AXIOM_POLY_VCF}" \
  --resource:mills,known=false,training=true,truth=true,prior=12 "${MILLS_INDELS_VCF}" \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2 "${DBSNP_VCF}" \
  -O "${INDEL_RECAL}" \
  --tranches-file "${INDEL_TRANCHES}" \
  > "${LOGDIR}/${SAMPLE_ID}_vqsr_indel_recal.log" 2>&1

log "Applying VQSR for INDELs..."
INDEL_FILTERED_VCF="${WORKDIR}/${SAMPLE_ID}.indels.filtered.vcf.gz"

gatk ApplyVQSR \
  -R "${REF_FASTA}" \
  -V "${RAW_VCF}" \
  --java-options "-XX:ParallelGCThreads=${THREADS}" \
  -mode INDEL \
  --recal-file "${INDEL_RECAL}" \
  --tranches-file "${INDEL_TRANCHES}" \
  --truth-sensitivity-filter-level 99.0 \
  -O "${INDEL_FILTERED_VCF}" \
  > "${LOGDIR}/${SAMPLE_ID}_applyvqsr_indel.log" 2>&1

log "Running VQSR for SNPs..."
SNP_RECAL="${WORKDIR}/${SAMPLE_ID}.snps.recal"
SNP_TRANCHES="${WORKDIR}/${SAMPLE_ID}.snps.tranches"

gatk VariantRecalibrator \
  -R "${REF_FASTA}" \
  -V "${INDEL_FILTERED_VCF}" \
  --java-options "-XX:ParallelGCThreads=${THREADS}" \
  -mode SNP \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
  --resource:hapmap,known=false,training=true,truth=true,prior=15 "${HAPMAP_VCF}" \
  --resource:omni,known=false,training=true,truth=true,prior=12 "${OMNI_VCF}" \
  --resource:1000G,known=false,training=true,truth=false,prior=10 "${PHASE1_SNPS_VCF}" \
  --resource:dbsnp,known=true,training=false,truth=false,prior=7 "${DBSNP_VCF}" \
  -O "${SNP_RECAL}" \
  --tranches-file "${SNP_TRANCHES}" \
  > "${LOGDIR}/${SAMPLE_ID}_vqsr_snp_recal.log" 2>&1

log "Applying VQSR for SNPs..."
FINAL_VCF="${OUTDIR}/${SAMPLE_ID}.final.vqsr.vcf.gz"

gatk ApplyVQSR \
  -R "${REF_FASTA}" \
  -V "${INDEL_FILTERED_VCF}" \
  --java-options "-XX:ParallelGCThreads=${THREADS}" \
  -mode SNP \
  --recal-file "${SNP_RECAL}" \
  --tranches-file "${SNP_TRANCHES}" \
  --truth-sensitivity-filter-level 99.0 \
  -O "${FINAL_VCF}" \
  > "${LOGDIR}/${SAMPLE_ID}_applyvqsr_snp.log" 2>&1

# Step 11: Quick peek at final VCF (first 10 non-header rows)
log "Saving first 10 variant records from final VCF..."
zgrep -v "^#" "${FINAL_VCF}" | head -10 > "${OUTDIR}/${SAMPLE_ID}_final_first10.txt"

# Step 12: Variant calling metrics
log "Collecting variant calling metrics..."
gatk CollectVariantCallingMetrics \
  -I "${FINAL_VCF}" \
  --DBSNP "${DBSNP_VCF}" \
  -O "${OUTDIR}/${SAMPLE_ID}_variant_calling_metrics" \
  > "${LOGDIR}/${SAMPLE_ID}_collect_variant_metrics.log" 2>&1

log "Pipeline complete. Outputs in: ${OUTDIR}"
