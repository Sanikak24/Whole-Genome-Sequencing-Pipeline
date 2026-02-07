# Whole-Genome-Sequencing-Pipeline
End-to-end WGS pipeline using SLURM, GATK, and VQSR. Covers FASTQ QC, trimming, alignment to GRCh38, duplicate marking, BQSR, variant calling, and filtering

# WGS Variant Calling Pipeline — SLURM + GATK + VQSR

End-to-end Whole Genome Sequencing (WGS) variant calling pipeline developed as part of graduate coursework in the MS Genome Bioinformatics program at the University of Pittsburgh.

This project demonstrates production-style HPC workflows for transforming paired-end FASTQ files into high-confidence, VQSR-filtered VCFs using GATK Best Practices.


##  Highlights

- Runs on an HPC cluster using SLURM
- Handles raw FASTQs → cleaned BAM → gVCF → final VCF
- Uses GRCh38 reference genome
- Implements GATK VQSR for SNPs and INDELs
- Captures QC metrics and alignment statistics
- Modular, portable, and reproducible


##  Workflow Overview

1. Raw FASTQ QC — FastQC  
2. Adapter trimming — Cutadapt  
3. Post-trim QC — FastQC  
4. Alignment — BWA-MEM (GRCh38)  
5. Sorting & indexing — SAMtools  
6. Duplicate marking — GATK MarkDuplicatesSpark  
7. Base Quality Score Recalibration (BQSR)  
8. Alignment QC — flagstat + mean depth  
9. Variant calling — HaplotypeCaller (gVCF)  
10. Joint genotyping — GenotypeGVCFs  
11. Variant Quality Score Recalibration (VQSR)  
12. Variant calling metrics



##  Tools Used

- FastQC
- Cutadapt
- BWA
- SAMtools
- GATK 



