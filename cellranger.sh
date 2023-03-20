#!/bin/bash
#SBATCH -J donor18_cellrngrcount
#SBATCH -e /scratch/ysun81_lab/Nora/CovidDonor18_starsolo.err
#SBATCH -o /scratch/ysun81_lab/Nora/CovidDonor18_starsolo.out
#SBATCH -t 72:00:00 
#SBATCH -c 8
#SBATCH --partition=standard
#SBATCH --mem=200G

module purge all
module load rnastar/2.7.8a

STAR --runThreadN 8 --soloType CB_UMI_Simple --soloFeatures Gene --genomeDir /gpfs/fs2/scratch/ysun81_lab/GRCh38_Covid19/star/ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outSAMattributes CR UR CY UY CB UB --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 9 --readFilesCommand "gzip -cd" --readFilesIn /gpfs/fs2/scratch/ysun81_lab/Nora/COVID_DONOR18_fastq/COVID_DONOR18_S1_L001_R2_001.fastq.gz,/gpfs/fs2/scratch/ysun81_lab/Nora/COVID_DONOR18_fastq/COVID_DONOR18_S1_L001_R1_001.fastq.gz