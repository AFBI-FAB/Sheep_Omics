#!/bin/bash

SECONDS=0

# Software Requirements
# FastQC v0.12.1 
# cutadapt version 1.15
# TrimmomaticPE
# HISAT2 version 2.1.0
# featureCounts v2.0.1
# StringTie 2.2.1

# Path to the working directory where fastq files are present. In this example, the samples are present in a sub-directory of sheep named 24_samples. Hence, the current working directory is set.
workdir="/mnt/sda1/00_fastq/Sheep/"
cd $workdir

# Create directories to store outputs and give full permission to those folders for the output files to be saved.
echo "Creating directories to store outputs!"
mkdir 1.fastqc
mkdir 2.trimmomatic
mkdir 3.fastqc.after.trimmomatic
mkdir 4.hisat2
mkdir 5.featurecounts
mkdir logs

# To get the total number of reads from all the fastq files, run the command below. For this, make sure you have installed seqkit. If its not installed, run 
# "conda install -c bioconda seqkit" or refer to https://github.com/shenwei356/seqkit

seqkit stats -To stats.tsv *.fastq.gz

# Step 1: Quality Control
# Run FASTQC (2 minutes per sample)

echo "Starting FastQC"
for file in 24_Samples/*R1_001.fastq.gz; do
    prefix="${file%R1_001.fastq,gz}"
    reverse="${file%R1_001.fastq.gz}R2_001.fastq.gz"
    fastqc -1 $file -2 $reverse -t 15 -o 1.fastqc
done
echo "FastQC finished succesfully!"


# Step 2: Adapter trimming
# Run Trimmomatic

echo "Starting Trimmomatic"
threads="8" #higher the number faster the speed; but make sure you have enough CPU support
for R1 in 24_Samples/*_R1_001.fastq.gz ; do
    R2="${R1%_R1_001.fastq.gz}_R2_001.fastq.gz"
    sample=$(echo $R1|sed 's/_R1_001.fastq.gz//'|sed 's/24_Samples\///'); 
    TrimmomaticPE -threads $threads "$R1" "$R2" 2.trimmomatic/"${sample}_R1_paired.fastq.gz" 2.trimmomatic/"${sample}_R1_unpaired.fastq.gz"  \
    2.trimmomatic/"${sample}_R2_paired.fastq.gz" 2.trimmomatic/"${sample}_R2_unpaired.fastq.gz" ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 
done
echo "Trimmomatic run finished succesfully!"

# Step3: Fastqc on adapter trimmed reads

echo "Starting FastQC on trimmed reads"
for file in 2.trimmomatic/*_R1_paired.fastq.gz; do
    prefix="${file%_R1_paired.fastq.gz}"
    reverse="${file%_R1_paired.fastq.gz}_R2_paired.fastq.gz"
    fastqc -1 $file -2 $reverse -t 15 -o 3.fastqc.after.trimmomatic
done
echo "FastQC on trimmed reads finished succesfully!"

# Step 4: Alignment to the reference genome using HISAT2
# 4.1: Build reference genome index. The genome should be in unzipped format or else it wont work

echo "Starting to build the genome index files"
refgenome="Reference_genome/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna"
hisat2-build  $refgenome $refgenome
chmod a+x /mnt/sda1/00_fastq/Sheep/Reference_genome/*.ht2

# 4.2 Now Align to the reference genome

echo "Index files generated. Now performing the alignment. This might take a while.."
threads=8
for R1 in 2.trimmomatic/*_R1_paired.fastq.gz; do
    R2=$(echo $R1| sed 's/_R1_/_R2_/'); 
    sample=$(echo $R1|sed 's/_R1_paired.fastq.gz//'|sed 's/2.trimmomatic\///'); 
    echo hisat2 -q --summary-file 4.hisat2/$sample.summary.txt --threads $threads \
    -x /mnt/sda1/00_fastq/Sheep/Reference_genome/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna \
    -1 $R1 -2 $R2 | samtools sort -O BAM | tee 4.hisat2/$sample.bam | samtools index - 4.hisat2/$sample.bam.bai;
done
echo "Alignment finished succesfully!" 

# To get mapped reads from the bam file,
# samtools view -b -F 4 7220.bam > 7220.mapped.bam

# Step 5 Generate read counts matrix using featureCounts
# When you want to analyze the data for differential gene expression analysis, it would be convenient to have counts for all samples in a single file (gene count matrix).

echo "Generating the read counts matrix.."
gtffile="/mnt/sda1/00_fastq/Sheep/Reference_genome/genomic.gtf"

featureCounts -T 8 -t 'gene' -g 'gene_id' -f -a $gtffile -o 5.featurecounts/Lambs.featurecounts.hisat2 4.hisat2/*.bam

# -T specfies the number of threads to be used
# -t specifies the feature type. By default it is 'exon', but in the above command I have specified 'gene'
# -g specifies the attribute type in GTF annotation. By default it is 'gene_id'
# -f specifies to perform the read counting at the feature level (eg. counting reads for exons rather than genes)
# -a specifies the annotation file
# -o specifies the name of the output file

## Since the featureCounts output has additional columns with information about genomic coordinates, gene length etc., 
## we can use the cut command to select only those columns that you are interested in. Columns 1 and sample wise counts columns

cut -f1,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32 5.featurecounts/Lambs.featurecounts.hisat2 > 5.featurecounts/Lambs.featurecounts.hisat2.Rmatrix
echo "Read count matrix generated!"
echo "Have a look at the counts matrix. You may change the column names if needed.."
echo "To ready this text file (count matrix) for the next step of differential gene expression analysis, you will need to clean it up further by removing the first header line, and modifying the column names (headers) to simpler, smaller sampleIDs.
echo "An example of tidying the file is as shown below:
## vim 5.featurecounts/Lambs.featurecounts.hisat2.Rmatrix
##dd
##gg
##:%s/4.hisat2\///g

# After these step, try to remove extra texts in the column names (eg..bam in the sample names)

# Step 6 Multiqc Report generation
echo "Generating Multiqc report.."
multiqc -o 9.multiqc .
echo "Report generated succesfully!"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
