# Basecall
dorado basecaller rna004_130bps_sup@v5.0.0   pod5/ -r > basecall.bam
samtools bam2fq basecall.bam>final.fastq

# Alignment
minimap2 --MD -t 64 -ax map-ont ecoli.fa final.fastq  | samtools view -hbS -F 260 - | samtools sort -@ 64 -o genomic.bam
samtools index genomic.bam
samtools depth -a genomic.bam>coverage.txt

# Giraffe commands
giraffe estimate --read final.fastq
giraffe observe --read final.fastq --ref ecoli.fa

# Seqkit
seqkit stat -a final.fastq

# Epinano
python Epinano_Variants.py -c 8 -r ecoli.fa -b genomic.bam