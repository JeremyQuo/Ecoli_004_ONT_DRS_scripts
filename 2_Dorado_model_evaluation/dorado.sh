# dorado
dorado basecaller rna004_130bps_sup@v5.1.0 --modified-bases pseU \
--reference ecoli.fa  pod5/ -r > basecall_pseU.bam
modkit pileup basecall_pseU.bam basecall_pseU.bed

dorado basecaller rna004_130bps_sup@v5.1.0 --modified-bases inosine_m6A \
--reference ecoli.fa  pod5/ -r > basecall_inosine_m6A.bam
modkit pileup basecall_inosine_m6A.bam basecall_inosine_m6A.bed

dorado basecaller rna004_130bps_sup@v5.1.0 --modified-bases pseU \
--reference ecoli.fa  pod5/ -r > basecall_m5C.bam
modkit pileup basecall_m5C.bam basecall_m5C.bed

