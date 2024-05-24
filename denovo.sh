# YAK
yak count -k31 -b37 -t16 -o $2.yak <fasta>

# Hifiasm
hifiasm -o $2 -t 16 $fastq

# GFA to FA
awk '/^S/{print ">"$2;print $3}' hap.gfa > $2.fa

# align FA
minimap2 -R "@RG\tID:$2\tSM:$2" -a -t 16 --MD -x map-pb -asm5 -Y -y $ref $2.fa | samtools view -Sbh - | samtools sort -m 4G -@16 - > $2.bam
samtools index $2.bam

# SVIM
svim-asm diploid $PWD <hap1> <hap2> $ref --sample $3

# Quast
python quast <hap1> <hap2>  -o $3.quast -t 16 -r $ref
