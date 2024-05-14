

## Create BAM file ref-to-ref /ref1=ref ref2=asm ##

minimap2 -ax asm5 <ref1> <ref2> | samtools view -Sbh - > tmp.ref1.ref2
samtools sort -m 4G -@1 tmp.ref1.ref2 > ref1_ref2.bam
samtools index ref1_ref2.bam


## Run TIDDIT coverage module to find gaps ##

tiddit=container.sif
singularity exec --bind /dataset $tiddit tiddit --cov --bam ref1_ref2.bam -o  ref1_ref2.cov -z 10000

## Remove known gap regions in template ##

bedtools intersect -v -a ref1_ref2.cov -b gapsRef1.bed  > ref1_ref2.novelgap.bed

## Collect all gap regions : coverage of 0 across the bin ##
grep -P '\t0\.0' ref1_ref2.novelgap.bed >> ref1_ref2.novelgap.bed
