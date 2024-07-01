# Align sr

bwa mem -p -t 16 <ref> <fastq>

# Align and SV calling of linked-read
longranger wgs –id <id> --reference <ref> --fastq <fastq> --vcmode freebayes
