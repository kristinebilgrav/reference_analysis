import sys
import pysam

#bamfile = pysam.AlignmentFile(sys.argv[2], "rb")
bamfile = pysam.AlignmentFile(sys.argv[2], "rc", reference_filename='/sw/data/uppnex/ToolBox/hg38bundle/Homo_sapiens_assembly38.fasta' )
#bamfile = pysam.AlignmentFile(sys.argv[2], "rc", reference_filename='/proj/sens2017106/nobackup/kristine/reference_stuff/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna')
regions=sys.argv[1]

def myregions(file):
	reg= {}
	for line in open(file):
		chr=line.split(':')[0]
		if chr == 'chrMT':
			continue
		start=int(line.split(':')[-1].split('-')[0])
		end=int(line.rstrip('\n').split(':')[-1].split('-')[1])
		if chr not in reg:
			reg[chr]={}

		if start not in reg[chr]:
			reg[chr][start]=end

	return reg




def findMQ(bam, regdict):
	mq=[]
	for chr in regdict:
		for start in regdict[chr]:
			end=regdict[chr][start]
			for read in bamfile.fetch(chr, start, end, until_eof=True):
				mapq= str(read.mapping_quality)
				mq.append(mapq)
				print('\t'.join([str(chr), str(start), str(end), str(mapq)]))
	return mq


regionsdict=myregions(regions)
mappingquality= findMQ(bamfile, regionsdict)
#print('\n'.join(mappingquality))
