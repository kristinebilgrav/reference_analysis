import sys
import pysam

#bamfile = pysam.AlignmentFile(sys.argv[2], "rb")
#bamfile = pysam.AlignmentFile(sys.argv[2], "rc", reference_filename='/sw/data/uppnex/ToolBox/hg38bundle/Homo_sapiens_assembly38.fasta' )
bamfile = pysam.AlignmentFile(sys.argv[2], "rc", reference_filename='/proj/sens2017106/nobackup/kristine/reference_stuff/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna')
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




def findMM(bam, regdict):
	mm={}
	flagmm={}
	myqual=[256,2048,2057, 2056, 2129, 2209, 2193]
	for chr in regdict:
		for start in regdict[chr]:
			end=regdict[chr][start]
			for read in bamfile.fetch(chr, start, end, until_eof=True):
				flag= int(read.flag)
				#print(flag)
				name=read.query_name
				if flag in myqual:
					if name not in flagmm:
						flagmm[name]=0
					flagmm[name]+=1
				if name not in mm:
					mm[name]=0
				mm[name]+=1
	return mm, flagmm


regionsdict=myregions(regions)
multimapping= findMM(bamfile, regionsdict)

total=0
for r in multimapping[1]:
	if multimapping[1][r] > 1:
		total += 1

perc1=len(multimapping[1])/len(multimapping[0])*100
print(str(bamfile), len(multimapping[0]), len(multimapping[1]), perc1, total )
