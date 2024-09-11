import sys
import statistics

"""
find Mbp of DRR sequence that have low-medium-high coverage
AND
find amount average coverage across gap regions
AND
identify variable regions

INPUT: cov file from TIDDIT
"""

mychr = list(range(0,23))
mychr = [str(i) for i in mychr]
mychr.append('X')
mychr.append('Y')

genome='T2T' 
#genome ='GRCh38'
#genome='GRCh37'
#genome='Chimpanzee' 
#genome='Bonobo'

def collect_cov(file): # collects info from file

	covdict = {}
	bindict = {}
	for line in open(file):
		tabline = line.rstrip('\n').split()
		chr=tabline[0]

		#if chr not in mychr:
		#	continue

		start = tabline[1]
		end = tabline[2]
		cov = float(tabline[3])

		if chr not in covdict:
			covdict[chr] = {}
			bindict[chr]={}
		if start not in covdict[chr]:
			covdict[chr][start] = []
			bindict[chr][start] = end

		covdict[chr][start].append(cov)
		lst=[chr, start, end, str(cov)]
		#print('\t'.join(lst))
	return covdict, bindict


def count_presence(covdict, bindict):

	for chr in covdict:
		for s in covdict[chr]:

			average = sum(covdict[chr][s])/len(covdict[chr][s])
			covlist = covdict[chr][s]

			#count presence
			lowcov = len([i for i in covlist if i <= 8])
			highcov = len([i for i in covlist if i >8 and i <100 ])

			if highcov >= 5:
				status = 'Common'

			else:
				status = 'Absent'


			lst = [chr, s, bindict[chr][s], str(average), genome, status]
			print('\t'.join(lst))
	return

def find_amount(file):
	amount=[]
	for line in open(file):
		am = int(line.rstrip('\n'))
		amount.append(am)
	mean = sum(amount)/len(amount)
	median = statistics.median(amount)

	meanMBP=mean*10000/1000000
	medianMBP=median*10000/1000000

	print('median', medianMBP)
	print('mean', meanMBP)

find_amount(sys.argv[1]) #find amount of DRRs
all_cov= collect_cov(sys.argv[1]) #coverage across all regions
count_presence(all_cov[0], all_cov[1])

