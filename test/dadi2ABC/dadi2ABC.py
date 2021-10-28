#!/usr/bin/python
from numpy.random import choice
import sys
import os

infileName = sys.argv[1] # dadi's input file
nA = int(sys.argv[2]) # number of sequences from population A
nB = int(sys.argv[3]) # number of sequences from population B
use_outgroup = int(sys.argv[4]) # 0 to fold; 1 to unfold
binpath = sys.argv[5] # path of DILS directory containing all binarys

# function converting a recorded gene into a ms like output
def dadi2ms(gene, use_outgroup):
	align = {}
	positions = []
	nSNPs = gene['nSNPs']
	cnt = 0
	for pos in range(nSNPs):
		if pos not in gene['SNPs_to_remove']:
			outgroup = gene['outgroup'][pos][1]
			allele1 = '0'
			allele2 = '1'
			if use_outgroup == 1 and outgroup != '-':
				if outgroup == gene['allele2'][pos]:
					allele1 = '1'
					allele2 = '0'
			align[cnt] = allele1*gene['nAllele1_spA'][pos] + allele2*gene['nAllele2_spA'][pos] + allele1*gene['nAllele1_spB'][pos] + allele2*gene['nAllele2_spB'][pos]
			positions.append(gene['position'][pos])
			cnt += 1

	positions = '\t'.join([ str(i/(1.0*max(positions))) for i in positions ])
	ms = '//\nsegsites: {0}\npositions: {1}\n'.format(nSNPs - gene['nSNPs_to_remove'], positions)
	nInd = len(align[0])
	for ind in range(nInd):
		for pos in align:
			ms += '{0}'.format(align[pos][ind])
		ms += '\n'
	ms += '\n'
	return(ms)


# parse the DaDi' input file
genes = {}
infile = open(infileName, 'r')
header = infile.readline().strip().split('\t')
for line in infile:
	line = line.strip().split('\t')
	gene_tmp = line[8]
	position_tmp = int(line[9])
	outgroup_tmp = line[1]
	
	# allele 1	
	allele1_tmp = line[2]
	nAllele1_spA_tmp = int(line[3])
	nAllele1_spB_tmp = int(line[4])

	# allele 2	
	allele2_tmp = line[5]
	nAllele2_spA_tmp = int(line[6])
	nAllele2_spB_tmp = int(line[7])
	
	if gene_tmp not in genes:
		genes[gene_tmp] = {}
		genes[gene_tmp]['nSNPs'] = 0
		genes[gene_tmp]['position'] = []
		genes[gene_tmp]['allele1'] = []
		genes[gene_tmp]['allele2'] = []
		genes[gene_tmp]['outgroup'] = []
		genes[gene_tmp]['nAllele1_spA'] = []
		genes[gene_tmp]['nAllele2_spA'] = []
		genes[gene_tmp]['nAllele1_spB'] = []
		genes[gene_tmp]['nAllele2_spB'] = []
	genes[gene_tmp]['nSNPs'] += 1
	genes[gene_tmp]['position'].append(position_tmp)
	genes[gene_tmp]['allele1'].append(allele1_tmp)
	genes[gene_tmp]['allele2'].append(allele2_tmp)
	genes[gene_tmp]['outgroup'].append(outgroup_tmp)
	genes[gene_tmp]['nAllele1_spA'].append(nAllele1_spA_tmp)
	genes[gene_tmp]['nAllele2_spA'].append(nAllele2_spA_tmp)
	genes[gene_tmp]['nAllele1_spB'].append(nAllele1_spB_tmp)
	genes[gene_tmp]['nAllele2_spB'].append(nAllele2_spB_tmp)
infile.close()


# standardize all SNPs according to nA and nB
for gene_tmp in genes:
	genes[gene_tmp]['nSNPs_to_remove'] = 0
	genes[gene_tmp]['SNPs_to_remove'] = []
	for snp in range(genes[gene_tmp]['nSNPs']):
		nA_tmp = genes[gene_tmp]['nAllele1_spA'][snp] + genes[gene_tmp]['nAllele2_spA'][snp] 
		nB_tmp = genes[gene_tmp]['nAllele1_spB'][snp] + genes[gene_tmp]['nAllele2_spB'][snp]
		if nA_tmp < nA or nB_tmp < nB: # if not enough individuals in at least one population
			genes[gene_tmp]['nSNPs_to_remove'] += 1
			genes[gene_tmp]['SNPs_to_remove'].append(snp)
		else:
			if nA_tmp > nA:
				nAlleles_tmp = list(choice(['allele1', 'allele2'], nA, p=[genes[gene_tmp]['nAllele1_spA'][snp]/(1.0*nA_tmp), genes[gene_tmp]['nAllele2_spA'][snp]/(1.0*nA_tmp)], replace=True))
				genes[gene_tmp]['nAllele1_spA'][snp] = nAlleles_tmp.count('allele1')
				genes[gene_tmp]['nAllele2_spA'][snp] = nAlleles_tmp.count('allele2')
			if nB_tmp > nB:
				nAlleles_tmp = list(choice(['allele1', 'allele2'], nB, p=[genes[gene_tmp]['nAllele1_spB'][snp]/(1.0*nB_tmp), genes[gene_tmp]['nAllele2_spB'][snp]/(1.0*nB_tmp)], replace=True))
				genes[gene_tmp]['nAllele1_spB'][snp] = nAlleles_tmp.count('allele1')
				genes[gene_tmp]['nAllele2_spB'][snp] = nAlleles_tmp.count('allele2')
			
			# remove monomorphic position created by subsampling, i.e, when there is a single allele in both species
			if genes[gene_tmp]['nAllele1_spA'][snp] + genes[gene_tmp]['nAllele1_spB'][snp] == 0 or genes[gene_tmp]['nAllele2_spA'][snp] + genes[gene_tmp]['nAllele2_spB'][snp] == 0:
				genes[gene_tmp]['nSNPs_to_remove'] += 1
				genes[gene_tmp]['SNPs_to_remove'].append(snp)


# produce a bpfile
retainedLoci = []
line1 = [] # nSNPs
line2 = [] # nsamA
line3 = [] # nsamB
for gene_tmp in genes:
	nSNPs = genes[gene_tmp]['nSNPs'] - genes[gene_tmp]['nSNPs_to_remove']
	# remove genes with zero SNP
	if nSNPs > 0:
		retainedLoci.append(gene_tmp)
		line1.append(nSNPs)
		line2.append(nA)
		line3.append(nB)

outfile = open('bpfile', 'w')
outfile.write('# nGenes={0}; nSNPs={1}; nsamA={2}; nsamB={3}; input file={4}\n'.format(len(line1), sum(line1), nA, nB, infileName))
outfile.write('\t'.join( [ str(i) for i in line1 ] ) + '\n')
outfile.write('\t'.join( [ str(i) for i in line2 ] ) + '\n')
outfile.write('\t'.join( [ str(i) for i in line3 ] ) + '\n')
outfile.close()

# produce an ms like file
outfile = open('obs.ms', 'w')
for gene_tmp in retainedLoci:
	tmp = dadi2ms(genes[gene_tmp], use_outgroup)
	outfile.write(tmp)
outfile.close()

# execute mscalc
command = 'cat obs.ms | pypy {0}/mscalc_2pop_SFS.py {1}; mv ABCstat.txt obs_stat.txt; mv ABCjsfs.txt obs_sfs.txt'.format(binpath, use_outgroup)
os.system(command)

