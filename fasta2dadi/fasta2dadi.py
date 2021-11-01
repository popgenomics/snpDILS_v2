#!/usr/bin/python
#################################################################################################################################
#################################################################################################################################
#####                                                                                                                       #####
#####    This file is part of Demographic Inferences with Linked Selection : DILS.                                          #####
#####                                                                                                                       #####   
#####    DILS is free software: you can redistribute it and/or modify                                                       #####
#####    it under the terms of the GNU General Public License as published by                                               #####
#####    the Free Software Foundation, either version 3 of the License, or                                                  #####
#####    (at your option) any later version.                                                                                #####
#####                                                                                                                       #####    
#####    DILS is distributed in the hope that it will be useful,                                                            #####
#####    but WITHOUT ANY WARRANTY; without even the implied warranty of                                                     #####
#####    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                      #####
#####    GNU General Public License for more details.                                                                       #####
#####                                                                                                                       #####    
#####    You should have received a copy of the GNU General Public License                                                  #####
#####    along with DILS.  If not, see <https://www.gnu.org/licenses/>.                                                     #####
#####                                                                                                                       #####    
#####    Please send bugreports with examples or suggestions to                                                             #####
#####    camille.roux@univ-lille.fr                                                                                         #####
#####                                                                                                                       #####    
#####    Or write a post on https://groups.google.com/forum/#!forum/dils---demographic-inferences-with-linked-selection     #####
#####                                                                                                                       #####
#################################################################################################################################
#################################################################################################################################

import sys
import os
from math import ceil
import random
import getopt

def main(argv):
	helpmessage = 'test.py -i <fastaFile> -o <outputfile> -s <outgroup> -r <region>\n'
	helpmessage += '    <fastaFile> name of the input file in fasta format\n'
	helpmessage += '    <outputfile> name of the written outputfile containing the snp data used for DaDi\n'
	helpmessage += '    <outgroup> species name of the outgroup\n'
	helpmessage += '    <region> coding or noncoding\n\n'
	helpmessage += 'example: pypy3 test.py -i mytilus_renamed.fas -o coding_noOut_mussels.out -r coding # only 3rd coding position, without outgroup\n'
	helpmessage += 'example: pypy3 test.py -i mytilus_renamed.fas -o coding_out_mussels.out -r coding -s Lighthouse # only 3rd coding position, using Lighthouse as outgroup\n'
	helpmessage += 'example: pypy3 test.py -i mytilus_renamed.fas -o nonCoding_noOut_mussels.out -r noncoding # all positions, without outgroup\n'
	helpmessage += 'example: pypy3 test.py -i mytilus_renamed.fas -o nonCoding_out_mussels.out -r noncoding -s Lighthouse # all positions, using Lighthouse as outgroup\n\n'
	helpmessage += 'python3 or pypy3 can be used without modifying the code, but pypy3 much faster.'
	
	fastaFile = ''
	outputfile = ''
	outgroup = ''
	region = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:s:r:",["ifile=", "ofile=", "outgroup=", "region="])
	except getopt.GetoptError:
		print(helpmessage)
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print(helpmessage)
			sys.exit()
		elif opt in ("-i", "--ifile"):
			fastaFile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
		elif opt in ("-s", "--outgroup"):
			outgroup = arg
		elif opt in ("-r", "--region"):
			region = arg
	arguments = {}
	arguments['fastaFile'] = fastaFile
	arguments['outputfile'] = outputfile 
	arguments['outgroup'] = outgroup 
	arguments['region'] = region
	return(arguments) 

if __name__ == "__main__":
	arguments = main(sys.argv[1:])
	fastaFile = arguments['fastaFile']
	outputfile = arguments['outputfile']
	outgroup = arguments['outgroup']
	region = arguments['region']
	print('Input file is {fastaFile}'.format(fastaFile=fastaFile))
	print('Output file is {outputfile}'.format(outputfile=outputfile))
	print('Region is {region}'.format(region=region))
	print('Outgroup {outgroup}\n'.format(outgroup=outgroup))

def coloredSeq(seq):
	# print sequences with the standard color code
	seq = seq.replace("A", '\x1b[5;30;41m' + 'A' + '\x1b[0m')
	seq = seq.replace("T", '\x1b[5;30;44m' + 'T' + '\x1b[0m')
	seq = seq.replace("G", '\x1b[6;30;43m' + 'G' + '\x1b[0m')
	seq = seq.replace("C", '\x1b[5;30;42m' + 'C' + '\x1b[0m')
	return(seq)

def trunc2triplets(size):
	# trunc a value to its closest and smaller multiple of 3
	size = int(size)
	for i in range(2):
		if size%3 != 0:
			size -= 1
	return(size)

# nN = number of non-synonymous sites in the codon i: example for CGG -> nN = 2/3 + 3/3 + 0/3
# nS = number of synonymous sites in the codon i: example for CGG -> n> = 1/3 + 0/3 + 3/3
codonTable = {'AAA': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAC': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAG': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAT': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ACA': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACC': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACG': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACT': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'AGA': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGC': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AGG': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGT': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ATA': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATC': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATG': {'aa': 'M', 'nN': 3.0, 'nS': 0.0}, 'ATT': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'CAA': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAC': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAG': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAT': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CCA': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCC': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCG': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCT': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CGA': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGC': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CGG': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGT': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CTA': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTC': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'CTG': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTT': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'GAA': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAC': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAG': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAT': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GCA': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCC': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCG': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCT': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GGA': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGC': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGG': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGT': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GTA': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTC': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTG': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTT': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'TAC': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TAT': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TCA': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCC': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCG': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCT': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TGC': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TGG': {'aa': 'W', 'nN': 3.0, 'nS': 0.0}, 'TGT': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTA': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTC': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTG': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTT': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}}

def fasta2list(fastaFile, outgroup, region):
	fasta = open(fastaFile).readlines()
	seqName = [x.split(" ")[0].rstrip().replace('>','') for x in fasta if x[0] == '>']
	if region=="coding":
		seq = ''.join([x.rstrip()[2::3] if x[0]!='>' else '@' for x in fasta])[1:].split('@')
	else:
		seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')
	# get outgroup (0) or ingroup (1)
	# if no outgroup: only 1 values
	groups = [ 1 if i.split('|')[1]==outgroup else 0 for i in seqName ]
	species = [ i.split('|')[1] for i in seqName ]
	loci = [ i.split('|')[0] for i in seqName ]
	
	# returns final object
	res = {}
	for i in range(len(loci)):
		loc_tmp = loci[i]
		if loc_tmp not in res:
			res[loc_tmp] = {}
			res[loc_tmp]['ingroup'] = {}
			res[loc_tmp]['outgroup'] = {}
			res[loc_tmp]['ingroup']['species'] = []
			res[loc_tmp]['ingroup']['seq'] = []
			res[loc_tmp]['outgroup']['species'] = []
			res[loc_tmp]['outgroup']['seq'] = []
		
		if species[i]==outgroup:
			group_tmp = 'outgroup'
		else:
			group_tmp = 'ingroup'
		res[loc_tmp][group_tmp]['species'].append(species[i])
		res[loc_tmp][group_tmp]['seq'].append(seq[i])
#	return({'seq': seq, 'id': seqName})
	return(res)

def getConsensus(align):
	# align[locus]['seq', 'id']
	# L[locus] = n nucleotides
	L_locus = len(align[0])
	consensus='' # consensus = sequence
	for pos in range(L_locus):
		position = []
		for sequence in align:
			base = sequence[pos]
			if base != 'N':
				position.append(base)
		if len(position) == 0:
			consensus += 'N'
		else:
			consensus += random.sample(position, 1)[0]
	return(consensus)

def getSNP(locus):
	# returns list of polymorphic positions within the ingroup
	bases = ['A', 'T', 'C', 'G']
	# locus = one locus found in the object returned by fasta2list
	L = len(locus['ingroup']['seq'][0])
	SNPs = []
	allele1 = []
	allele2 = []
	for pos_tmp in range(L):
		align = []
		
		for ind_tmp in locus['ingroup']['seq']:
			base_tmp = ind_tmp[pos_tmp]
			if base_tmp in bases:
				align.append(base_tmp)
		alleles = list(set(align))
		
		if len(alleles) == 2:
			SNPs.append(pos_tmp)
			allele1.append(alleles[0])
			allele2.append(alleles[1])
	res = {}
	res['positions'] = SNPs
	res['allele1'] = allele1 
	res['allele2'] = allele2 
	return(res)

def getSpecies(align):
	list_species = []
	for locus in align:
		list_species_tmp = list(set(align[locus]['ingroup']['species']))
		for i in list_species_tmp:
			if i not in list_species:
				list_species.append(i)
	list_species.sort()
	return(list_species)

def countAlleles(locus, list_species, consensus, outgroup, region, locus_name):
	# scalar
	if region=="coding":
		scalar1 = 3
		scalar2 = 2
	else:
		scalar1 = 1
		scalar2 = 0
	
	# get outgroup
	if outgroup != '':
		consensus_outgroup = getConsensus(locus['outgroup']['seq'])
	else:
		consensus_outgroup = '-'*len(consensus)
	
	# get SNPs
	SNPs = getSNP(locus=locus)
	nSNP = len(SNPs['positions'])
	
	if nSNP==0:
		return('')
	else:
		# ingroup
		ingroup = {}
		for i in range( len(locus['ingroup']['seq'])):
			if locus['ingroup']['species'][i] not in ingroup:
				ingroup[locus['ingroup']['species'][i]] = []
			ingroup[ locus['ingroup']['species'][i] ].append( locus['ingroup']['seq'][i] )

		# loop over the SNPs
		res = ''
		for snp_tmp in range(nSNP):
			pos = SNPs['positions'][snp_tmp]
			Allele1 = consensus[pos]
			
			### TMP ###
#			print(coloredSeq("".join([ i[pos] for i in locus['ingroup']['seq'] ])))
			###########
			
			if SNPs['allele1'][snp_tmp] == Allele1:
				Allele2 = SNPs['allele2'][snp_tmp]
			else:
				Allele2 = SNPs['allele1'][snp_tmp]
			
			if region=="coding":
				Ing = '-{0}-'.format(Allele1)
				Out = '-{0}-'.format(consensus_outgroup[pos])
			else:
				if pos>0:
					pos0 = consensus[pos-1]
					pos0_outgroup = consensus_outgroup[pos-1]
				else:
					pos0 = '-'
					pos0_outgroup = '-'
				
				if pos<(len(consensus)-1):
					pos2 = consensus[pos+1]
					pos2_outgroup = consensus_outgroup[pos+1]
				else:
					pos2 = '-'
					pos2_outgroup = '-'
				
				Ing = '{pos0}{pos1}{pos2}'.format(pos0=pos0, pos1=Allele1, pos2=pos2)
				Out = '{pos0}{pos1}{pos2}'.format(pos0=pos0_outgroup, pos1=consensus_outgroup[pos], pos2=pos2_outgroup)
			
			# Allele 1
			res += "{Ing}\t{Out}\t{Allele1}".format(Ing=Ing, Out=Out, Allele1=Allele1)
			for species in list_species:
				if species in ingroup:
					res += "\t{nAllele1}".format(nAllele1="".join([ i[pos] for i in ingroup[species] ]).count(Allele1))
				else:
					res += "\t0"
			
			# Allele 2
			res += "\t{Allele2}".format(Allele2=Allele2)
			for species in list_species:
				if species in ingroup:
					res += "\t{nAllele2}".format(nAllele2="".join([ i[pos] for i in ingroup[species] ]).count(Allele2))
				else:
					res += "\t0"
			res += "\t{Gene}\t{Position}\n".format(Gene=locus_name, Position=pos*scalar1 + scalar2)
	return(res)

# get the data from fasta file
align = fasta2list(fastaFile=fastaFile, outgroup=outgroup, region=region)

# get the names of species located in the file
list_species = getSpecies(align)

# format the header
header = 'Ing\tOut\tAllele1\t'
for species in list_species:
	header += "{species}\t".format(species=species)
header += "Allele2"
for species in list_species:
	header += "\t{species}".format(species=species)
header += "\tGene\tPosition\n"

# produces the outfile
outfile = open(outputfile, "w")
outfile.write(header)

for locus_name in align:
	consensus_ingroup = getConsensus(align=align[locus_name]['ingroup']['seq'])
#	countAlleles(locus=align[locus], list_species=list_species, consensus=consensus_ingroup)
	res = countAlleles(locus=align[locus_name], list_species=list_species, consensus=consensus_ingroup, outgroup=outgroup, region=region, locus_name=locus_name)
	outfile.write(res)
outfile.close()


