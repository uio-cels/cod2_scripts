#!/usr/bin/env python
"""
Python 2.7+ is needed to run this. Splits an assembly based on alignment of flanking
sequence of SNPs in a linkage map

"""


import numpy as np
#import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet

import sys, argparse, os, re
import math
from operator import itemgetter, attrgetter


def read_lg(lg_file):
	"""
	Parses the linkage map. Originally called CodMap_InfoFeb2014_til_Lex.txt with this
	format:
	Name	Comment_FINAL	LG	Order	Female	Male	Average		Meioses	
	
	In the first case, I'm only interested in the columns Name, LG and Order.
	"""
	
	lg_snps = {}
	
	#Throw away the header
	lg_file.readline()
	
	for line in lg_file:
		elements = line.split()
		lg_snps[elements[0]] = [elements[2], elements[3]]
	
	return lg_snps
	
def parse_psl_file(psl_file, head):
	"""
	Parses a psl created by BLAT and returns a dictionary with all the SNPs.
	SNPs with multiple alignments will have multiple values (a list of alignments), while
	SNPs with a unique alignment will just have one entry in the list.
	
	It also creates a score and percentage identity for each entry. And a field for whether
	or not it is the best alignment
	
	"""
	snps = {}
	
	#Throw away first 5 lines: Not needed when using -noHead with BLAT
	if not head:
		for i in range(0, 5):
			psl_file.readline()
	
	identity = 0.0
	score = 0
	best_align = False
	#Parses the psl file and adds the entries to the dictionary called snps		
	for line in psl_file:
		elements = line.split()
		score = calc_score(int(elements[0]), int(elements[2]), int(elements[1]), int(elements[4]), int(elements[6]))
		identity = 100.0 - calc_millibad(int(elements[12]), int(elements[11]), int(elements[16]), int(elements[15]), int(elements[0]), int(elements[2]), int(elements[1]), int(elements[4])) * 0.1
		elements.append(score)
		elements.append(identity)
		elements.append(best_align)
		#a bit of a hack, but I want to sort on the alignment against the longest target seq
		#and need that as an int. Not sure how else I'd do it.
		elements[14] = int(elements[14])

		if elements[9] in snps:
			snps[elements[9]].append(elements)
		else:
			snps[elements[9]] = [elements]

			
	return snps
	
	
def find_lgs_in_sequences(best_score_snps, lg_snps):
	"""
	Creates a data structure with sequence as key and a list of dictionaries with the 
	start and stop positions of the different linkage groups.
	"""
	
	sequences_in_lgs = {}
	sequences_with_snps_not_in_lg = set()
	
	for snp in best_score_snps:
		#Go through the equal scoring alignments
		for ali in best_score_snps[snp]:
			seq = ali[13]
			#If SNP in linkage map
			if ali[9] in lg_snps:
				#If sequence name already in sequences_in_lgs
				lg = lg_snps[ali[9]][0]
				if seq in sequences_in_lgs:
					#If linkage group exists already
					
					if lg in sequences_in_lgs[seq]:
						sequences_in_lgs[seq][lg][1].append({snp: [int(ali[15]), int(ali[16])]})
						if int(ali[15]) < int(sequences_in_lgs[seq][lg][0][0]):
							sequences_in_lgs[seq][lg][0][0] = int(ali[15]) 
						#End pos
						if int(ali[16]) > int(sequences_in_lgs[ali[13]][lg][0][1]):
							sequences_in_lgs[seq][lg][0][1] = int(ali[16])

					else:
						sequences_in_lgs[seq][lg] = [[int(ali[15]), int(ali[16])], [{snp: [int(ali[15]), int(ali[16])]}]]
						#sequences_in_lgs[seq][lg] = [[int(ali[15]), int(ali[16])], [snp]]
				#Create an entry for this sequence
				else:
					sequences_in_lgs[seq] = {lg: [[int(ali[15]), int(ali[16])], [{snp: [int(ali[15]), int(ali[16])]}]]}
					#sequences_in_lgs[seq] = {lg: [[int(ali[15]), int(ali[16])], [snp]]}
			#SNP is not in linkage map
			else:
				#Just add it to a set. Don't care that much about start and stop etc now
				sequences_with_snps_not_in_lg.add(seq)
					
	#Need to go through the sequences_with_snps_not_in_lg properly
	for seq in sequences_in_lgs:
		if seq in sequences_with_snps_not_in_lg:
			sequences_with_snps_not_in_lg.remove(seq)
				
	return sequences_in_lgs, sequences_with_snps_not_in_lg



#Rewritten from biopython's implementation: http://biopython.org/DIST/docs/api/Bio.SearchIO.BlatIO-pysrc.html
def calc_millibad(qend, qstart, tend, tstart, matches, repmatches, mismatches, qnuminsert): 
	# calculates millibad 
	# adapted from http://genome.ucsc.edu/FAQ/FAQblat.html#blat4 
	size_mul = 1 
	millibad = 0.0
	
	qali_size = size_mul * (qend - qstart) 
	tali_size = tend - tstart
	ali_size = min(qali_size, tali_size) 
	if ali_size <= 0: 
		return 0 

	size_dif = qali_size - tali_size 
	size_dif = 0 if size_dif < 0 else size_dif 
 
	total = size_mul * (matches + repmatches + mismatches) 
	if total != 0: 
		millibad = (1000 * (mismatches * size_mul + qnuminsert + 
			round(3 * math.log(1 + size_dif)))) / total 
	
	return millibad 
	
	
def calc_score(matches, repmatches, mismatches, qnuminsert, tnuminsert): 
	# calculates score
	# adapted from http://genome.ucsc.edu/FAQ/FAQblat.html#blat4 
	size_mul = 1 

	return (size_mul * (matches + (repmatches >> 1)) - size_mul * mismatches - qnuminsert - tnuminsert)


#Taken from Lex's /projects/cees/bin/scaffoldgap2bed.py
def get_gap_def(min_gap_len = 1):
    """
    Sets up the regular expression for splitting sequences into contigs and gaps
    example:
        get_gap_def(20)
    returns
         re.compile(r'(N{20,})')
    """
    return re.compile('(N{%i,})' % min_gap_len)
    
#Also taken from Lex's /projects/cees/bin/scaffoldgap2bed.py    
def get_coord(seq, gap_def):
    """
    Returns, for a given sequence,
    a list with on each line the coordinates of a gap or a contig
    when splitting on gaps of minimum lengths min_gap_len.
    """
    pos = 0
    coord = []
    # test for empty sequence
    if len(seq) < 1:
        return coord
    for subseq in split_seq(seq, gap_def):
        start = pos
        end = pos + len(subseq)
        coordline = [str(start), str(end), seq_type(subseq)]
        coord.append(coordline)
        pos = end
    return coord

def split_seq(seq, gap_def):
    """
    Splits a DNA sequence into contigs and gaps
    by splitting on min_gap_length.
    Gaps consisting of fewer N's than min_gap_len are ignored
    """
    for part in gap_def.split(seq.upper()):
        yield part


def seq_type(seq):
    """
    Determines whether a sequence consists of 'N's only 
    (i.e., represents a gap)
    """
    return 'gap' if set(seq.upper()) == {'N'} else 'bases'


def parse_fasta_file(fasta_file):
	
	fasta_entries = {}
		
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		fasta_entries[seq_record.id] = seq_record
		#fasta_entries[seq_record.id] = np.zeros(len(seq_record), dtype=np.int)
		
		
	return fasta_entries

def print_sequence(lg, seq, sequence):

	pwd = os.getcwd()
	
	cur = os.path.join(pwd, lg, "")
	#print cur
	#Create the directory if not existing
	if not os.path.isdir(cur):
		os.mkdir(cur)
		
	#print cur + prefix + "_" + lg + ".fasta"
		
	with open(cur + prefix + "_" + lg + ".fasta", "a") as output_fasta:
		SeqIO.write(sequence, output_fasta, "fasta")
		
	with open(cur + prefix + "_" + lg + ".bed", "a") as output_record:
		output_record.write("%s\t0\t%d\n" % (seq, len(sequence)))
		
	with open(cur + prefix + "_" + lg + "_gaps.bed", "a") as output_gaps:
		#output the gaps in the sequence also. 
		for line in get_coord(str(sequence), gap_def):
			if line[2] =='gap':
				# add sequence identifier at the beginning
				line.insert(0, seq)
				output_gaps.write('\t'.join(line))
				output_gaps.write('\n')


def output_sequences(sequences_in_lgs, sequences_with_snps_not_in_lg, prefix, fasta_entries, lg_snps):

	for seq in fasta_entries:
		if seq in sequences_in_lgs:
			lgs = sequences_in_lgs[seq].keys()
			#Easy, just one LG connected with the sequence
			if len(lgs) == 1:
				lg = lgs[0]
				#print lg
				#print fasta_entries[seq]
				print_sequence(lg, seq, fasta_entries[seq])
			else:
				#Create something that has the uninterrupted interval of a particular LG in it
				#partition_seq = {}
				
				#Need a sorted list of the SNPs in some way
				
				all_snps = []
				#Create a list of all SNPs connected with the sequence:
				for lg in lgs:
					#print(sequences_in_lgs[seq][lg][1])
					for snp in sequences_in_lgs[seq][lg][1]:
						all_snps.append(snp)
						
				#print all_snps	
				sorted_snps = sorted(all_snps, key=lambda x: x.values())
				
				#print sorted_snps
				#current_lg is set to the LG of the first SNP
				start_pos = 0
				end_pos = 0
				#Inelegant
				first_snp = True
				lg = ''
				for snp_dict in sorted_snps:
					snp = snp_dict.keys()[0]
					lg = lg_snps[snp][0]
					if first_snp:
						current_lg = lg
						first_snp = False
					if lg == current_lg:
						#-1 because python is 0 based, but these alignments are 1 based
						end_pos = snp_dict[snp][1] - 1
					#Change of LG
					else:
						#Have to create a SeqRecord
						name = seq + '_%d-%d' % (start_pos, end_pos)
						print_sequence(lg, name, SeqRecord(Seq(str(fasta_entries[seq].seq[start_pos:end_pos]), Alphabet.single_letter_alphabet), id=name, name=name, description=name))
						#The stuff between the alignments of the SNPs flanking sequences is 
						#defined as trash
						start_pos = snp_dict[snp][0] - 1
						name = seq + '_trash_%d-%d' % (end_pos+1, start_pos-1)
						print_sequence('trash', name, SeqRecord(Seq(str(fasta_entries[seq].seq[end_pos+1:start_pos-1]), Alphabet.single_letter_alphabet), id=name, name=name, description=name))
						end_pos = snp_dict[snp][1] - 1
						current_lg = lg
				#I guess the rest now is to the end of the sequence
				end_pos = len(fasta_entries[seq]) - 1
				name = seq + '_%d-%d' % (start_pos, end_pos)
				print_sequence(lg, name, SeqRecord(Seq(str(fasta_entries[seq].seq[start_pos:end_pos]), Alphabet.single_letter_alphabet), id=name, name=name, description=name))
				
								
		#Sequences not in LG, but has SNPs
		elif seq in sequences_with_snps_not_in_lg:
			print_sequence('no_lg_snps', seq, fasta_entries[seq])
		
		#No SNPs connected with sequence
		else:
			print_sequence('no_lg_no_snps', seq, fasta_entries[seq])
			
					#print pos

if __name__ == '__main__':
	
	
	parser = argparse.ArgumentParser(description=
	'Splits an assembly into linkage groups, sequence with unplaced SNPs and sequence without SNPs. Also, trash sequence.')
	parser.add_argument('-f', '--asm', action='store', help='fasta file for assembly', required=True, type=argparse.FileType('r'))
	parser.add_argument('-p', '--psl', action='store', help='BLAT psl file with SNPs mapped to assembly', required=True, type=argparse.FileType('r'))
	parser.add_argument('-l', '--lg',  action='store', help='linkage map file', required=False, type=argparse.FileType('r'))
	parser.add_argument('-s', '--score',  action='store', help='fraction score', required=False, type=float, default=1.00)
	parser.add_argument('-i', '--ident',  action='store', help='only work with alignments with higher than this identity', required=False, type=float, default=95.0)
	parser.add_argument('-g', '--length',  action='store', help='alignment length has to be more than this fraction of total flank sequence length', required=False, type=float, default=0.5)
	parser.add_argument('-e', '--ead',  action='store_false', help='BLAT ran with -noHead', required=False, default=True)
	parser.add_argument('-m', '--multiple', action='store_true', help='do not use SNPs mapping multiple times for splitting', required=False, default=False)
			
	args = parser.parse_args()
	
	#Define the gap sizes for later
	gap_def = get_gap_def(1)
	
	#Define the prefix
	prefix = '.'.join(args.asm.name.split('.')[0:-1])
	
	#Read in the fasta file
	fasta_entries = parse_fasta_file(args.asm)
	
	#Read in the linkage map, 9355 SNPs
	lg_snps = read_lg(args.lg)	
	
	#Get the mapped SNP flanks from the BLAT mapping
	snps = parse_psl_file(args.psl, args.ead)
	
	
	#Only interested in the best mapping flanking sequence of the SNPs. If there are several, save those too
	best_score_snps = {}
	not_multiple_mapping_snps = {}
	
	#Best score for all, plus all other alignments within the same score
	for snp in snps:
		#Sorting by score, then identity and last length of target sequence
		sorted_score_list = sorted(snps[snp], key=itemgetter(21, 22, 14), reverse=True)
		#Always add the best scoring alignment, and the one against the longest sequence
		#sorted_score_list[0][23] = True
		#Alignment length has to be more than 50 % of the flanking sequences
		if ((sorted_score_list[0][22] >= args.ident) and ((float(sorted_score_list[0][0])/float(sorted_score_list[0][10])) >= args.length)):
			best_score_snps[snp] = [sorted_score_list[0]]
			not_multiple_mapping_snps[snp] = [sorted_score_list[0]]

		for j in sorted_score_list[1:]:
			#By default equal, but has to be above args.score
			if (((float(j[21])/sorted_score_list[0][21]) >= args.score) and (j[22] >= args.ident) and ((float(j[0])/float(j[10])) >= args.length)):
				best_score_snps[snp].append(j)
				#If multiple mapping SNP, remove it
				not_multiple_mapping_snps.pop(snp, None)
	
	
	#Then the hard part, split the sequences based on LG.
	#I want a data structure with sequence name as key, and a list of LG: pos1-pos2 as values.
	
	if args.multiple:
		sequences_in_lgs, sequences_with_snps_not_in_lg = find_lgs_in_sequences(not_multiple_mapping_snps, lg_snps)
	else:	
		sequences_in_lgs, sequences_with_snps_not_in_lg = find_lgs_in_sequences(best_score_snps, lg_snps)
	
	#print (len(sequences_in_lgs))
	#print (len(sequences_with_snps_not_in_lg))
	
	output_sequences(sequences_in_lgs, sequences_with_snps_not_in_lg, prefix, fasta_entries, lg_snps)
	

	