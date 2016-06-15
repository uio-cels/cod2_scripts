#!/usr/bin/env python
"""
Python 2.7+ is needed to run this.

parses blat .psl files and orders and orients scaffolds based on alignment of flanking 
sequences of SNPs.

Example of input linkage map file:
Name	Comment_FINAL	LG	Order		Female	Male	Average	Meioses
Gdist:574061_1677	SNP-0	gm-01	0		0	0	0	956
Gdist:483922_6360	SNP-0	gm-01	1		0	0.2	0.5	464
Gdist:483921_286	SNP-0	gm-01	2		0	0.2	0.5	1067


"""
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet


import sys, argparse
import math
from operator import itemgetter, attrgetter

translate_to_canadian_lgs = {
'gm-01': 'LG07',
'gm-02': 'LG12',
'gm-03': 'LG01',
'gm-04': 'LG16',
'gm-05': 'LG21',
'gm-06': 'LG03',
'gm-07': 'LG06',
'gm-08': 'LG11',
'gm-09': 'LG05',
'gm-10': 'LG08',
'gm-11': 'LG02',
'gm-12': 'LG15',
'gm-13': 'LG20',
'gm-14': 'LG14',
'gm-15': 'LG09',
'gm-16': 'LG10',
'gm-17': 'LG13',
'gm-18': 'LG04',
'gm-19': 'LG23',
'gm-20': 'LG18',
'gm-21': 'LG22',
'gm-22': 'LG19',
'gm-23': 'LG17'
}



def parse_fasta_file(fasta_file):
	
	fasta_entries = {}
		
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		fasta_entries[seq_record.id] = seq_record
				
	return fasta_entries


def read_lg(lg_file, exclude_snps):
	"""
	Parses the linkage map. In this
	format:
	Name	Comment_FINAL	LG	Order	Female	Male	Average		Meioses	
	
	In the first case, I'm only interested in the columns Name, LG and Order.
	"""
	
	lg_snps = {}
	
	#Throw away the header
	lg_file.readline()
	
	for line in lg_file:
		elements = line.split()
		if elements[0] not in exclude_snps:
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
	
	identity = 0
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
		elements[14]= int(elements[14])
		#elements.append(int(elements[14]))
		#print (elements[0], elements[21])
		if elements[9] in snps:
			snps[elements[9]].append(elements)
		else:
			snps[elements[9]] = [elements]
			
	return snps

#Creating the data structure sequences should be separate from this
def extract_longest_hit(snps, lg_snps):
	#Getting longest alignment for each SNP
	
	equal_length_ali = []
	#Sequences: seq_name is key, a list as values with length of sequence,
	#a list of SNPs and a dictionary of linkage groups with count as values
	sequences = {}
	unique = 0
	unique_lg = 0
	unique_other = 0
	multiple_best_hits = {}
	
	for snp in snps:
		
		#Get longest alignment, 9251 of them
		sorted_hit_list = sorted(snps[snp], key = lambda ele: int(ele[0]), reverse=True)
		#print (sorted_hit_list)
		#print(len(sorted_hit_list))
		longest_ali = sorted_hit_list[0]
		#Add all with equal alignment length to a list
		for ali in sorted_hit_list:
			if (int(ali[0]) == int(longest_ali[0])):
				equal_length_ali.append(ali)
	
		#Go through that list and create some datastructures and so
		for best in equal_length_ali:
			#If sequence name already in sequences
			if best[13] in sequences:
				sequences[best[13]][1].append(best[9])
				#if linkage group already in datastructure
				#  best[9] 
				if best[9] in lg_snps:
					#print(lg_snps[best[9]][0])
					#print(sequences[best[13]][2])
					if lg_snps[best[9]][0] in sequences[best[13]][2]:
						#add count of one SNP
						sequences[best[13]][2][lg_snps[best[9]][0]] += 1	
					else:
						sequences[best[13]][2][lg_snps[best[9]][0]] = 1
			
			else:
				#Don't care if SNP is not in linkage map
				if best[9] in lg_snps:
					#print ({lg_snps[best[9]][0]: 1})
					sequences[best[13]] = [best[14], [best[9]], {lg_snps[best[9]][0]: 1}]
	
			if (len(equal_length_ali) == 1):
					unique += 1
					if (best[9] in lg_snps):
						unique_lg += 1
					else:
						unique_other += 1
			else:
				if best[9] in multiple_best_hits:
					multiple_best_hits[best[9]].append(best[13])
				else:
					multiple_best_hits[best[9]] = [best[13]]
				
				
		equal_length_ali = []	
			#for equal in equal_length_ali.append:	
				
			#print(longest_ali)
	return 	sequences, unique, unique_lg, unique_other, multiple_best_hits
	

def summary_counts(sequences):

	#seq would not exist without being connected to a lg
	count_seq_lg = len(sequences)
	bases_seq_lg = 0
	count_seq_multi_lg = 0
	bases_seq_multi_lg = 0	
	conflict_count_seq_lg = 0
	conflict_bases_seq_lg = 0
	
	linked_to_multiple_lgs = {}

	sum_multi = 0

	for seq in sequences:
		#print(sequences[seq])
		#if sequences[seq]
		bases_seq_lg += int(sequences[seq][0])
		#if (len(sequences[seq][1]) > 1):
		for lg in sequences[seq][2]:
			sum_multi += sequences[seq][2][lg]
		if (sum_multi > 1):
			#print (sequences[seq])
			count_seq_multi_lg += 1
			bases_seq_multi_lg += int(sequences[seq][0])
		sum_multi = 0		
		if (len(sequences[seq][2]) > 1):
			conflict_count_seq_lg += 1
			conflict_bases_seq_lg += int(sequences[seq][0])
			linked_to_multiple_lgs[seq] = sequences[seq]
	
	return count_seq_lg, bases_seq_lg, count_seq_multi_lg, bases_seq_multi_lg, conflict_count_seq_lg, conflict_bases_seq_lg, linked_to_multiple_lgs 
		



def extract_hits(file_prefix, snps, lg_snps, ident, length):
	"""
	Creates some data structures that is more easily parseable than before. 
	Connects SNPs to sequences (scaffolds).
	"""
	#Just get all alignments of SNPs and filter before this point
	
	#equal_length_ali = []
	#Sequences: seq_name is key, a list as values with length of sequence,
	#a list of SNPs and a dictionary of linkage groups with count as values
	sequences = {}
	unique = 0
	unique_lg = 0
	unique_other = 0
	multiple_ali_hits = {}
	

	for snp in snps:	
		#sys.stderr.write("%s\n" % snp)		
		#Go through that list and create some datastructures and so
		for ali in snps[snp]:
			#sys.stderr.write("ali\n" )
			#If sequence name already in sequences
			#sys.stderr.write("%s\n" %ali)
			if ali[13] in sequences:
				
				if ali[9] in lg_snps:
					sequences[ali[13]][1].append([ali[9], int(ali[15]), int(ali[16]), int(ali[21]), int(ali[22]), ali[23]])
					#if linkage group already in datastructure
					if lg_snps[ali[9]][0] in sequences[ali[13]][2]:
						#add count of one SNP
						sequences[ali[13]][2][lg_snps[ali[9]][0]] += 1	
					else:
						sequences[ali[13]][2][lg_snps[ali[9]][0]] = 1
			
			else:
				#Don't care if SNP is not in linkage map
				if ali[9] in lg_snps:
					#print ({lg_snps[ali[9]][0]: 1})
					sequences[ali[13]] = [int(ali[14]), [[ali[9], int(ali[15]), int(ali[16]), int(ali[21]), int(ali[22]), ali[23]]], {lg_snps[ali[9]][0]: 1}]
	
			if (len(snps[snp]) == 1):
					unique += 1
					if (ali[9] in lg_snps):
						unique_lg += 1
					else:
						unique_other += 1
			else:
				if ali[9] in multiple_ali_hits:
					multiple_ali_hits[ali[9]].append(ali[13])
				else:	
					multiple_ali_hits[ali[9]] = [ali[13]]

	return 	sequences, unique, unique_lg, unique_other, multiple_ali_hits

def oles_method(file_prefix, snps, mapped_snps, lg_snps, ident, length, fasta_entries):
	

	sequences, unique, unique_lg, unique_other, multiple_best_hits = extract_hits(file_prefix, mapped_snps, lg_snps, ident, length)

	multi_lg = 0
	multi_no_lg = 0
	for snp in multiple_best_hits:
		if snp in lg_snps:
			multi_lg += 1
		else:
			multi_no_lg += 1

	#print sequences
	
	
	#Order the sequences according to linkage map below here. Presume they are split
	
	
	#LG is key, with SNP and order as values	
	lgs = {}
	
	seqs_linked_lg = {}
	#Can just add SeqRecords to each other, with 100N in between
		
	for snp in lg_snps:
		if lg_snps[snp][0] in lgs:
			lgs[lg_snps[snp][0]].append([snp, int(lg_snps[snp][1])])
		else:
			lgs[lg_snps[snp][0]] = [[snp, int(lg_snps[snp][1])]]
	
	#print("file_prefix: %s" % file_prefix)
	
	output_handle = open(file_prefix + "_linked.fasta", "w")

	sequences_output = set()

	for lg in sorted(lgs):
		#print lg
		#print (lgs[lg])
		sorted_snps_order = sorted(lgs[lg], key = lambda ele: ele[1])
		#print sorted_snps_order
		#I want LG, SNP, order, seq, seq len
		
		#create an empty SeqRecord
		lg_seq = SeqRecord('', id=lg)
		
		#100 Ns between the placed scaffolds:
		#"In the chromosome agp file, we prefer to use 100 to represent gaps of unknown length." NCBI
		#We linked the scaffolds onto chromosomes with a string of 100 Ns representing the gap between 2 adjacent scaffolds on the basis of a high-resolution genetic map.: http://www.nature.com/ng/journal/v46/n11/full/ng.3098.html#methods Common carp
			
		gap = SeqRecord('N'*100)
		
		#print (gap)
		
		print ("length sorted_snps_order: %d" % len(sorted_snps_order))
		
		for snp, order in sorted_snps_order:
			#first, does the sequence contain just 1 SNP? Just write it then.
			#second, if more than 1 SNP, is the orientation correct?
			
			print ("lg: %s, snp: %s, order: %d" % (lg, snp, order))	
			

			#Sequences: seq_name is key, a list as values with length of sequence,
			#a list of SNPs and a dictionary of linkage groups with count as values		
			if ((snp in mapped_snps) and (mapped_snps[snp][0][22] >= ident) and ((float(mapped_snps[snp][0][0])/float(mapped_snps[snp][0][10])) >= length)):
				
				seq = mapped_snps[snp][0][13]
			
				print ("sequence is %s" % seq)
				
				print("len(sequences[seq][1]): %d" % len(sequences[seq][1]))		
				
				if len(sequences[seq][1]) == 1:
					#only one SNP associated
					print ("Only one SNP: %s " % snp)
					if len(lg_seq) > 0:
						lg_seq = lg_seq + gap + fasta_entries[seq]
					else:
						lg_seq = fasta_entries[seq]		
				else:
					#if this sequence has not already been output
					if seq not in sequences_output:
					
						#more than one SNP associated, need to figure out orientation
					
						#first time we enter here, the current SNP should have the lowest value, could use that
						#start alignment pos of current snp is: mapped_snps[snp][0][15] 
						rev_comp = False
						for snip in sequences[seq][1]:
							print("snip: ", snip)
							print("snip[0]: ", snip[0])
							print("mapped_snps[snip[0]][0][15]: ", mapped_snps[snip[0]][0][15])
							print("mapped_snps[snp][0][15]: ", mapped_snps[snp][0][15])
							if mapped_snps[snip[0]][0][15] < mapped_snps[snp][0][15]:
								rev_comp = True
							
						if rev_comp:
							sequence_to_add = fasta_entries[seq].reverse_complement()
						else:
							sequence_to_add = fasta_entries[seq]
					
						if len(lg_seq) > 0:
					
							lg_seq = lg_seq + gap + sequence_to_add
						else:
							lg_seq = sequence_to_add
						
						sequences_output.add(seq)	
					#pass
					
				#print (lg_seq)
		
		#print(lg_seq)
		lg_seq.name = lg
		lg_seq.id = lg
		lg_seq.description = ''
		
		SeqIO.write(lg_seq, output_handle, "fasta")		

	#Print all non-linkage mapped linked sequences
	for seq in fasta_entries:
		if seq not in sequences:
			SeqIO.write(fasta_entries[seq], output_handle, "fasta")		

	output_handle.close()

#Rewritten from biopython's implementation: http://biopython.org/DIST/docs/api/Bio.SearchIO.BlatIO-pysrc.html
def calc_millibad(qend, qstart, tend, tstart, matches, repmatches, mismatches, qnuminsert): 
	# calculates millibad 
	# adapted from http://genome.ucsc.edu/FAQ/FAQblat.html#blat4 
	size_mul = 1 
	millibad = 0 
	
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

if __name__ == '__main__':
	
	
	parser = argparse.ArgumentParser(description=
	'The blat .psl file.')
	parser.add_argument('-p', '--psl', action='store', help='blat psl file', required=True, type=argparse.FileType('r'))
	parser.add_argument('-l', '--lg',  action='store', help='linkage map file', required=False, type=argparse.FileType('r'))
	parser.add_argument('-s', '--score',  action='store', help='fraction score', required=False, type=float, default=0.99)	
	parser.add_argument('-i', '--ident',  action='store', help='only output SNPs with higher than this percent identity', required=False, type=float, default=95.0)	
	parser.add_argument('-e', '--ead',  action='store_false', help='BLAT ran with -noHead', required=False, default=True)
	parser.add_argument('-g', '--length',  action='store', help='alignment length has to be more than this fraction of total flank sequence length', required=False, type=float, default=0.5)
	parser.add_argument('-x', '--exclude',  action='store', help='exclude these SNPs from LG analysis', required=False, type=argparse.FileType('r'))
	parser.add_argument('-m', '--multiple', action='store_true', help='do not use SNPs mapping multiple times for ordering and orienting sequences, recommended', required=False, default=False)
	parser.add_argument('-f', '--asm', action='store', help='fasta file for assembly', required=True, type=argparse.FileType('r'))
	parser.add_argument('-t', '--mito', action='store', help='name of the fasta entry corresponding to the mitochondrial genome', required=False)
	parser.add_argument('-n', '--ns', action='store', help='number of Ns between linkage mapped scaffolds', required=False, type=int, default=100)
	
	args = parser.parse_args()
	"""
	TODO:
	- some logging 
	- bed/gff file?
	
	"""

	snps = parse_psl_file(args.psl, args.ead)
	
	prefix = '.'.join(args.psl.name.split('.')[0:-1])
	
	#print("prefix: %s" % prefix)
	
	exclude_snps = []
	if args.exclude:
		prefix += "_exclude"
		for line in args.exclude:
			exclude_snps.append(line.strip())	
		
	
	#9355 SNPs in the linkage map - excluded snps
	lg_snps = read_lg(args.lg, exclude_snps)
	
	#Read in the fasta file
	fasta_entries = parse_fasta_file(args.asm)
	
	uniquely_mapping_snps = {}
	ambigously_mapping_snps = {}

	uniquely_mapping_snps_in_lg = {}
	longest_alignment_for_all_snps = {}
	
	#print (lg_snps)
	for snp in snps:
		#print (snp)
		#Maps only once
		if (len(snps[snp]) == 1):
			uniquely_mapping_snps[snp] = snps[snp]
			if (snp in lg_snps):
				#print ("Here")
				uniquely_mapping_snps_in_lg[snp] = snps[snp]
		else:
			ambigously_mapping_snps[snp] = snps[snp]
			
			
		
	#Filter based on Lex's rule (0.95 alignment length compared to flank sequence length)
	total_num_lex_snp_maps = 0
	lex_snps = {}
	for snp in snps:
		for i in snps[snp]:
			#0.95 alignment length compared to flank sequence length
			if ((float(i[0])/(int(i[10])-1)) >= 0.95):
				total_num_lex_snp_maps += 1
				if i[9] in lex_snps:
					lex_snps[i[9]].append(i)
				else:
					lex_snps[i[9]] = [i]
			
		
	for snp in snps:
		#sorted based on the fraction of alignment length divided by flanking sequence length
		sorted_hit_list = sorted(snps[snp], key = lambda ele: float(ele[0])/(int(ele[10])-1), reverse=True)
		#print (sorted_hit_list[0])
		longest_alignment_for_all_snps[snp] = [sorted_hit_list[0]]
		
	ole_snps = {}
	ole_snps_ident = {}
	best_ali_snps = {}
	#
	not_multiple_mapping_snps = {}

	
	best_alignments = 0
	secondary_alignments = 0
	
	#I want a separate file with all the "best" alignments by score
	#best_ali = open(prefix + ".best_alignments.tsv", "w")
	
	#My way, best score for all, plus all other alignments within 99 % of that score
	for snp in snps:
		#Sort by score first, then identity, then length of target sequence. Should prioritize hits 
		#against scaffolds instead of contigs
		sorted_score_list = sorted(snps[snp], key=itemgetter(21,22,14), reverse=True)
		#Always add the best scoring alignment, and the one against the longest sequence
		sorted_score_list[0][23] = True
		ole_snps[snp] = [sorted_score_list[0]]
		if ((sorted_score_list[0][22] >= args.ident) and ((float(sorted_score_list[0][0])/float(sorted_score_list[0][10])) >= args.length)):
			ole_snps_ident[snp] = [sorted_score_list[0]]
			not_multiple_mapping_snps[snp] = [sorted_score_list[0]]
		for j in sorted_score_list[1:]:
			#Not equal, has to be above args.score
			if ((float(j[21])/sorted_score_list[0][21]) >= args.score):
				#print("j")
				#print j
				ole_snps[snp].append(j)
				secondary_alignments += 1
			if (((float(j[21])/sorted_score_list[0][21]) >= args.score) and (j[22] >= args.ident) and ((float(j[0])/float(j[10])) >= args.length)):
				not_multiple_mapping_snps.pop(snp, None)
				
				#This should not happen, but apparently does
				if snp not in ole_snps_ident:						
					ole_snps_ident[snp] = [j]
				else:					
					ole_snps_ident[snp].append(j)
	
	
					
	

	if args.multiple:
		oles_method(prefix, snps, not_multiple_mapping_snps, lg_snps, args.ident, args.length, fasta_entries)	
	else:
		oles_method(prefix, snps, ole_snps_ident, lg_snps, args.ident, args.length, fasta_entries)	
	



"""	

1. matches - Number of matching bases that aren't repeats.
2. misMatches - Number of bases that don't match.
3. repMatches - Number of matching bases that are part of repeats.
4. nCount - Number of 'N' bases.
5. qNumInsert - Number of inserts in query.
6. qBaseInsert - Number of bases inserted into query.
7. tNumInsert - Number of inserts in target.
8. tBaseInsert - Number of bases inserted into target.
9. strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
10. qName - Query sequence name.
11. qSize - Query sequence size.
12. qStart - Alignment start position in query.
13. qEnd - Alignment end position in query.
14. tName - Target sequence name.
15. tSize - Target sequence size.
16. tStart - Alignment start position in query.
17. tEnd - Alignment end position in query.
18. blockCount - Number of blocks in the alignment.
19. blockSizes - Comma-separated list of sizes of each block.
20. qStarts - Comma-separated list of start position of each block in query.
21. tStarts - Comma-separated list of start position of each block in target.

Get score and % identity: http://genome.ucsc.edu/FAQ/FAQblat.html#blat4

And for the data structure in this program
22. Score 
23. Identity
24. Best alignment (true or false)


"""
