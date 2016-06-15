#!/usr/bin/env python
"""

Python 2.7+ is needed to run this.

"""

import sys, argparse
#from Bio import SeqIO
#import numpy as np
import pysam
import re

if __name__ == '__main__':
        
        parser = argparse.ArgumentParser(description=
        'Extract the position of start and end alignment of non-soft or hard-clipped reads/contigs. Output is bed format (points)')
        parser.add_argument('-b', '--bam', action='store', help='bam file with reads aligned to the assembly', required=True)

        args = parser.parse_args()
        
        #samfile = pysam.Samfile(args.bam, "r" )
        samfile = pysam.Samfile(args.bam, "rb" )

     	#for simplicity, likely only works with BWA output and single reads/contigs, not paired
        for read in samfile.fetch():
        	#first check if mapped or not
        	if not read.is_unmapped and read.mapq > 3 and read.mapq != 255:
        		#deal with hard clip in some way. Only 0x100 and 0x800 can have hard clip if mapped by BWA
        		if not read.is_secondary and not read.is_supplementary:
        			#can only have soft clip
        			#if zero, no softclip here:

        			if read.query_alignment_start == 0:
        				print ("%s\t%d\t%d\t%s" % (samfile.getrname(read.reference_id), read.get_reference_positions()[0], read.get_reference_positions()[0]+1, read.query_name))
        			#if these match, no soft clip
        			if read.query_alignment_end == read.query_length:
        				print ("%s\t%d\t%d\t%s" % (samfile.getrname(read.reference_id), read.get_reference_positions()[-1], read.get_reference_positions()[-1]+1, read.query_name))
        				#print (read.get_reference_positions()[-1])
        		else:
        			#likely hard clipped. If any of the sides match, get that.
        			m = re.findall("[A-Z]", read.cigarstring)
        			#print(m)
        			if m[0] == "M":
        				#first occurence is match, print pos
        				print ("%s\t%d\t%d\t%s" % (samfile.getrname(read.reference_id), read.get_reference_positions()[0], read.get_reference_positions()[0]+1, read.query_name))
        			elif m[-1] == "M":
        				#last occurence is match, print pos
        				print ("%s\t%d\t%d\t%s" % (samfile.getrname(read.reference_id), read.get_reference_positions()[-1], read.get_reference_positions()[-1]+1, read.query_name))
        			
