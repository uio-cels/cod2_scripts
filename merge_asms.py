#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
module load python2/2.7.6
"""

from __future__ import division

import sys
import argparse
import string
import copy
#import bx.align.maf

#To get a reverse enumerate
from itertools import izip

#from heapq import heappush, heappop

import heapq

import psyco_full

from bx.align import maf
from bx.align import core
from operator import itemgetter, attrgetter
from collections import defaultdict
from collections import deque
from bisect import bisect_left

debug = False

def getSequencePath(path, seq, asm_output, asm_output_gap_fraction):					
	sequence = ''
	
	for node in path:
		if asm_output:
			node_seq, node_src, amount = node.get_sequence_from(asm_output, seq, asm_output_gap_fraction)
		else:
			node_seq, node_src, amount = node.get_sequence_with_fewest_gaps(seq)	
		#node_seq, node_src, amount = node.get_sequence_with_fewest_gaps(seq)
		printDebug("in getSequencePath")
		printDebug(node.print_without_seq())
		#printDebug("Seq with fewest gaps: %s, amount %d " % (node_src, amount))
		#if len(sequence) > 0:
		#	printDebug("Fraction length of node/length of sequence: %.2g" % (float(node.size)/len(sequence)))
		#printDebug("Rev_comp on node %s" % node.rev_comp)				
		sequence = sequence + node_seq
	return sequence

class SeqPath(object):
	def __init__(self, id, path):
		self.id = id
		self.path = path
		#Number supporting links, length of path, connecting path, target sequence, index of starting node, index of target node, target path
		self.head_head = []
		self.head_tail = []
		self.tail_head = []
		self.tail_tail = []
		self.connections = []
		self.marked = False
		#Index of first node with max alignments/components
		self.head = None
		#Index of last node with max alignments/components
		self.tail = None
		
	def setTailAndHead(self, head, tail):	
		self.head = head
		self.tail = tail
		
	def add_head_head(self, head_head):     
		self.head_head.append(head_head)
		self.connections.append(head_head)
	
	def add_head_tail(self, head_tail):	
		self.head_tail.append(head_tail)
		self.connections.append(head_tail)
		
	def add_tail_head(self, tail_head):	
		self.tail_head.append(tail_head)
		self.connections.append(tail_head)
	
	def add_tail_tail(self, tail_tail):	
		self.tail_tail.append(tail_tail)
		self.connections.append(tail_tail)
		
	def info(self):
		s = "I contain %s, of length %d \n" % (self.id, len(self.path))
		self.head_head.sort()
		self.head_tail.sort()
		self.tail_head.sort()
		self.tail_tail.sort()
		
		if self.head_head:
			for hh in self.head_head:
				s = s + "head to head path to %s, %d common sequences, of length %d, starting at index %d and ending at index %d (of %d nodes)\n" % (hh[3], hh[0], hh[1], hh[4], hh[5], len(hh[6].path))
		
		if self.head_tail:
			for ht in self.head_tail:
				s = s + "head to tail path to %s, %d common sequences, of length %d, starting at index %d and ending at index %d (of %d nodes)\n" % (ht[3], ht[0], ht[1], ht[4], ht[5], len(ht[6].path))

		if self.tail_head:
			for th in self.tail_head:
				s = s + "tail to head path to %s, %d common sequences, of length %d, starting at index %d and ending at index %d (of %d nodes)\n" % (th[3], th[0], th[1], th[4], th[5], len(th[6].path))

		if self.tail_tail:
			for tt in self.tail_tail:
				s = s + "tail to tail path to %s, %d common sequences, of length %d, starting at index %d and ending at index %d (of %d nodes)\n" % (tt[3], tt[0], tt[1], tt[4], tt[5], len(tt[6].path))	
									
		return s 


	def getPath(self):
		"""
		Returns a path.
		"""

		printDebug("In getPath")
		printDebug("self.marked: %s" % self.marked)
		printDebug("len(self.connections): %d" %  len(self.connections))
		if self.marked and len(self.connections) == 0:
			printDebug("Returning nothing")
			return None, None
		elif len(self.connections) == 0:
			printDebug("Returning something")
			self.marked = True
			return self.id, self.path
		else:
			printDebug("Default")
			return self.id, self.path
		
		
	def getPathSequence(self, asm_output, entered, asm_output_gap_fraction):
		"""
		Returns the sequence of a path.
		"""
		printDebug('')		
		printDebug("In getPathSequence")
		printDebug("self.marked: %s" % self.marked)
		#printDebug("self.id: %s" % self.id)
		printDebug(self.info(), endline='')
		if self.head:
			printDebug("self.head: %d" % self.head)
		if self.tail:
			printDebug("self.tail: %d" % self.tail)			
		
		#printDebug("len(self.connections): %d" %  len(self.connections))
		#printDebug("len(self.head_head): %d" %  len(self.head_head))
		#printDebug("len(self.head_tail): %d" %  len(self.head_tail))
		#printDebug("len(self.tail_head): %d" %  len(self.tail_head))
		#printDebug("len(self.tail_tail): %d" %  len(self.tail_tail))
		if self.marked and len(self.connections) == 0:
			printDebug("Returning nothing")
			return [], None
		elif self.marked:
			#Been used already
			printDebug("Been marked")
			return [], None	
		elif len(self.connections) == 0:
			printDebug("Returning self")
			self.marked = True
			sequence = getSequencePath(self.path, self.id, asm_output, asm_output_gap_fraction)
			return [self.id], sequence
		elif not entered:
			printDebug("In not entered")
			sequence = ''
			seqs = [self.id]
			headExtended = False
			tailExtended = False
			headSequence = ''
			tailSequence = ''
			self.marked = True
			#this sequence has connections to other sequences, so we have to get those
			
			"""
			#Sort all entries here at first
			#Sort by path length first
			self.head_tail.sort(key=itemgetter(1))
			#Then get the one with the most connections
			self.head_tail.sort(key=itemgetter(0), reverse=True)
			
			self.head_head.sort(key=itemgetter(1))
			self.head_head.sort(key=itemgetter(0), reverse=True)	
			
			self.tail_tail.sort(key=itemgetter(1))
			self.tail_tail.sort(key=itemgetter(0), reverse=True)
			
			self.tail_head.sort(key=itemgetter(1))
			self.tail_head.sort(key=itemgetter(0), reverse=True)	
			"""
			
			#Sort by path length first
			self.head_tail.sort(key=itemgetter(1), reverse=True)
			self.head_head.sort(key=itemgetter(1), reverse=True)	
			self.tail_tail.sort(key=itemgetter(1), reverse=True)
			self.tail_head.sort(key=itemgetter(1), reverse=True)
			#self.head_tail.sort(reverse=True)
			#self.head_head.sort(reverse=True)	
			#self.tail_tail.sort(reverse=True)
			#self.tail_head.sort(reverse=True)										
			lengthHeadConnectionPath = -1
			lengthTailConnectionPath = -1
			
			workWith = ''
				
			if self.head_tail and self.head_head:
				#Should sort with most connections first, then the longest.
				#self.head_tail.sort(reverse=True)
				#self.head_head.sort(reverse=True)
				"""
				#If equal number of connections				
				if self.head_head[0][0] == self.head_tail[0][0]:
					if self.head_head[0][1] <= self.head_tail[0][1]:
						#work with head_tail if equal or longer length
						workWith = 'head_tail'
					else:
						workWith = 'head_head'
				#more connections with head_tail		
				elif self.head_head[0][0] < self.head_tail[0][0]:		
					workWith = 'head_tail'
				else:
					workWith = 'head_head'
				
				#If number of connections is equal or more than 2		
				if self.head_head[0][0] >= 2 and self.head_tail[0][0] >= 2:
					if self.head_head[0][1] <= self.head_tail[0][1]:
						#work with head_tail if equal or longer length
						workWith = 'head_tail'
					else:
						workWith = 'head_head'
				#more connections with head_tail		
				elif self.head_head[0][0] < self.head_tail[0][0]:		
					workWith = 'head_tail'
				else:
					workWith = 'head_head'				
				"""	
				#If both head_head and head_tail point to the same scaffold, take the shortest path
				if self.head_head[0][3] == self.head_tail[0][3]:
					if self.head_head[0][1] < self.head_tail[0][1]:
						workWith = 'head_head'
					else:
						workWith = 'head_tail'
				elif self.head_head[0][1] <= self.head_tail[0][1]:
					workWith = 'head_tail'	
				else:
					workWith = 'head_head'
				
				
				#if self.head_head[0][1] <= self.head_tail[0][1]:
				if 	workWith == 'head_tail':	
					printDebug("self.head_tail")
					print("self.head_tail")
					print self.head_tail
					
					connectingPathLength = 0
					for node in self.head_tail[0][2]:
						connectingPathLength = connectingPathLength + node.size
					
					printDebug("Path connecting %s and %s, consists of %d nodes and has a size of %d." % (self.id, self.head_tail[0][3], self.head_tail[0][1], connectingPathLength))
					"""
					combinedPathLength = 0
					combinedPath = self.path[:self.head] + self.head_tail[0][6].path[self.head_tail[0][5]:]
					for node in combinedPath:
						combinedPathLength = combinedPathLength + node.size		
						
					printDebug("Path consisting of start of %s and end of %s, consists of %d nodes and has a size of %d." % (self.id, self.head_tail[0][3], len(combinedPath), combinedPathLength))
					"""
					
					#Then get that path
					printDebug("Getting head_sequences")
					head_sequences, head_sequence = self.head_tail[0][6].getPathSequence(asm_output, 'tail', asm_output_gap_fraction)
				
					#If that exists, add the connection sequence, and merge
					if head_sequences and head_sequence:
						#This path goes from head of this sequence to tail of another, I think I should reverse it to be correct
						printDebug("Get connectionSeq")
						lengthHeadConnectionPath = self.head_tail[0][1]
						connectionSeq = getSequencePath(self.head_tail[0][2], self.id, asm_output, asm_output_gap_fraction)
						#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
						seqs = head_sequences + seqs
						headSequence = head_sequence + connectionSeq	
				else:
					#If head_head has more connections, we go with that
					printDebug("self.head_head")
					print("self.head_head")
					#Find the one furthest away
					print self.head_head
					connectingPathLength = 0
					
					for node in self.head_head[0][2]:
						connectingPathLength = connectingPathLength + node.size
					
					printDebug("Path connecting %s and %s, consists of %d nodes and has a size of %d." % (self.id, self.head_head[0][3], self.head_head[0][1], connectingPathLength))				
					
					
					#Then get that path
					printDebug("Get head_sequences")
					head_sequences, head_sequence = self.head_head[0][6].getPathSequence(asm_output, 'head', asm_output_gap_fraction)
				
					#If that exists, add the connection sequence, and merge
					if head_sequences and head_sequence:
						#This path goes from head of this sequence to head of another, I think I should reverse it to be correct
						printDebug("Get connectionSeq")
						lengthHeadConnectionPath = self.head_head[0][1]
						connectionSeq = getSequencePath(self.head_head[0][2], self.id, asm_output, asm_output_gap_fraction)
						#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
						head_sequences.reverse()
						seqs = head_sequences + seqs
						headSequence = revcomp(head_sequence) + connectionSeq	
																
						
			#The head of this has a connection to the tail of another.
			elif self.head_tail:
				printDebug("self.head_tail")
				print("self.head_tail")
				#Find the one with most connections
				#self.head_tail.sort(reverse=True)
				print self.head_tail
				connectingPathLength = 0
				for node in self.head_tail[0][2]:
					connectingPathLength = connectingPathLength + node.size
				
				printDebug("Path connecting %s and %s, consists of %d nodes and has a size of %d." % (self.id, self.head_tail[0][3], self.head_tail[0][1], connectingPathLength))				
				#Then get that path
				printDebug("Get head_sequences")
				head_sequences, head_sequence = self.head_tail[0][6].getPathSequence(asm_output, 'tail', asm_output_gap_fraction)
				
				#If that exists, add the connection sequence, and merge
				if head_sequences and head_sequence:
					#This path goes from head of this sequence to tail of another, I think I should reverse it to be correct
					printDebug("Get connectionSeq")
					lengthHeadConnectionPath = self.head_tail[0][1]
					connectionSeq = getSequencePath(self.head_tail[0][2], self.id, asm_output, asm_output_gap_fraction)

					#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
					seqs = head_sequences + seqs					
					headSequence = head_sequence + connectionSeq
					
			elif self.head_head:
				#If head_head is longer, then we go with that
				printDebug("self.head_head")
				print("self.head_head")
				#Find the one furthest away
				print self.head_head
				connectingPathLength = 0
				for node in self.head_head[0][2]:
					connectingPathLength = connectingPathLength + node.size
				
				printDebug("Path connecting %s and %s, consists of %d nodes and has a size of %d." % (self.id, self.head_head[0][3], self.head_head[0][1], connectingPathLength))			
				#Then get that path
				printDebug("Get head_sequences")
				head_sequences, head_sequence = self.head_head[0][6].getPathSequence(asm_output, 'head', asm_output_gap_fraction)
			
				#If that exists, add the connection sequence, and merge
				if head_sequences and head_sequence:
					#This path goes from head of this sequence to head of another, I think I should reverse it to be correct
					printDebug("Get connectionSeq")
					lengthHeadConnectionPath = self.head_head[0][1]
					connectionSeq = getSequencePath(self.head_head[0][2], self.id, asm_output, asm_output_gap_fraction)					
					#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
					head_sequences.reverse()
					seqs = head_sequences + seqs
					headSequence = revcomp(head_sequence) + connectionSeq					
			
			#At this point, a path from head would have been found. We want to close down
			#all other potential paths, for instance a shorter path to a short sequence
			printDebug("At head connections in entered")
			if self.head_head:
				for connection in self.head_head:
					#Only mark paths shorter than the head connection path
					if len(connection[6].path) < lengthHeadConnectionPath:
						printDebug(connection[3])
						printDebug(connection[6].marked)
						connection[6].marked = True
			if self.head_tail:
				for connection in self.head_tail:
					#Only mark paths shorter than the head connection path
					if len(connection[6].path) < lengthHeadConnectionPath:				
						printDebug(connection[3])
						printDebug(connection[6].marked)
						connection[6].marked = True			
			
				
				#print reducedPaths[p.head_tail[0][2]]
				#for node in reducedPaths[p.head_tail[0][2]][0:p.head_tail[0][4]]:				
				
			workWith = ''
			
			#Then get longest path from tail	
			if self.tail_tail and self.tail_head:			
				#self.tail_tail.sort(reverse=True)
				#self.tail_head.sort(reverse=True)
				"""
				#If equal number of connections				
				#if self.tail_head[0][0] == self.tail_tail[0][0]:
				#If both are equal to or more than 2
				if self.tail_head[0][0] >= 2 and self.tail_tail[0][0] >= 2:
					if self.tail_head[0][1] >= self.tail_tail[0][1]:
						#work with tail_head if equal or longer length
						workWith = 'tail_head'
					else:
						workWith = 'tail_tail'
				#more connections with head_tail		
				elif self.tail_head[0][0] > self.tail_tail[0][0]:		
					workWith = 'tail_head'
				else:
					workWith = 'tail_tail'
				"""
				if self.tail_head[0][3] == self.tail_tail[0][3]:
					if self.tail_head[0][1] > self.tail_tail[0][1]:
						workWith = 'tail_tail'
					else:
						workWith = 'tail_head'
				elif self.tail_head[0][1] >= self.tail_tail[0][1]:
					workWith = 'tail_head'
				else:
					workWith = 'tail_tail'
				
				
				#if self.head_head[0][1] <= self.head_tail[0][1]:
				if 	workWith == 'tail_head':				
				#if self.tail_head[0][1] >= self.tail_tail[0][1]:				
				
					#Go with tail_head if that has equal or more connections
					printDebug("self.tail_head")
					print("self.tail_head")
					#self.tail_head.sort(reverse=True)
					print self.tail_head
					connectingPathLength = 0
					for node in self.tail_head[0][2]:
						connectingPathLength = connectingPathLength + node.size
					
					printDebug("Path connecting %s and %s, consists of %d nodes and has a size of %d." % (self.id, self.tail_head[0][3], self.tail_head[0][1], connectingPathLength))				
					#the sequences with tail first
					#print p.head_tail
					#print reducedPaths[p.head_tail[0][2]]
					#for node in reducedPaths[p.head_tail[0][2]][0:p.head_tail[0][4]]:				
				
					printDebug("Get tail_sequences")
					tail_sequences, tail_sequence = self.tail_head[0][6].getPathSequence(asm_output, 'head', asm_output_gap_fraction)
					#If that exists, add the connection sequence, and merge
					if tail_sequences and tail_sequence:
						printDebug("Get connectionSeq")
						lengthTailConnectionPath = self.tail_head[0][1]
						connectionSeq = getSequencePath(self.tail_head[0][2], self.id, asm_output, asm_output_gap_fraction)						
						#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
						seqs = seqs + tail_sequences
						tailSequence = connectionSeq + tail_sequence
					
				else:
					#If tail_tail is longer, then we go with that
					printDebug("self.tail_tail")
					print("self.tail_tail")
					#Find the one furthest away
					print self.tail_tail
					connectingPathLength = 0
					for node in self.tail_tail[0][2]:
						connectingPathLength = connectingPathLength + node.size
					
					printDebug("Path connecting %s and %s, consists of %d (%d) nodes and has a size of %d." % (self.id, self.tail_tail[0][3], self.tail_tail[0][1], len(self.tail_tail[0][2]), connectingPathLength))							
					#Then get that path
					printDebug("Get tail_sequences")
					tail_sequences, tail_sequence = self.tail_tail[0][6].getPathSequence(asm_output, 'tail', asm_output_gap_fraction)
				
					#If that exists, add the connection sequence, and merge
					if tail_sequences and tail_sequence:
						#This path goes from tail of this sequence to tail of another, I think I should reverse it to be correct
						printDebug("Get connectionSeq")
						lengthTailConnectionPath = self.tail_tail[0][1]
						connectionSeq = getSequencePath(self.tail_tail[0][2], self.id, asm_output, asm_output_gap_fraction)
						
						#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
						tail_sequences.reverse()
						seqs = seqs + tail_sequences
						tailSequence = connectionSeq + revcomp(tail_sequence)		
			#The tail of this is connected to the head of another.
			elif self.tail_head:	
				printDebug("self.tail_head")
				print("self.tail_head")
				#self.tail_head.sort(reverse=True)
				print self.tail_head

				connectingPathLength = 0
				for node in self.tail_head[0][2]:
					connectingPathLength = connectingPathLength + node.size
				
				printDebug("Path connecting %s and %s, consists of %d nodes and has a size of %d." % (self.id, self.tail_head[0][3], self.tail_head[0][1], connectingPathLength))							
				
				printDebug("Get tail_sequences")
				tail_sequences, tail_sequence = self.tail_head[0][6].getPathSequence(asm_output, 'head', asm_output_gap_fraction)
				#If that exists, add the connection sequence, and merge
				if tail_sequences and tail_sequence:
					printDebug("Get connectionSeq")
					lengthTailConnectionPath = self.tail_head[0][1]
					connectionSeq = getSequencePath(self.tail_head[0][2], self.id, asm_output, asm_output_gap_fraction)
						
					#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
					seqs = seqs + tail_sequences
					tailSequence = connectionSeq + tail_sequence
			
			elif self.tail_tail:
				printDebug("self.tail_tail")
				print("self.tail_tail")
				#Find the one furthest away
				print self.tail_tail

				connectingPathLength = 0
				for node in self.tail_tail[0][2]:
					connectingPathLength = connectingPathLength + node.size
				
				printDebug("Path connecting %s and %s, consists of %d (%d) nodes and has a size of %d." % (self.id, self.tail_tail[0][3], self.tail_tail[0][1], len(self.tail_tail[0][2]), connectingPathLength))							
			
				#Then get that path
				printDebug("Get tail_sequences")
				tail_sequences, tail_sequence = self.tail_tail[0][6].getPathSequence(asm_output, 'tail', asm_output_gap_fraction)
			
				#If that exists, add the connection sequence, and merge
				if tail_sequences and tail_sequence:
					#This path goes from tail of this sequence to tail of another, I think I should reverse it to be correct	
					printDebug("Get connectionSeq")		
					lengthTailConnectionPath = self.tail_tail[0][1]	
					connectionSeq = getSequencePath(self.tail_tail[0][2], self.id, asm_output, asm_output_gap_fraction)
					#printDebug(seqs)
					#printDebug(tail_sequences)				
					#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
					tail_sequences.reverse()
					seqs = seqs + tail_sequences
					tailSequence = connectionSeq + revcomp(tail_sequence)	

			#At this point, a path from tail would have been found. We want to close down
			#all other potential paths, for instance a shorter path to a short sequence
			printDebug("At tail connections in entered")
			if self.tail_head:
				for connection in self.tail_head:
					if len(connection[6].path) < lengthTailConnectionPath:	
						printDebug(connection[3])
						printDebug(connection[6].marked)
						connection[6].marked = True
						
			if self.tail_tail:
				for connection in self.tail_tail:
					if len(connection[6].path) < lengthTailConnectionPath:	
						printDebug(connection[3])
						printDebug(connection[6].marked)
						connection[6].marked = True	
									
			if headSequence and tailSequence:
				selfSequence = getSequencePath(self.path[self.head:self.tail+1], self.id, asm_output, asm_output_gap_fraction)
					
				return seqs, headSequence + selfSequence + tailSequence
				
			elif headSequence:
				#printDebug("elif headSequence")
				#printDebug(headSequence)
				selfSequence = getSequencePath(self.path[self.head:], self.id, asm_output, asm_output_gap_fraction)				
				return seqs, headSequence + selfSequence
				
			elif tailSequence:
				selfSequence = getSequencePath(self.path[:self.tail+1], self.id, asm_output, asm_output_gap_fraction)
					
				return seqs, selfSequence + tailSequence
				
			else:
				selfSequence = getSequencePath(self.path, self.id, asm_output, asm_output_gap_fraction)

				return seqs, selfSequence		
			
		elif entered == 'head':
			#Entered this sequence via head, check if connections via tail, else return self
			self.marked = True
			sequence = ''
			seqs = [self.id]
			headExtended = False
			tailExtended = False
			tailSequence = ''
			lengthHeadConnectionPath = -1
			lengthTailConnectionPath = -1			
			"""		
			self.tail_tail.sort(key=itemgetter(1))
			self.tail_tail.sort(key=itemgetter(0), reverse=True)
			
			self.tail_head.sort(key=itemgetter(1))
			self.tail_head.sort(key=itemgetter(0), reverse=True)			
			"""	
			#self.tail_tail.sort(reverse=True)
			#self.tail_head.sort(reverse=True)					
	
			self.tail_tail.sort(key=itemgetter(1), reverse=True)
			self.tail_head.sort(key=itemgetter(1), reverse=True)			
			printDebug("entered via head")
			#printDebug("printing connections")
			#print("entered via head")
			#print("printing connections")				
			#print self.connections
			
			workWith = ''
			
			#Then get longest path from tail	
			if self.tail_tail and self.tail_head:			
				#self.tail_tail.sort(reverse=True)
				#self.tail_head.sort(reverse=True)
				"""
				#If equal number of connections				
				#if self.tail_head[0][0] == self.tail_tail[0][0]:
				if self.tail_head[0][0] >= 2 and self.tail_tail[0][0] >= 2:
					if self.tail_head[0][1] >= self.tail_tail[0][1]:
						#work with tail_head if equal or longer length
						workWith = 'tail_head'
					else:
						workWith = 'tail_tail'
				#more connections with head_tail		
				elif self.tail_head[0][0] > self.tail_tail[0][0]:		
					workWith = 'tail_head'
				else:
					workWith = 'tail_tail'
				"""

									
				#If both head_head and head_tail point to the same scaffold, take the shortest path
				if self.tail_head[0][3] == self.tail_tail[0][3]:
					if self.tail_head[0][1] > self.tail_tail[0][1]:
						workWith = 'tail_tail'
					else:
						workWith = 'tail_head'
				elif self.tail_head[0][1] >= self.tail_tail[0][1]:
					workWith = 'tail_head'
				else:
					workWith = 'tail_tail'
				
				
				#if self.head_head[0][1] <= self.head_tail[0][1]:
				if 	workWith == 'tail_head':					
				#if self.tail_head[0][1] >= self.tail_tail[0][1]:				
					
					printDebug("self.tail_head")

					connectingPathLength = 0
					for node in self.tail_head[0][2]:
						connectingPathLength = connectingPathLength + node.size
				
					printDebug("Path connecting %s and %s, consists of %d nodes and has a size of %d." % (self.id, self.tail_head[0][3], self.tail_head[0][1], connectingPathLength))							
					
					printDebug("Get tail_sequences")			
					tail_sequences, tail_sequence = self.tail_head[0][6].getPathSequence(asm_output, 'head', asm_output_gap_fraction)
					#If that exists, add the connection sequence, and merge
					if tail_sequences and tail_sequence:
						printDebug("Get connectionSeq")
						lengthTailConnectionPath = self.tail_head[0][1]	
						
						connectionSeq = getSequencePath(self.tail_head[0][2], self.id, asm_output, asm_output_gap_fraction)						
						#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
						seqs = seqs + tail_sequences
						tailSequence = connectionSeq + tail_sequence
					
				else:
					#If tail_tail is longer, then we go with that
					printDebug("self.tail_tail")

					connectingPathLength = 0
					for node in self.tail_tail[0][2]:
						connectingPathLength = connectingPathLength + node.size
				
					printDebug("Path connecting %s and %s, consists of %d (%d) nodes and has a size of %d." % (self.id, self.tail_tail[0][3], self.tail_tail[0][1], len(self.tail_tail[0][2]), connectingPathLength))							
								
					#Then get that path
					printDebug("Get tail_sequences")
					tail_sequences, tail_sequence = self.tail_tail[0][6].getPathSequence(asm_output, 'tail', asm_output_gap_fraction)
				
					#If that exists, add the connection sequence, and merge
					if tail_sequences and tail_sequence:
						#This path goes from tail of this sequence to tail of another, I think I should reverse it to be correct
						printDebug("Get connectionSeq")
						lengthTailConnectionPath = self.tail_tail[0][1]	
						connectionSeq = getSequencePath(self.tail_tail[0][2], self.id, asm_output, asm_output_gap_fraction)
						
						#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
						tail_sequences.reverse()
						seqs = seqs + tail_sequences
						tailSequence = connectionSeq + revcomp(tail_sequence)		
			#The tail of this is connected to the head of another.
			elif self.tail_head:	
				printDebug("self.tail_head")
				print("self.tail_head")
				#self.tail_head.sort(reverse=True)
				print self.tail_head
				
				#the sequences with tail first
				#print p.head_tail
				#print reducedPaths[p.head_tail[0][2]]
				#for node in reducedPaths[p.head_tail[0][2]][0:p.head_tail[0][4]]:				
				connectingPathLength = 0
				for node in self.tail_head[0][2]:
					connectingPathLength = connectingPathLength + node.size
			
				printDebug("Path connecting %s and %s, consists of %d (%d) nodes and has a size of %d." % (self.id, self.tail_head[0][3], self.tail_head[0][1], len(self.tail_head[0][2]), connectingPathLength))							
				
				
				printDebug("Get tail_sequences")
				tail_sequences, tail_sequence = self.tail_head[0][6].getPathSequence(asm_output, 'head', asm_output_gap_fraction)
				#If that exists, add the connection sequence, and merge
				if tail_sequences and tail_sequence:
					printDebug("Get connectionSeq")
					lengthTailConnectionPath = self.tail_head[0][1]	
					connectionSeq = getSequencePath(self.tail_head[0][2], self.id, asm_output, asm_output_gap_fraction)
						
					#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
					seqs = seqs + tail_sequences
					tailSequence = connectionSeq + tail_sequence
			
			elif self.tail_tail:
				printDebug("self.tail_tail")
				print("self.tail_tail")
				#Find the one furthest away
				print self.tail_tail

				connectingPathLength = 0
				for node in self.tail_tail[0][2]:
					connectingPathLength = connectingPathLength + node.size
			
				printDebug("Path connecting %s and %s, consists of %d (%d) nodes and has a size of %d." % (self.id, self.tail_tail[0][3], self.tail_tail[0][1], len(self.tail_tail[0][2]), connectingPathLength))							
			
				#Then get that path
				printDebug("Get tail_sequences")
				tail_sequences, tail_sequence = self.tail_tail[0][6].getPathSequence(asm_output, 'tail', asm_output_gap_fraction)
	
			
				#If that exists, add the connection sequence, and merge
				if tail_sequences and tail_sequence:
					#This path goes from tail of this sequence to tail of another, I think I should reverse it to be correct	
					printDebug("Get connectionSeq")		
					lengthTailConnectionPath = self.tail_tail[0][1]		
					connectionSeq = getSequencePath(reversed(self.tail_tail[0][2]), self.id, asm_output, asm_output_gap_fraction)
				
					#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
					tail_sequences.reverse()
					seqs = seqs + tail_sequences
					tailSequence = connectionSeq + revcomp(tail_sequence)
											
			#At this point, a path from tail would have been found. We want to close down
			#all other potential paths, for instance a shorter path to a short sequence
			printDebug("At tail connections in entered via head")
			if self.tail_head:
				for connection in self.tail_head:
					if len(connection[6].path) < lengthTailConnectionPath:					
						printDebug(connection[3])
						printDebug(connection[6].marked)
						connection[6].marked = True
			if self.tail_tail:
				for connection in self.tail_tail:
					if len(connection[6].path) < lengthTailConnectionPath:	
						printDebug(connection[3])
						printDebug(connection[6].marked)
						connection[6].marked = True	
			
			if tailSequence:
				#This should hopefully work 
				selfSequence = getSequencePath(self.path[self.head:self.tail+1], self.id, asm_output, asm_output_gap_fraction)
					
				return seqs, selfSequence + tailSequence	
														
			else:
				#If no other connections has been found, we end up here
				selfSequence = getSequencePath(self.path[self.head:], self.id, asm_output, asm_output_gap_fraction)

				return seqs, selfSequence						
								
		elif entered == 'tail':	
			self.marked = True
			sequence = ''
			seqs = [self.id]
			headExtended = False
			tailExtended = False
			headSequence = ''
			lengthHeadConnectionPath = -1
			lengthTailConnectionPath = -1			
			"""
			#Sort all entries here at first
			#Sort by path length first
			self.head_tail.sort(key=itemgetter(1))
			#Then get the one with the most connections
			self.head_tail.sort(key=itemgetter(0), reverse=True)
			
			self.head_head.sort(key=itemgetter(1))
			self.head_head.sort(key=itemgetter(0), reverse=True)	
			"""
			#self.head_tail.sort(reverse=True)
			#self.head_head.sort(reverse=True)	
			
			self.head_tail.sort(key=itemgetter(1), reverse=True)
			self.head_head.sort(key=itemgetter(1), reverse=True)	
							
			printDebug("entered via tail")

			workWith = ''
				
			if self.head_tail and self.head_head:
				#Should sort with most connections first, then the longest.
				#self.head_tail.sort(reverse=True)
				#self.head_head.sort(reverse=True)
				"""
				#If equal number of connections				
				#if self.head_head[0][0] == self.head_tail[0][0]:
				if self.head_head[0][0] >= 2 and self.head_tail[0][0] >= 2:
					if self.head_head[0][1] <= self.head_tail[0][1]:
						#work with head_tail if equal or longer length
						workWith = 'head_tail'
					else:
						workWith = 'head_head'
				#more connections with head_tail		
				elif self.head_head[0][0] < self.head_tail[0][0]:		
					workWith = 'head_tail'
				else:
					workWith = 'head_head'
				"""
				if self.head_head[0][3] == self.head_tail[0][3]:
					if self.head_head[0][1] < self.head_tail[0][1]:
						workWith = 'head_head'
					else:
						workWith = 'head_tail'
				elif self.head_head[0][1] <= self.head_tail[0][1]:
					workWith = 'head_tail'	
				else:
					workWith = 'head_head'
				
				
				#if self.head_head[0][1] <= self.head_tail[0][1]:
				if workWith == 'head_tail':	

					#Go with head_tail if that has equal or more connections. Don't care about length now
					printDebug("self.head_tail")

					connectingPathLength = 0
					for node in self.head_tail[0][2]:
						connectingPathLength = connectingPathLength + node.size
				
					printDebug("Path connecting %s and %s, consists of %d nodes and has a size of %d." % (self.id, self.head_tail[0][3], self.head_tail[0][1], connectingPathLength))							
					

				
					#Then get that path
					printDebug("Get head_sequences")
					head_sequences, head_sequence = self.head_tail[0][6].getPathSequence(asm_output, 'tail', asm_output_gap_fraction)
				
					#If that exists, add the connection sequence, and merge
					if head_sequences and head_sequence:
						printDebug("Get connectionSeq")
						lengthHeadConnectionPath = self.head_tail[0][1]	
						connectionSeq = getSequencePath(self.head_tail[0][2], self.id, asm_output, asm_output_gap_fraction)
						seqs = head_sequences + seqs
						headSequence = head_sequence + connectionSeq	
				else:
					#If head_head is longer, then we go with that
					printDebug("self.head_head")
					print("self.head_head")
					#Find the one furthest away
					print self.head_head

					connectingPathLength = 0
					for node in self.head_head[0][2]:
						connectingPathLength = connectingPathLength + node.size
				
					printDebug("Path connecting %s and %s, consists of %d nodes and has a size of %d." % (self.id, self.head_head[0][3], self.head_head[0][1], connectingPathLength))							
					

				
					#Then get that path
					printDebug("Get head_sequences")
					head_sequences, head_sequence = self.head_head[0][6].getPathSequence(asm_output, 'head', asm_output_gap_fraction)
				
					#If that exists, add the connection sequence, and merge
					if head_sequences and head_sequence:
						#This path goes from head of this sequence to head of another, I think I should reverse it to be correct
						printDebug("Get connectionSeq")
						lengthHeadConnectionPath = self.head_head[0][1]	
						connectionSeq = getSequencePath(self.head_head[0][2], self.id, asm_output, asm_output_gap_fraction)
						#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
						head_sequences.reverse()
						seqs = head_sequences + seqs
						headSequence = revcomp(head_sequence) + connectionSeq											
						
			#The head of this has a connection to the tail of another.
			elif self.head_tail:
				printDebug("self.head_tail")
				print("self.head_tail")
				#Find the one with most connections
				#self.head_tail.sort(reverse=True)
				#printDebug(self.head_tail)

				connectingPathLength = 0
				for node in self.head_tail[0][2]:
					connectingPathLength = connectingPathLength + node.size
			
				printDebug("Path connecting %s and %s, consists of %d nodes and has a size of %d." % (self.id, self.head_tail[0][3], self.head_tail[0][1], connectingPathLength))							
				

				
				#Then get that path
				printDebug("Get head_sequences")
				head_sequences, head_sequence = self.head_tail[0][6].getPathSequence(asm_output, 'tail', asm_output_gap_fraction)
				
				#If that exists, add the connection sequence, and merge
				if head_sequences and head_sequence:
					#This path goes from head of this sequence to tail of another, I think I should reverse it to be correct
					printDebug("Get connectionSeq")
					lengthHeadConnectionPath = self.head_tail[0][1]	
					connectionSeq = getSequencePath(self.head_tail[0][2], self.id, asm_output, asm_output_gap_fraction)

					#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
					seqs = head_sequences + seqs
					headSequence = head_sequence + connectionSeq
					
			elif self.head_head:
				#If head_head is longer, then we go with that
				printDebug("self.head_head")
				print("self.head_head")
				#Find the one furthest away
				print self.head_head

				connectingPathLength = 0
				for node in self.head_head[0][2]:
					connectingPathLength = connectingPathLength + node.size
				
				printDebug("Path connecting %s and %s, consists of %d nodes and has a size of %d." % (self.id, self.head_head[0][3], self.head_head[0][1], connectingPathLength))							
					

			
				#Then get that path
				printDebug("Get head_sequences")
				head_sequences, head_sequence = self.head_head[0][6].getPathSequence(asm_output, 'head', asm_output_gap_fraction)
			
				#If that exists, add the connection sequence, and merge
				if head_sequences and head_sequence:
					#This path goes from head of this sequence to head of another, I think I should reverse it to be correct
					printDebug("Get connectionSeq")
					lengthHeadConnectionPath = self.head_head[0][1]	
					connectionSeq = getSequencePath(self.head_head[0][2], self.id, asm_output, asm_output_gap_fraction)					
					#Probably need to check if I need to reverse complement connectionSeq, but dunno yet.
					head_sequences.reverse()
					seqs = head_sequences + seqs
					headSequence = revcomp(head_sequence) + connectionSeq	
						
			#At this point, a path from head would have been found. We want to close down
			#all other potential paths, for instance a shorter path to a short sequence
			printDebug("At head connections in entered via tail")
			if self.head_head:
				for connection in self.head_head:
					if len(connection[6].path) < lengthHeadConnectionPath:					
						printDebug(connection[3])
						printDebug(connection[6].marked)
						connection[6].marked = True
			if self.head_tail:
				for connection in self.head_tail:
					if len(connection[6].path) < lengthHeadConnectionPath:	
						printDebug(connection[3])
						printDebug(connection[6].marked)
						connection[6].marked = True	
								
			if headSequence:	
				selfSequence = getSequencePath(self.path[self.head:self.tail+1], self.id, asm_output, asm_output_gap_fraction)				
				return seqs, headSequence + selfSequence

			else:
				selfSequence = getSequencePath(self.path[:self.tail+1], self.id, asm_output, asm_output_gap_fraction)
			
				return seqs, selfSequence
								
		else:
			printDebug("Default")
			sequence = getSequencePath(self.path, self.id, asm_output, asm_output_gap_fraction)
			
			return [self.id], sequence			

class ExtComponent(core.Component):
	def __init__( self, src='', start=0, size=0, strand=None, src_size=None, text='', num_n_bases=0, plus_strand_start=None):
		"""
		Copied from /lib/bx/align/core.py in bxpython. Couldn't find a license.
		"""
		
		self._alignment = None
		self.src = src
		self.start = start          # Nota Bene:  start,size,strand are as they
		self.size = size            # .. appear in a MAF file-- origin-zero, end
		self.strand = strand        # .. excluded, and minus strand counts from
		self._src_size = src_size   # .. end of sequence
		self.text = text
		self.quality = None
		# Optional fields to keep track of synteny status (only makes sense
		# when the alignment is part of an ordered set)
		self.synteny_left = None
		self.synteny_right = None
		self.synteny_empty = None
		# If true, this component actually represents a non-aligning region,
		# and has no text.
		self.empty = False
		# Index maps a coordinate (distance along + strand from + start) to alignment column
		self.index = None
		self.num_n_bases = num_n_bases		#Added to extend the class
		self.plus_strand_start = plus_strand_start


	def calc_n_bases(self):
		n = 0
		for i in self.text:
			if i.lower() == 'n':
				n += 1
		self.num_n_bases = n
		
	def calc_gap_bases(self):
		"""
		Calculates the amount of both N based and - bases
		"""
		n = 0
		for i in self.text:
			if i.lower() == 'n' or i == '-':
				n += 1
		self.num_gap_bases = n
	
	def gap_less(self):
		
		gap_less = ''
		for i in self.text:
			if i != '-':
				gap_less += i
				
		self.gap_less = gap_less

	def print_without_seq(self):
		return "s %s %d %d %d %s %d" % ( self.src, self.plus_strand_start, self.start,
                                             self.size, self.strand,
                                             self.src_size ) 		
		
		
	def __str__( self ):
		return "s %s %d %d %d %s %d %s" % (self.src, self.plus_strand_start, 
											self.start, self.size, self.strand, 
											self.src_size, self.text )
	
	def __len__(self):
		return len(self.text)

	def __deepcopy__( self, memo ):
		new = ExtComponent( src=self.src, start=self.start, size=self.size, strand=self.strand, src_size=self._src_size, text=self.text, num_n_bases=self.num_n_bases, plus_strand_start=self.plus_strand_start )
		new._alignment = self._alignment
		new.quality = self.quality
		new.synteny_left = self.synteny_left
		new.synteny_right = self.synteny_right
		new.synteny_empty = self.synteny_empty
		new.empty = self.empty
		new.index = self.index 
		new.gap_less = self.gap_less
		new.num_gap_bases = self.num_gap_bases
		new.num_n_bases = self.num_n_bases		    
		return new

        
class AlignmentBlock(core.Alignment):
	#def __init__( self, score=0, attributes={}, species_to_lengths=None ):
	def __init__(self, ali):
		#self.sortby = sortby
		self.score = ali.score
		self.text_size = ali.text_size
		self.attributes = ali.attributes
		self.species_to_lengths = ali.species_to_lengths
		self.components = []
		#a list of all the assemblies connected with this node
		self.sources = []
		#A list of all sequences in this node
		self.seqs = []
		#I'll expose this for sorting purposes
		self.src_start = []
		#Whether sequences connected with this node should be reverse complemented or not
		self.rev_comp = False
		self.lg = ''			#Can be used as a help when traversing graph
		self.marked = False
		self.temp_marked = False
		self.bestNeighbor = None
		self.bestForwardNeighbor = None
		self.bestBackwardNeighbor = None
		self.forward_neighbors = set()
		self.backward_neighbors = set()
		self.neighbors = set()
		for c in ali.components:
			comp = ExtComponent()
			comp.src = c.src
			self.seqs.append(c.src)
			#Hack, should work on all mine assemblies
			#print c.src
			asm, seq = c.src.split(".")
			self.sources.append(asm)
			comp.start = c.start
			comp.size = c.size
			comp.strand = c.strand
			if c.strand == '-':
				comp.plus_strand_start = c.src_size - c.start - c.size
			else:
				 comp.plus_strand_start = c.start 
			comp.src_size = c.src_size
			comp.text = c.text
			comp.calc_n_bases()
			comp.calc_gap_bases()
			comp.gap_less()
			self.components.append(comp)
			self.src_start.append([asm, seq, c.start])
			
		#Since all components, the aligned sequences, are the same length I can do this. 
		#more sophisticated later on.
		#print len(self.components[0].text)
		self.size = len(self.components[0].text)
			
	def __str__( self ):
		s = "a score=" + str( self.score )
		for key in self.attributes:
			s += " %s=%s" % ( key, self.attributes[key] )
		s += "\n"
		# Components
		for c in self.components:
			s += str( c )
			s += "\n"
		return s
		
	def print_without_seq(self):
		s = "a score=" + str( self.score )
		for key in self.attributes:
			s += " %s=%s" % ( key, self.attributes[key] )
		s += "\n"
		# Components
		for c in self.components:
			s += c.print_without_seq()
			s += "\n"
		return s		
		
	def addToNeighborsSet(self, newNeighbors):
		"""
		Accepts a tuple of neighbors, and adds forward and backward neighbors if they exist
		Forward means that one of the sequence in the alignment had higher coordinates
		than one of the sequences in this alignment
		"""
		if newNeighbors[0]:
			self.backward_neighbors.add(newNeighbors[0])
			self.neighbors.add(newNeighbors[0])
		if newNeighbors[1]:
			self.forward_neighbors.add(newNeighbors[1])
			self.neighbors.add(newNeighbors[1])

	def getArbitraryPath(self):
		"""
		Returns a path.
		"""
		if self.marked and len(self.neighbors) == 0:
			return []
		elif len(self.neighbors) == 0:
			self.marked = True
			return [self]
		
		elif self.marked:
			neighbor = self.neighbors.pop()
			path = neighbor.getArbitraryPath()
			return path
		
		else:
			self.marked = True
			neighbor = self.neighbors.pop()
			path = neighbor.getArbitraryPath()
			return [self] + path		


	def get_sequence_with_fewest_gaps(self, seq):
		"""
		TODO: record the one chosen, and the ones not chosen
		"""
		#index_of_lowest_fraction = 0
		#lowest_fraction = 1.0
		index_of_lowest_amount = 0
		lowest_amount = 10000000
		strand = ''
		for i, c in enumerate(self.components):
			amount = c.num_gap_bases
			if amount < lowest_amount:
				index_of_lowest_amount = i
				lowest_amount = amount
			if c.src == seq and c.strand == '-':
				#printDebug("set rev_comp on this node, in get_sequence_with_fewest_gaps")
				#printDebug(self.print_without_seq()) 
				self.rev_comp = True
	
		if self.rev_comp:
			#Gap_less in this case is without the '-' the aligner puts in
			return (revcomp(self.components[index_of_lowest_amount].gap_less), self.components[index_of_lowest_amount].src, lowest_amount)
		else:
			return (self.components[index_of_lowest_amount].gap_less, self.components[index_of_lowest_amount].src, lowest_amount)
		
		
	def get_lowest_fraction_of_gaps(self):
		"""
		Returns the lowest fraction of gaps of all components. Does not return the component
		itself.
		"""
		index_of_lowest_fraction = 0
		lowest_fraction = 1.0
		for i, c in enumerate(self.components):
			fraction = float(c.num_gap_bases)/c.size
			if fraction < lowest_fraction:
				lowest_fraction = fraction
		return lowest_fraction
	
	def get_lowest_amount_of_gaps(self):
		"""
		Returns the lowest amount of gaps of all components. Does not return the component
		itself.
		"""
		index_of_lowest_amount = 0
		lowest_amount = 10000000
		for i, c in enumerate(self.components):
			amount = c.num_gap_bases
			if amount < lowest_amount:
				lowest_amount = amount
		return lowest_amount


	def get_sequence_from(self, asm_output, seq, asm_output_gap_fraction):
		"""
		TODO: record the one chosen, and the ones not chosen
		"""
		
		index_of_lowest_fraction = -1
		index_of_asm_output = -1
		lowest_fraction = 1.0
		index_of_lowest_amount = -1
		lowest_amount = 10000000
		strand = ''
		asm_amount = 10000000
		for i, c in enumerate(self.components):
			#fraction = float(c.num_gap_bases)/c.size
			amount = c.num_gap_bases
			printDebug("c.src: %s, num_gap_bases: %d" % (c.src,  c.num_gap_bases))
			if amount < lowest_amount:
				lowest_amount = amount
				index_of_lowest_amount = i	
			#print "Sequence: %s, num gaps: %d\n" % (c.src, c.num_gap_bases)
			#if fraction <= lowest_fraction:
			#	index_of_lowest_fraction = i
			#	lowest_fraction = fraction
			if asm_output in c.src:
				#print "asm_output (%s) in c.src (%s)" % (asm_output, c.src)
				index_of_asm_output = i
				asm_amount = amount
			if c.src == seq and c.strand == '-': 
				printDebug("set rev_comp on this node, in get_sequence_from")
				printDebug(self.print_without_seq()) 			
				self.rev_comp = True			
			#if c.src == seq:
				#print "seq: %s, strand: %s" % (seq, c.strand)
			#	strand = c.strand	
		
				
		if index_of_asm_output != -1:
			#if strand == '-':
			#asm_output_gap_fraction has to be set, and has to be less than asm_amount/lowest_amount and lowest_amount has to be more than 500.
			if asm_output_gap_fraction and ((asm_amount/(float(lowest_amount) + 0.01)) > asm_output_gap_fraction) and asm_amount > 500:
				printDebug("asm_output_gap_fraction, %g, is less than %g" % (asm_output_gap_fraction, (asm_amount/(float(lowest_amount) + 0.01))))
				if self.rev_comp:
					return (revcomp(self.components[index_of_lowest_amount].gap_less), self.components[index_of_lowest_amount].src, lowest_amount)
				else:
					return (self.components[index_of_lowest_amount].gap_less, self.components[index_of_lowest_amount].src, lowest_amount)
			else:
				if self.rev_comp:
					#Gap_less in this case is without the '-' the aligner puts in
					return (revcomp(self.components[index_of_asm_output].gap_less), self.components[index_of_asm_output].src, lowest_fraction)
				else:
					#Gap_less in this case is without the '-' the aligner puts in
					return (self.components[index_of_asm_output].gap_less, self.components[index_of_asm_output].src, lowest_fraction)
			#elif index_of_lowest_fraction != -1:
			#Gap_less in this case is without the '-' the aligner puts in
		#	if self.components[index_of_lowest_fraction].strand == '-':
		#		return (revcomp(self.components[index_of_lowest_fraction].gap_less), self.components[index_of_lowest_fraction].src, lowest_fraction)
		#	else:
		#		return (self.components[index_of_lowest_fraction].gap_less, self.components[index_of_lowest_fraction].src, lowest_fraction)
		elif index_of_lowest_amount != -1:
			#Gap_less in this case is without the '-' the aligner puts in
			#if strand == '-':
			if self.rev_comp:
				return (revcomp(self.components[index_of_lowest_amount].gap_less), self.components[index_of_lowest_amount].src, lowest_amount)
			else:
				return (self.components[index_of_lowest_amount].gap_less, self.components[index_of_lowest_amount].src, lowest_amount)
		else:
			print "Shoot, none of these components were good enough:"
			for i, c in enumerate(self.components):
				print "Index: %d" % i
				print c.print_without_seq()				

	def calc_n_bases(self):
		"""
		For all components, calculates the number of N bases, that is, gap bases.
		"""
		for c in self.components:
			c.calc_n_bases()

	def get_n_bases(self):
		"""
		Returns the number of N bases, gap bases, for each component.
		"""
		s = "a score=" + str( self.score )
		for key in self.attributes:
			s += " %s=%s" % ( key, self.attributes[key] )
		s += "\n"
		for c in self.components:
			s += "Num N: %d " % c.num_n_bases
			s += str( c )
			s += "\n"
		return s
		

	def get_gap_less(self):
		"""
		Returns the gapless version of the alignment.
		"""
		s = "a score=" + str( self.score )
		for key in self.attributes:
			s += " %s=%s" % ( key, self.attributes[key] )
		s += "\n"
		for c in self.components:
			s += "s %s %d %d %s %d %s" % ( c.src, c.start,
                                             c.size, c.strand,
                                             c.src_size, c.gap_less )
			s += "\n"
		
		return s
		
	def get_component( self, src ):
		for c in self.components:
			if c.src == src: return c
		return None

	def __deepcopy__( self, memo ):
		from copy import deepcopy
		new = AlignmentBlock(self)

		new.sources = self.sources
		new.seqs = self.seqs
		new.src_start = self.src_start
		new.rev_comp = self.rev_comp 
		new.marked = self.marked
		new.temp_marked = self.temp_marked 
		new.bestNeighbor = deepcopy(self.bestNeighbor)
		new.bestForwardNeighbor = deepcopy(self.bestForwardNeighbor)
		new.bestBackwardNeighbor = deepcopy(self.bestBackwardNeighbor)
		new.forward_neighbors = self.forward_neighbors
		new.backward_neighbors = self.backward_neighbors
		new.neighbors = self.neighbors
		new.size = self.size
		new.components = [] 
		for component in self.components:
			new.add_component( deepcopy( component ) )  
		return new

# doesn't handle "U" in RNA sequences
def revcomp(seq):
	#From last-417, maf-swap.py	
	complement = string.maketrans('ACGTNSWRYKMBDHVacgtnswrykmbdhv', 'TGCANSWYRMKVHDBtgcanswyrmkvhdb')                          
	
	return seq[::-1].translate(complement)

def getSortedVersionBasedOnPositionInSuppliedSeq(nodeList, seq):
	
	tempList = []
	tempList = [(node.get_component(seq).plus_strand_start, node) for node in nodeList]

	tempList.sort()
	
	return [node for start, node in tempList]
	
def filterNodeList(nodeList, seq):
	tempList = []
		
	#print "Filtering %s" % seq	
	prev_start = nodeList[0].get_component(seq).plus_strand_start
	prev_size = nodeList[0].get_component(seq).size
	for index, node in enumerate(nodeList):
		#print node.print_without_seq()
		#print "prev_start: %d" % prev_start
		#print "prev_size: %d" % prev_size
		#print node.get_component(seq).plus_strand_start
		#Skipping the node if it overlaps with previous and not the first element of the list
		if node.get_component(seq).plus_strand_start != (prev_start + prev_size) and index != 0:
			#skip the node
			#print "Skipping node of index %d:" % index
			#print node.print_without_seq()
			continue
		prev_start = node.get_component(seq).plus_strand_start
		prev_size = node.get_component(seq).size
		tempList.append(node)
		
	return tempList
				
	
#def getPath(sortedNodeListCollection, start, end, seqToFollow, size_of_node, nodesUsedAndConsidered):
def getPath(sortedNodeListCollection, start, end, seqToFollow, size_of_node, nodesUsed, asm_output):
	
	nodesConsidered = set()
	
	#size of node is used to get an approximately distance of the path to find
	
	#First, find which sequences common to both start and end
	commonSeq = set(start.seqs).intersection(end.seqs)
	printDebug("commonSeq:")
	printDebug(commonSeq)
	#sortedNodeListCollection[c.src], c.plus_strand_start, c.src
	
	startNodeStrand = start.get_component(seqToFollow).strand
	endNodeStrand = end.get_component(seqToFollow).strand
	
	printDebug("start node: ")
	printDebug(start.print_without_seq())
	
	printDebug("end node: ")
	printDebug(end.print_without_seq())
	printDebug("Size of gap node: %d" % size_of_node)
		
	path = []
	preferredPath = None
	preferredLength = 0
	preferredCost = 100000	
	cheapestPath = []
	cheapestLength = 0 
	cheapestCost = 1000000
	
	#Get the size of the sequence we're following, so we can avoid nodes with larger sequence
	seqToFollowLength = start.get_component(seqToFollow).src_size
	asmToFollow = seqToFollow.split('.')[0]
	foundLongerSequence = False
	foundNodeInLongPath = False
	falsePulling = False
	usedNode = False
	
	#Find the cheapest path
	for seq in commonSeq:
		#print "Seq: %s" % seq
		printDebug("Seq: %s" % seq)
		#print start.get_component(seq).plus_strand_start
		#print end.get_component(seq).plus_strand_start
		start_pos = binary_search(sortedNodeListCollection[seq], start.get_component(seq).plus_strand_start, seq)
		printDebug("Start pos: %d" % start_pos)
		end_pos = binary_search(sortedNodeListCollection[seq], end.get_component(seq).plus_strand_start, seq)
		printDebug("End pos: %d" % end_pos)
		
		#Skipping paths more than 10 nodes (20 nodes?)
		if abs(start_pos - end_pos) > 10:
			printDebug("path involving common seq %s is longer than 10 nodes, skipping" % seq)
			continue
		
		if start_pos < end_pos:
			path = sortedNodeListCollection[seq][start_pos+1:end_pos]
			
			#print "path:"		
			#Checking length of path here, and skipping if too long 	
			length = 0	
			for node in path:
				length = length + node.size
				#print node.print_without_seq()
				if node in nodesUsed:
					printDebug("node already used. Skipping this path...")
					#print node.print_without_seq()	
					usedNode = True
					
				#Do not touch nodes connected with longer sequence from same assembly. Safety measure I hope
				for comp in node.components:
					#Skip if a node is part of another path with more than 5 nodes
					if asmToFollow in comp.src and len(sortedNodeListCollection[comp.src]) > 5 and comp.src != seqToFollow:
						printDebug("touching nodes from a path longer than 5 nodes (length: %d): " % len(sortedNodeListCollection[comp.src]))
						#print node.print_without_seq()
						foundNodeInLongPath = True
					
					#I get a problem with these paths pulling in other parts of seqToFollow and messing up. Need to stop it	
					if seq != seqToFollow and seqToFollow in comp.src:
						printDebug("falsePulling")
						falsePulling = True	
							
					#if asmToFollow in comp.src and comp.src_size > seqToFollowLength:
					#	print "touching nodes with sequence larger than the one we are following (%d): " % seqToFollowLength
					#	print node.print_without_seq()
					#	foundLongerSequence = True
			
			if (length >= size_of_node*10):
				printDebug("Length is more than 10 times the size of gap node")
				usedNode = False
				falsePulling = False
				foundNodeInLongPath = False
				continue
			
			#Can't skip the main loop from inside the loops above, but can here.			
			#if foundLongerSequence:
			#	foundLongerSequence = False
			#	continue
			if foundNodeInLongPath:
				printDebug("foundNodeInLongPath")
				usedNode = False
				falsePulling = False
				foundNodeInLongPath = False
				continue
						 
			if usedNode:
				usedNode = False
				falsePulling = False
				foundNodeInLongPath = False
				printDebug("usedNode")
				continue
				
			if falsePulling:
				usedNode = False
				falsePulling = False
				foundNodeInLongPath = False
				printDebug("falsePulling")
				continue	
			
			for node in path:
				#this is a bit tricky. When start is less than end, this sequence is on the 
				#same strand as seqToFollow. However, if one node and its sequence is on minus
				#strand, we add another component from seqToFollow to be sure it is reverse
				#complemented
				if node.get_component(seq).strand == '-':
					node.rev_comp = True
					#printDebug("made this node rev comp, in getPath:")
					#printDebug(node.print_without_seq())

		else:
			#Shit, I think these sequences need to be reverse complemented...
			path = list(reversed(sortedNodeListCollection[seq][end_pos+1:start_pos]))
			#if start.get_component(seqToFollow).strand == start.get_component(seq).strand:

			#print "path:"

			#Checking length of path here, and skipping if too long
			length = 0	
			for node in path:
				length = length + node.size
				#print node.print_without_seq()
				#Do not touch nodes already used (in a longer sequence)
				if node in nodesUsed:
					printDebug("node already used. Skipping this path...")
					#print node.print_without_seq()	
					usedNode = True
								
				#Do not touch nodes connected with longer sequence from same assembly. Safety measure I hope. 
				for comp in node.components:
					#Skip if a node is from the same assembly, part of another path with more than 5 nodes and not part of the sequence we're following
					if asmToFollow in comp.src and len(sortedNodeListCollection[comp.src]) > 5 and comp.src != seqToFollow:
						printDebug("touching nodes from a path longer than 5 nodes (length: %d): " % len(sortedNodeListCollection[comp.src]))
						#print node.print_without_seq()
						foundNodeInLongPath = True

					#I get a problem with these paths pulling in other parts if seqToFollow and messing up. Need to stop it	
					if seq != seqToFollow and seqToFollow in comp.src:
						printDebug("falsePulling")
						falsePulling = True	
													
					#if asmToFollow in comp.src and comp.src_size > seqToFollowLength:
					#	print "touching nodes with sequence larger than the one we are following (%d): " % seqToFollowLength
					#	print node.print_without_seq()
					#	foundLongerSequence = True
			
			if (length >= size_of_node*10):
				printDebug("Length is more than 10 times the size of gap node")
				usedNode = False
				falsePulling = False
				foundNodeInLongPath = False				
				continue
			
			#Can't skip the main loop from inside the loops above, but can here.			
			#if foundLongerSequence:
			#	foundLongerSequence = False
			#	continue
			if foundNodeInLongPath:
				usedNode = False
				falsePulling = False
				foundNodeInLongPath = False
				printDebug("foundNodeInLongPath")
				continue
				
			#Skip path if node used
			if usedNode:
				usedNode = False
				falsePulling = False
				foundNodeInLongPath = False
				printDebug("usedNode")
				continue
		
			if falsePulling:
				usedNode = False
				falsePulling = False
				foundNodeInLongPath = False
				printDebug("falsePulling")
				continue				
				
			for node in path:
				#this is also a bit tricky. When end is less than start, this sequence is on the 
				#opposite strand as seqToFollow. If one node and its sequence is on plus
				#strand, we add another component from seqToFollow to be sure it is reverse
				#complemented. We deal with this later
				if node.get_component(seq).strand == '+':
					node.rev_comp = True
					#printDebug("made this node rev comp, in getPath:")
					#printDebug(node.print_without_seq())				

		cost = 0
		length = 0			
		for node in path:
			print node.print_without_seq()
			nodesConsidered.add(node)
			length = length + node.size
			cost = cost + node.get_lowest_amount_of_gaps()
			
		#Special for Newbler. Gaps of 25 is when Newbler can't place the contigs end properly
		if cost == 25 and length == 25:
			cost = cost + 1000
		#Special for Celera. Gaps of 20 is when Celera can't place the contigs end properly
		if cost == 20 and length == 20:
			cost = cost + 1000					
		printDebug("cost: %d, length: %d, size of potential gap %d, fraction length/size %.2g " % (cost, length, size_of_node,  length/float(size_of_node+0.01)))

		if asm_output: 
			if asm_output in seq:
				preferredCost = cost
				preferredPath = path
				preferredLength = length	
	
		#Want the path with the least amount of gaps, but only if it is within 0.1 - 10 times the size of 
		#potential gap. It seems, at least with PacBio-based assembly, that a lot of gaps are almost 
		#overfilled, constricting the accepted sizes might help.
		#if (cost < cheapestCost):
		#if (cost < cheapestCost) and (length > size_of_node*0.1) and (length < size_of_node*10) :
		if (cost < cheapestCost):	
			cheapestCost = cost
			cheapestPath = path
			cheapestLength = length
	

		
	if preferredPath != None:
		printDebug("preferred cost %d, length %d, fraction length/size %.2g" % (preferredCost, preferredLength, preferredLength/float(size_of_node+0.01)))
		for node in preferredPath:
			#print node.print_without_seq()	
			nodesUsed.add(node)	
	elif cheapestPath != None:
		printDebug("cheapest cost %d, length %d, fraction length/size %.2g" % (cheapestCost, cheapestLength,  cheapestLength/float(size_of_node+0.01)))
		for node in cheapestPath:
			#print node.print_without_seq()	
			nodesUsed.add(node)	
	else:
		printDebug("couldn't find a proper path")
		
		
	#return cheapestPath, nodesUsedAndConsidered
	#if preferredPath != None:
	#	return preferredPath, nodesUsedAndConsidered
	#else:
	#	return cheapestPath, nodesUsedAndConsidered
			
	#return cheapestPath, nodesUsedAndConsidered
	if preferredPath != None:
		return preferredPath, nodesUsed, nodesConsidered
	else:
		return cheapestPath, nodesUsed, nodesConsidered		
		
def checkPath(path, nodesUsed, asmToFollow, seqToFollow, sortedNodeListCollection, seq):

	foundLongerSequence = False
	foundNodeInLongPath = False
	falsePulling = False
	usedNode = False

	#Checking length of path here, and skipping if too long
	for node in path:
		#print node.print_without_seq()
		#Do not touch nodes already used (in a longer sequence)
		#if node in nodesUsed:
			#print "node already used. Skipping this path..."
			#print node.print_without_seq()	
			#usedNode = True
			
		#Do not touch nodes connected with longer sequence from same assembly. Safety measure I hope. 
		for comp in node.components:
			#Skip if a node is from the same assembly, part of another path with more than 5 nodes and not part of the sequence we're following
			if asmToFollow in comp.src and len(sortedNodeListCollection[comp.src]) > 10 and comp.src != seqToFollow:
				printDebug("touching nodes from a path longer than 10 nodes (length: %d)" % len(sortedNodeListCollection[comp.src]))
				#print node.print_without_seq()
				foundNodeInLongPath = True

			#I get a problem with these paths pulling in other parts of seqToFollow and messing up. Need to stop it	
			if seq != seqToFollow and seqToFollow in comp.src:
				printDebug("falsePulling")
				falsePulling = True	

	return foundNodeInLongPath or usedNode or falsePulling
	
		
def getConnectionPath(sortedNodeListCollection, commonSeq, start, end, seqToFollow, whichEnd, firstIsLast, asm_output):

	#Several operations might touch the nodes in a connectionPath, so I need to deepcopy them
	#before doing anything. Should work


	nodesUsedAndConsidered = set()
	nodesUsed = set()
	
	print "trying getConnectionPath"
	printDebug("trying getConnectionPath")
	#First, find which sequences common to both start and end
	#commonSeq = set(start.seqs).intersection(end.seqs)
	
	printDebug("commonSeq:")
	printDebug(commonSeq)
	
	printDebug("start node: ")
	printDebug(start.print_without_seq())
	
	printDebug("end node: ")
	printDebug(end.print_without_seq())
	
	printDebug("seqToFollow: %s" % seqToFollow)
	printDebug("whichEnd is %s" % whichEnd)
		
	path = []
	preferredPath = None
	preferredLength = 0
	preferredCost = 100000
	cheapestPath = None
	cheapestLength = 0 
	cheapestCost = 1000000
	
	#Get the size of the sequence we're following, so we can avoid nodes with larger sequence
	
	seqToFollowLength = start.get_component(seqToFollow).src_size
	asmToFollow = seqToFollow.split('.')[0]
	foundLongerSequence = False
	foundNodeInLongPath = False
	falsePulling = False
	usedNode = False
	
	direction = ''
	
	#Find the cheapest path
	for seq in commonSeq:
		printDebug("Seq: %s" % seq)
		#print start.get_component(seq).plus_strand_start
		#print end.get_component(seq).plus_strand_start
		start_pos = binary_search(sortedNodeListCollection[seq], start.get_component(seq).plus_strand_start, seq)
		printDebug("Start pos: %d" % start_pos)
		end_pos = binary_search(sortedNodeListCollection[seq], end.get_component(seq).plus_strand_start, seq)
		printDebug("End pos: %d" % end_pos)
		printDebug("Difference: %d" % abs(start_pos - end_pos))
		printDebug("firstIsLast: %s" % firstIsLast) 
		
		#Skipping paths more than 10 nodes (20 nodes?)
		if abs(start_pos - end_pos) > 20:
			printDebug("path involving common seq %s is longer than 20 nodes, skipping" % seq)
			continue
		
		#Find the direction
		if start.get_component(seqToFollow).strand == start.get_component(seq).strand and whichEnd == 'tail':
			#If same strand and going off tail end, I think this most be rightward direction
			if start_pos < end_pos:
				direction = 'right'
				printDebug("whichEnd is tail and direction is right")
			elif firstIsLast:
				direction = 'left'
				printDebug("whichEnd is tail and direction is left, firstIsLast")
			else:
				printDebug("whichEnd is tail and direction is left, not valid in this case")	
				continue
								
		elif start.get_component(seqToFollow).strand != start.get_component(seq).strand and whichEnd == 'tail':
			#when not same strand, the direction needs to be oposite to be correct direction
			if start_pos > end_pos:
				direction = 'right'
				printDebug("whichEnd is tail and direction is right")
			elif firstIsLast:
				direction = 'left'
				printDebug("whichEnd is tail and direction is left, firstIsLast")
			else:
				printDebug("whichEnd is tail and direction is left, not valid in this case")
				continue					
								
		elif start.get_component(seqToFollow).strand == start.get_component(seq).strand and whichEnd == 'head':
			#If same strand and going off head end, start pos most be higher than end pos
			if start_pos > end_pos:
				direction = 'left'
				printDebug("whichEnd is head and direction is left")
			elif firstIsLast:
				direction = 'right'
				printDebug("whichEnd is head and direction is right, firstIsLast")	
			else:
				printDebug("whichEnd is head and direction is right, not valid in this case")
				continue				
								
		elif start.get_component(seqToFollow).strand != start.get_component(seq).strand and whichEnd == 'head':
			#
			if start_pos < end_pos:
				direction = 'left'
				printDebug("whichEnd is head and direction is left")
			elif firstIsLast:
				direction = 'right'
				printDebug("whichEnd is head and direction is right, firstIsLast")
			else:
				printDebug("whichEnd is head and direction is right, not valid in this case")
				continue									
		else: 
			printDebug("Direction seems to be wrong way compared to whichEnd")
			continue
					
	
		if whichEnd == 'tail' and direction == 'right':
			#printDebug("whichEnd is tail")

			if start_pos < end_pos:
				printDebug("start_pos is less than end_pos")
				path = copy.deepcopy(sortedNodeListCollection[seq][start_pos+1:end_pos])
				
				if checkPath(path, nodesUsed, asmToFollow, seqToFollow, sortedNodeListCollection, seq):
					printDebug("checkPath returned true, skipping")
					continue
	
				for node in path:
					#this is a bit tricky. When start is less than end, this sequence is on the 
					#same strand as seqToFollow. However, if one node and its sequence is on minus
					#strand, we add another component from seqToFollow to be sure it is reverse
					#complemented
					if node.get_component(seq).strand == '-':
						node.rev_comp = True
						printDebug("made this node rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq())
					elif node.get_component(seq).strand == '+':
						node.rev_comp = False
						printDebug("made this node not rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq()) 						
						
			elif start_pos > end_pos:
				printDebug("start_pos is more than end_pos")
				#Shit, I think these sequences need to be reverse complemented...
				path = copy.deepcopy(list(reversed(sortedNodeListCollection[seq][end_pos+1:start_pos])))
				#if start.get_component(seqToFollow).strand == start.get_component(seq).strand:


				if checkPath(path, nodesUsed, asmToFollow, seqToFollow, sortedNodeListCollection, seq):
					printDebug("checkPath returned true, skipping")
					continue			
		
				for node in path:
					#this is also a bit tricky. When end is less than start, this sequence is on the 
					#opposite strand as seqToFollow. If one node and its sequence is on plus
					#strand, we add another component from seqToFollow to be sure it is reverse
					#complemented. We deal with this later
					if node.get_component(seq).strand == '+':
						node.rev_comp = True
						printDebug("made this node rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq())
					elif node.get_component(seq).strand == '-':
						node.rev_comp = False
						printDebug("made this node not rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq()) 						
												
								
		elif whichEnd == 'head' and direction == 'left':
			#Should only enter here when index of first is same as last
			#printDebug("whichEnd is head")

			if start_pos < end_pos:
				printDebug("start_pos is less than end_pos")
				#need to reverse complement these
				path = copy.deepcopy(list(reversed(sortedNodeListCollection[seq][start_pos+1:end_pos])))

				if checkPath(path, nodesUsed, asmToFollow, seqToFollow, sortedNodeListCollection, seq):
					continue
				
				for node in path:
					#this is also a bit tricky. When end is less than start, this sequence is on the 
					#opposite strand as seqToFollow. If one node and its sequence is on plus
					#strand, we add another component from seqToFollow to be sure it is reverse
					#complemented. We deal with this later
					if node.get_component(seq).strand == '+':
						node.rev_comp = True
						printDebug("made this node rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq())
					elif node.get_component(seq).strand == '-':
						node.rev_comp = False
						printDebug("made this node not rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq()) 																	

						
			elif start_pos > end_pos:
				printDebug("start_pos is more than end_pos")
				path = copy.deepcopy(sortedNodeListCollection[seq][end_pos+1:start_pos])

				if checkPath(path, nodesUsed, asmToFollow, seqToFollow, sortedNodeListCollection, seq):
					printDebug("checkPath returned true, skipping")
					continue
	
				for node in path:
					#this is a bit tricky. When start is less than end, this sequence is on the 
					#same strand as seqToFollow. However, if one node and its sequence is on minus
					#strand, we add another component from seqToFollow to be sure it is reverse
					#complemented
					if node.get_component(seq).strand == '-':
						node.rev_comp = True
						printDebug("made this node rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq())
					elif node.get_component(seq).strand == '+':
						node.rev_comp = False
						printDebug("made this node not rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq()) 							
		
		
		#The stuff below only happens when first is last				
		elif whichEnd == 'head' and direction == 'right':
				
			#printDebug("whichEnd is tail")

			if start_pos < end_pos:
				printDebug("start_pos is less than end_pos")
				#need to reverse complement these
				path = copy.deepcopy(list(reversed(sortedNodeListCollection[seq][start_pos+1:end_pos])))


				if checkPath(path, nodesUsed, asmToFollow, seqToFollow, sortedNodeListCollection, seq):
					printDebug("checkPath returned true, skipping")
					continue

				for node in path:
					#this is also a bit tricky. When end is less than start, this sequence is on the 
					#opposite strand as seqToFollow. If one node and its sequence is on plus
					#strand, we add another component from seqToFollow to be sure it is reverse
					#complemented. We deal with this later
					if node.get_component(seq).strand == '+':
						node.rev_comp = True
						printDebug("made this node rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq())											
					elif node.get_component(seq).strand == '-':
						node.rev_comp = False
						printDebug("made this node not rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq()) 

			elif start_pos > end_pos:
				printDebug("start_pos is more than end_pos")
				path = copy.deepcopy(sortedNodeListCollection[seq][end_pos+1:start_pos])

				if checkPath(path, nodesUsed, asmToFollow, seqToFollow, sortedNodeListCollection, seq):
					printDebug("checkPath returned true, skipping")
					continue
	
				for node in path:
					#this is a bit tricky. When start is less than end, this sequence is on the 
					#same strand as seqToFollow. However, if one node and its sequence is on minus
					#strand, we add another component from seqToFollow to be sure it is reverse
					#complemented
					if node.get_component(seq).strand == '-':
						node.rev_comp = True
						printDebug("made this node rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq())
					elif node.get_component(seq).strand == '+':
						node.rev_comp = False
						printDebug("made this node not rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq()) 												
								
		elif whichEnd == 'tail' and direction == 'left':
			#printDebug("whichEnd is head")

			if start_pos < end_pos:
				printDebug("start_pos is less than end_pos")
				path = copy.deepcopy(sortedNodeListCollection[seq][start_pos+1:end_pos])
				
				if checkPath(path, nodesUsed, asmToFollow, seqToFollow, sortedNodeListCollection, seq):
					printDebug("checkPath returned true, skipping")
					continue
	
				for node in path:
					#this is a bit tricky. When start is less than end, this sequence is on the 
					#same strand as seqToFollow. However, if one node and its sequence is on minus
					#strand, we add another component from seqToFollow to be sure it is reverse
					#complemented
					if node.get_component(seq).strand == '-':
						node.rev_comp = True
						printDebug("made this node rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq())
					elif node.get_component(seq).strand == '+':
						node.rev_comp = False
						printDebug("made this node not rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq()) 						
						
			elif start_pos > end_pos:
				printDebug("start_pos is more than end_pos")
				#Shit, I think these sequences need to be reverse complemented...
				path = copy.deepcopy(list(reversed(sortedNodeListCollection[seq][end_pos+1:start_pos])))
				#if start.get_component(seqToFollow).strand == start.get_component(seq).strand:


				if checkPath(path, nodesUsed, asmToFollow, seqToFollow, sortedNodeListCollection, seq):
					continue			
		
				for node in path:
					#this is also a bit tricky. When end is less than start, this sequence is on the 
					#opposite strand as seqToFollow. If one node and its sequence is on plus
					#strand, we add another component from seqToFollow to be sure it is reverse
					#complemented. We deal with this later
					if node.get_component(seq).strand == '+':
						node.rev_comp = True
						printDebug("made this node rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq())
					elif node.get_component(seq).strand == '-':
						node.rev_comp = False
						printDebug("made this node not rev comp, in getConnectionPath:")
						printDebug(node.print_without_seq()) 													
									
		cost = 0
		length = 0			
		for node in path:
			printDebug(node.print_without_seq())
			nodesUsedAndConsidered.add(node)
			length = length + node.size
			cost = cost + node.get_lowest_amount_of_gaps()
			
		#Special for Newbler. Gaps of 25 is when Newbler can't place the contigs end properly
		if cost == 25 and length == 25:
			cost = cost + 1000
		#Special for Celera. Gaps of 20 is when Celera can't place the contigs end properly
		if cost == 20 and length == 20:
			cost = cost + 1000		
		printDebug( "cost: %d, length: %d " % (cost, length))
		
		if asm_output: 
			if asm_output in seq:
				preferredCost = cost
				preferredPath = path
				preferredLength = length
	
		if (cost < cheapestCost):
		#if (cost < cheapestCost) and (length > size_of_node*0.1) and (length < size_of_node*10) :
		#if (cost < cheapestCost and path):	
			cheapestCost = cost
			cheapestPath = path
			cheapestLength = length
		
		direction = ''	

	if preferredPath != None:
		printDebug("preferred cost %d, length %d" % (preferredCost, preferredLength))
		for node in preferredPath:
			#print node.print_without_seq()	
			nodesUsed.add(node)	
	elif cheapestPath != None:
		printDebug("cheapest cost %d, length %d" % (cheapestCost, cheapestLength))
		for node in cheapestPath:
			#print node.print_without_seq()	
			nodesUsed.add(node)
	else:
		printDebug("couldn't find a proper path")
		
	if preferredPath != None:
		printDebug("preferred cost %d, length %d" % (preferredCost, preferredLength))
			
	#return cheapestPath, nodesUsedAndConsidered
	if preferredPath != None:
		return preferredPath, nodesUsed
	else:
		return cheapestPath, nodesUsed
		
	
def binary_search(a, x, seq, lo=0, hi=None):
	#Modified from http://hg.python.org/cpython/file/2.7/Lib/bisect.py
    """Return the index where to insert item x in list a, assuming a is sorted.

    The return value i is such that all e in a[:i] have e < x, and all e in
    a[i:] have e >= x.  So if x already appears in the list, a.insert(x) will
    insert just before the leftmost x already there.
    	
    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    
    a is a list of nodes, x is the position (plus_strand_start), seq is the scaffold we are interested in
    """
    
    if lo < 0:
    	raise ValueError("lo must be non-negative")
    hi = hi if hi is not None else len(a)
    #print "len(a): %d" % len(a)
    #print "hi: %d" % hi
    
    while lo < hi:
        mid = (lo+hi)//2
        #sys.stderr.write("lo: %d, " %  lo) 
        #sys.stderr.write("hi: %d, " % hi)
        #sys.stderr.write("mid: %d \n" % mid)
        #Only special case
        #sys.stderr.write("a[mid].get_component(seq).plus_strand_start: %d, " % a[mid].get_component(seq).plus_strand_start)
        #sys.stderr.write("a[mid].get_component(seq).src_size: %d\n" % a[mid].get_component(seq).src_size)
        if a[mid].get_component(seq).plus_strand_start < x: lo = mid+1
        else: hi = mid
    return lo

def findProximalNodesByBinarySearch(sortedNodeList, sequencePosition, seq):
	
	#print sortedNodeList, sequencePosition, seq
	
	#newNeighbors = []
	
	#print "seq: %s" % seq
	#print "sequencePosition: %d" % sequencePosition
	#print "len(sortedNodeList): %d" % len(sortedNodeList)
	#sys.stderr.write("seq: %s " % seq)
	#sys.stderr.write("sequencePosition: %d " % sequencePosition)
	#sys.stderr.write("len(sortedNodeList): %d \n" % len(sortedNodeList))
	#Find position in the list
	pos = binary_search(sortedNodeList, sequencePosition, seq)
	
	
	#print "pos: %d" % pos
	#sys.stderr.write("pos: %d \n" % pos)	
	
	#Because I filter the lists, remove overlapping alignments, the two cases below can happen.
	#They happen because the node/alignment is gone from the list.
	if pos >= len(sortedNodeList):
		#print "Position exceeds or is equal to length of sortedNodeList, pos: %d" % pos
		return (None, None)
	elif sortedNodeList[pos].get_component(seq).plus_strand_start != sequencePosition:
		#print "Did not find exact match, plus_strand_start: %d, seqPos: %d" % (sortedNodeList[pos].get_component(seq).plus_strand_start, sequencePosition)
		#print sortedNodeList[pos].print_without_seq()
		return (None, None)
		
	if (pos == 0 and len(sortedNodeList) == 1):
		#Got no neigbors
		return (None, None)
	elif pos == 0:
		#First element
		return (None, sortedNodeList[pos + 1])
	elif pos == (len(sortedNodeList) - 1):
		#Last element
		return (sortedNodeList[pos - 1], None)
	else:
		return (sortedNodeList[pos - 1], sortedNodeList[pos + 1])

	#return newNeigbors

def printSeq(s, n, prefix, header):
	"""
	Prints to prefix.fasta in blocks n bases long.
	"""
	with open(prefix + ".fasta", "a") as fasta_file:
		#First the header
		fasta_file.write("%s\n" % header)
					
		if (n <= 0):
			fasta_file.write("%s\n" % s)

		else:
			p = 0
			while p < len( s ):
				fasta_file.write("%s\n" % s[p:min(p+n,len(s))])
				p += n        
	
def printDebug(text, endline='\n'):
	"""
	Prints to stderr.
	"""

	sys.stderr.write("DEBUG: %s%s" % (text, endline))


def connectScaffolds(sortedNodeListCollection, asm_path, cov_fraction, prefix, links, numAsms, asm_output, extend, asm_output_gap_fraction):

	#Start by reducing the search space
	reducedPaths = closeGaps(sortedNodeListCollection, asm_path, cov_fraction, numAsms, asm_output, extend)
			
	#sortedBySizeNodeList = [(len(nodeList), seq, nodeList) for seq, nodeList in sortedNodeListCollection.items()]
	
	#sortedBySizeNodeList.sort(reverse=True)
	
	#Hardcoding accepted length of connection path 
	#acceptedLengthPath = 10		
			
	#
	#startAndEndNodes = {}
	startAndEndNodes = []
	connectionGraph = {}
	#tempList = [(node.get_component(seq).plus_strand_start, node) for node in nodeList]

	#for nodeListLen, seq, nodeList in sortedBySizeNodeList:	
	for seq, nodeList in reducedPaths.items():
		
		if asm_path in seq:
			seqPath = SeqPath(seq, reducedPaths[seq])
			connectionGraph[seq] = seqPath
			nodeListLen = len(nodeList)
			printDebug("seq: %s, length in nodes: %d" % (seq, nodeListLen))
			first = [0, 0]
			last = [0, 0]
			start_node = None
			start_index = -1
			end_node = None
			end_index = -1
			
			#I need to rethink this a bit. I can't really have identical indexes anyhow.
			#First, first has to be in first half to be valid.					
			for index, node in enumerate(nodeList):
				#Get the first occurrence of the max number of alignments plus seq has to be present 
				#(if a unplaced sequence has filled a gap, this might not be the case)
				if len(node.components) > first[0] and node.get_component(seq):
					first = [len(node.components), index]
				if index >= (len(nodeList)/2) - 1:
					break
					


			#Also, last can't be in first half of the sequence
			#http://christophe-simonis-at-tiny.blogspot.no/2008/08/python-reverse-enumerate.html
			reverse_enumerate = lambda l: izip(range(len(l)-1, -1, -1), reversed(l))
			for index, node in reverse_enumerate(nodeList):
				if len(node.components) > last[0] and node.get_component(seq):
					last = [len(node.components), index]
				if index <= (len(nodeList)/2):
					break					
				"""	
				#For sequences with several nodes, sometimes only one node have all assemblies in alignment,
				#but if there is one with one less assembly, I want to include that one instead
				if index != first[1]:
					#If the index is not equal, then it's OK
					last = [len(node.components), index]
				elif len(node.components)-1 != last[0]:
					#If it's even less than one less assembly, then it's OK
					last = [len(node.components), index]
				else:
					printDebug("Keeping the previously found entry of %d assemblies at index %d, instead of the last found at %d of %d assemblies" %
					(last[0], last[1], len(node.components), index))
				"""	
			if first[0] == 0:
				#didn't find a good node in first half. Just set it to node 0.
				first[0] = len(nodeList[0].components)
				first[1] = 0
				printDebug("set first to index 0")
			else:
				printDebug("found first at index %d" % first[1])
				
			if last[0] == 0:
				#didn't find a good node in last half. Just set it to the last node. 
				last[0] = len(nodeList[len(nodeList)-1].components)
				last[1] = len(nodeList)-1
				printDebug("set last to index %d" % len(nodeList)-1)
			else:				
				printDebug("found last at index %d" % last[1])	
			
			
			
			
			#startAndEndNodes[seq] = {'start': (set(nodeList[first[1]].seqs), nodeList[first[1]], first[1]), 'end': (set(nodeList[last[1]].seqs), nodeList[last[1]], last[1], 'length': nodeListLen)}
			#Sequence id, length of list of nodes, list consisting of set of sequences, the node affected and the index of that node, and a second list of the same, but for the last position
			startAndEndNodes.append([seq, nodeListLen, [set(nodeList[first[1]].seqs), nodeList[first[1]], first[1]], [set(nodeList[last[1]].seqs), nodeList[last[1]], last[1]]])
			connectionGraph[seq].setTailAndHead(first[1], last[1])
											
	firstIsLast = False
	for s_seq, s_len, s_first, s_last in startAndEndNodes:
		printDebug("s_seq: %s" % s_seq)
		if s_first[2] == s_last[2]:
			printDebug("the index of first node for s_seq is identical to the last node, firstIsLast")
			firstIsLast = True	
		for t_seq, t_len, t_first, t_last in startAndEndNodes:
			if s_seq != t_seq:
				#printDebug("s_seq: %s, t_seq: %s" %(s_seq, t_seq))

				#printDebug("s_first:")
				#printDebug(s_first)
				#printDebug("t_first:")
				#printDebug(t_first)		
				#printDebug("s_last:")
				#printDebug(s_last)
				#printDebug("t_last:")
				#printDebug(t_last)								
				s_s_commonSeq = s_first[0].intersection(t_first[0])
				s_e_commonSeq = s_first[0].intersection(t_last[0])
				e_s_commonSeq = s_last[0].intersection(t_first[0])
				e_e_commonSeq = s_last[0].intersection(t_last[0])
				if len(s_s_commonSeq) >= links:
					printDebug("s_seq: %s, t_seq: %s" %(s_seq, t_seq))
					
					print "%s and %s share these sequences (start - start):" % (s_seq, t_seq)
					print s_s_commonSeq	
					print "node one:" 
					print s_first[1].print_without_seq()
					print "node two:" 
					print t_first[1].print_without_seq()													
					s_s_commonPath, s_s_nodesUsed = getConnectionPath(sortedNodeListCollection, s_s_commonSeq, s_first[1], t_first[1], s_seq, 'head', firstIsLast, asm_output)	
					#if len(s_s_commonPath) <= acceptedLengthPath:
					if s_s_commonPath != None:
						connectionGraph[s_seq].add_head_head([len(s_s_commonSeq), len(s_s_commonPath), s_s_commonPath, t_seq, s_first[2], t_first[2], connectionGraph[t_seq]])
					
				if len(s_e_commonSeq) >= links:
					printDebug("s_seq: %s, t_seq: %s" %(s_seq, t_seq))
					print "%s and %s share these sequences (start - end):" % (s_seq, t_seq)
					print s_e_commonSeq
					print "node one:" 
					print s_first[1].print_without_seq()
					print "node two:" 
					print t_last[1].print_without_seq()														
					s_e_commonPath, s_e_nodesUsed = getConnectionPath(sortedNodeListCollection, s_e_commonSeq, s_first[1], t_last[1], s_seq, 'head', firstIsLast, asm_output)						
					#connectionGraph[s_seq].append(['se', s_e_commonPath, s_seq, t_seq, True, True])
					#if len(s_e_commonPath) <= acceptedLengthPath:
					if s_e_commonPath != None:
						connectionGraph[s_seq].add_head_tail([len(s_e_commonSeq), len(s_e_commonPath), s_e_commonPath, t_seq, s_first[2], t_last[2], connectionGraph[t_seq]])
					
				if len(e_s_commonSeq) >= links:
					printDebug("s_seq: %s, t_seq: %s" %(s_seq, t_seq))
					print "%s and %s share these sequences (end - start):" % (s_seq, t_seq)
					print e_s_commonSeq
					print "node one:" 
					print s_last[1].print_without_seq()
					print "node two:" 
					print t_first[1].print_without_seq()													
					e_s_commonPath, e_s_nodesUsed = getConnectionPath(sortedNodeListCollection, e_s_commonSeq, s_last[1], t_first[1], s_seq, 'tail', firstIsLast, asm_output)												
					#connectionGraph[s_seq].append(['es', e_s_commonPath, s_seq, t_seq, False, False])
					#if len(e_s_commonPath) <= acceptedLengthPath:	
					if e_s_commonPath != None:
						connectionGraph[s_seq].add_tail_head([len(e_s_commonSeq), len(e_s_commonPath), e_s_commonPath, t_seq,  s_last[2], t_first[2], connectionGraph[t_seq]])
					
				if len(e_e_commonSeq) >= links:
					printDebug("s_seq: %s, t_seq: %s" %(s_seq, t_seq))
					print "%s and %s share these sequences (end - end):" % (s_seq, t_seq)
					print e_e_commonSeq
					print "node one:" 
					print s_last[1].print_without_seq()
					print "node two:" 
					print t_last[1].print_without_seq()													
					e_e_commonPath, e_e_nodesUsed = getConnectionPath(sortedNodeListCollection, e_e_commonSeq, s_last[1], t_last[1], s_seq, 'tail', firstIsLast, asm_output)	
					#connectionGraph[s_seq].append(['ee', e_e_commonPath, s_seq, t_seq, False, True])	
					#if len(e_e_commonPath) <= acceptedLengthPath:
					if e_e_commonPath != None:
						connectionGraph[s_seq].add_tail_tail([len(e_e_commonSeq), len(e_e_commonPath), e_e_commonPath, t_seq,  s_last[2], t_last[2], connectionGraph[t_seq]])
				
		firstIsLast = False
	
	sortedConnectionGraph = [(len(p.path), seq, p) for seq, p in connectionGraph.items()]
	sortedConnectionGraph.sort(reverse=True)
					
	connectionPaths = {}
	connectedPaths = {}
	
	for path_length, seq, p in sortedConnectionGraph:
		print seq
		print p
		#seqs, connectedPath = p.getPath()
		seqs, pathSequence = p.getPathSequence(asm_output, '', asm_output_gap_fraction)
		if len(seqs) > 1:
			seq_name = '_'.join(seqs)
		elif len(seqs) == 1:
			seq_name = seqs[0]
		print seqs
		print "PathSequence"
		#print pathSequence
		#print connectedPath
		#connectionPaths[seqs] = connectedPath
		if seqs and pathSequence:
			connectedPaths[seq_name] = pathSequence
		
	#sortedConnectionGraph = [(len(p.path), seq, p) for seq, p in connectionGraph.items()]
		#sortedBySizeNodeList = [(len(nodeList), seq, nodeList) for seq, nodeList in sortedNodeListCollection.items()]
	
	#sortedConnectionGraph.sort(reverse=True)
						
	with open(prefix + ".fasta", "w") as fasta_file_out:
		fasta_file_out.write('')
		
	print "Number of sequences: %d " % len(connectedPaths)
	sortedPathsByLength = [(len(path), seq, path) for seq, path in connectedPaths.items()]

	sortedPathsByLength.sort(reverse=True)
	
	sequence = ''
	path_num = 0
	path_length = 0

	for length, seq, path in sortedPathsByLength:
		#print path
		printDebug("path_%d original sequence: %s actual_length: %d" % (path_num, seq, length))

		printSeq(path, 60, prefix, ">path_%d original sequence: %s actual_length: %d" % (path_num, seq, length))
					
		path_num = path_num + 1							
	
def leftExtend(start_node, asmToFollow, seqToFollow, asm_output, sortedNodeListCollection):
	nodesUsedAndConsidered = set()
	nodesUsed = set()
	#size of node is used to get an approximately distance of the path to find
	
	#First, find which sequences common to both start and end
	seqs = set(start_node.seqs)
	#sortedNodeListCollection[c.src], c.plus_strand_start, c.src
	
	
	printDebug("leftExtend, start node: ")
	printDebug(start_node.print_without_seq())
		
	path = []
	preferredPath = None
	preferredLength = 0
	preferredCost = 100000	
	cheapestPath = []
	cheapestLength = 0 
	cheapestCost = 1000000
	increasingIndex = True
	nodesWithoutFullSupport = 3
	lastIndexOfFullSupport = None
	#extendPath = True
	doBreak = False
	
	for seq in seqs:
		#Find the position
		printDebug("Seq: %s" % seq)
		#print start.get_component(seq).plus_strand_start
		#print end.get_component(seq).plus_strand_start
		start_pos = binary_search(sortedNodeListCollection[seq], start_node.get_component(seq).plus_strand_start, seq)
		seqToFollow_pos = binary_search(sortedNodeListCollection[seqToFollow], start_node.get_component(seqToFollow).plus_strand_start, seqToFollow)
		printDebug("Start pos: %d" % start_pos)
		printDebug("Length of sequence in nodes: %d" % len(sortedNodeListCollection[seq]))
					
		#if start_pos >= 0 and start_pos < len(sortedNodeListCollection[seq]) - 1:

		
		#Get direction
		if start_node.get_component(seqToFollow).strand == start_node.get_component(seq).strand:
			#Decreasing index goes to the left
			increasingIndex = False
			printDebug("Not increasing index")
		elif start_node.get_component(seqToFollow).strand != start_node.get_component(seq).strand:
			#Increasing index goes to the left
			increasingIndex = True
			printDebug("Increasing index")
			
		index = start_pos	
		lastIndexOfFullSupport = start_pos
		while True:

			#change to new node	
			if increasingIndex:
				index = index + 1
			else:
				index = index - 1
						
			printDebug("Index: %d" % index)
			if index >= 0 and index < len(sortedNodeListCollection[seq]) - 1:
				printDebug(sortedNodeListCollection[seq][index].print_without_seq())
				
			#break if first node in the sequence		
			if index <= 0:
				printDebug("Index <= 0")
				break
				
			#break if last node in the sequence	
			if index >= len(sortedNodeListCollection[seq]):
				printDebug("index >= len(sortedNodeListCollection[seq])")
				break 

			for comp in sortedNodeListCollection[seq][index].components:
				if asmToFollow in comp.src and len(sortedNodeListCollection[comp.src]) > 5 and comp.src != seqToFollow:
					#Not acceptable to include sequences that are from asmToFollow that are longer than 5 nodes and that are not seqToFollow
					printDebug("Tried to include sequences from asmToFollow longer than 5 nodes")
					doBreak = True
				if comp.src == seqToFollow:
					printDebug(comp.src)
					printDebug(seqToFollow)
					printDebug("length of seqToFollow: %d" % len(sortedNodeListCollection[seqToFollow]))
					pos = binary_search(sortedNodeListCollection[seqToFollow], sortedNodeListCollection[seq][index].get_component(seqToFollow).plus_strand_start, seqToFollow)	
					if pos > seqToFollow_pos:
						#We're moving the wrong way!
						doBreak = True
										
			if doBreak:
				doBreak = False
				break						
						
			#Check if node has full support (more than 2 of the starting sequences)	
			printDebug("length intersection %d" % len(seqs.intersection(sortedNodeListCollection[seq][index].seqs)))		
			if len(seqs.intersection(sortedNodeListCollection[seq][index].seqs)) >= 2:
			#if len(sortedNodeListCollection[seq][index].seqs) >= 2:
				lastIndexOfFullSupport = index
				nodesWithoutFullSupport = 3
			elif nodesWithoutFullSupport == 0:
				printDebug("Have had 3 nodes without full support, breaking")
				nodesWithoutFullSupport = 3
				break
			else:
				nodesWithoutFullSupport = nodesWithoutFullSupport - 1
			
		printDebug("lastIndexOfFullSupport: %d, difference start_pos-lastIndexOfFullSupport %d" % (lastIndexOfFullSupport, abs(start_pos-lastIndexOfFullSupport)))	
			
		#Let's see if this works
		if lastIndexOfFullSupport <	start_pos:
			#Same strand as seqToFollow
			path = copy.deepcopy(sortedNodeListCollection[seq][lastIndexOfFullSupport:start_pos])
			for node in path:
				if node.get_component(seq).strand == '-':
					node.rev_comp = True
					#printDebug("made this node rev comp, in leftExtend:")
					#printDebug(node.print_without_seq())
				elif node.get_component(seq).strand == '+':
					node.rev_comp = False
					#printDebug("made this node not rev comp, in leftExtend:")
					#printDebug(node.print_without_seq()) 			
		else:
			path = copy.deepcopy(list(reversed(sortedNodeListCollection[seq][start_pos+1:lastIndexOfFullSupport+1])))
			for node in path:
				if node.get_component(seq).strand == '+':
					node.rev_comp = True
					#printDebug("made this node rev comp, in leftExtend:")
					#printDebug(node.print_without_seq())
				elif node.get_component(seq).strand == '-':
					node.rev_comp = False
					#printDebug("made this node not rev comp, in leftExtend:")
					#printDebug(node.print_without_seq()) 			
		
		printDebug("For %s, path is %d nodes long" % (seq, len(path)))
			
		if path:
			cost = 0
			length = 0			
			for node in path:
				#print node.print_without_seq()
				nodesUsedAndConsidered.add(node)
				length = length + node.size
				cost = cost + node.get_lowest_amount_of_gaps()
			printDebug("leftExtend cost: %d, length: %d" % (cost, length))

			if asm_output: 
				if asm_output in seq:
					preferredCost = cost
					preferredPath = path
					preferredLength = length	
	
			#Want the path with the least amount of gaps, but only if it is within 0.1 - 10 times the size of 
			#potential gap. It seems, at least with PacBio-based assembly, that a lot of gaps are almost 
			#overfilled, constricting the accepted sizes might help.
			#if (cost < cheapestCost):
			#if (cost < cheapestCost) and (length > size_of_node*0.1) and (length < size_of_node*10) :
			if (cost < cheapestCost):	
				cheapestCost = cost
				cheapestPath = path
				cheapestLength = length

	if preferredPath != None:
		printDebug("preferred cost %d, length %d" % (preferredCost, preferredLength))
		for node in preferredPath:
			#print node.print_without_seq()	
			nodesUsed.add(node)	
	else:
		for node in cheapestPath:
			#print node.print_without_seq()	
			nodesUsed.add(node)			
		
	#return cheapestPath, nodesUsedAndConsidered
	if preferredPath != None:
		return preferredPath, nodesUsed
	else:
		return cheapestPath, nodesUsed			
	
def rightExtend(start_node, asmToFollow, seqToFollow, asm_output, sortedNodeListCollection):
	nodesUsedAndConsidered = set()
	nodesUsed = set()
	#size of node is used to get an approximately distance of the path to find
	
	#First, find which sequences common to both start and end
	seqs = set(start_node.seqs)
	#sortedNodeListCollection[c.src], c.plus_strand_start, c.src
	
	
	printDebug("rightExtend, start node: ")
	printDebug(start_node.print_without_seq())
		
	path = []
	preferredPath = None
	preferredLength = 0
	preferredCost = 100000	
	cheapestPath = []
	cheapestLength = 0 
	cheapestCost = 1000000
	increasingIndex = True
	nodesWithoutFullSupport = 3
	lastIndexOfFullSupport = None
	#extendPath = True
	doBreak = False
	
	for seq in seqs:
		#Find the position
		printDebug("Seq: %s" % seq)
		#print start.get_component(seq).plus_strand_start
		#print end.get_component(seq).plus_strand_start
		start_pos = binary_search(sortedNodeListCollection[seq], start_node.get_component(seq).plus_strand_start, seq)
		seqToFollow_pos = binary_search(sortedNodeListCollection[seqToFollow], start_node.get_component(seqToFollow).plus_strand_start, seqToFollow)

		printDebug("Start pos: %d" % start_pos)
		printDebug("Length of sequence in nodes: %d" % len(sortedNodeListCollection[seq]))
		
		#Get direction
		if start_node.get_component(seqToFollow).strand == start_node.get_component(seq).strand:
			#Increasing index goes to the right
			increasingIndex = True
			printDebug("Not increasing index")
		elif start_node.get_component(seqToFollow).strand != start_node.get_component(seq).strand:
			#Decreasing index goes to the right
			increasingIndex = False
			printDebug("Increasing index")
			
		index = start_pos	
		lastIndexOfFullSupport = start_pos
		while True:

			#change to new node	
			if increasingIndex:
				index = index + 1
			else:
				index = index - 1
						
			printDebug("Index: %d" % index)
			if index >= 0 and index < len(sortedNodeListCollection[seq]) - 1:
				printDebug(sortedNodeListCollection[seq][index].print_without_seq())
				
			#break if first node in the sequence		
			if index <= 0:
				printDebug("Index <= 0")
				break
				
			#break if last node in the sequence	
			if index >= len(sortedNodeListCollection[seq]):
				printDebug("index >= len(sortedNodeListCollection[seq]) ")
				break 

			for comp in sortedNodeListCollection[seq][index].components:
				if asmToFollow in comp.src and len(sortedNodeListCollection[comp.src]) > 5 and comp.src != seqToFollow:
					#Not acceptable to include sequences that are from asmToFollow that are longer than 5 nodes and that are not seqToFollow
					printDebug("Tried to include sequences from asmToFollow longer than 5 nodes")
					doBreak = True
				if comp.src == seqToFollow:
					printDebug(comp.src)
					printDebug(seqToFollow)
					printDebug("length of seqToFollow: %d" % len(sortedNodeListCollection[seqToFollow]))				
					pos = binary_search(sortedNodeListCollection[seqToFollow], sortedNodeListCollection[seq][index].get_component(seqToFollow).plus_strand_start, seqToFollow)	
					if pos < seqToFollow_pos:
						#We're moving the wrong way!
						doBreak = True					
					
			if doBreak:
				doBreak = False
				break						
						
			#Check if node has full support (more than 2 of the starting sequences)	
			printDebug("length intersection %d" % len(seqs.intersection(sortedNodeListCollection[seq][index].seqs)))		
			if len(seqs.intersection(sortedNodeListCollection[seq][index].seqs)) >= 2:
			#if len(sortedNodeListCollection[seq][index].seqs) >= 2:
				lastIndexOfFullSupport = index
				nodesWithoutFullSupport = 3
			elif nodesWithoutFullSupport == 0:
				printDebug("Have had 3 nodes without full support, breaking")
				nodesWithoutFullSupport = 3
				break
			else:
				nodesWithoutFullSupport = nodesWithoutFullSupport - 1
			
		printDebug("lastIndexOfFullSupport: %d, difference start_pos-lastIndexOfFullSupport: %d" % (lastIndexOfFullSupport, abs(start_pos-lastIndexOfFullSupport)))	
			
		#Let's see if this works
		if lastIndexOfFullSupport >	start_pos:
			#Same strand as seqToFollow
			path = copy.deepcopy(sortedNodeListCollection[seq][start_pos+1:lastIndexOfFullSupport+1])
			for node in path:
				if node.get_component(seq).strand == '-':
					node.rev_comp = True
					#printDebug("made this node rev comp, in rightExtend:")
					#printDebug(node.print_without_seq())
				elif node.get_component(seq).strand == '+':
					node.rev_comp = False
					#printDebug("made this node not rev comp, in rightExtend:")
					#printDebug(node.print_without_seq()) 			
		else:
			path = copy.deepcopy(list(reversed(sortedNodeListCollection[seq][lastIndexOfFullSupport:start_pos])))
			for node in path:
				if node.get_component(seq).strand == '+':
					node.rev_comp = True
					#printDebug("made this node rev comp, in rightExtend:")
					#printDebug(node.print_without_seq())
				elif node.get_component(seq).strand == '-':
					node.rev_comp = False
					#printDebug("made this node not rev comp, in rightExtend:")
					#printDebug(node.print_without_seq()) 			
		
		printDebug("For %s, path is %d nodes long" % (seq, len(path)))
			
		if path:
			cost = 0
			length = 0			
			for node in path:
				#print node.print_without_seq()
				nodesUsedAndConsidered.add(node)
				length = length + node.size
				cost = cost + node.get_lowest_amount_of_gaps()
			printDebug("rightExtend cost: %d, length: %d" % (cost, length))

			if asm_output: 
				if asm_output in seq:
					preferredCost = cost
					preferredPath = path
					preferredLength = length	
	
			#Want the path with the least amount of gaps, but only if it is within 0.1 - 10 times the size of 
			#potential gap. It seems, at least with PacBio-based assembly, that a lot of gaps are almost 
			#overfilled, constricting the accepted sizes might help.
			#if (cost < cheapestCost):
			#if (cost < cheapestCost) and (length > size_of_node*0.1) and (length < size_of_node*10) :
			if (cost < cheapestCost):	
				cheapestCost = cost
				cheapestPath = path
				cheapestLength = length

	if preferredPath != None:
		printDebug("preferred cost %d, length %d" % (preferredCost, preferredLength))
		for node in preferredPath:
			#print node.print_without_seq()	
			nodesUsed.add(node)	
	else:
		for node in cheapestPath:
			#print node.print_without_seq()	
			nodesUsed.add(node)			
		
	#return cheapestPath, nodesUsedAndConsidered
	if preferredPath != None:
		return preferredPath, nodesUsed
	else:
		return cheapestPath, nodesUsed		
	
def closeGaps(sortedNodeListCollection, asm_path, cov_fraction, numAsms, asm_output, extend):
	#nodesUsedAndConsidered = set()
	nodesUsed = set()
	nodesConsidered = set()
	usedNodes = set()
	paths = {}
	leftPath = []
	rightPath = []
	last = [0, 0]
	first = [0, 0]

	#Sort sortedNodeListCollection after longest sequence
	sortedBySizeNodeList = [(len(nodeList), seq, nodeList) for seq, nodeList in sortedNodeListCollection.items()]
	
	sortedBySizeNodeList.sort(reverse=True)
	
	for nodeListLen, seq, nodeList in sortedBySizeNodeList:	
		#Follow the sequences from one assembly				
		if asm_path in seq:
			printDebug("Trying to close gaps in sequence %s" % seq)
			
			inNodeListAndUsed = set(nodeList).intersection(nodesUsed)
			#inNodeListAndConsidered = set(nodeList).intersection(nodesUsedAndConsidered)
			#if inNodeListAndConsidered:
			if inNodeListAndUsed:
				nodesUsedSize = 0
				nodesSize = 0
				for n in inNodeListAndUsed:
					nodesUsedSize = nodesUsedSize + n.size
				for n in nodeList:	
					nodesSize = nodesSize + n.size
				printDebug("Some or all of nodes from sequence %s, and of length %d and size %d has been used. These number %d and are of size %d. " % (seq, len(nodeList), nodesSize, len(inNodeListAndUsed), nodesUsedSize))
				#for node in inNodeListAndUsed:
					#print node.print_without_seq()
				#We would expect some dual use of nodes because of repeats for instance. But, if more than
				#20 % of a path/nodeList is used, then skip it. Only small scaffolds/contigs should be affected.
				#if  len(inNodeListAndUsed)/float(len(nodeList)) > 0.25:
				if  float(nodesUsedSize)/nodesSize >= cov_fraction:
					printDebug("Skipping %s ..." % seq)
					continue				
			
			path = []
			#The first element must be OK.
			#prev_start = nodeList[0].get_component(seq).plus_strand_start
			#prev_size = nodeList[0].get_component(seq).size
			
			#First, see if we can extend the path to the left of the sequence
			if extend:
				printDebug("Trying to extend path to the left")
				start_node = None
				start_index = -1
				#We want to try to extend from the first node with the max number of assemblies
				#involved. Let's get it:
				for index, node in enumerate(nodeList):
					#Get the first occurrence of the 3 alignments
					#(if a unplaced sequence has filled a gap, this might not be the case)
					#if len(node.components) >= first[0] and node.get_component(seq):
					if len(node.components) >= 3 and node.get_component(seq):
						first = [len(node.components), index]
						printDebug("len(node.components): %d, first index: %d" % (len(node.components), index))
						break
			
				start_node = nodeList[first[1]]			
			
				leftPath, usedNodes = leftExtend(start_node, asm_path, seq, asm_output, sortedNodeListCollection)	
				printDebug("Length of leftPath: %d" % len(leftPath))
				printDebug("Nodes extending to left")
			
				for node in leftPath:
					printDebug(node.print_without_seq())
			
			for index, node in enumerate(nodeList):	
				
				#print "at node: "
				#print node.print_without_seq()
				#If the alignment block only has one component, it is likely a gap
				if (len(node.components) == 1) and (index != 0) and (index != len(nodeList)-1):

					#a_path, usedAndConsidered = getPath(sortedNodeListCollection, nodeList[index-1], nodeList[index+1], seq, node.size, nodesUsed)
					a_path, used, considered = getPath(sortedNodeListCollection, nodeList[index-1], nodeList[index+1], seq, node.size, nodesUsed, asm_output)
					#nodesUsedAndConsidered.update(usedAndConsidered)
					
					if a_path and used:
						printDebug("found a path!")
						#for n in s_path:
						#	print n.print_without_seq()	
						path = path + a_path
						nodesUsed.update(used)
						nodesConsidered.update(considered)
					elif used:
						printDebug("No path, but some nodes used")
						#a_path can be empty, if shortest path is without the gap node
						#if so, just add the nodes to used and don't do anything with the path
						nodesUsed.update(used)
						nodesConsidered.update(considered)
					elif extend:
						#no paths were gotten. Let's try to extend into the gap. Experimental!
						printDebug("Trying to extend into gap!")							
						gapRightPath, gapRUsedNodes = rightExtend(nodeList[index-1], asm_path, seq, asm_output, sortedNodeListCollection)	
						gapLeftPath, gapLUsedNodes = leftExtend(nodeList[index+1], asm_path, seq, asm_output, sortedNodeListCollection)
					else:	
						printDebug("Just adding the node")
						path = path + [node]
						nodesUsed.add(node)
				else:
					path = path + [node]
					nodesUsed.add(node)
			"""
			index = 0
			while index < len(nodeList):
				node = nodeList[index]

				#If the alignment block only has one component, it is likely a gap
				if (len(node.components) < numAsms) and (index != 0) and (index != len(nodeList)-1):
				
					firstNodeIndex = index - 1 
					printDebug("firstNodeIndex: %d" % firstNodeIndex)
					#index - 1 should have the same number of alignments as numAsms. We need to find the index of the next node 
					#with alignments equal to numAsms
					lastNodeIndex = index + 1
					while len(nodeList[lastNodeIndex].components) != numAsms and lastNodeIndex < len(nodeList) - 1:
						lastNodeIndex = lastNodeIndex + 1
					
					printDebug("lastNodeIndex: %d" % lastNodeIndex)
					
					
					a_path, used = getPath(sortedNodeListCollection, nodeList[firstNodeIndex], nodeList[lastNodeIndex], seq, node.size, nodesUsed)
					
					if a_path and used:
						#for n in s_path:
						#	print n.print_without_seq()	
						path = path + a_path
						nodesUsed.update(used)
						index = lastNodeIndex
					elif used:
						#a_path can be empty, if shortest path is without the gap node
						#if so, just add the nodes to used and don't do anything with the path
						nodesUsed.update(used)
						index = lastNodeIndex
					else:	
						path = path + [node]
						nodesUsed.add(node)
						index = index + 1
				else:
					path = path + [node]
					nodesUsed.add(node)
					index = index + 1
			"""		
			if extend:
				#Then we want to see if we can extend to the right of the scaffold
								
				end_node = None
				end_index = -1
				if len(nodeList) != 1:
					#http://christophe-simonis-at-tiny.blogspot.no/2008/08/python-reverse-enumerate.html
					reverse_enumerate = lambda l: izip(range(len(l)-1, -1, -1), reversed(l))
					for index, node in reverse_enumerate(nodeList):
						if len(node.components) >= 3 and node.get_component(seq):
							last = [len(node.components), index]
							printDebug("len(node.components): %d, last index: %d" % (len(node.components), index))
							break
			
					printDebug("found last at index %d" % last[1])	
				else:
					last = [len(nodeList[0].components), 0]
			
				start_node = nodeList[last[1]]			
			
				rightPath, usedNodes = rightExtend(start_node, asm_path, seq, asm_output, sortedNodeListCollection)	
				printDebug("Length of rightPath: %d" % len(rightPath))
				printDebug("Nodes extending to right")
			
				for node in rightPath:
					printDebug(node.print_without_seq())						
			
				printDebug("len(leftPath): %d, len(path[:first[1]]): %d" % (len(leftPath), len(path[:first[1]])))
				printDebug("len(rightPath): %d, len(path[last[1]:]): %d" % (len(rightPath), len(path[last[1]:])))
			
			if extend:
				if len(leftPath) > len(path[:first[1]]) and len(rightPath) > len(path[last[1]:]):
					paths[seq] = leftPath + path + rightPath
				elif len(leftPath) > len(path[:first[1]]):
					paths[seq] = leftPath + path
				elif len(rightPath) > len(path[last[1]:]):
					paths[seq] = path + rightPath
				else:	
					paths[seq] = path
			else:
				paths[seq] = path
				
	return paths
	
def getStats(sortedNodeListCollection):
	"""
	Prints some statistics about the alignment going into the process, after
	filtering. Return the number of assemblies being compared
	"""
	stats = defaultdict(list)
	for seq, nodeList in sortedNodeListCollection.items():
	
		asm = seq.split('.')[0]
		num_nodes = len(nodeList)
		alignment_length = 0
		sequence_length = 0
		num_asms_in_nodes = defaultdict(int)
		
		#Get the alignment length, sequence length and how many nodes with how many sequences
		#for each list of nodes.
		for node in nodeList:
			alignment_length = alignment_length + node.size
			sequence_length = sequence_length + len(node.get_component(seq).gap_less)
			num_asms_in_nodes[len(node.components)] = num_asms_in_nodes[len(node.components)] + 1
		
		#Update the assembly stats with the addition of these stats	
		if stats[asm]:	
			stats[asm][0] = stats[asm][0] + num_nodes
			stats[asm][1] = stats[asm][1] + alignment_length
			stats[asm][2] = stats[asm][2] + sequence_length
			for l in num_asms_in_nodes:
				if l in stats[asm][3]:
					stats[asm][3][l] = stats[asm][3][l] + num_asms_in_nodes[l]
				else:
					stats[asm][3][l] = num_asms_in_nodes[l]
		
		#Or create the stats			
		else:
			stats[asm] = [num_nodes, alignment_length, sequence_length, dict()]
			for l in num_asms_in_nodes:
				if l in stats[asm][3]:
					stats[asm][3][l] = stats[asm][3][l] + num_asms_in_nodes[l]
				else:
					stats[asm][3][l] = num_asms_in_nodes[l]			
			#for l in num_asms_in_nodes:
			#	stats[asm][3][l] = num_asms_in_nodes[l]	
		#print num_asms_in_nodes		
	
	for asm in stats:
		print "Assembly %s consists of %d nodes, with alignment length of %d and sequence length of %d" % (asm, stats[asm][0], stats[asm][1], stats[asm][2])
		for l in stats[asm][3]:
			print "%d nodes have alignment of %d assemblies" % (stats[asm][3][l], l)
	
	#Return the number of assemblies		
	return len(stats)		
							
	
#def parse_maf(maf_file, asm_path, asm_output, connect, close_gaps, arbi, cov_fraction):
def parse_maf(maf_file, asm_path, asm_output, connect, close_gaps, cov_fraction, prefix, extend, asm_output_gap_fraction):
	
	AllNodes = []
	sortedNodeListCollection = defaultdict(list)
	
	#Returns an Alignment class, described here: https://bitbucket.org/james_taylor/bx-python/src/babadb4d4bf2d71e50a6f3569c10691ec9f3bc81/lib/bx/align/core.py?at=default
	#This function is found here: https://bitbucket.org/james_taylor/bx-python/src/babadb4d4bf2d71e50a6f3569c10691ec9f3bc81/lib/bx/align/maf.py?at=default
	print "Reading the MAF file..."
	for m in maf.Reader(maf_file):
		a = AlignmentBlock(m)
		AllNodes.append(a)
	
	print "Creating a dictionary of sequences..."
	#Create a dictionary of sequences, with a list of all the nodes with that sequence as value
	for node in AllNodes:
		for c in node.components:
			sortedNodeListCollection[c.src].append(node)
	"""
	with open("sortedNodeListCollection", "w") as output_sorted_lists:
			output_sorted_lists.write("")
	"""
	print "Creating a sorted list for the sequences..."
	#For all the sequences, create a sorted list of based on position in sequence
	for seq, nodeList in sortedNodeListCollection.items():
		sortedNodeListCollection[seq] = getSortedVersionBasedOnPositionInSuppliedSeq(nodeList, seq)
		#print "%s is of length %d" % (seq, len(sortedNodeListCollection[seq]))
		#sortedNodeListCollection[seq] = getCleanSortedVersionBasedOnPositionInSuppliedSeq(nodeList, seq)
		
	print "Filtering the sorted lists for overlapping alignments..."
	for seq, nodeList in sortedNodeListCollection.items():
		sortedNodeListCollection[seq] = filterNodeList(nodeList, seq)
		
		#print "%s is of length %d" % (seq, len(sortedNodeListCollection[seq]))
	
	#I want number of nucleotides/gaps, alignment length, number of nodes, number of nodes of how many assemblies	
	print "Getting some statistics of the alignments..."
	numAsms = getStats(sortedNodeListCollection)
	
	if debug:
	#if False:
		with open(prefix + ".sortedNodeListCollection", "w") as output_sorted_lists:
				output_sorted_lists.write("")		
				
		for seq, nodeList in sortedNodeListCollection.items():			
			with open(prefix + ".sortedNodeListCollection", "a") as output_sorted_lists:
				output_sorted_lists.write("Sequence: %s, num nodes: %d\n" % (seq, len(sortedNodeListCollection[seq])))
			
				for i, n in enumerate(sortedNodeListCollection[seq]):
					output_sorted_lists.write("Node number: %d\n" % i)
					output_sorted_lists.write(n.print_without_seq())
	
	
	print "Adding neighbors to all nodes..."
	sys.stdout.flush()
	for node in AllNodes:
		for c in node.components:
   			newNeighbors = findProximalNodesByBinarySearch(sortedNodeListCollection[c.src], c.plus_strand_start, c.src)
   			#print newNeighbors
   			node.addToNeighborsSet(newNeighbors)


	#print sortedNodeListCollection
	"""
	with open("neighbors", "w") as neighbors:
		neighbors.write("")
	
	for seq in sortedNodeListCollection:
		with open("neighbors", "a") as neighbors:
			neighbors.write("\nSelf: ")
			neighbors.write(seq)
			
			for node in sortedNodeListCollection[seq]:
				neighbors.write("\nNode:")
				neighbors.write(node.print_without_seq())
				
				neighbors.write("\nForward Neighbors:\n")
				for neighbor in node.forward_neighbors:
					neighbors.write(neighbor.print_without_seq())
					neighbors.write("\n")

				neighbors.write("\nBackward Neighbors:\n")
				for neighbor in node.backward_neighbors:
					neighbors.write(neighbor.print_without_seq())
					neighbors.write("\n")
	
	"""
	paths = {}	
	
	if asm_path:
		#Connect implies close gaps
		if connect:
			print "Trying to connect scaffolds..."
			#connect has to be at least one less than number of assemblies
			if (numAsms-1) >= connect:
				connectedScaffolds = connectScaffolds(sortedNodeListCollection, asm_path, cov_fraction, prefix, connect, numAsms, asm_output, extend, asm_output_gap_fraction)
			else:
				print "The number of connections required, %d, is equal or larged than the number of assemblies, %d, and can't work." % (connect, numAsms) 
																						
		elif close_gaps:
			print "Trying to close some gaps..."
			sys.stdout.flush()
			paths = closeGaps(sortedNodeListCollection, asm_path, cov_fraction, numAsms, asm_output, extend)
			
		else:
			for seq, nodeList in sortedNodeListCollection.items():
				if asm_path in seq:
					paths[seq] = nodeList
	
	else:
		print "Unsupported mode (for now at least)"
		#for node in AllNodes:
			#this does not work, probably
		#	path = node.getImprovedPathIter()
		#	paths.append( path )

  	if not connect:
		#Print the sequences	
		path_length = 0
		path_num = 0
		sequence = ''
	
		with open(prefix + ".fasta", "w") as fasta_file_out:
			fasta_file_out.write("")
	
		#print paths
		print "Number of sequences: %d " % len(paths)
	
		#Sorting by number of nodes in each path, will output the one with the most nodes first
		sortedPathsByNumberNodes = [(len(path), seq, path) for seq, path in paths.items()]
	
		sortedPathsByNumberNodes.sort(reverse=True)
	
		for num_nodes, seq, path in sortedPathsByNumberNodes:		
			for node in path:
				print node.print_without_seq()
				print node.size
				path_length = path_length + node.size
				if asm_output:
					node_seq, node_src, amount = node.get_sequence_from(asm_output, seq, asm_output_gap_fraction)
				else:
					node_seq, node_src, amount = node.get_sequence_with_fewest_gaps(seq)
				printDebug("Seq with fewest gaps: %s, amount %d " % (node_src, amount))
				if len(sequence) > 0:
					printDebug("Fraction length of node/length of sequence: %.2g" % (float(node.size)/len(sequence)))
				printDebug("Rev_comp on node %s" % node.rev_comp)
				sequence = sequence + node_seq
			#print "Length of alignment path: %d" % path_length

			printDebug("path_%d num_elements: %d path_length: %d actual_length: %d" % (path_num, len(path), path_length, len(sequence)))

			printSeq(sequence, 60, prefix, ">path_%d num_elements: %d path_length: %d actual_length: %d" 
							% (path_num, len(path), path_length, len(sequence)))
						
			path_num = path_num + 1
			path_length = 0
			sequence = ''	
				

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description=
	'')
	parser.add_argument('-m', '--maf', action='store', help='MAF file', required=True, type=argparse.FileType('r'))	
	parser.add_argument('-p', '--path', action='store', help='Path along this assembly', required=False, default=None)
	parser.add_argument('-r', '--prefix', action='store', help='Prefix of output files', required=False, default="prefix")	
	parser.add_argument('-o', '--output', action='store', help='Preferably output sequences from this assembly. If it does not exist, it is replaced with the one with least gaps', required=False, default=None)
	parser.add_argument('-c', '--connect', action='store', help='If following a path along one assembly, try to connect the scaffolds in that assembly. Specify how many links (information from other assemblies) are needed to merge scaffolds.', default=0, type=int)
	parser.add_argument('-g', '--close_gaps', action='store_true', help='If following a path along one assembly, try to close gaps by following alternative routes', default=False)
	parser.add_argument('-d', '--debug', action='store_true', help='Print lots of debug information', default=False)
	#parser.add_argument('-a', '--arbi', action='store_true', help='Arbitary path', default=False)
	parser.add_argument('-f', '--fraction', action='store', help='If a sequence is covered to this degree by another, it is thrown out', default=0.75, type=float)
	parser.add_argument('-n', '--output_gap_fraction', action='store', help='If the output sequence in an alignment has more than this fraction of gaps compared to the one with the least gaps, chose the least gappy', default=0.0, type=float)
	parser.add_argument('-e', '--extend', action='store_true', help='Try to extend scaffolds if supported by other sequence. Also try to extend into gaps if gapclosing fails', default=False)


	args = parser.parse_args()
	
	debug = args.debug
		
	#parse_maf(args.maf, args.path, args.output, args.connect, args.close_gaps, args.arbi, args.fraction)
	parse_maf(args.maf, args.path, args.output, args.connect, args.close_gaps, args.fraction, args.prefix, args.extend, args.output_gap_fraction)
	
