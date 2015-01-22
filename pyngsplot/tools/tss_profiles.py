#!/usr/bin/python

########################################################################
# 28 July 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, re, sys
import argparse
import HTSeq
import numpy
from matplotlib import pyplot
import pysam
import ConfigParser
from multiprocessing import Pool, Manager
import itertools

def read_counts(count):
	with open(count) as f:
		header= next(f)				
		header = header.rstrip()
		word = header.split("\t")
		gapdh = int(word[1])
	return gapdh

def sam_size(ibam):
	results=  {}
	for bam in ibam:
		size = reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(bam) ])
		results[bam] = size
	return results

def read_bam_tss(bam, norm, positions, size_dict, halfwinwidth, ucsc, return_dict):
	if norm:
		gapdh = read_counts(norm[bam])
		constant = 1000/float(gapdh)
	else:
		constant = 1000000/float(size_dict[bam])
	bamfile = HTSeq.BAM_Reader(bam)
	profile = numpy.zeros( 2*halfwinwidth, dtype="f" )
	c = 0 #Consider adding diviser to make them compariable between plots
	fragmentsize = 150 #Ajust later
	for p in positions: #Problem if ensembl vs UCSC
		if ucsc == True:
			if p.chrom == "MT":
				chrom = "chrM"
			else:
				chrom = "chr" + p.chrom
		else:
			chrom = p.chrom
		start = p.pos - halfwinwidth - fragmentsize 
		end = p.pos + halfwinwidth + fragmentsize
		if start < 0:
			start = 0
		window = HTSeq.GenomicInterval( chrom, start, end, "." )
		for almnt in bamfile[ window ]:	
			almnt.iv.length = fragmentsize
			if p.strand == "+":
				start_in_window = almnt.iv.start - p.pos + halfwinwidth 
				end_in_window   = almnt.iv.end - p.pos + halfwinwidth 
			else:
				start_in_window = p.pos + halfwinwidth - almnt.iv.end
				end_in_window   = p.pos + halfwinwidth - almnt.iv.start
			start_in_window = max( start_in_window, 0 )
			end_in_window = min( end_in_window, 2*halfwinwidth )
			if start_in_window >= 2*halfwinwidth or end_in_window < 0:
				continue
			profile[ start_in_window : end_in_window ] += constant
		c += 1
	profile = profile/float(c) #Average over the number of TSS regions
	return_dict[bam] = profile

def read_bam_tss_function(args):
	return read_bam_tss(*args)

def plot_tss_profile(ibams, gff, total_mapped, flank_size, label, gene_filter, ucsc=False, counts=None):
	gtffile = HTSeq.GFF_Reader( gff )
	halfwinwidth = int(flank_size)
	tsspos = set()
	#Add a filter for genes which are in the input file
	if gene_filter:
		g_filter = {}
		with open(gene_filter) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				g_filter[word[0]] = 1
		for feature in gtffile:
		#	feature.name is gene name, useful for delving into GFF files again!
			if feature.name in g_filter:
				if feature.type == "exon" and feature.attr["exon_number"] == "1":
					tsspos.add( feature.iv.start_d_as_pos )
	else:
		#Just use every gene from the GTF file
		for feature in gtffile:
			if feature.type == "exon" and feature.attr["exon_number"] == "1":
				tsspos.add( feature.iv.start_d_as_pos )
	#Manager dictionary to take the results from the pooled results
	manager = Manager()
	return_dict = manager.dict()
	pool = Pool(8)
	if ucsc:
		pool.map(read_bam_tss_function, itertools.izip(list(ibams.keys()), itertools.repeat(counts), itertools.repeat(tsspos), itertools.repeat(total_mapped), itertools.repeat(halfwinwidth), 
			itertools.repeat(True), itertools.repeat(return_dict))) ##Running annotation in parallel
	else:
		pool.map(read_bam_tss_function, itertools.izip(list(ibams.keys()), itertools.repeat(counts), itertools.repeat(tsspos), itertools.repeat(total_mapped), itertools.repeat(halfwinwidth), 
			itertools.repeat(False), itertools.repeat(return_dict))) ##Running annotation in parallel
	pool.close()
	pool.join()	
	#Simple line plot using the boundaries given and legend
	for key in return_dict.keys():
		pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), return_dict[key], label=ibams[key])  
	pyplot.legend(prop={'size':8})
	pyplot.savefig(label+".pdf")