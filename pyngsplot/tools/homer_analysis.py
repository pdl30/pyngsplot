#!/usr/bin/python

########################################################################
# 2 September 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, re, sys
import argparse
import numpy as np
import HTSeq
import subprocess

#This just works with peaks
def peak_homer_analysis(peaks, bam_file=None, tag=None):
	#Bam file has to be sorted and indexed!
	if bam_file:
		name = re.sub(".bam", "", bam_file)
		name = re.sub("_sort", "", name)
		command = ["makeTagDirectory", "{}_homer".format(name), bam_file]
		subprocess.call(command)
		command2 = "annotatePeaks.pl {} mm10 -size 4000 -hist 10 -ghist -d {}_homer".format(peaks, name)
		out = name + "_homer_output.txt"
		output = open(out, "w")
		subprocess.call(command2.split(), stdout=output)
		output.close()
	elif tag:
		name = re.sub("_homer", "", tag)
		command2 = "annotatePeaks.pl {} mm10 -size 4000 -hist 10 -ghist -d {}_homer".format(peaks, name)
		out = name + "_homer_output.txt"
		print command2
		output = open(out, "w")
		subprocess.call(command2.split(), stdout=output)
		output.close()
	return out

#TSS heatmap, create a list of lists and then convert it to numpy array, problem is you need values for individual positions which are tricky!
#For now just use current approach
def tss_homer_analysis(gtf, gene_filter, ucsc, bam_file=None, tag=None):
	#Bam file has to be sorted and indexed!
	peaks = "tmp.bed"
	output = open("tmp.bed", "w")
	tsspos = set()
	gtffile = HTSeq.GFF_Reader( gtf )
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
	for p in tsspos: #Problem if ensembl vs UCSC
		if ucsc == True:
			if p.chrom == "MT":
				chrom = "chrM"
			else:
				chrom = "chr" + p.chrom
		else:
			chrom = p.chrom
		start = p.pos - 200 
		end = p.pos + 200
		if start < 0:
			start = 0
		output.write("{}\t{}\t{}\n".format(chrom, start, end)),
	output.close()
	if bam_file:
		name = re.sub(".bam", "", bam_file)
		name = re.sub("_sort", "", name)
		command = ["makeTagDirectory", "{}_homer".format(name), bam_file]
		subprocess.call(command)
		command2 = "annotatePeaks.pl {} mm10 -size 10000 -hist 10 -ghist -d {}_homer".format(peaks, name)
		out = name + "_homer_output.txt"
		output = open(out, "w")
		subprocess.call(command2.split(), stdout=output)
		output.close()
	elif tag:
		name = re.sub("_homer", "", tag)
		command2 = "annotatePeaks.pl {} mm10 -size 10000 -hist 10 -ghist -d {}_homer".format(peaks, name)
		out = name + "_homer_output.txt"
		output = open(out, "w")
		subprocess.call(command2.split(), stdout=output)
		output.close()
	return out

