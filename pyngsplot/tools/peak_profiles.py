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
import matplotlib 
matplotlib.use('Agg')
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
	
def read_bam(bam, positions, halfwinwidth, return_dict, norm, size_dict):
	bamfile = HTSeq.BAM_Reader(bam)
	profile = numpy.zeros( 2*halfwinwidth, dtype="f" )
	if norm:
		gapdh = read_counts(norm[bam])
		constant = 1000/float(gapdh)
	elif size_dict:
		constant = 1000000/float(size_dict[bam])
	else:
		constant = 1
	c = 0
	for p in positions:
		mid = (p.end - p.start) /2
		mid = mid+p.start
		window = HTSeq.GenomicInterval( p.chrom, mid - halfwinwidth, mid + halfwinwidth, "." )
		for almnt in bamfile[ window ]:
			start_in_window = almnt.iv.start - mid + halfwinwidth
			end_in_window   = almnt.iv.end   - mid +  halfwinwidth
		#	print almnt.iv.start, almnt.iv.end, p.start, p.end, start_in_window, end_in_window
			start_in_window = max( start_in_window, 0 ) #A more elegant solution for bed files under 0
			end_in_window = min( end_in_window, 2*halfwinwidth )
			if start_in_window >= 2*halfwinwidth or end_in_window < 0:
				continue
			profile[ start_in_window : end_in_window ] += constant
		c += 1
	profile = profile/float(c) #Average over the number of TSS regions
	return_dict[bam] = profile

def read_bam_function(args):
	return read_bam(*args)

#Should make this more simple, do I use size/gapdh/none as normalisation?
def parallel_peak_file_plot(ibams, bed_file, size_dict, halfwinwidth, norm, controls):
	positions = set()
	for line in open(bed_file):
		fields = line.split( "\t" )
		name = re.sub("chr", "", fields[0])
		window = HTSeq.GenomicInterval( name, int(fields[1]), int(fields[2]), "." )
		positions.add(window)
	if controls == None:
		manager = Manager()
		return_dict = manager.dict()
		pool = Pool(8)
		if norm: #Normalisation provided
			pool.map(read_bam_function, itertools.izip(list(ibams.keys()), itertools.repeat(positions), itertools.repeat(halfwinwidth), 
				itertools.repeat(return_dict), itertools.repeat(norm), itertools.repeat(None))) ##Running annotation in parallel
		else:
			pool.map(read_bam_function, itertools.izip(list(ibams.keys()), itertools.repeat(positions), itertools.repeat(halfwinwidth), 
				itertools.repeat(return_dict), itertools.repeat(None), itertools.repeat(size_dict)))
		pool.close()
		pool.join()		
		for key in return_dict.keys():
			pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), return_dict[key], label=ibams[key])  
		pyplot.legend(prop={'size':8})
		pyplot.savefig("Average_peak_profile.pdf".format(ibams[key]))
	else:
		manager = Manager()
		return_dict = manager.dict()
		pool = Pool(8)
		pool.map(read_bam_function, itertools.izip(list(ibams.keys()), itertools.repeat(positions), itertools.repeat(halfwinwidth), 
				itertools.repeat(return_dict), itertools.repeat(None), itertools.repeat(None)))
		control_dict = manager.dict()
		control_bam = []
		for key in controls:
			control_bam.append(controls[key])
		control_sizes = sam_size(control_bam)
		pool = Pool(8)
		pool.map(read_bam_function, itertools.izip((control_bam), itertools.repeat(positions), itertools.repeat(halfwinwidth), 
			itertools.repeat(control_dict),  itertools.repeat(None),  itertools.repeat(None))) 
		pool.close()
		pool.join()	
	#	colors = ["b", "g", "r", "y", "k"] #To make it more robust, just use default colors
		c = 0
		for key in return_dict.keys():
			control = controls[key]
			new_profile = return_dict[key] - control_dict[control] #Unsure if working properly, maybe make it more intelligent?
			gapdh = read_counts(norm[key])
			constant = 1000/float(gapdh)
			new_profile = new_profile*constant
			pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), new_profile, label=ibams[key])#, color=colors[c])
		pyplot.legend(prop={'size':8})
		pyplot.savefig("Average_peak_profile.pdf".format(ibams[key]))
