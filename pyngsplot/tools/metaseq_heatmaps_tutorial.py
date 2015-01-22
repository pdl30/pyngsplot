#!/usr/bin/python 

########################################################################
# 19 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, re, os
import ConfigParser
import itertools
import argparse
from collections import defaultdict
import gffutils
import pybedtools
import metaseq
import multiprocessing
import numpy as np
from matplotlib import pyplot as plt
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval

def tss_generator(gtf):
	"""
	Generator function to yield TSS +/- 1kb of each annotated transcript
	"""
	for transcript in db.features_of_type('transcript'):
		yield TSS(asinterval(transcript), upstream=1000, downstream=1000)

def change_bed_size(bed, obed):
	output = open(obed, "w")
	with open(bed) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			chrom = word[0].strip("chr")
			start = int(word[1])-1000
			end = int(word[2])+1000
			if start < 0:
				start = 0
			output.write("{}\t{}\t{}\n".format(chrom, start, end)),

def metaseq_heatmap(ip_filename, bed, outname):
	change_bed_size(bed, "tmp.bed")
	peak = pybedtools.BedTool("tmp.bed")
#	db = gffutils.create_db(gtf, dbfn='test.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
	processes = multiprocessing.cpu_count()
	ip_signal = genomic_signal(ip_filename, 'bam')
	# Create arrays in parallel
	ip_array = ip_signal.array(peaks, bins=100, processes=processes)
	# Normalize to library size
	ip_array /= ip_signal.mapped_read_count() / 1e6
	# Sanity-checks
	x = np.linspace(-2000, 2000, 100)
	plt.rcParams['font.family'] = 'Arial'
	plt.rcParams['font.size'] = 10

	fig = metaseq.plotutils.imshow(ip_array, x=x, figsize=(3, 7),
	vmin=5, vmax=99,  percentile=True,
	line_kwargs=dict(color='k', label='All'),
	fill_kwargs=dict(color='k', alpha=0.3),
	sort_by=ip_array.mean(axis=1)
	)
	fig.line_axes.set_ylabel('Average enrichment');
	fig.line_axes.set_xlabel('Distance from Center (bp)');
	fig.array_axes.set_ylabel('Peaks')
	fig.array_axes.set_xticklabels([])
	fig.array_axes.axvline(0, linestyle=':', color='k')
	fig.line_axes.axvline(0, linestyle=':', color='k')
	fig.savefig('{}.png'.format(outname))