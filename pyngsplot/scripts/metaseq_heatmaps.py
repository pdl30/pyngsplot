#!/usr/bin/python

########################################################################
# 28 July 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, re, sys
import HTSeq
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import multiprocessing
import metaseq
import pybedtools
import numpy as np

#Changes bed window to custom size
def change_bed_size(bed, obed, size):
	output = open(obed, "w")
	with open(bed) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			mid = int(word[2]) - int(word[1])
			mid = mid/2
			mid = mid+int(word[1])
			chrom = word[0].strip("chr")
			if chrom == "M":
				chrom = "MT"
			to_add = float(size)/2
			start = int(mid)-to_add
			end = int(mid)+to_add
			if start < 0:
				start = 0
			output.write("{}\t{}\t{}\n".format(chrom, int(start), int(end))),

def read_counts(count):
	with open(count) as f:
		header= next(f)				
		header = header.rstrip()
		word = header.split("\t")
		gapdh = int(word[1])
	return gapdh

def metaseq_heatmap(conditions, bed, counts, window, controls, threads):
	#Must figure out how to work with window
#	db = gffutils.create_db(gtf, dbfn='test.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
	threads = int(threads)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	for key in sorted(conditions):
		ip_signal = metaseq.genomic_signal(key, 'bam')
		# Create arrays in parallel
		ip_array = ip_signal.array(bed, bins=100, processes=threads)
		# Normalize to library size
		if counts:
			gapdh = read_counts(counts[key])
			#print key, gapdh, float(gapdh)/1000,  1000/float(gapdh)
			ip_array *= 1000 / float(gapdh)
		else:
			ip_array /= ip_signal.mapped_read_count() / 1e6
		
		if controls:
			input_signal = metaseq.genomic_signal(controls[key], 'bam')
			input_array = input_signal.array(bed, bins=100, processes=threads)
			if counts:
				gapdh = read_counts(counts[key])
				input_array *= 1000/ float(gapdh)  #Test!!!
			else:
				input_array /= input_array.mapped_read_count() / 1e6
				
		x = np.linspace(-1000, 1000, 100)
		ax.plot(x, ip_array.mean(axis=0), label=conditions[key])
		# Add a vertical line at the TSS
		
		if controls:
			ip_array = ip_array - input_array #needs testing, not sure if working

		x = np.linspace(-1000, 1000, 100)
		fig2 = metaseq.plotutils.imshow(ip_array, x=x, figsize=(7, 10),
		vmin=5, vmax=99,  percentile=True,
		line_kwargs=dict(color='k', label='All'),
		fill_kwargs=dict(color='k', alpha=0.3),
		sort_by=ip_array.mean(axis=1))

		fig2.line_axes.set_ylabel('Average enrichment');
		fig2.line_axes.set_xlabel('Distance from Center (bp)');
		fig2.array_axes.set_ylabel('Peaks')
		#fig.array_axes.set_xticklabels([])
		fig2.array_axes.axvline(0, linestyle=':', color='k')
		fig2.line_axes.axvline(0, linestyle=':', color='k')
		fig2.savefig('{}_heatmap.png'.format(conditions[key]))
		plt.close(fig2)
	ax.axvline(0, linestyle=':', color='k')
	ax.set_xlabel('Distance from Center (bp)')
	ax.set_ylabel('Average read coverage (per million mapped reads)')
	ax.legend(loc=1, fancybox=True, framealpha=0.5, prop={'size':7})
	fig.savefig('Average_profile.png')
	plt.close(fig)

def preprocess_gtf(gtf, gene_filter):
	gtffile = HTSeq.GFF_Reader( gtf )
	output =  open("tmp_meta.bed", "w") 
	g_filter = {}
	tsspos = set()
	if gene_filter:
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
		for feature in gtffile:
			if feature.type == "exon" and feature.attr["exon_number"] == "1":
				tsspos.add( feature.iv.start_d_as_pos )
	for p in tsspos: #Problem if ensembl vs UCSC
		start = p.pos - 200 
		end = p.pos + 200
		if start < 0:
			start = 0
		output.write("{}\t{}\t{}\n".format(p.chrom, start, end)),
	output.close()
	tss = pybedtools.BedTool("tmp_meta.bed")
	return tss