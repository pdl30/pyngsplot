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
import ConfigParser
import HTSeq
import matplotlib 
import subprocess
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import multiprocessing
import metaseq
import pybedtools
import numpy as np
import pkg_resources
import pysam 

def ConfigSectionMap(section):
	dict1 = {}
	options = Config.options(section)
	for option in options:
		try:
			dict1[option] = Config.get(section, option)
			if dict1[option] == -1:
				DebugPrint("skip: %s" % option)
		except:
			print("exception on %s!" % option)
			dict1[option] = None
	return dict1

def index_bam(conditions):
	for key in conditions:
		command = "samtools index {}".format(key)

def metaseq_heatmap(conditions, bed, glen, window, threads):
	#Must figure out how to work with window
	threads = int(threads)
	fig = plt.figure()
	ax = fig.add_subplot(111)

	for key in sorted(conditions):
		ip_signal = metaseq.genomic_signal(key, 'bam')

		# Create arrays in parallel
		ip_array = ip_signal.array(bed, bins=100, processes=threads)
		ip_array = ip_array*1e9

		# Normalize to library size
		
		diviser = glen * ip_signal.mapped_read_count() 
		ip_array = np.transpose(ip_array)
		ip_array = ip_array/diviser
		ip_array = np.transpose(ip_array)
		x = np.linspace(-window, 0, 100)
		ax.plot(x, ip_array.mean(axis=0), label=conditions[key])
		# Add a vertical line at the TSS
		
		x = np.linspace(-window, 0, 100)
		fig2 = metaseq.plotutils.imshow(ip_array, x=x, figsize=(7, 10),
		vmin=5, vmax=99,  percentile=True,
		line_kwargs=dict(color='k', label='All'),
		fill_kwargs=dict(color='k', alpha=0.3),
		sort_by=ip_array.mean(axis=1))

		fig2.line_axes.set_ylabel('Average read coverage');
		fig2.line_axes.set_xlabel('Distance from TTS (bp)');
		#fig.array_axes.set_xticklabels([])
		fig2.array_axes.axvline(0, linestyle=':', color='k')
		fig2.line_axes.axvline(0, linestyle=':', color='k')
		fig2.savefig('{}_heatmap.png'.format(conditions[key]))
		plt.close(fig2)
	ax.axvline(0, linestyle=':', color='k')
	ax.set_xlabel('Distance from TTS (bp)')
	ax.set_ylabel('Average read coverage (per million mapped reads)')
	ax.legend(loc=1, fancybox=True, framealpha=0.5, prop={'size':7})
	fig.savefig('Average_profile.png')
	plt.close(fig)

def sam_size(ibam):
	results=  {}
	for bam in ibam:
		size = reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(bam) ])
		results[bam] = size
	return results

def using_htseq(conditions, bed, glen, w_size, threads, size_dict):
	fragmentsize = 200
	for key in sorted(conditions):
		#constant = 1000000/float(size_dict[key])
		bamfile = HTSeq.BAM_Reader( key )
		coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )
		profile = np.zeros( w_size, dtype='f')   
		for i, p in enumerate(bed):
			diviser = float(size_dict[key])*glen[i]
			constant = 1e9/float(diviser)
			#print i, diviser, glen[i], size_dict[key], constant
			#print glen[i], size_dict[key], diviser, 1e-9
			window = HTSeq.GenomicInterval( p[0], int(p[1]), int(p[2]), "." )
			for almnt in bamfile[ window ]:
				almnt.iv.length = fragmentsize
				start_in_window = almnt.iv.start - int(p[1])
				end_in_window   = almnt.iv.end - int(p[1]) 

				start_in_window = max( start_in_window, 0 )
				end_in_window = min( end_in_window, w_size )
				
				if start_in_window >= w_size or end_in_window < 0:
					continue
			#	print constant
				profile[ start_in_window : end_in_window ] += constant

		plt.plot( np.arange( -w_size, 0), profile, label=conditions[key])
		plt.legend(prop={'size':6})
	plt.savefig("total_profile.pdf")

def process_bed(ibed, size, ens):
	bed = ""
	glen = []
	with open(ibed) as f:
		for line in f:
			line = line.rstrip()
			word = line.split()
			length = int(word[2]) - int(word[1])
			if length >= 3000:
				coord = int(word[2]) - int(size)
				if ens:
					if is_int(word[0]):
						chrom = "chr"+word[0]
						coords = "{}\t{}\t{}\n".format(chrom, coord, word[2])
						bed = bed+coords
						ilen = int(word[2]) - int(word[1])
						glen.append(ilen)
					elif word[0] == "MT":
						chrom = "chrM"
						coords = "{}\t{}\t{}\n".format(chrom, coord, word[2])
						bed = bed+coords
						ilen = int(word[2]) - int(word[1])
						glen.append(ilen)
				else:
					coords = "{}\t{}\t{}\n".format(word[0], coord, word[2])
					bed = bed+coords
					ilen = int(word[2]) - int(word[1])
					glen.append(ilen)
					
	tss = pybedtools.BedTool(bed, from_string=True)
	glen = np.array(glen)
	return tss, glen #Checked and working perfectly

def is_int(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

def convert_ens_ucsc(gtf, out):
	output = open(out, "w")
	with open(gtf) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if is_int(word[0]):
				chrom = "chr"+word[0]
				output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, word[1], word[2], word[3], word[4], word[5], word[6], word[7], word[8])),
			elif word[0] == "MT":
				chrom = "chrM"
				output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, word[1], word[2], word[3], word[4], word[5], word[6], word[7], word[8])),
	output.close()

def main():
#if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Plots 3 prime biases in sequencing long genes.')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	metaseq_parser = subparsers.add_parser('meta', help="Metaseq plots")
	htseq_parser = subparsers.add_parser('htseq', help="HTseq plots")

	metaseq_parser.add_argument('-c', '--config', help='''ConfigParser input file with [Conditions]. This contains the output names for the plots. See examples''', required=True)
	metaseq_parser.add_argument('-b', '--bed', help='Input BED file with gene start and ends', required=True)
	metaseq_parser.add_argument('-s', '--size', help='Size of downstream region, default=2000', default=2000, required=False)
	metaseq_parser.add_argument('-e', action='store_true', help='Convert ensembl bed file to ucsc format', required=False)
	metaseq_parser.add_argument('-i', action='store_true', help='Will index all input bam files.', required=False)
	metaseq_parser.add_argument('-t', '--threads', help='Threads to use, default=8', default=8, required=False)

	htseq_parser.add_argument('-c', '--config', help='''ConfigParser input file with [Conditions]. This contains the output names for the plots. See examples''', required=True)
	htseq_parser.add_argument('-b', '--bed', help='Input BED file with gene start and ends', required=True)
	htseq_parser.add_argument('-s', '--size', help='Size of downstream region, default=2000', default=2000, required=False)
	htseq_parser.add_argument('-e', action='store_true', help='Convert ensembl bed file to ucsc format', required=False)
	htseq_parser.add_argument('-i', action='store_true', help='Will index all input bam files.', required=False)
	htseq_parser.add_argument('-t', '--threads', help='Threads to use, default=8', default=8, required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions")

	if args["subparser_name"] == "meta":
		if args["i"]:
			index_bam(conditions)
		gtf, glen = process_bed(args["bed"], int(args["size"]), args["e"])
		metaseq_heatmap(conditions, gtf, glen, int(args["size"]), args["threads"])
	elif args["subparser_name"] == "htseq":
		if args["i"]:
			index_bam(conditions)
		mapped_reads = sam_size(conditions)
		gtf, glen = process_bed(args["bed"], int(args["size"]), args["e"])
		using_htseq(conditions, gtf, glen, int(args["size"]), args["threads"], mapped_reads)