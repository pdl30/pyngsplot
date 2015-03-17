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
from pyngsplot.tools import peak_profiles, tss_profiles, metaseq_heatmaps, homer_analysis
import pybedtools

def ConfigSectionMap(section, Config):
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

def convert_ucsc_ens(peak, out):
	output = open(out, "w")
	with open(peak) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			chrom = word[0].strip("chr")
			output.write("{}\t{}\t{}\n".format(chrom, word[1], word[2])),
	output.close()

def main():
	parser = argparse.ArgumentParser(description='Similar to NGS plot, uses RPM normalisation by default. Plots both heatmaps and average profiles. Do not run in parallel\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")

	peak_parser = subparsers.add_parser('peak', help="Plot Peak profiles and heatmap")
	tss_parser = subparsers.add_parser('tss', help="Plot TSS's from GTF")
	
	peak_parser.add_argument('-c', '--config', help='ConfigParser input file with [Conditions] and [Counts] for normalisation if required. Assumes BAM files are aligned to Ensembl. Conditions contains the output names for the plots. See examples', required=False)
	peak_parser.add_argument('-i', '--bed', help='Input peak file', required=True)
	peak_parser.add_argument('-s', '--size', help='Size of flanking region', default=2000, required=False)
	peak_parser.add_argument('-a', help='Use [Counts] from config', action="store_true", required=False)
	peak_parser.add_argument('-b', help='Use [Controls] from config', action="store_true", required=False)
	peak_parser.add_argument('-e', help='Use if BAM file is annotated to Ensembl', action="store_true", required=False)
	peak_parser.add_argument('-t', '--threads', help='Threads to use, default=8', default=8, required=False)
	peak_parser.add_argument('-n', '--name', help='Prefix to add to end of all plots', required=True)

	tss_parser.add_argument('-c', '--config', help='ConfigParser input file with [Conditions] and [Counts], [Controls] for normalisation if required. Assumes BAM files are aligned to Ensembl. Conditions contains the output names for the plots. See examples', required=False)
	tss_parser.add_argument('-i', '--gtf', help='Input GTF file in Ensembl format', required=True)
	tss_parser.add_argument('-s', '--size', help='Size of flanking region', default=2000, required=False)
	tss_parser.add_argument('-f', '--filter', help='Input list of genes to use for a filter, IDs must be be on first column', default=None, required=False)
	tss_parser.add_argument('-a', help='Use [Counts] from config', action="store_true", required=False)
	tss_parser.add_argument('-b', help='Use [Controls] from config', action="store_true", required=False)
	tss_parser.add_argument('-e', help='Use if BAM file is annotated to Ensembl', action="store_true", required=False)
	tss_parser.add_argument('-t', '--threads', help='Threads to use, default=8', default=8, required=False)
	tss_parser.add_argument('-n', '--name', help='Prefix to add to end of all plots', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)
	if args["a"]:
		ConfigSectionMap("Counts")
		counts = ConfigSectionMap("Counts", Config)
	else:
		counts = None
	if args["b"]:
		ConfigSectionMap("Controls", Config)
		controls = ConfigSectionMap("Controls", Config)
	else:
		controls = None

	if args["subparser_name"] == "peak":
		if args["e"]:
			bed_file = metaseq_heatmaps.change_bed_size(args["bed"], args["size"], True)
		else:
			bed_file = metaseq_heatmaps.change_bed_size(args["bed"], args["size"], False)
		peak = pybedtools.BedTool(bed_file)
		metaseq_heatmaps.metaseq_heatmap(conditions, peak, counts, args["size"], controls, args["threads"], args["name"])

	elif args["subparser_name"] == "tss":

		gtf = metaseq_heatmaps.preprocess_gtf(args["gtf"], args["filter"])
		metaseq_heatmaps.metaseq_heatmap(conditions, gtf, counts, args["size"], controls, args["threads"], args["name"])