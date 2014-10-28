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
from pyngsplot.scripts import peak_profiles, tss_profiles, metaseq_heatmaps, homer_analysis
import pybedtools

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

def convert_ucsc_ens(peak, out):
	output = open(out, "w")
	with open(peak) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			chrom = word[0].strip("chr")
			output.write("{}\t{}\t{}\n".format(chrom, word[1], word[2])),
	output.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Similar to NGS plot, uses RPM normalisation by default. Plots both heatmaps and average profiles. Do not run in parallel\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")

	peak_parser = subparsers.add_parser('peak', help="Plot Peak profiles and heatmap")
	tss_parser = subparsers.add_parser('tss', help="Plot TSS's from GTF")
	
	peak_parser.add_argument('-c', '--config', help='ConfigParser input file with [Conditions] and [Counts] for normalisation if required. Assumes BAM files are aligned to Ensembl. Conditions contains the output names for the plots. See examples', required=False)
	peak_parser.add_argument('-i', '--bed', help='Input peak file', required=True)
	peak_parser.add_argument('-s', '--size', help='Size of flanking region', default=2000, required=False)
	#peak_parser.add_argument('-meth', help='What method to use: Options are meta (default) and custom. Custom is still a work in progress and should not be used', default="meta", required=False)
	peak_parser.add_argument('-a', '--counts', help='Use [Counts] from config', action="store_true", required=False)
	peak_parser.add_argument('-b', '--cont', help='Use [Controls] from config', action="store_true", required=False)
	peak_parser.add_argument('-t', '--threads', help='Threads to use, default=8', default=8, required=False)

	tss_parser.add_argument('-c', '--config', help='ConfigParser input file with [Conditions] and [Counts], [Controls] for normalisation if required. Assumes BAM files are aligned to Ensembl. Conditions contains the output names for the plots. See examples', required=False)
	tss_parser.add_argument('-i', '--gtf', help='Input GTF file', required=False)
	tss_parser.add_argument('-s', '--size', help='Size of flanking region', default=2000, required=False)
	tss_parser.add_argument('-f', '--filter', help='Input list of genes to use for a filter, IDs must be be on first column', default=None, required=False)
	#tss_parser.add_argument('-meth', help='What method to use: Options are meta (default) and custom. Custom is still a work in progress and should not be used', required=False, default="meta")
	tss_parser.add_argument('-a', '--counts', help='Use [Counts] from config', action="store_true", required=False)
	tss_parser.add_argument('-b', '--cont', help='Use [Controls] from config', action="store_true", required=False)
	tss_parser.add_argument('-t', '--threads', help='Threads to use, default=8', default=8, required=False)

	args = vars(parser.parse_args())
	
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions")
	if args["counts"]:
		ConfigSectionMap("Counts")
		counts = ConfigSectionMap("Counts")
	else:
		counts = None
	if args["cont"]:
		ConfigSectionMap("Controls")
		controls = ConfigSectionMap("Controls")
	else:
		controls = None

	if args["subparser_name"] == "peak":
#		if args["meth"] == "meta":
			#Needs profiles but maybe this works
		metaseq_heatmaps.change_bed_size(args["bed"], "tmp_meta.bed", args["size"])
		peak = pybedtools.BedTool("tmp_meta.bed")
		metaseq_heatmaps.metaseq_heatmap(conditions, peak, counts, args["size"], controls, args["threads"])

#		elif args["meth"] == "custom":
			#Multiple problems to fix
			#Create profiles first
#			mapped_reads = peak_profiles.sam_size(conditions)
#			peak_profiles.parallel_peak_file_plot(conditions, args["ibed"], mapped_reads, args["size"], counts, controls) #Custom Average Profile
			#Now create heatmaps. First check if tag directory exists:
#			convert_ucsc_ens(args["ibed"], "tmp.bed")
			
			#Now have to include controls. No need for normalisation. Or is it library normalised by homer so be careful
#			for bam_file in sorted(conditions):
#				name = re.sub(".bam", "", bam_file)
#				name = re.sub("_sort", "", name)
#				if os.path.isdir("{}_homer".format(name)):
#					homer_file = homer_analysis.peak_homer_analysis("tmp.bed", tag="{}_homer".format(name)) 
#				else:
#					homer_file = homer_analysis.peak_homer_analysis("tmp.bed", bam_file=bam_file) #OK this is sorted
#				custom_heatmaps.plot_matrix(homer_file, conditions[bam_file]) #Size currently fixed to 4kb

	elif args["subparser_name"] == "tss":
#		if args["meth"] == "meta":
		gtf = metaseq_heatmaps.preprocess_gtf(args["gtf"], args["filter"])
		metaseq_heatmaps.metaseq_heatmap(conditions, gtf, counts, args["size"], controls, args["threads"])
#		elif args["meth"] == "custom":
#			#Same here
#			mapped_reads = peak_profiles.sam_size(conditions)
#			tss_profiles.plot_tss_profile(conditions, args["igtf"], mapped_reads, args["size"], args["filter"], counts, controls)
#			#Now for the heatmaps
#			for bam_file in sorted(conditions):
#				name = re.sub(".bam", "", bam_file)
#				name = re.sub("_sort", "", name)
#				if os.path.isdir("{}_homer".format(name)):
#					homer_file = homer_analysis.peak_homer_analysis("tmp.bed", tag="{}_homer".format(name)) 
#				else:
#					homer_file = homer_analysis.peak_homer_analysis("tmp.bed", bam_file=bam_file) #OK this is sorted
#				custom_heatmaps.plot_matrix(homer_file, conditions[bam_file]) #Size currently fixed to 4kb