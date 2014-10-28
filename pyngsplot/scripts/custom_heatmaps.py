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
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm as CM
import pandas as pd
import subprocess
from plot_ngs_data.scripts import homer_heatmaps
import metaseq
import pybedtools
import multiprocessing
import ConfigParser

def plot_matrix(homer_file, label):
	data = list()
	with open(homer_file) as f:
		header = next(f)
		head = header.split()
		del head[0]
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			del word[0]
			data.append(word)
	b = np.array(data, dtype="f" )
	c = b[b[:,1].argsort()] #Sorts from lowest to highest!
	d = np.log2(c+1) #Not sure if correct with +1
	plt.locator_params(axis = 'x', nbins = 10) #Took ages to figure this out, do not forget! This sets the number of bins and must match the arange length!
	cmap = CM.get_cmap('Reds', 100) #Sets color scheme
	plt.imshow(d, interpolation='none', cmap=cmap, aspect="auto")
	plt.xticks(np.arange(0, 401, 50), np.arange(-2000, 2001, 500), rotation='vertical') #Key is to make them both exactly the same length!!
	plt.colorbar()
	plt.savefig(label+"_custom_heatmap.png") #Doesnt work with PDF for some reason, stick to PNG
	plt.close()

def pandas_plot_matrix(homer_file, label):
	data = list()
	index = []
	with open(homer_file) as f:
		header = next(f)
		head = header.split()
		del head[0]
		#for pandas, this is the columns
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			index.append(word[0])
			del word[0]
			data.append(word)
	b = np.array(data, dtype="f" ) #This is fine
	c = b[b[:,1].argsort()] #Sorts from lowest to highest!
	d = np.log2(c+1) #Not sure if correct with +1
	#Now convert this to pandas dataframe
	df = pd.DataFrame(d, index=index, columns=head)
	print df
	cmap = CM.get_cmap('Reds', 100) #Sets color scheme
	plt.imshow(df, interpolation='none', cmap=cmap, aspect="auto")
#	plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
	#plt.xticks(np.arange(0, len(df.columns), 500), df.columns)
	plt.xticks(np.arange(0, 1001, 100), np.arange(-5000, 5001, 1000), rotation='vertical') 
	plt.colorbar()
	plt.savefig(label+".png") #Doesnt work with PDF for some reason, stick to PNG
	plt.close()