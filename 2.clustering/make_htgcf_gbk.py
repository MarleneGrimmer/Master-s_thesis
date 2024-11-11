#!/usr/bin/python

import sys, os, glob
from Bio import SeqIO

for f in glob.glob("gbk/*.gbk"):          #find all files in gbk/ that end in .gbk
	fname = f.split("/")[-1].strip()  #take those files and split them at / -> take only filename (remove directory)
	fname = fname.replace(".gbk", "") #remove the .gbk from all files
	with open(f, "r") as infile:      #open all those files for reading
		for record in SeqIO.parse(infile, "genbank"):   #input sequence files in genbank format
			record.name = fname
			with open("htgcf.gbk", "a") as outfile: #open file for writing
				SeqIO.write(record, outfile, "genbank") #add the genbank sequences to the outfile
