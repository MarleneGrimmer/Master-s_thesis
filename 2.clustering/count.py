#!/usr/bin/python

import sys, os, glob
from Bio import SeqIO

counter = 0
with open ("htgcf_mibig.gbk", "r") as infile:
	for record in SeqIO.parse(infile, "genbank"):
		counter += 1

print(counter)
