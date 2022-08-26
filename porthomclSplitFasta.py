#!/usr/bin/env python
import sys, os
from optparse import OptionParser

options = None



if __name__ == '__main__':

	usage = "usage: %prog arg"
	parser = OptionParser(usage)
	parser.add_option("-i", "--input", dest="input", help="input fasta")
	parser.add_option("-o", "--outputFolder", dest="outputFolder", help="input fasta")
	# parser.add_option("-n", "--minlen", dest="minlen", help="min acceptable length <= (default=0)", default=0, type=int)
	# parser.add_option("-x", "--maxlen", dest="maxlen", help="max acceptable length >= (default=-1: all)", default=-1, type=int)
	parser.add_option("", "--numberFileNames", dest="numberFileNames", help="Generate numbered output files instead of TAXON.fasta", default=False, action='store_true')
	#parser.add_option("-s", "--filesize", dest="filesize", help="size of output files in number of seqs", type=int)
	
	
	(options, args) = parser.parse_args()

	if len(args) != 0 or not options.input or not options.outputFolder:
		parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")
	


	file_count = 0
	
	output_handle = None 
	output_taxon = None 

	infile = open(options.input)
	for line in infile:
		if line[0]=='#':
			continue
		if line[0]=='>':
			record_id = line[1:-1]
			taxon = record_id.split('|')[0]
			
			if output_taxon != taxon:
				file_count += 1
				output_taxon = taxon
				if output_handle:
					output_handle.close()
				if options.numberFileNames:
					output_handle = open(os.path.join(options.outputFolder, str(file_count) + '.fasta'), "w")
				else:
					output_handle = open(os.path.join(options.outputFolder, output_taxon + '.fasta'), "w")
				
				print file_count, output_taxon

		output_handle.write(line)

	output_handle.close()
	print 'Done.'
