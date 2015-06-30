#!/usr/bin/python
import sys, os
from optparse import OptionParser
from Bio import SeqIO
import Bio



options = None



if __name__ == '__main__':

	usage = "usage: %prog arg"
	parser = OptionParser(usage)
	parser.add_option("-i", "--input", dest="input", help="input fasta")
	parser.add_option("-o", "--outputFolder", dest="outputFolder", help="input fasta")
	
	parser.add_option("-n", "--minlen", dest="minlen", help="min acceptable length <= (default=0)", default=0, type=int)
	parser.add_option("-x", "--maxlen", dest="maxlen", help="max acceptable length >= (default=-1: all)", default=-1, type=int)
	parser.add_option("", "--numberFileNames", dest="numberFileNames", help="Generate numbered output files instead of TAXON.fasta", default=False, action='store_true')

	
	#parser.add_option("-s", "--filesize", dest="filesize", help="size of output files in number of seqs", type=int)
	
	
	(options, args) = parser.parse_args()

	if len(args) != 0 or not options.input or not options.outputFolder:
		#print len(args)
		parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")
	
	#print options.filesize

	file_count = 0
	
	output_handle = None 
	output_taxon = None 

	for record in SeqIO.parse(options.input, "fasta") :
		seq_size = len(record.seq)
		if seq_size >=options.minlen and (options.maxlen<0 or seq_size<=options.maxlen):

			taxon = record.id.split('|')[0]


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


			SeqIO.write([record], output_handle, 'fasta')
			
#	input_file.close()
	output_handle.close()
	print 'Done.'
