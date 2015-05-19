#!/usr/bin/python
import sys
from optparse import OptionParser
from Bio import SeqIO
import Bio



options = None



if __name__ == '__main__':

	usage = "usage: %prog arg"
	parser = OptionParser(usage)
	parser.add_option("-i", "--input", dest="input", help="input fasta")
	
	parser.add_option("-n", "--minlen", dest="minlen", help="min acceptable length <= (default=0)", default=0, type=int)
	parser.add_option("-x", "--maxlen", dest="maxlen", help="max acceptable length >= (default=-1: all)", default=-1, type=int)
	
	parser.add_option("-s", "--filesize", dest="filesize", help="size of output files in number of seqs", type=int)
	
	
	(options, args) = parser.parse_args()

	if len(args) != 0 or not options.input or not options.filesize:
		#print len(args)
		parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")
	
	print options.filesize

	file_count = 1
	count = 0
	
	output_handle = open(options.input + '.' + str(file_count), "w")

	for record in SeqIO.parse(options.input, "fasta") :
		seq_size = len(record.seq)
		if seq_size >=options.minlen and (options.maxlen<0 or seq_size<=options.maxlen):
			count += 1
			SeqIO.write([record], output_handle, 'fasta')
			

			if count > options.filesize:
				file_count += 1
				output_handle.close()
				print file_count,
				sys.stdout.flush()
				output_handle = open(options.input + '.' + str(file_count), "w")
				count = 0

#	input_file.close()
	output_handle.close()
