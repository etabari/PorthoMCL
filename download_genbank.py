#!/usr/bin/python

from ftplib import FTP
from optparse import OptionParser
import os

options = None




if __name__ == '__main__':

	usage = "usage: %prog [options] LOCALFOLDER"
	parser = OptionParser(usage)
	# parser.add_option("-f", "--localfolder", dest="localfolder", help="The folder of The Genome containing the fna files will be copied to")
	parser.add_option("-h", "--host", dest="host", help="ncbi FTP Host Address (default=ftp.ncbi.nlm.nih.gov)", default='ftp.ncbi.nlm.nih.gov')
	parser.add_option("-d", "--directory", dest="directory", help="ncbi FTP Host directory (default=/genomes/genbank/bacteria/)", default='/genomes/genbank/bacteria/')

	
	(options, args) = parser.parse_args()

	if len(args) != 1:
		parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")

	localfolder = args[1]


	print('Connecting to FTP server: ', options.host)
	ftp = FTP(options.host)
	print('login annonynously...')
	ftp.login()

	print('change directory to ', options.directory)
	ftp.cwd(options.directory)

	genomes = ftp.nlst()
	total = len(genomes)

	print 'Total genomes in this ftp folder: ', total

	errors = dict()
	index = 0


	for genome in genomes:
		try:
			ftp.cwd(directory+genome+'/latest_assembly_versions')
			index += 1
			print index,'/',total, ': ', genome ,

			# try:
			# 	os.makedirs(os.path.join(localfolder, genome))
			# except:
			# 	if not errors.has_key(genome):
			# 		errors[genome] = ['FOLDER']

			files = ftp.nlst()

			print ' (contains ', len(files), ' assemblies)'
			continue

			for anyfile in files:
				if anyfile.endswith('faa') or anyfile.endswith('fna') or anyfile.endswith('ptt'):
					try:
						local_filename = os.path.join(localfolder, genome ,anyfile)
						print 'getting ', anyfile, '->', local_filename
						local_file = open(local_filename, 'wb')
						ftp.retrbinary('RETR ' + anyfile, local_file.write)
						local_file.close()
					except:
						if not errors.has_key(genome):
							errors[genome] = []
						errors[genome] += [anyfile]

			
		except Exception as detail:
			index -= 1
			print detail


	ftp.close()

	print 'ERROR'
	for v,k in errors.iteritems():
		print v,k
	print
	print 'STATS'