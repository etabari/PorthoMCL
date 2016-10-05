#!/usr/bin/python

from ftplib import FTP
from optparse import OptionParser
import os, re

options = None




if __name__ == '__main__':

	usage = "usage: %prog [options] LOCALFOLDER\n\nExample: %prog -g Escherichia_coli downloaded\n         will download all assemblies for 'E. coli' into folder 'downloaded'"
	parser = OptionParser(usage)
	# parser.add_option("-f", "--localfolder", dest="localfolder", help="The folder of The Genome containing the fna files will be copied to")
	parser.add_option("-s", "--server", dest="host", help="ncbi FTP Host Address (default=ftp.ncbi.nlm.nih.gov)", default='ftp.ncbi.nlm.nih.gov')
	parser.add_option("-d", "--directory", dest="directory", help="ncbi FTP Host directory (default=/genomes/genbank/bacteria/)", default='/genomes/genbank/bacteria/')
	parser.add_option("-g", "--genome", dest='genome', help='a regular expression for the name of the genome (default=.*)', default='.*')
	parser.add_option("-a", "--all", dest='getall', help='download all files for assembly (default=FALSE, only FAA)', action="store_true", default=False)
	parser.add_option("-n", "--noprotein", dest='noprotein', help='download even no protein data available (default=FALSE)', action="store_true", default=False)

	
	(options, args) = parser.parse_args()

	if len(args) != 1:
		parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")

	localfolder = args[0]


	print('Connecting to FTP server: ', options.host)
	ftp = FTP(options.host)
	print('login annonynously...')
	ftp.login()

	print 'change directory to ', options.directory
	ftp.cwd(options.directory)

	genomes = ftp.nlst()
	total = len(genomes)

	print 'Total genomes in this ftp folder: ', total

	reg = re.compile(options.genome)
	genomes = filter(reg.match, genomes)

	total = len(genomes)
	print 'Total genomes in match criteria: ', total


	errors = dict()
	index = 0


	for genome in genomes:
		try:
			if index>10:
				break
			ftp.cwd(options.directory+genome+'/latest_assembly_versions')
			index += 1
			print index,'/',total, ': ', genome ,

			# try:
			# 	os.makedirs(os.path.join(localfolder, genome))
			# except:
			# 	if not errors.has_key(genome):
			# 		errors[genome] = ['FOLDER']

			assemblies = ftp.nlst()

			print ' contains ', len(assemblies), ' assemblies...'

			for assembly in assemblies:

				ftp.cwd(options.directory+genome+'/latest_assembly_versions/'+assembly)
				files = ftp.nlst()

				if not options.noprotein and not assembly+'_protein.faa.gz' in files:
					continue

				os.makedirs(os.path.join(localfolder, genome+'_'+assembly))

				for anyfile in files:

					if options.getall or anyfile.endswith('report.txt') or  anyfile.endswith('md5checksums.txt') or anyfile.endswith('faa.gz'):
						try:
							local_filename = os.path.join(localfolder, genome+'_'+assembly,anyfile)
							local_file = open(local_filename, 'wb')
							ftp.retrbinary('RETR ' + anyfile, local_file.write)
							local_file.close()
							print '\t\t', anyfile
						except Exception as detail:
							if not errors.has_key(genome):
								errors[genome] = []
							errors[genome] += [anyfile]
							print '\t\t', anyfile,'\t','[FAILED]\n\t\t', detail

		except Exception as detail:
			index -= 1
			print detail


	ftp.close()
