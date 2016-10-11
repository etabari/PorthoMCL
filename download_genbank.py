#!/usr/bin/python

from ftplib import FTP
from optparse import OptionParser
import os, re
import hashlib

options = None


def hashfile(afile, hasher, blocksize=65536):
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    return hasher.hexdigest()

def download_anyfile(ftp, genome, assembly, anyfile, localfolder, options, md5data=None):
	try:
		if options.flatten:
			if not assembly in anyfile:
				local_filename = os.path.join(localfolder, genome+'_'+assembly+'_'+anyfile)
			else:
				local_filename = os.path.join(localfolder, genome+'_'+anyfile)
		else:
			local_filename = os.path.join(localfolder, genome+'_'+assembly,anyfile)

		local_file = open(local_filename, 'wb')
		ftp.retrbinary('RETR ' + anyfile, local_file.write)
		local_file.close()

		HashCheckResult = None
		if options.checkmd5 and md5data is not None:
			if anyfile in md5data:
				md5hash = hashfile(open(local_filename, 'rb'), hashlib.md5())
				HashCheckResult = md5hash == md5data[anyfile]
		print '\t\t', anyfile, HashCheckResult
	except Exception as detail:
		if not errors.has_key(genome):
			errors[genome] = []
		errors[genome] += [anyfile]
		print '\t\t', anyfile,'\t','[FAILED]\n\t\t', detail
	return local_filename


if __name__ == '__main__':

	usage = "usage: %prog [options] LOCALFOLDER\n\nExample: %prog -g Escherichia_coli downloaded\n         will download all assemblies for 'E. coli' into folder 'downloaded'"
	parser = OptionParser(usage)
	# parser.add_option("-f", "--localfolder", dest="localfolder", help="The folder of The Genome containing the fna files will be copied to")
	parser.add_option("-s", "--server", dest="host", help="ncbi FTP Host Address (default=ftp.ncbi.nlm.nih.gov)", default='ftp.ncbi.nlm.nih.gov')
	parser.add_option("-d", "--directory", dest="directory", help="ncbi FTP Host directory (default=/genomes/genbank/bacteria/)", default='/genomes/genbank/bacteria/')
	parser.add_option("-g", "--genome", dest='genome', help='a regular expression for the name of the genome (default=.*)', default='.*')
	parser.add_option("-a", "--all", dest='getall', help='download all files for assembly (default=FALSE, only FAA)', action="store_true", default=False)
	parser.add_option("-n", "--noprotein", dest='noprotein', help='download even no protein data available (default=FALSE)', action="store_true", default=False)
	parser.add_option("-f", "--flatten", dest='flatten', help="Download all in the same folder (don't make forlder for each assembly", action="store_true", default=False)
	parser.add_option("-m", "--md5", dest='checkmd5', help="check MD5 hash for downloaded files", action="store_true", default=False)

	
	(options, args) = parser.parse_args()

	if len(args) != 1:
		parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")

	localfolder = args[0]


	print 'Connecting to FTP server: ', options.host
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

	if options.flatten:
		try:
			os.makedirs(localfolder)
		except:
			print 'WARNING: output folder already exists'


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

				if not options.flatten:
					try:
						os.makedirs(os.path.join(localfolder, genome+'_'+assembly))
					except:
						print '\t\tWARNING: folder already exists'

				md5data = None
				if 'md5checksums.txt' in files:
					local_filename = download_anyfile(ftp, genome, assembly, 'md5checksums.txt',localfolder, options)
					if options.checkmd5:
						md5data = {}
						with open(local_filename) as md5file:
							for md5line in md5file:
								md5value = md5line.split('  ')
								md5data[md5value[1][2:-1]] = md5value[0]

				for anyfile in files:
					if anyfile != 'md5checksums.txt' and (options.getall or anyfile.endswith('report.txt') or anyfile.endswith('faa.gz')):
						download_anyfile(ftp, genome, assembly, anyfile, localfolder, options, md5data)

				exit(0)

		except Exception as detail:
			index -= 1
			print '>>', detail


	ftp.close()
