#!/usr/bin/python
import time
from datetime import datetime
import sys, os
import gzip
import random, math
from optparse import OptionParser


options = None
DEBUG = True


options = None


def readParalogTemp(taxon, file_name):


	return paralog_temp




def readTaxonList(filename):

	taxon_list = []
	taxon_list_file = open(filename)
	for line in taxon_list_file:
		line = line.strip()
		if line:
			taxon_list += [line]
	taxon_list_file.close()
	return taxon_list

def memory_usage_resource():
	import resource
	rusage_denom = 1024.
	if sys.platform == 'darwin':
		# ... it seems that in OSX the output is different units ...
		rusage_denom = rusage_denom * rusage_denom
	mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
	return round(mem,0)


def log(s):
	global options
	print >> sys.stderr, s
	if options.logfile:
		l = open(options.logfile, 'a')
		l.write(s+'\n')
		l.close()



if __name__ == '__main__':
	usage = "This is STEP 5.3 of PorthoMcl.\n\nusage: %prog arg\n"
	parser = OptionParser(usage)

	parser.add_option("-t", "--taxonlist", dest="taxonlistfile", help="A single column file containing the list of taxon to work with")
	parser.add_option("-x", "--index", dest="index", help="an integer number identifying which taxon to work on" , type='int')
	parser.add_option("-l", "--logfile", dest="logfile", help="logfile")

	parser.add_option('-o', '--inOrthologGeneFolder', dest='inOrthologGeneFolder', help='folder that stores TaxonID.og.tsv files (Optional)')
	parser.add_option('-q', '--inInParalogTempFolder', dest='inInParalogTempFolder', help='folder thast stores TempParalog unnormalized score  (TaxonID.pt.tsv)')
	
	parser.add_option('-p', '--outInParalogFolder', dest='outInParalogFolder', help='folder that will stores TaxonID.par.tsv files')
	
	#	#
	
	(options, args) = parser.parse_args()

	#log('-----')

	#print best_hit

	log('{2} | InParalogs | {0} | {1} | {3} | {4} MB | {5}'.format(1 , 'reading taxon list', options.index, '' ,memory_usage_resource(), datetime.now() ))

	taxon_list = readTaxonList(options.taxonlistfile)

	if options.index <= 0 or options.index > len(taxon_list):
		log('{2} | InParalogs | {0} | {1} | {3} | {4} MB | {5}'.format('ERROR' , 'Error in index', options.index, '', memory_usage_resource(), datetime.now() ))
		exit()


	taxon1s = taxon_list[options.index - 1]

	# if not options.OverwiteOutput and os.path.exists(os.path.join(options.outputfolder, taxon1s + '.ort.tsv')):
	# 	exit(0)


	log('{2} | InParalogs | {0} | {1} | {3} | {4} MB | {5}'.format(2 , 'reading Ortholog Gene file (og file)', options.index, taxon1s , memory_usage_resource(), datetime.now() ))

	OrthologUniqueId = []
	if options.inOrthologGeneFolder: 
		OrthologUniqueId = readTaxonList(options.inOrthologGeneFolder)


	log('{2} | InParalogs | {0} | {1} | {3} | {4} MB | {5}'.format(3 , 'reading TempParalog file (pt file)', options.index, taxon1s , memory_usage_resource(), datetime.now() ))

	paralog_temp = {}
	InParalogTaxonAvg = 0
	InplgOrthTaxonAvg = 0
	InplgOrthTaxonAvg_Count = 0
	with open(file_name) as best_hit_file:
		for line in best_hit_file:
			cols = line.strip().split('\t')
			
			query_id = taxon + '|' + cols[0]

			subject_id =  taxon + '|' + cols[1]

			unnormalized_score = float(cols[3])
			paralog_temp[(query_id,subject_id)] = unnormalized_score

			InParalogTaxonAvg += unnormalized_score

			if cols[0] in OrthologUniqueId or cols[1] in OrthologUniqueId:
				InplgOrthTaxonAvg += unnormalized_score
				InplgOrthTaxonAvg_Count += 1

	
	if InplgOrthTaxonAvg_Count > 0
		InplgOrthTaxonAvg /= InplgOrthTaxonAvg_Count
	else:
		InplgOrthTaxonAvg = 'NA'

	InParalogTaxonAvg /= len(paralog_temp)


	for (query_id,subject_id) in paralog_temp:
		paralog_temp[(query_id,subject_id)] = paralog_temp[(query_id,subject_id)] / InplgOrthTaxonAvg




	log('{2} | Orthology | {0} | {1} | {3} | {4} MB | {5}'.format(5, 'Finished' , options.index, taxon2s , memory_usage_resource(), datetime.now() ))


