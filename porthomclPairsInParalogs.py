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
	parser.add_option("-x", "--index", dest="index", help="An integer number identifying which taxon to work on [1-size_of_taxon_list]" , type='int')

	parser.add_option('-q', '--inInParalogTempFolder', dest='inInParalogTempFolder', help='Folder that contains Temporary Paralog pairs (pt files)')
	parser.add_option('-o', '--inOrthologGeneFolder', dest='inOrthologGeneFolder', help='Folder that contains List of genes that have orthologous relationship (og files)')
	
	parser.add_option('-p', '--outInParalogFolder', dest='outInParalogFolder', help='folder that will stores InParalogs (par files)')
	parser.add_option("-l", "--logfile", dest="logfile", help="log file (optional, if not supplied STDERR will be used)")
	
	parser.add_option('', '--KeepUnnormalizedScore', dest='KeepUnnormalizedScore', help='Write the un-normalize scores for paralog pairs (default=False)', default=False, action="store_true")
	
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
		OrthologUniqueId = readTaxonList(os.path.join(options.inOrthologGeneFolder, taxon1s + '.og.tsv' ))


	log('{2} | InParalogs | {0} | {1} | {3} | {4} MB | {5}'.format(3 , 'reading TempParalog file (pt file)', options.index, taxon1s , memory_usage_resource(), datetime.now() ))


	file_name = os.path.join(options.inInParalogTempFolder, taxon1s + '.pt.tsv' )
	paralog_temp = {}
	InParalogTaxonAvg = 0
	InplgOrthTaxonAvg = 0
	InplgOrthTaxonAvg_Count = 0
	with open(file_name) as paralog_temp_file:
		for line in paralog_temp_file:
			cols = line.strip().split('\t')
			
			query_id = cols[0]

			subject_id =  cols[1]

			unnormalized_score = float(cols[2])
			paralog_temp[(query_id,subject_id)] = unnormalized_score

			InParalogTaxonAvg += unnormalized_score

			if query_id in OrthologUniqueId or subject_id in OrthologUniqueId:
				InplgOrthTaxonAvg += unnormalized_score
				InplgOrthTaxonAvg_Count += 1


	if len(paralog_temp)>0:
		InParalogTaxonAvg /= len(paralog_temp)
	else:
		InParalogTaxonAvg = 'NA'
	
	if InplgOrthTaxonAvg_Count > 0:
		InplgOrthTaxonAvg /= InplgOrthTaxonAvg_Count
	else:
		InplgOrthTaxonAvg = 'NA'

	log('{2} | InParalogs | {0} | {1} | {3} | {4} MB | {5}'.format(4 , 'InplgOrthTaxonAvg: '+ str(InplgOrthTaxonAvg) +' : InParalogTaxonAvg:' + str(InParalogTaxonAvg) , options.index, taxon1s , memory_usage_resource(), datetime.now() ))

	if InplgOrthTaxonAvg == 'NA':
		InplgOrthTaxonAvg = InParalogTaxonAvg



	log('{2} | InParalogs | {0} | {1} | {3} | {4} MB | {5}'.format(5, 'writing the paralog file' , options.index, taxon1s , memory_usage_resource(), datetime.now() ))

	out_f = open (os.path.join(options.outInParalogFolder , taxon1s + '.par.tsv'), 'w')
	
	if options.KeepUnnormalizedScore:
		out_f.write('query_id\tsubject_id\tunnormalized_score\tnormalized_score\n')
	for (query_id,subject_id) in paralog_temp:
		norm = paralog_temp[(query_id,subject_id)] / InplgOrthTaxonAvg

		out_f.write(taxon1s +'|' + query_id + '\t')
		out_f.write(taxon1s +'|' + subject_id + '\t')
		if options.KeepUnnormalizedScore:
			out_f.write(str(paralog_temp[(query_id,subject_id)]) + '\t')
		else:
			norm = int( norm * 1000 + .5) / 1000.0 

		out_f.write(str(norm) + '\n')
	out_f.close()

