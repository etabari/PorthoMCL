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

class SimilarSequenceLine:
	def __init__(self, query_taxon, query_id, subject_taxon,  subject_id, evalue_mant, evalue_exp, percent_ident, percent_match):
		self.query_taxon = query_taxon
		self.query_id = query_id
		self.subject_taxon = subject_taxon
		self.subject_id = subject_id
		self.evalue_mant = evalue_mant
		self.evalue_exp = evalue_exp
		self.percent_ident = percent_ident
		self.percent_match = float(percent_match)


def readBestHit(file_name):
	best_hits = {}
	with open(file_name) as best_hit_file:
		for line in best_hit_file:
			cols = line.strip().split('\t')
			
			
			query_id = cols[0].split('|')[1]

			(subject_taxon, subject_id) = cols[1].split('|')

			try: 
				a = best_hits[subject_taxon]
			except:
				best_hits[subject_taxon] = {}

			try:
				best_hits[subject_taxon][query_id] += [(subject_id, int(cols[2]), float(cols[3]))]
			except:
				best_hits[subject_taxon][query_id] = [(subject_id, int(cols[2]), float(cols[3]))]
	return best_hits




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

	parser.add_option('-s', '--inSimSeq', dest='inSimSeq', help='folder that stores Similar Sequence files (TaxonID.ss.tsv) ')
	parser.add_option('-b', '--inBestHitFolder', dest='inBestHitFolder', help='folder that stores Best Hit files (TaxonID.bh.tsv)')
	parser.add_option('-q', '--inQueryTaxonScoreFolder', dest='inQueryTaxonScoreFolder', help='folder thast stores best query-taxon evalue score  (TaxonID.q-t.tsv)')
	
	parser.add_option('-p', '--outInParagogFolder', dest='outInParagogFolder', help='folder that will stores TaxonID.ort.tsv files')
	
	parser.add_option('', '--evalueExponentCutoff', dest='evalueExponentCutoff', help='evalue Exponent Cutoff (a nebative value, default=-5)', default=-5, type='int')
	parser.add_option('', '--percentMatchCutoff', dest='percentMatchCutoff', help='percent Match Cutoff (integer value, default=50)', default=50, type='int')
	
	parser.add_option('', '--cacheInputFile', dest='cacheInputFile', help='Cache input file or read it again. (Only use if I/O is very slow)', default=False, action="store_true")
	
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


	log('{2} | InParalogs | {0} | {1} | {3} | {4} MB | {5}'.format(2 , 'reading query-taxon evalue score (q-t file)', options.index, taxon1s , memory_usage_resource(), datetime.now() ))

	BestInterTaxonScore = {}
	taxon1_qt_filename =  os.path.join(options.inQueryTaxonScoreFolder , taxon1s + '.q-t.tsv')
	with open(taxon1_qt_filename) as taxon1_qt_file:
		for line in taxon1_qt_file:
			cols = line.strip().split()
			query_id = cols[0]
			taxon2s = cols[1]
			ev_exp = int(cols[2])
			ev_mant = float(cols[3])

			try:
				(min_exp, mants) = BestInterTaxonScore[query_id]
				if ev_exp < min_exp:
					BestInterTaxonScore[query_id] = (ev_exp, [ev_mant])
				elif ev_exp == min_exp:
					BestInterTaxonScore[query_id] = (ev_exp, mants+[ev_mant])
			except:
				BestInterTaxonScore[query_id] = (ev_exp, [ev_mant])


	for query_id in BestInterTaxonScore.keys():
		(ev_exp, mants) = BestInterTaxonScore[query_id]
		BestInterTaxonScore[query_id] = (ev_exp, min(mants))

		if DEBUG:
			print query_id, ev_exp, mants


	log('{2} | Orthology | {0} | {1} | {3} | {4} MB | {5}'.format(5, 'Finished' , options.index, taxon2s , memory_usage_resource(), datetime.now() ))
