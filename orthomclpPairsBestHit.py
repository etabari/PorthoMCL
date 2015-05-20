#!/usr/bin/python
import time
from datetime import datetime
import sys, os
import gzip
import random
from optparse import OptionParser



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
	usage = "This is STEP 4 of orthomclp.\n\nusage: %prog arg\n"
	parser = OptionParser(usage)

	parser.add_option("-t", "--taxonlist", dest="taxonlistfile", help="A single column file containing the list of taxon to work with")

	parser.add_option('-i', '--inputfolder', dest='inputfolder', help='folder that stores TaxonID.ss.tsv files (Split SimilarSeuence.tsv) ')
	parser.add_option('-o', '--outputfolder', dest='outputfolder', help='folder that will stores Best Hit files (If not set, current folder)')

	parser.add_option("-x", "--index", dest="index", help="an integer number identifying which taxon to work on [1-size_of_taxon_list]" , type='int')
	parser.add_option("-l", "--logfile", dest="logfile", help="log file")

	
	parser.add_option('', '--evalueExponentCutoff', dest='evalueExponentCutoff', help='evalue Exponent Cutoff (a nebative value, default=-5)', default=-5, type='int')
	parser.add_option('', '--percentMatchCutoff', dest='percentMatchCutoff', help='percent Match Cutoff (integer value, default=50)', default=50, type='int')
	
	parser.add_option('', '--cacheInputFile', dest='cacheInputFile', help='Cache input file or read it again. (Only use if I/O is very slow)', default=False, action="store_true")
	parser.add_option('', '--bestQueryTaxonScore', dest='bestQueryTaxonScore', help='folder to generate best query-taxon evalue score (optional)')
	
	#
	
	(options, args) = parser.parse_args()


	log('1   | reading taxon list | * |'+  str(memory_usage_resource())+ 'MB | '+ str(datetime.now()) )
	taxon_list = readTaxonList(options.taxonlistfile)
	log('1   | reading taxon list |  loaded ' + str(len(taxon_list)) +' taxons | '+  str(memory_usage_resource())+ 'MB | '+ str(datetime.now()) )

	if options.index <= 0 or options.index > len(taxon_list):
		log('ERROR | Error in index 1 to ' + str(len(taxon_list)))
		exit()

	taxon1s = taxon_list[options.index - 1]


	if options.cacheInputFile:
		log('OPTIONS   | cacheInputFile | Cache Enabled | ')


	log('2   | Reading similar sequences (ss file) for ' + taxon1s +' | * | '+  str(memory_usage_resource())+ 'MB | '+ str(datetime.now()) )

	best_query_taxon_score = {}

	input_file_cache = []

	with open(os.path.join(options.inputfolder, taxon1s+'.ss.tsv')) as input_file:
		for line in input_file:

			column = line.strip().split('\t')

			query = column[0].split('|')
			query_taxon = query[0]
			query_id = query[1]

			subject = column[1].split('|')
			subject_taxon = subject[0]

			if query_taxon != subject_taxon:

				evalue_mant = float(column[2])
				evalue_exp = int(column[3])

				if options.cacheInputFile:
					input_file_cache += [SimilarSequenceLine(query_taxon, query_id, subject_taxon, subject[1], evalue_mant, evalue_exp, column[4], column[5])]

				try:
					best_query_taxon_score[(query_id, subject_taxon)] += [(evalue_mant, evalue_exp)]
				except:
					best_query_taxon_score[(query_id, subject_taxon)] = [(evalue_mant, evalue_exp)]


	for (query_id,subject_taxon) in best_query_taxon_score:

		evalues = best_query_taxon_score[(query_id, subject_taxon)]
		
		min_exp = sys.maxint  #min(evalues, key = lambda t: t[1])

		min_mants = []

		for (evalue_mant, evalue_exp) in evalues:
			if evalue_exp < min_exp:
				min_exp = evalue_exp
				min_mants += [evalue_mant]

			if evalue_mant == 0 and evalue_exp == 0:
				min_mants += [evalue_mant]

		best_query_taxon_score[(query_id,subject_taxon)] = (min_exp, min(min_mants))


	if options.bestQueryTaxonScore:

		log('3   | Creating bestQueryTaxonScore (q-t file) for ' + taxon1s +' | * | '+  str(memory_usage_resource())+ 'MB | '+ str(datetime.now()) )

		with open(os.path.join(options.bestQueryTaxonScore, taxon1s+'.q-t.tsv'), 'w') as out_file:

			for (query_id,subject_taxon) in sorted(best_query_taxon_score):

				evalue = best_query_taxon_score[(query_id,subject_taxon)]
				out_file.write('{0}\t{1}\t{2}\t{3}\n'.format(query_id, subject_taxon, evalue[0], evalue[1]))


	log('4   | Creating BestHit (bh file) for ' + taxon1s +' | * | '+  str(memory_usage_resource())+ 'MB | '+ str(datetime.now()) )

	if not options.outputfolder:
		options.outputfolder = '.'
	out_file = open(os.path.join(options.outputfolder, taxon1s+'.bh.tsv') ,'w')

	if not options.cacheInputFile:
		with open(os.path.join(options.inputfolder, taxon1s+'.ss.tsv')) as input_file:
			for line in input_file:

				column = line.strip().split('\t')

				
				query = column[0].split('|')
				query_taxon = query[0]
				query_id = query[1]

				subject = column[1].split('|')
				subject_taxon = subject[0]

				if query_taxon != subject_taxon:

					try:
						cutoff = best_query_taxon_score[(query_id,subject_taxon)]
					except:
						continue



					evalue_mant = float(column[2])
					evalue_exp = int(column[3])
					percent_match = float(column[5])

					if evalue_exp < options.evalueExponentCutoff and percent_match > options.percentMatchCutoff and (evalue_mant < 0.01 or evalue_exp==cutoff[0] and evalue_mant==cutoff[1]):
						out_file.write('{0}|{1}\t{2}|{3}\t{4}\t{5}\n'.format(query_taxon, query_id, subject_taxon, subject[1], evalue_exp, evalue_mant))

	else:

		for s in input_file_cache:

			try:
				cutoff = best_query_taxon_score[(s.query_id, s.subject_taxon)]
			except:
				continue

			if s.evalue_exp < options.evalueExponentCutoff and s.percent_match > options.percentMatchCutoff and (s.evalue_mant < 0.01 or s.evalue_exp==cutoff[0] and s.evalue_mant==cutoff[1]):
				out_file.write('{0}|{1}\t{2}|{3}\t{4}\t{5}\n'.format(s.query_taxon, s.query_id, s.subject_taxon, s.subject_id, s.evalue_exp, s.evalue_mant))

	out_file.close()

	log('5   | Done ' + taxon1s +' | * | '+  str(memory_usage_resource())+ 'MB | '+ str(datetime.now()) )







	



