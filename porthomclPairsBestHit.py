#!/usr/bin/python
import time
from datetime import datetime
import sys, os
import gzip
import random
from optparse import OptionParser



options = None
best_query_taxon_score = {}
BestInterTaxonScore = {}


class SimilarSequenceLine:
	def __init__(self, line):
		column = line.strip().split('\t')

		self.query_id = column[0]
		(self.query_taxon, self.query_seq) = column[0].split('|')

		self.subject_id = column[1]
		(self.subject_taxon,self.subject_seq)  = column[1].split('|')

		self.evalue_mant = float(column[2])
		self.evalue_exp = int(column[3])

		self.percent_ident = column[4]
		self.percent_match = float(column[5])



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


def writeStoOutputFiles(s, out_bh_file, out_br_file):
	global best_query_taxon_score, BestInterTaxonScore, options
	try:
		(cutoff_exp, cutoff_mant) = best_query_taxon_score[(s.query_id, s.subject_taxon)]

		if ( 
			s.query_taxon != s.subject_taxon and
			s.evalue_exp < options.evalueExponentCutoff and 
			s.percent_match > options.percentMatchCutoff and 
			(s.evalue_mant < 0.01 or s.evalue_exp==cutoff_exp and s.evalue_mant==cutoff_mant)
		   ):
			out_bh_file.write('{0}\t{1}\t{2}\t{3}\n'.format(s.query_seq, s.subject_id, s.evalue_exp, s.evalue_mant))

	except KeyError:
		pass
		
	try:
		(cutoff_exp, cutoff_mant) = BestInterTaxonScore[s.query_id]

		if (
			s.query_taxon == s.subject_taxon and 
			s.query_id != s.subject_id and 
			s.evalue_exp <= options.evalueExponentCutoff and 
			s.percent_match >= options.percentMatchCutoff and 
			(s.evalue_mant < 0.01 or s.evalue_exp<cutoff_exp or (s.evalue_exp == cutoff_exp and s.evalue_mant<=cutoff_mant))
		   ):
			out_br_file.write('{0}\t{1}\t{2}\t{3}\n'.format(s.query_seq, s.subject_seq, s.evalue_exp, s.evalue_mant))

	except KeyError:
		# Include the ones with
		if (
			s.query_taxon == s.subject_taxon and 
			(options.keepOrthoMCLBug or s.query_id != s.subject_id) and  #### THIS IS an OrthoMCL bug
			s.evalue_exp <= options.evalueExponentCutoff and 
			s.percent_match >= options.percentMatchCutoff
		   ):
			out_br_file.write('{0}\t{1}\t{2}\t{3}\n'.format(s.query_seq, s.subject_seq, s.evalue_exp, s.evalue_mant))


if __name__ == '__main__':
	usage = "This is STEP 5.1 of PorthoMCL.\n\nusage: %prog arg\n"
	parser = OptionParser(usage)

	parser.add_option("-t", "--taxonlist", dest="taxonlistfile", help="A single column file containing the list of taxon to work with")

	parser.add_option('-s', '--inSimSeq', dest='inSimSeq', help='folder that stores TaxonID.ss.tsv files (Split SimilarSequence.tsv) ')
	parser.add_option('-b', '--outBestHitFolder', dest='outBestHitFolder', help='folder that will stores Best Hit files (If not set, current folder)')
	parser.add_option('-q', '--outQueryTaxonScoreFolder', dest='outQueryTaxonScoreFolder', help='folder to generate best query-taxon evalue scores (q-t, q+t files) (required for Paralogs)')

	parser.add_option("-x", "--index", dest="index", help="an integer number identifying which taxon to work on [1..size_of_taxon_list]" , type='int')
	parser.add_option("-l", "--logfile", dest="logfile", help="log file")

	
	parser.add_option('', '--evalueExponentCutoff', dest='evalueExponentCutoff', help='evalue Exponent Cutoff (a nebative value, default=-5)', default=-5, type='int')
	parser.add_option('', '--percentMatchCutoff', dest='percentMatchCutoff', help='percent Match Cutoff (integer value, default=50)', default=50, type='int')
	
	parser.add_option('', '--cacheInputFile', dest='cacheInputFile', help='Cache input file or read it again. (Only use if I/O is very slow)', default=False, action="store_true")
	parser.add_option('', '--keepOrthoMCLBug', dest='keepOrthoMCLBug', help='Keep the OrthoMCL bug in creating BetterHit files (br) where self hits are included', default=False, action="store_true")
	
	#
	
	(options, args) = parser.parse_args()


	log('{2} | Best Hit | {0} | {1} | {3} | {4} MB | {5}'.format(1 , 'reading taxon list', options.index, '', memory_usage_resource(), datetime.now() ))
	taxon_list = readTaxonList(options.taxonlistfile)

	if options.index <= 0 or options.index > len(taxon_list):
		log('{2} | Best Hit | {0} | {1} | {3} | {4} MB | {5}'.format('ERROR' , 'Error in index', options.index,'', memory_usage_resource(), datetime.now() ))
		exit()

	taxon1s = taxon_list[options.index - 1]


	if options.cacheInputFile:
		log('{2} | Best Hit | {0} | {1} | {3} | {4} MB | {5}'.format('OPTION' , 'Caching Input files', options.index, taxon1s, memory_usage_resource(), datetime.now() ))


	log('{2} | Best Hit | {0} | {1} | {3} | {4} MB | {5}'.format(2 , 'Reading similar sequences (ss file)', options.index, taxon1s, memory_usage_resource(), datetime.now() ))


	input_file_cache = []

	with open(os.path.join(options.inSimSeq, taxon1s+'.ss.tsv')) as input_file:
		for line in input_file:

			ss = SimilarSequenceLine(line)


			if options.cacheInputFile:
				input_file_cache += [ss]


			if ss.query_taxon != ss.subject_taxon:

				try:
					best_query_taxon_score[(ss.query_id, ss.subject_taxon)] += [(ss.evalue_mant, ss.evalue_exp)]
				except:
					best_query_taxon_score[(ss.query_id, ss.subject_taxon)] = [(ss.evalue_mant, ss.evalue_exp)]


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


	if options.outQueryTaxonScoreFolder:

		# log('{2} | Best Hit | {0} | {1} | * | {3} MB | {4}'.format(3 , 'Creating bestQueryTaxonScore (q-t file)', options.index, memory_usage_resource(), datetime.now() ))
		# with open(os.path.join(options.outQueryTaxonScoreFolder, taxon1s+'.q-t.tsv'), 'w') as out_file:
		# 	for (query_id,subject_taxon) in sorted(best_query_taxon_score):
		# 		(ev_exp, ev_mant) = best_query_taxon_score[(query_id,subject_taxon)]
		# 		out_file.write('{0}\t{1}\t{2}\t{3}\n'.format(query_id, subject_taxon, ev_exp, ev_mant))


		log('{2} | Best Hit | {0} | {1} | {3} | {4} MB | {5}'.format(3 , 'Creating BestInterTaxonScore (q-t file)', options.index,taxon1s, memory_usage_resource(), datetime.now() ))

		
		for (query_id,subject_taxon) in best_query_taxon_score:

			(ev_exp, ev_mant) = best_query_taxon_score[(query_id,subject_taxon)]
			try:
				(min_exp, mants) = BestInterTaxonScore[query_id]
				if ev_exp < min_exp:
					BestInterTaxonScore[query_id] = (ev_exp, [ev_mant])
				elif ev_exp == min_exp:
					BestInterTaxonScore[query_id] = (ev_exp, mants+[ev_mant])
			except:
				BestInterTaxonScore[query_id] = (ev_exp, [ev_mant])



		with open(os.path.join(options.outQueryTaxonScoreFolder, taxon1s+'.q-t.tsv'), 'w') as out_file:

			for query_id in sorted(BestInterTaxonScore):

				(ev_exp, ev_mants) = BestInterTaxonScore[query_id]
				out_file.write('{0}\t{1}\t{2}\n'.format(query_id, ev_exp, min(ev_mants)))



	log('{2} | Best Hit | {0} | {1} | {3} | {4} MB | {5}'.format(4 , 'Creating BestHit/BetterHit (bh/br file)', options.index, taxon1s, memory_usage_resource(), datetime.now() ))

	BestHit = {}

	if not options.outBestHitFolder:
		options.outBestHitFolder = '.'
	out_bh_file = open(os.path.join(options.outBestHitFolder, taxon1s+'.bh.tsv') ,'w')
	out_br_file = open(os.path.join(options.outBestHitFolder, taxon1s+'.br.tsv') ,'w')

	if not options.cacheInputFile:
		with open(os.path.join(options.inSimSeq, taxon1s+'.ss.tsv')) as input_file:
			for line in input_file:
				s = SimilarSequenceLine(line)
				writeStoOutputFiles(s, out_bh_file, out_br_file)


	else:

		for s in input_file_cache:
			writeStoOutputFiles(s, out_bh_file, out_br_file)


	out_bh_file.close()
	out_br_file.close()
	log('{2} | Best Hit | {0} | {1} | {3} | {4} MB | {5}'.format(5 , 'Done', options.index, taxon1s, memory_usage_resource(), datetime.now() ))







	



