#!/usr/bin/env python
import time
from datetime import datetime
import sys, os
import gzip
import random, math
from optparse import OptionParser



options = None

## User for Orthology
best_query_taxon_score = {}
## Used for the Paralogy
BestInterTaxonScore = {}
BetterHit = {}


# class SimilarSequenceLine:
# 	def __init__(self, line):
# 		column = line.strip().split('\t')

# 		self.query_id = column[0]
# 		(self.query_taxon, self.query_seq) = column[0].split('|')

# 		self.subject_id = column[1]
# 		(self.subject_taxon,self.subject_seq)  = column[1].split('|')

# 		self.evalue_mant = float(column[2])
# 		self.evalue_exp = int(column[3])

# 		#self.percent_ident = column[4]
# 		self.percent_match = float(column[4])



class SimilarSequenceLine(namedtuple('SimilarSequenceLine', 'query_id,query_taxon,query_seq,subject_id,subject_taxon,subj_seq,evalue_mant,evalue_exp,percemt_match')):
	__slots__ = ()
	@classmethod
	def _fromLine(cls, line, new=tuple.__new__, len=len):
		'Make a new SimilarSequenceLine object from a sequence or iterable'
		column = line.strip().split('\\t')
		(query_taxon, query_seq) = column[0].split('|')
		(subject_taxon, subject_seq)  = column[1].split('|')
		iterable = (column[0], query_taxon, query_seq, column[1], subject_taxon, subject_seq, float(column[2]), int(column[3]), float(column[4]))
		result = new(cls, iterable)
		if len(result) != 9:
			raise TypeError('Expected 9 arguments, got %d' % len(result))
		return result


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


def writeStoOutputFiles(s, out_bh_file):
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
		
	if options.outInParalogTempFolder:
		try:
			(cutoff_exp, cutoff_mant) = BestInterTaxonScore[s.query_id]

			if (
				s.query_taxon == s.subject_taxon and 
				s.query_id != s.subject_id and 
				s.evalue_exp <= options.evalueExponentCutoff and 
				s.percent_match >= options.percentMatchCutoff and 
				(s.evalue_mant < 0.01 or s.evalue_exp<cutoff_exp or (s.evalue_exp == cutoff_exp and s.evalue_mant<=cutoff_mant))
			   ):

				# try:
				# 	BetterHit[(s.query_seq, s.subject_seq)] += [(s.evalue_exp, s.evalue_mant)]
				# except KeyError:
				BetterHit[(s.query_seq, s.subject_seq)] = (s.evalue_exp, s.evalue_mant)

		except KeyError:
			# Include the ones with
			if (
				s.query_taxon == s.subject_taxon and 
				(options.keepOrthoMCLBug or s.query_id != s.subject_id) and  #### THIS IS an OrthoMCL bug
				s.evalue_exp <= options.evalueExponentCutoff and 
				s.percent_match >= options.percentMatchCutoff
			   ):
				# try:
				# 	BetterHit[(s.query_seq, s.subject_seq)] += [(s.evalue_exp, s.evalue_mant)]
				# except KeyError:
				BetterHit[(s.query_seq, s.subject_seq)] = (s.evalue_exp, s.evalue_mant)


if __name__ == '__main__':
	usage = "This is STEP 5.1 of PorthoMCL.\n\nusage: %prog options\n"
	parser = OptionParser(usage)

	parser.add_option("-t", "--taxonlist", dest="taxonlistfile", help="A single column file containing the list of taxon to work with")
	parser.add_option("-x", "--index", dest="index", help="An integer number identifying which taxon to work on [1-size_of_taxon_list]" , type='int')

	parser.add_option('-s', '--inSimSeq', dest='inSimSeq', help='Input folder that contains split similar sequences files (ss files)')
	parser.add_option('-p', '--inInParalogFolder', dest='inInParalogFolder', help='Input folder that contains InParalogs (par files)')
	parser.add_option('-o', '--inOrthologFolder', dest='inOrthologFolder', help='Input folder that contains Orthologous pairs (ort files)')

	parser.add_option('-c', '--outCoOrthFolder', dest='outCoOrthFolder', help='folder that will stores CoOrtholog files (cor files)')
	parser.add_option("-l", "--logfile", dest="logfile", help="log file (optional, if not supplied STDERR will be used)")

	
	parser.add_option('', '--evalueExponentCutoff', dest='evalueExponentCutoff', help='evalue Exponent Cutoff (a nebative value, default=-5)', default=-5, type='int')
	parser.add_option('', '--percentMatchCutoff', dest='percentMatchCutoff', help='percent Match Cutoff (integer value, default=50)', default=50, type='int')
	
	#
	
	(options, args) = parser.parse_args()


	if len(args) != 0 or not options.taxonlistfile  or not options.index or not options.inSimSeq or not options.inInParalogFolder or not options.inOrthologFolder or not options.outCoOrthFolder:
		parser.error("incorrect arguments.\n\t\tUse -h to get more information or refer to the MANUAL.md")


	log('{2} | CoOrtholog | {0} | {1} | {3} | {4} MB | {5}'.format(1 , 'reading taxon list', options.index, '', memory_usage_resource(), datetime.now() ))
	taxon_list = readTaxonList(options.taxonlistfile)

	if options.index <= 0 or options.index > len(taxon_list):
		log('{2} | Best Hit | {0} | {1} | {3} | {4} MB | {5}'.format('ERROR' , 'Error in index', options.index,'', memory_usage_resource(), datetime.now() ))
		exit()

	taxon1s = taxon_list[options.index - 1]



	log('{2} | CoOrtholog | {0} | {1} | {3} | {4} MB | {5}'.format(2 , 'Reading similar sequences (ss file)', options.index, taxon1s, memory_usage_resource(), datetime.now() ))

	with open(os.path.join(options.inSimSeq, taxon1s+'.ss.tsv')) as input_file:
		for line in input_file:

			ss = SimilarSequenceLine._fromLine(line)


	if options.outInParalogTempFolder:
		log('{2} | CoOrtholog | {0} | {1} | {3} | {4} MB | {5}'.format(5 , 'Creating InParalogTemp file needed for InParalogs (pt file)', options.index, taxon1s, memory_usage_resource(), datetime.now() ))

		out_pt_file = open(os.path.join(options.outInParalogTempFolder, taxon1s+'.pt.tsv') ,'w')

		for (seq1, seq2) in BetterHit:
			
			if seq1 < seq2:

				(bh1_evalue_exp, bh1_evalue_mant) = BetterHit[(seq1, seq2)]

				try:
					(bh2_evalue_exp, bh2_evalue_mant) =  BetterHit[(seq2, seq1)]
				except:
					continue

				if bh1_evalue_mant < 0.01 or bh2_evalue_mant < 0.01:
					unnormalized_score = (bh1_evalue_exp + bh2_evalue_exp) / -2 
				else:
					unnormalized_score = (math.log10(bh1_evalue_mant * bh2_evalue_mant) + bh1_evalue_exp + bh2_evalue_exp) / -2


				out_pt_file.write('{0}\t{1}\t{2}\n'.format(seq1, seq2, unnormalized_score))


		out_pt_file.close()








	log('{2} | CoOrtholog | {0} | {1} | {3} | {4} MB | {5}'.format(6 , 'Done', options.index, taxon1s, memory_usage_resource(), datetime.now() ))







	



