#!/usr/bin/python
import time
from datetime import datetime
import sys, os
import gzip
import random, math
from optparse import OptionParser


options = None
DEBUG = False

def readTaxonList(filename):

	taxon_list = []
	taxon_list_file = open(filename)
	for line in taxon_list_file:
		line = line.strip()
		if line:
			taxon_list += [line]
	taxon_list_file.close()
	return taxon_list


def readBestHit(query_taxon, file_name):
	best_hits = {}
	with open(file_name) as best_hit_file:
		for line in best_hit_file:
			cols = line.strip().split('\t')
			
			
			query_id = query_taxon+'|'+cols[0]

			subject_id = cols[1]

			subject_taxon = subject_id.split('|')[0]

			try: 
				a = best_hits[subject_taxon]
			except:
				best_hits[subject_taxon] = {}

			try:
				best_hits[subject_taxon][query_id] += [(subject_id, int(cols[2]), float(cols[3]))]
			except:
				best_hits[subject_taxon][query_id] = [(subject_id, int(cols[2]), float(cols[3]))]
	return best_hits




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
	usage = "This is STEP 5.2 of PorthoMcl.\n\nusage: %prog options\n"
	parser = OptionParser(usage)

	parser.add_option("-t", "--taxonlist", dest="taxonlistfile", help="A single column file containing the list of taxon to work with")
	parser.add_option("-x", "--index", dest="index", help="An integer number identifying which taxon to work on [1-size_of_taxon_list]" , type='int')

	parser.add_option('-b', '--inBestHitFolder', dest='inBestHitFolder', help='Input folder that contains Best Hit files (bh files)')

	parser.add_option('-o', '--outOrthologFolder', dest='outOrthologFolder', help='Folder that will store Orthologous pairs (ort files) in')
	parser.add_option("-l", "--logfile", dest="logfile", help="log file (optional, if not supplied STDERR will be used)")

	parser.add_option('', '--OverwiteOutput', dest='OverwiteOutput', help='If the output file exists, overwrite it. (default=process terminates)', default=False, action="store_true")
	parser.add_option('', '--KeepUnnormalizedScore', dest='KeepUnnormalizedScore', help='Write the un-normalize scores for ortholog pairs (default=False)', default=False, action="store_true")
	

	#
	
	(options, args) = parser.parse_args()


	if len(args) != 0 or not options.taxonlistfile or not options.inBestHitFolder or not options.index or not options.outOrthologFolder:
		parser.error("incorrect arguments.\n\t\tUse -h to get more information or refer to the MANUAL.md")


	#print best_hit

	log('{2} | Orthology | {0} | {1} | {3} | {4} MB | {5}'.format(1 , 'reading taxon list', options.index, '' ,memory_usage_resource(), datetime.now() ))

	taxon_list = readTaxonList(options.taxonlistfile)

	if options.index <= 0 or options.index > len(taxon_list):
		log('{2} | Orthology | {0} | {1} | {3} | {4} MB | {5}'.format('ERROR' , 'Error in index', options.index, '', memory_usage_resource(), datetime.now() ))
		exit()


	taxon1s = taxon_list[options.index - 1]

	# if not options.OverwiteOutput and os.path.exists(os.path.join(options.outputfolder, taxon1s + '.ort.tsv')):
	# 	exit(0)


	log('{2} | Orthology | {0} | {1} | {3} | {4} MB | {5}'.format(2 , 'Reading Best hit (bh file)', options.index, taxon1s , memory_usage_resource(), datetime.now() ))

	taxon1_filename =  os.path.join(options.inBestHitFolder , taxon1s + '.bh.tsv')
	taxon1_dic = readBestHit(taxon1s, taxon1_filename)

	# for taxon2s in sorted(taxon1_dic.keys()):
	# 	print taxon2s
	# 	for query_id in taxon1_dic[taxon2s]:
	# 		print '\t\t', query_id, ':   ',
	# 		for (s_id, ev_exp, ev_mnt) in taxon1_dic[taxon2s][query_id]:
	# 			print s_id, ',',
	# 		print 

	# exit()

	log('{2} | Orthology | {0} | {1} | {3} | {4} MB | {5}'.format(2.1 , 'Size: ' + str(len(taxon1_dic)) + ' taxons.', options.index, taxon1s , memory_usage_resource(), datetime.now() ))

	orthologs = []
	orthologs_index = 0

	taxon2_index = 0

	for taxon2s in sorted(taxon1_dic.keys()):

		if taxon1s < taxon2s:
			taxon2_filename =  os.path.join(options.inBestHitFolder , taxon2s + '.bh.tsv')
		
			taxon2_index += 1
			
			if not os.path.exists(taxon2_filename):
				log('{2} | Orthology | {0} | {1} | {3} | {4} MB | {5}'.format('3.' + str(taxon2_index) , 'ERROR: Targeted Best hit ('+taxon2s+'.bh.tsv) Does not exists', options.index, taxon1s , memory_usage_resource(), datetime.now() ))
				continue			

			taxon2_dic = readBestHit(taxon2s, taxon2_filename)

			log('{2} | Orthology | {0} | {1} | {3} | {4} MB | {5}'.format('3.' + str(taxon2_index) , 'Reading Targeted Best hit ('+taxon2s+'.bh.tsv)', options.index, taxon1s , memory_usage_resource(), datetime.now() ))
			
			taxon1_taxon2_score_sum = 0
			taxon1_taxon2_score_count = 0
			taxon1_taxon2_notfound_count = 0

			for seq1 in taxon1_dic[taxon2s]:

				for (seq2, evalue_exp1, evalue_mant1) in taxon1_dic[taxon2s][seq1]:

					try:
						for (seq0, evalue_exp2, evalue_mant2) in taxon2_dic[taxon1s][seq2]:

							unnormalized_score = None
							if evalue_mant1 < 0.01 or evalue_mant2 < 0.01:
								unnormalized_score = (evalue_exp1 + evalue_exp2) / -2 
							else:
								unnormalized_score = (math.log10(evalue_mant1 * evalue_mant2) + evalue_exp2 + evalue_exp1) / -2

							orthologs += [ [seq1, seq2, unnormalized_score, None] ]

							taxon1_taxon2_score_sum +=	unnormalized_score	
							taxon1_taxon2_score_count += 1
							
							if DEBUG:
								print '\t' , seq1 ,'(',len(taxon1_dic[taxon2s][seq1]),')', '<----->' , seq2 , str(unnormalized_score)

					except KeyError:
						taxon1_taxon2_notfound_count += 1
						if DEBUG:
							print '\t', seq1, '(',len(taxon1_dic[taxon2s][seq1]),')', '<----> NOTHING'
							

			log('{2} | Orthology | {0} | {1} | {3} | {4} MB | {5}'.format('4.' + str(taxon2_index) , 'Orthologs: ' + str( taxon1_taxon2_score_count) , options.index, taxon1s , memory_usage_resource(), datetime.now() ))
			log('{2} | Orthology | {0} | {1} | {3} | {4} MB | {5}'.format('4.' + str(taxon2_index), 'No hits: ' + str( taxon1_taxon2_notfound_count) , options.index, taxon1s , memory_usage_resource(), datetime.now() ))
			
			average = 1
			if taxon1_taxon2_score_count>0:
				average = taxon1_taxon2_score_sum/taxon1_taxon2_score_count

			if DEBUG:
				print  taxon1s,'<->',taxon2s, 'Count:',taxon1_taxon2_score_count ,'Average:' ,average
				print

			
			for i in xrange(taxon1_taxon2_score_count):
				if average>0:
					orthologs[orthologs_index+i][3] = orthologs[orthologs_index+i][2] / average 
				else:
					orthologs[orthologs_index+i][3] = 1

			orthologs_index += taxon1_taxon2_score_count

	out_f = open (os.path.join(options.outOrthologFolder , taxon1s + '.ort.tsv'), 'w')
	if options.KeepUnnormalizedScore:
		out_f.write('query_id\tsubject_id\tunnormalized_score\tnormalized_score\n')
	for ortholog in orthologs:
		out_f.write(ortholog[0] + '\t')
		out_f.write(ortholog[1] + '\t')
		if options.KeepUnnormalizedScore:
			out_f.write(str(ortholog[2]) +'\t')
		else:
			ortholog[3] = int( ortholog[3] * 1000 + .5) / 1000.0 
		out_f.write(str(ortholog[3]) + '\n')

	out_f.close()
	log('{2} | Orthology | {0} | {1} | {3} | {4} MB | {5}'.format(5, 'Finished' , options.index, taxon1s , memory_usage_resource(), datetime.now() ))
