####################################################
##      				 PorthoMCL
## -------------------------------------------------
## This is a sample PBS Script to run BLASTs for 
## multiple input files against a single database
## 
## You need to update all the variable in enclosed 
## in [] for your PBS script
##
## Once this script is complete, you must be able to 
## submit the job using 
##		qsub porthomclRunPBS.sh
## 
## Ehsan Tabari
#################################################### 

## #PBS lines are not comments. PBS uses these as commands. 
## read more about it on openpbs.org 
##

## This is the Queue name variable. Where to submit the job to? 
## 
#PBS -q [QUEUE_NAME]

## This will indicate the indices of the input to run. 1-100 maybe?
##
#PBS -t 1-12

## Job Name: give it a good name so that you recognize it among the jobs running in the cluster.
## example: #PBS -N 3_porthomcl_run_blast
##
#PBS -N x_porthomcl

## Your requirements: I prefer to run each BLAST multithreaded as well. (-num_threads 8) so I ask for 8 cores.
##
#PBS -l nodes=1:ppn=8
#PBS -l walltime=100:00:00

## redirect errors and outputs
## 
#PBS -e /dev/null
#PBS -o /dev/null


## This is the folder where you keep your data 
## This is going to be big.
## Don't forget the / at the end of it
## example: ROOTFOLDER=/scratch/user/porthomcl/data
##
ROOTFOLDER=sample

SHORT_JOBID=`echo $PBS_JOBID |cut -d. -f1`
1NLOGO_JOBNAME=`echo $PBS_JOBNAME |cut -d- -f1`
exec 1>$PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.out  2>$PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.err


cd $PBS_O_WORKDIR

TAXON_FILE="`sed -n $PBS_ARRAYID\p $ROOTFOLDER/taxon_list`"


# Create a file to mark the start of this BLAST
echo "$PBS_ARRAYID|$TAXON_FILE|$(date)" > $PBS_JOBNAME-$SHORT_JOBID.str


############################################
### STEP 3: run the blast
###
### https://github.com/etabari/OrthoMCLP#step-4-parse-blast-results
### note that -num_threads 8 is passed, adjust it accordingly to your nodes.
###
#  blastp -query $ROOTFOLDER/3.blastquery/goodProteins.fasta.$PBS_ARRAYID  -db $ROOTFOLDER/3.blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out $ROOTFOLDER/3.blastres/blastres.$PBS_ARRAYID.tab

############################################
### STEP 5: Finding Best Hits
###
### https://github.com/etabari/OrthoMCLP#finding-best-hits
###
# porthomclPairsBestHit.py -t $ROOTFOLDER/taxon_list -i $ROOTFOLDER/5.input/ -b $ROOTFOLDER/5.besthit -q $ROOTFOLDER/5.bestquerytaxon  -l $PBS_JOBNAME-$SHORT_JOBID.log -x $PBS_ARRAYID


############################################
### STEP 5.0.1: Split the input file
###
### https://github.com/etabari/PorthoMCL#finding-orthologs
###
# awk -F'[|\t]' '{print $1"|"$2"\t"$3"|"$4"\t"$7"\t"$8"\t"$9"\t"$10 >> ("$ROOTFOLDER/5.splitSimSeq/"$1".ss.tsv")}' $ROOTFOLDER/4.parsedblast/similarSequences.txt



############################################
### STEP 5.1: Finding Orthologs
###
### https://github.com/etabari/PorthoMCL#finding-orthologs
###
# porthomclPairsOrthologs.py -t $ROOTFOLDER/taxon_list -b $ROOTFOLDER/5.besthit -o $ROOTFOLDER/5.orthologs -l $PBS_JOBNAME-$SHORT_JOBID.log -x $PBS_ARRAYID


############################################
### STEP 5.2: Finding Paralogs
###
### https://github.com/etabari/PorthoMCL#finding-paralogs
###
#########
###
### Step 5.2.1: Find the genes:
###
# awk -F'[|\t]' '{print $4 >> ("$ROOTFOLDER/5.ogenes/"$3".og.tsv")}' $ROOTFOLDER/5.orthologs/$TAXON_FILE.ort.tsv 
#
# awk -F'[|\t]' '{print $2 >> ("$ROOTFOLDER/5.ogenes/"$1".og.tsv")}' $ROOTFOLDER/5.orthologs/$TAXON_FILE.ort.tsv 
#


#########
###
### Step 5.2.2: sort unique files:
###
# sort -u $ROOTFOLDER/5.ogenes/$TAXON_FILE.og.tsv > $ROOTFOLDER/5.ogenes/$TAXON_FILE.og.tsv
#
#########
###
### Step 5.2.3: Run the paralog
###



# Create a file to mark the end of this BLAST
echo "$PBS_ARRAYID|$TAXON_FILE|$(date)" >  $PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.fin

