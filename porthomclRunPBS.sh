####################################################
##      				 OrthoMCLP
## -------------------------------------------------
## This is a sample PBS Script to run BLASTs for 
## multiple input files against a single database
## 
## You need to update all the variable in enclosed 
## in [] for your PBS script
##
## Once this script is complete, you must be able to 
## submit the job using 
##		qsub orthomclpRunBlasts.sh
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
#PBS -t 1-4

## Job Name: give it a good name so that you recognize it among the jobs running in the cluster.
## example: #PBS -N 3_orthomclp_run_blast
##
#PBS -N [JOB_NAME]

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
## example: ROOTFOLDER=/scratch/user/orthomclp/data/
##
ROOTFOLDER=[FOLDER_CONTAINING_INPUT_AND_DB_AND_OUTPUT]

SHORT_JOBID=`echo $PBS_JOBID |cut -d. -f1`
1NLOGO_JOBNAME=`echo $PBS_JOBNAME |cut -d- -f1`
exec 1>$PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.out  2>$PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.err


cd $PBS_O_WORKDIR

# Create a file to mark the start of this BLAST
touch $PBS_JOBNAME-$SHORT_JOBID.start


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
# orthomclpPairsBestHit.py -t sample/taxon_list -i sample/5.input/ -b sample/5.besthit -q sample/5.bestquerytaxon  -l sample/orthomclp.log -x $PBS_ARRAYID


# Create a file to mark the end of this BLAST
touch $PBS_O_WORKDIR/$PBS_JOBNAME-$SHORT_JOBID.end
