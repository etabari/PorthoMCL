# PorthoMCL
Parallel implementation of OrthoMCL


We have reimplemented the sections of OrthoMCL that rely on databases. This way, OrthoMCL could be ran in parallel for a large number of genomes.

Ehsan Tabari, May 19, 2015

The initial steps of PorthoMCL are exactly the same as OrthoMCL.
We have reimplemented steps 9 and 10. (main processing steps)

# Requirements

There are very few requirements for PorthoMCL. Here are the list of the things needed to run PorthoMCL

- Perl 5.x
- Python 2.x
- [BioPython](http://biopython.org/wiki/Download)
- [NCBI Blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [MCL](http://www.micans.org/mcl/sec_description1.html)

Perl and Python come preinstalled on most Linux and Unix (OS X). You need to install them on Windows. 

This implementation *removes* the need for a database server.

# Steps 

The sample folder contains execution of all these steps. Each folder created at each step is numbered by the step number. 
The starting data is located in sample/0.input_faa

## Step 1: Prepare the input sequences.

This is EXACTLY like the 5th step of OrhtoMCL. In this step fasta file will be generated to have sequence and protein identifiers.
 
You must run `orthomclAdjustFasta` on your input fasta files. orthomclAdjustFasta creates output in the same folder it's ran. 

The input arguments to orthomclAdjustFasta are:
- String to label each sequence
- Input Fasta file 
- A number identifying the field containing protein identification on the fasta header

```
orthomclAdjustFasta NC_000913 sample/0.input_faa/NC_000913.faa 4
```

If you have downloaded faa files from NCBI (for example: ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/all.faa.tar.gz)
We have supplied a bash script to run `orthomclAdjustFasta` on all faa files and produce the proper fasta files.

```shell
cd compliantFasta
orthomclAdjustFastaAll.sh <input_folder>
```

In the sample run, the results of this step is copied to samlple/1.compliantFasta

## Step 2: Filter the input 

This is EXACTLY like the 6th step of OrhtoMCL. 
This step produces a single goodProteins.fasta file to run BLAST on.  
It filters away poor-quality sequences (placing them in poorProteins.fasta).  
The filter is based on length and percent stop codons.

The input arguments to `orthomclFilterFasta` are:
- input_dir:               the directory containing a set of .fasta files
- min_length:              minimum allowed length of proteins.  (suggested: 10)
- max_percent_stop:        maximum percent stop codons.  (suggested 20)
- good_proteins_file:      _optional_  By default goodProteins.fasta in the current dir.
- poor_proteins_file:      _optional_  By default poorProteins.fasta in the current dir.

```shell
orthomclFilterFasta samlple/1.compliantFasta 10 20 
```

In the sample run, the results of this step is copied to samlple/2.filteredFasta



## Step 3: All-v-all BLAST

We have developed and tested this with NCBI BLAST 2.2.29+. 

At this stage, instead of running one huge blast that will never ends,
we will split the input files to smaller size and blast them against the same huge database. 
This way, you will have many blast runs that can be ran in parallel on a computing cluster. 
We have supplied a torque script for PBS based clusters. To do that we need to split the input file first.


#### 3.1 Create BLAST database

```
makeblastdb -in samlple/2.filteredFasta/goodProteins.fasta  -dbtype prot
```

The output files of the makeblastdb is copied to samlple/3.blastdb

#### 3.2 Split the input file 

Split goodProteins.fasta using porthomclSplitFasta.py so that each file has relatively a small size for each blast.
10K sequences in each query batch makes each blast run finish faster. 


The input arguments to `porthomclSplitFasta.py` are:
- input-fasta: The large fasta file to be split into mnultiple files 
- size: number of sequences to keep in each split file


```shell
porthomclSplitFasta.py -i samlple/2.filteredFasta/goodProteins.fasta  -s 10000
```

The output files of the makeblastdb is copied to samlple/3.blastquery

#### 3.3 Run blasts 

In this step the required blast arguments are as follows:
- query 3.blastquery/goodProteins.fasta  (BLAST Query)
- db 3.blastdb/goodProteins.fasta  (BLAST Database)
- seg yes  (Filter query sequence with [SEG](http://www.ncbi.nlm.nih.gov/pubmed/8743706) to remove low complexity regions)
- dbsize 100000000  (*constant* database size for adding other genomes and keeping e-values comparable)
- evalue 1e-5  (minimum e-value to report back)
- outfmt 6 (create tabluar output files)
- num_threads 8 (run it multi-threaded, only if you have a capable machine to run it on)
- out blastres/blastres.tab (tabulat output for the input)


Sample run for one of the input files:
```
cd sample
mkdir 3.blastres 
blastp -query 3.blastquery/goodProteins.fasta.1  -db 3.blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out 3.blastres/blastres.1.tab
 ```


#### 3.3.1 PBS Cluster BLAST

If you have access to a PBS base computing cluster, we have included a torque script to be ran on a given computing cluster. 
It initiates an array of jobs (for the sample 1 to 4) that each job array will execute one BLASTP commands on the input indetified by 
job array index. 

For example, job 4 will run a multithreaded blastp command with query goodProteins.fasta.4 and outputs blastres.4.tab

The script is `porthomclRunPBS.sh`

You have to customize the `#PBS` variables so that it matches the requirements of your cluster. 


#### 3.4 Put all blastp results together

You have parallelized the blast so far, you need to put the results back together

```shell
cd sample
cat 3.blastres/* > 3.blastresmerge/blastres.tab
```

_It's better to run this as a job on the cluster rather than running it interactively._


## Step 4: Parse BLAST results
This is a little different from the 8th step of OrhtoMCL. 
This step parses NCBI BLAST tabular output format into the format that can be loaded into the orthomcl database. 
It additionally prints out the minimum value of evalue_exp (A negative number) in STDERR.

This will serve as an input to the next step. 

```
cd sample
porthomclBlastParser 3.blastresmerge/blastres.tab 1.compliantFasta >> 4.parsedblast/similarSequences.txt
```

_It's better to run this as a job on the cluster rather than running it interactively._


## Step 5: Find Pair

This is the most time consuming step of OrthoMCL that would not finish using a MySQL server for large amount of data.
This is an embarrassingly parallel problem and can be easily splitted into multiple input files to be executed in parallel.

So the first step is to split the similarSequences.txt file into multiple files:

#### 5.0.1 Split input file

`awk` provides a lot of capabilities for such use, we remove the taxon columns from the input file.  

```shell
cd sample
mkdir 5.splitSimSeq
awk -F'[|\t]' '{print $1"|"$2"\t"$3"|"$4"\t"$7"\t"$8"\t"$9"\t"$10 >> ("5.splitSimSeq/"$1".ss.tsv")}' 4.parsedblast/similarSequences.txt
```

The files creared in this step from the sample run is in sample/5.splitSimSeq.

_It's better to run this as a job on the cluster rather than running it interactively._


#### 5.0.2 Taxon List file

Although we could have listed the compliantFasta folder like the original orthomcl, but we find it very inefficient. 
Also, having a list of taxons will help to easily exclude a taxon from the whole process.

This is achived by listing the folder and removing the extension, sending it to a file.

```shell
cd sample
ls -1 1.compliantFasta/ | sed -e 's/\..*$//'  > taxon_list
```


#### 5.1 Finding Best Hits

This step must be done for all the split input files before you could move forward to calculate Orthologs, inParalogs and coOrthologs.
OrthoMCL has a bug. read more about it [here](release_notes.md#1-orthomcl-bug): 


In this step actually the paralogs are found, and an unnormalized score is assigned to them.
Step 5.3 will normalize this score so that it be comparable among different genomes. 


The input parameters are: 

- Inputs:
	- **-t** | (--taxonlist)  | A single column file containing the list of taxon to work with
	- **-x** | (--index)	  |	an integer number identifying which taxon to work on [1-size_of_taxon_list]
	- **-s** | (--inSimSeq)	  | Input folder that stores *split similar sequences* files (ss files)
- Outputs:
	-  -b  (--outBestHitFolder)		folder that will stores Best Hit files (If not set, current folder)
	-  -q  (--outInParalogTempFolde)	folder to generate best InParalogTemp evalue scores (pt files) (required only for Paralogs)
	-  -l  (--logfile)	log file
- Options:
	-  --evalueExponentCutoff evalue Exponent Cutoff (a negative value, default=-5)
	-  --percentMatchCutoff percent Match Cutoff (integer value, default=50)
	-  --cacheInputFile      Cache input file or read it again. (Only use if I/O is very slow)
	-  **--keepOrthoMCLBug**     Keep the OrthoMCL bug in creating BetterHit files (br) where self hits are included

This is the example to run for the FIRST sample file (-x 1). To perform this for all the sample files, you must to run the same command for **-x 1** to **-x 12**. 


```shell
mkdir sample/5.paralogTemp
mkdir sample/5.besthit

porthomclPairsBestHit.py -t sample/taxon_list -s sample/5.splitSimSeq -b sample/5.besthit -p sample/5.paralogTemp -x 1

```

You can run this code in parallel with different values for -x. An example of such execution is included in the `porthomclRunPBS.sh` script.


#### 5.2 Finding Orthologs

The output of this step is all the ortholog genes. 


The input parameters are:
- **-t** (--taxonlist) A single column file containing the list of taxon to work with
- **-x**  (--index) an integer number identifying which taxon to work on [1-size_of_taxon_list]
- **-b**  (--inBestHitFolder) folder that stores Best Hit files
- **-o** (--outOrthologFolder) folder that will stores orthologs (TaxonID.ort.tsv files)
- --KeepUnnormalizedScore  Write the un-normalize scores also. (default=False)
This is the example to run for the FIRST sample file (-x 1). To perform this for all the sample files, you must to run the same command for -x 1 to -x 12. 

```shell
mkdir sample/5.orthologs

porthomclPairsOrthologs.py -t sample/taxon_list -b sample/5.besthit -o sample/5.orthologs -x 1

```
You can run this code in parallel with different values for -x. An example of such execution is included in the `porthomclRunPBS.sh` script.


#### 5.3 Finding Paralogs

This step normalizes the scores calculated in step 5.1. OrthoMCL uses the average paralog score of the genes that have orthologous relationships to 
normalized the score. 
Therefore a list of all the genes that have an orthology is essential. 

To find the paralogs, we need to extract every gene that has an orthologous relationship in every genome. 
This can be achieved using a set of bash commands as follows. 


```shell
cd sample
mkdir 5.ogenes

# genes in the second column
awk -F'[|\t]' '{print $4 >> ("5.ogenes/"$3".og.tsv")}' 5.orthologs/*.ort.tsv

# genes in the first column
awk -F'[|\t]' '{print $2 >> ("5.ogenes/"$1".og.tsv")}' 5.orthologs/*.ort.tsv

# keep unique ones
find 5.ogenes/ -maxdepth 1 -type f -exec sort -u -o {} {} \;
```

Parallelized version of preceding commands have been included in the `porthomclRunPBS.sh` for convenience.

Having these files (og files) in hand we can normalize the inparalog relationships:

The output of this step is all the paralogs genes with normalized scores. 


The input parameters are:
- **-t** (--taxonlist) A single column file containing the list of taxon to work with
- **-x**  (--index) an integer number identifying which taxon to work on [1-size_of_taxon_list]

- **-q**  (--inInParalogTempFolder) folder that stores Temp Paralogs with un-normalize scores (pt files).
- **-o**  (--inOrthologGeneFolder) folder that stores list of genes with orthology relationships (og files).

- **-p** (--outInParalogFolder) folder that will stores InParalogs (par files)
- --KeepUnnormalizedScore  Write the un-normalize scores also. (default=False)

This is the example to run for the FIRST sample file (-x 1). To perform this for all the sample files, you must to run the same command for -x 1 to -x 12. 

```shell
mkdir sample/5.paralogs

porthomclPairsInParalogs.py -t sample/taxon_list -q sample/5.paralogTemp -o sample/5.ogenes -p sample/5.paralogs -x 1

```
You can run this code in parallel with different values for -x. An example of such execution is included in the `porthomclRunPBS.sh` script.


.

.

.

.

THIS IS TO COTINUE
