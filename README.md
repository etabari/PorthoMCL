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

The detail manual is [here](MANUAL.md):


# Steps 

## Step 1: Prepare the input sequences.

## Step 2: Filter the input 

## Step 3: All-v-all BLAST

#### 3.1 Create BLAST database

#### 3.2 Split the input file 

#### 3.3 Run blasts 

#### 3.3.1 PBS Cluster BLAST

#### 3.4 Put all blastp results together

## Step 4: Parse BLAST results

## Step 5: Find Pair

#### 5.0.1 Split input file

#### 5.0.2 Taxon List file

#### 5.1 Finding Best Hits

#### 5.2 Finding Orthologs

#### 5.3 Finding Paralogs

#### 5.4 Finding CoOrthologs