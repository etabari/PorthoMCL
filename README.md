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

An application of PorthoMCL on 2,758 bacterial genome is available at [here](http://bioinfo.uncc.edu/ehsan.tabari/porthomcl/).

# Quick and Easy Run

It's possible to run porthomcl without reading the manual or tweeking the parameters. To this end we have provided `porthomcl.sh`.

All you need to do is to put all your protein fasta files inside a specificly named folder `0.input_faa` which is inside a container folder, just like the sample folder provided in this repository.

If you have a `container/0.input_faa` set up with all your fasta files inside the `0.input_faa`, all you need to do is to run:

```bash
orthomcl.sh container
```
## Requirements and options

you have to install all the requirements listed above.

The syntaxt to run PorthoMCL is as follows: 	
```bash
orthomcl.sh [OPTIONS] container_folder
```
Options can be: 


|  Options|                			|    Description                                                              
|--------|--------------------------|-----------------------------------------------------------------------------
|   -h   | --help 	 			    | Prints this help.              
|   -t   | --num_threads	  		| An integer number identifying the number of processes/threads to be used (default=4)
|   -s   | --startat				| Step to start at. in cace you stopped after a specific step, you can restart by runing from a step (values: 1-8. default=1)         
|   -l 	| --lib  					| The location of PorthoMCL files (you can alternatively add that location to your $PATH)
|    | --wait	| wait for a key press after each step

Example:

Run porthoMCL for the provided sample, from the repository folder

```bash
./orthomcl.sh -t 8 -l . sample
```
