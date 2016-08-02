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

All you need to do is to put all you protein fasta files inside a specificly named folder `0.input_faa` which is inside a container folder, just like the sample folder provided in this repository.

If you have a `container/0.input_faa` set up with all your fasta files inside the `0.input_faa`, all you need to do is to run:

```bash
orthomcl.sh container
```

