#!/bin/bash

################################################################################################
################################################################################################
################################################################################################
#######    THIS IS A STAND ALONE PORTHOMCL Executer. 
#######    ----------------------------------------
#######    All you need is a set of FAA files (amino acid fasta files)
#######          one file per chromosome/genome
#######          We use the file name as chromosome/genome identifier
#######
#######	   Put them all in a folder inside a CONTAINER folder 
#######         The container folder is where we produce other temp files or results files
#######         Put your files in a folder named 0.input_faa inside the container folder
#######     
#######    Run this file
#######       
#######         porthomcl.sh CONTAINER
#######
#######    for other options check porthomcl.sh --help
#######
################################################################################################
####### Ehsan Tabari (https://github.com/etabari/PorthoMCL)
################################################################################################
################################################################################################

NUM_THREADS=4
START=1
END=8

while [[ $# -gt 1 ]]
do
	key="$1"

	case $key in
		-h|--help)
	    SHOWHELP=YES
	    # shift # past argument
	    ;;
	    -t|--num_threads)
	    NUM_THREADS="$2"
	    shift # past argument
	    ;;
	    -s|--startat)
	    START="$2"
	    shift # past argument
	    ;;
	    -e|--endafter)
	    END="$2"
	    shift # past argument
	    ;;
	    -l|--lib)
	    LIBPATH=`cd "$2"; pwd`
	    shift # past argument
	    ;;
	    --wait)
	    WAITFORKEY=YES
	    ;;

	    *)
	            # unknown option
	    ;;
	esac
	shift # past argument or value
done


container_org="$1"
container=`cd "$1"; pwd`
LOGFILE="$container/porthomcl.log"

echo 'Container folder:' $container


if [ -z "$SHOWHELP" ] && ! [[ -d "$container/0.input_faa" ]]
then
	echo "ERROR: CONTAINER_FOLDER does not contain 1.compliantFasta"
	SHOWHELP="YES"
fi

if [ -z "$LIBPATH" ] 
then
	if (! hash orthomclAdjustFasta 2>/dev/null)
	then
		echo "ERROR: provided PorthoMCL does not exists or PorthoMCL was not found in \$PATH"
		SHOWHELP="YES"
	fi
else 
	LIBPATH="$LIBPATH"/
	if [ ! -f $LIBPATH/orthomclAdjustFasta ]
	then
		echo "ERROR: PorthoMCL was not found in $LIBPATH"
		SHOWHELP="YES"
	fi
fi


if (! hash makeblastdb 2>/dev/null)
then
	echo "ERROR: blast was not found in \$PATH"
	SHOWHELP="YES"
fi

if (! hash mcl 2>/dev/null)
then
	echo "ERROR: mcl was not found in \$PATH"
	SHOWHELP="YES"
fi

if [ "$NUM_THREADS" -lt 1 ] 
then 
	NUM_THREADS=1
fi

if [ ! -z "$SHOWHELP" ]
then
	echo "usage:"
	echo "   porthomcl.sh [OPTIONS] CONTAINER_FOLDER"
	echo 
	echo "   options"
	echo "     -l, --lib          path to PorthoMCL, not required if it's installed in the \$PATH"
	echo "     -t, --num_threads  number of threads/processes to run (default=4)"
	echo "     -s, --startat      an integer indicating the first Step to run (default=1)"
	echo "     -e, --endafter     an integer indicating the last Step to run (default=8)"
	echo "     --wait             wait for keypress at every step"
	echo 
	echo "   The CONTAINER_FOLDER must contain a folder named 0.input_faa"
    echo "   please put your fasta (amino acid fasta) files in container/0.input_faa"
    echo "   be aware that we will process any file in that folder. "
	exit
fi


# #####################################
# ### STEP 1: Prepare the input sequences.

if [ $START -le 1 ]; then
	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi
	
	echo "$1|1|orthomclAdjustFasta|Start|$(date)" 
	echo "$1|1|orthomclAdjustFasta|Start|$(date)" >> $LOGFILE

	mkdir $container/1.compliantFasta
	cd $container/1.compliantFasta
	for faa in $container/0.input_faa/*
	do
		

		strand=$(basename "$faa")
		strand="${strand%.*}" # remove .faa
		echo '    |processing|'$strand

		`"$LIBPATH"orthomclAdjustFasta $strand $faa 4 || exit 1`
	done
	echo "$1|1|orthomclAdjustFasta|End|$(date)"
	echo "$1|1|orthomclAdjustFasta|End|$(date)" >> $LOGFILE


	# #####################################
	# ### STEP 1.1: Taxon List file


	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi
	
	echo "$1|1.1|Create taxon list file|Start|$(date)"
	echo "$1|1.1|Create taxon list file|Start|$(date)" >> $LOGFILE

	cd $container
	ls -1 1.compliantFasta/ | sed -e 's/\..*$//'  > taxon_list

	echo "$1|1.1|Create taxon list file|End|$(date)"
	echo "$1|1.1|Create taxon list file|End|$(date)" >> $LOGFILE
fi

if [ $END -le 1 ]; then
	exit
fi
# #####################################
# ### STEP 2: Filter the input
if [ $START -le 2 ]; then
	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi
	
	echo "$1|2|orthomclFilterFasta|Start|$(date)" 
	echo "$1|2|orthomclFilterFasta|Start|$(date)" >> $LOGFILE

	`"$LIBPATH"orthomclFilterFasta $container/1.compliantFasta 10 20`
	mkdir $container/2.filteredFasta
	mv goodProteins.fasta $container/2.filteredFasta/
	mv poorProteins.fasta $container/2.filteredFasta/ 

	echo "$1|2|orthomclFilterFasta|End|$(date)" 
	echo "$1|2|orthomclFilterFasta|End|$(date)" >> $LOGFILE
fi
if [ $END -le 2 ]; then
	exit
fi
# #####################################
# ### STEP 3.1:  Create BLAST database
if [ $START -le 3 ]; then
	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi
	
	echo "$1|3.1|makeblastdb|Start|$(date)"
	echo "$1|3.1|makeblastdb|Start|$(date)" >> $LOGFILE

	makeblastdb -in $container/2.filteredFasta/goodProteins.fasta  -dbtype prot

	mkdir $container/3.blastdb

	mv $container/2.filteredFasta/goodProteins.fasta.* $container/3.blastdb/

	echo "$1|3.1|makeblastdb|End|$(date)" 
	echo "$1|3.1|makeblastdb|End|$(date)" >> $LOGFILE


	# #####################################
	# ### 3.2 Split the input file
	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi
	

	echo "$1|3.2|porthomclSplitFasta|Start|$(date)"
	echo "$1|3.2|porthomclSplitFasta|Start|$(date)" >> $LOGFILE

	mkdir $container/3.blastquery 

	`"$LIBPATH"porthomclSplitFasta.py -i $container/2.filteredFasta/goodProteins.fasta  -o $container/3.blastquery`

	echo "$1|3.2|porthomclSplitFasta|End|$(date)" 
	echo "$1|3.2|porthomclSplitFasta|End|$(date)" >> $LOGFILE

	# #####################################
	### 3.3 Run blasts
	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi
	

	echo "$1|3.3|Blasts|Start|$(date)" 
	echo "$1|3.3|Blasts|Start|$(date)" >> $LOGFILE

	cd $container
	mkdir 3.blastres 
	for query in 3.blastquery/*
	do
		echo "   $query"
	    strand=$(basename "$query")
	    strand="${strand%.*}" # remove .fasta
	    blastp -query 3.blastquery/$strand.fasta  -db 3.blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads $NUM_THREADS -out 3.blastres/$strand.tab
	done

	echo "$1|3.3|Blasts|End|$(date)" 
	echo "$1|3.3|Blasts|End|$(date)" >> $LOGFILE
fi
if [ $END -le 3 ]; then
	exit
fi
# #####################################
# ### 4 Parse blasts
if [ $START -le 4 ]; then
	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi
	
	echo "$1|4|Parse Blast Results|Start|$(date)"
	echo "$1|4|Parse Blast Results|Start|$(date)" >> $LOGFILE

	mkdir $container/4.splitSimSeq

	cd $container
	for blres in $container/3.blastres/*
	do
	    strand=$(basename "$blres")
	    strand="${strand%.*}" # remove .fasta
	    `"$LIBPATH"porthomclBlastParser 3.blastres/$strand.tab 1.compliantFasta >> 4.splitSimSeq/$strand.ss.tsv`
	done

	echo "$1|4|Parse Blast Results|End|$(date)" 
	echo "$1|4|Parse Blast Results|End|$(date)" >> $LOGFILE
fi
if [ $END -le 4 ]; then
	exit
fi
# #####################################
# ### 5 Finding Best Hits
if [ $START -le 5 ]; then
	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi
	
	echo "$1|5|porthomclPairsBestHit|Start|$(date)" 
	echo "$1|5|porthomclPairsBestHit|Start|$(date)" >> $LOGFILE

	cd $container
	taxon_count=$(sed -n '$=' taxon_list)

	mkdir $container/5.paralogTemp
	mkdir $container/5.besthit


	loop_start=1
	loop_end=$((loop_start + NUM_THREADS - 1))
	loop_end=$(($loop_end>$taxon_count?$taxon_count:$loop_end))

	while [ $loop_start -lt $taxon_count ]
	do
		for i in $(seq $loop_start $loop_end)
		do
			`"$LIBPATH"porthomclPairsBestHit.py -t $container/taxon_list -s $container/4.splitSimSeq -b $container/5.besthit -q $container/5.paralogTemp -x $i  -l $LOGFILE &`
		done
		wait
		loop_start=$((loop_end + 1))
		loop_end=$((loop_start + NUM_THREADS - 1))
		loop_end=($(($loop_end>$taxon_count?$taxon_count:$loop_end)))
		
	done

	echo "$1|5|porthomclPairsBestHit|End|$(date)" 
	echo "$1|5|porthomclPairsBestHit|End|$(date)" >> $LOGFILE
fi
if [ $END -le 5 ]; then
	exit
fi
# #####################################
# ### Step 6: Finding Orthologs
if [ $START -le 6 ]; then

	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi
	

	echo "$1|6|Finding orthologs (porthomclPairsOrthologs)|Start|$(date)"
	echo "$1|6|Finding orthologs (porthomclPairsOrthologs)|Start|$(date)" >> $LOGFILE

	cd $container
	taxon_count=$(sed -n '$=' taxon_list)

	mkdir $container/6.orthologs

	loop_start=1
	loop_end=$((loop_start + NUM_THREADS - 1))
	loop_end=$(($loop_end>$taxon_count?$taxon_count:$loop_end))

	while [ $loop_start -lt $taxon_count ]
	do
		for i in $(seq $loop_start $loop_end)
		do
			`"$LIBPATH"porthomclPairsOrthologs.py -t $container/taxon_list -b $container/5.besthit -o $container/6.orthologs -x $i  -l $LOGFILE &`
		done
		wait
		loop_start=$((loop_end + 1))
		loop_end=$((loop_start + NUM_THREADS - 1))
		loop_end=($(($loop_end>$taxon_count?$taxon_count:$loop_end)))
		
	done

	echo "$1|6|Finding orthologs (porthomclPairsOrthologs)|End|$(date)" 
	echo "$1|6|Finding orthologs (porthomclPairsOrthologs)|End|$(date)" >> $LOGFILE
fi
if [ $END -le 6 ]; then
	exit
fi
# #####################################
# ### Step 7:  Finding Paralogs
if [ $START -le 7 ]; then
	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi


	echo "$1|7.1|Finding Paralogs (awk)|Start|$(date)" 
	echo "$1|7.1|Finding Paralogs (awk)|Start|$(date)" >> $LOGFILE

	cd $container
	taxon_count=$(sed -n '$=' taxon_list)

	mkdir 7.ogenes

	# # genes in the second column
# 
	loop_start=1
	loop_end=$((loop_start + NUM_THREADS - 1))
	loop_end=$(($loop_end>$taxon_count?$taxon_count:$loop_end))

	while [ $loop_start -lt $taxon_count ]
	do
		echo $loop_start $loop_end
		for i in $(seq $loop_start $loop_end)
		do
			{ TAXON_FILE="`sed -n $i\p taxon_list`"; awk -F'[|\t]' '{print $4 >> ("7.ogenes/"$3".og.tsv")}' 6.orthologs/$TAXON_FILE.ort.tsv;}&
		done
		wait
		loop_start=$((loop_end + 1))
		loop_end=$((loop_start + NUM_THREADS - 1))
		loop_end=($(($loop_end>$taxon_count?$taxon_count:$loop_end)))
		
	done


	# genes in the first column
	awk -F'[|\t]' '{print $2 >> ("7.ogenes/"$1".og.tsv")}' 6.orthologs/*.ort.tsv

	loop_start=1
	loop_end=$((loop_start + NUM_THREADS - 1))
	loop_end=$(($loop_end>$taxon_count?$taxon_count:$loop_end))

	while [ $loop_start -lt $taxon_count ]
	do
		for i in $(seq $loop_start $loop_end)
		do
			{ TAXON_FILE="`sed -n $i\p taxon_list`";awk -F'[|\t]' '{print $2 >> ("7.ogenes/"$1".og.tsv")}' 6.orthologs/$TAXON_FILE.ort.tsv;} & 
		done
		wait
		loop_start=$((loop_end + 1))
		loop_end=$((loop_start + NUM_THREADS - 1))
		loop_end=($(($loop_end>$taxon_count?$taxon_count:$loop_end)))
		
	done

	echo "$1|7.1|Finding Paralogs (awk)|End|$(date)" 
	echo "$1|7.1|Finding Paralogs (awk)|End|$(date)" >> $LOGFILE
	###########
	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi


	echo "$1|7.2|Finding Paralogs (porthomclPairsOrthologs)|Start|$(date)" 
	echo "$1|7.2|Finding Paralogs (porthomclPairsOrthologs)|Start|$(date)" >> $LOGFILE

	mkdir $container/7.paralogs

	loop_start=1
	loop_end=$((loop_start + NUM_THREADS - 1))
	loop_end=$(($loop_end>$taxon_count?$taxon_count:$loop_end))

	while [ $loop_start -lt $taxon_count ]
	do
		for i in $(seq $loop_start $loop_end)
		do
			`"$LIBPATH"porthomclPairsInParalogs.py -t $container/taxon_list -q $container/5.paralogTemp -o $container/7.ogenes -p $container/7.paralogs -x $i -l $LOGFILE &`
		done
		wait
		loop_start=$((loop_end + 1))
		loop_end=$((loop_start + NUM_THREADS - 1))
		loop_end=($(($loop_end>$taxon_count?$taxon_count:$loop_end)))
		
	done


	echo "$1|7.2|Finding Paralogs (porthomclPairsOrthologs)|Start|$(date)" 
	echo "$1|7.2|Finding Paralogs (porthomclPairsOrthologs)|Start|$(date)" >> $LOGFILE

fi
if [ $END -le 7 ]; then
	exit
fi
# #####################################
# ### Step 8:  Running MCL
if [ $START -le 8 ]; then

	if [ ! -z "$WAITFORKEY" ]; then
		read -p "Press any key to continue... " -n1 -s
		echo
	fi
	

	echo "$1|8|Running MCL|Start|$(date)" 
	echo "$1|8|Running MCL|Start|$(date)" >> $LOGFILE


	cat $container/6.orthologs/*.tsv >> $container/8.all.ort.tsv

	mcl $container/8.all.ort.tsv  --abc -I 1.5 -t $NUM_THREADS -o $container/8.all.ort.group

	cat $container/7.paralogs/*.tsv >> $container/8.all.par.tsv

	mcl $container/8.all.par.tsv  --abc -I 1.5 -t $NUM_THREADS -o $container/8.all.par.group


	echo "$1|8|Running MCL|End|$(date)" 
	echo "$1|8|Running MCL|End|$(date)" >> $LOGFILE
fi