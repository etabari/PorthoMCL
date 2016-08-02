#!/bin/bash


NUM_THREADS=4


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

	    -l|--lib)
	    LIBPATH=`cd "$2"; pwd`
	    shift # past argument
	    ;;
	    # --default)
	    # DEFAULT=YES
	    # ;;
	    *)
	            # unknown option
	    ;;
	esac
	shift # past argument or value
done

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
	echo 
	echo "   The container folder must contain a folder named 0.input_faa"
    echo "   please put your fasta (amino acid fasta) files in container/0.input_faa"
    echo "   be aware that we will process any file in that folder. "
	exit
fi


# #####################################
# ### STEP 1: Prepare the input sequences.
read -p "Press any key to continue... " -n1 -s

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
read -p "Press any key to continue... " -n1 -s

echo "$1|1.1|Create taxon list file|Start|$(date)"
echo "$1|1.1|Create taxon list file|Start|$(date)" >> $LOGFILE

cd $container
ls -1 1.compliantFasta/ | sed -e 's/\..*$//'  > taxon_list

echo "$1|1.1|Create taxon list file|End|$(date)"
echo "$1|1.1|Create taxon list file|End|$(date)" >> $LOGFILE

# #####################################
# ### STEP 2: Filter the input
read -p "Press any key to continue... " -n1 -s


echo "$1|2|orthomclFilterFasta|Start|$(date)" 
echo "$1|2|orthomclFilterFasta|Start|$(date)" >> $LOGFILE

`"$LIBPATH"orthomclFilterFasta $container/1.compliantFasta 10 20`
mkdir $container/2.filteredFasta
mv goodProteins.fasta $container/2.filteredFasta/
mv poorProteins.fasta $container/2.filteredFasta/ 

echo "$1|2|orthomclFilterFasta|End|$(date)" 
echo "$1|2|orthomclFilterFasta|End|$(date)" >> $LOGFILE

# #####################################
# ### STEP 3.1:  Create BLAST database
read -p "Press any key to continue... " -n1 -s


echo "$1|3.1|makeblastdb|Start|$(date)"
echo "$1|3.1|makeblastdb|Start|$(date)" >> $LOGFILE

makeblastdb -in $container/2.filteredFasta/goodProteins.fasta  -dbtype prot

mkdir $container/3.blastdb

mv $container/2.filteredFasta/goodProteins.fasta.* $container/3.blastdb/

echo "$1|3.1|makeblastdb|End|$(date)" 
echo "$1|3.1|makeblastdb|End|$(date)" >> $LOGFILE


# #####################################
# ### 3.2 Split the input file
read -p "Press any key to continue... " -n1 -s


echo "$1|3.2|porthomclSplitFasta|Start|$(date)"
echo "$1|3.2|porthomclSplitFasta|Start|$(date)" >> $LOGFILE

mkdir $container/3.blastquery 

`"$LIBPATH"porthomclSplitFasta.py -i $container/2.filteredFasta/goodProteins.fasta  -o $container/3.blastquery`

echo "$1|3.2|porthomclSplitFasta|End|$(date)" 
echo "$1|3.2|porthomclSplitFasta|End|$(date)" >> $LOGFILE

# #####################################
# ### 3.3 Run blasts
read -p "Press any key to continue... " -n1 -s


echo "$1|3.3|Blasts|Start|$(date)" 
echo "$1|3.3|Blasts|Start|$(date)" >> $LOGFILE

cd $container
mkdir 3.blastres 
for query in 3.blastquery/*
do
    strand=$(basename "$query")
    strand="${strand%.*}" # remove .fasta
    blastp -query 3.blastquery/$strand.fasta  -db 3.blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads $NUM_THREADS -out 3.blastres/$strand.tab
done

echo "$1|3.3|Blasts|End|$(date)" 
echo "$1|3.3|Blasts|End|$(date)" >> $LOGFILE

# #####################################
# ### 4 Parse blasts
read -p "Press any key to continue... " -n1 -s


echo "$container|4|Parse Blast Results|Start|$(date)"
echo "$container|4|Parse Blast Results|Start|$(date)" >> $LOGFILE

mkdir $container/4.splitSimSeq

cd $container
for blres in $container/*
do
    strand=$(basename "$blres")
    strand="${strand%.*}" # remove .fasta
    `"$LIBPATH"porthomclBlastParser 3.blastres/$strand.tab 1.compliantFasta >> 4.splitSimSeq/$strand.ss.tsv`
done

echo "$container|4|Parse Blast Results|End|$(date)" 
echo "$container|4|Parse Blast Results|End|$(date)" >> $LOGFILE

# #####################################
# ### 5 Finding Best Hits
read -p "Press any key to continue... " -n1 -s


echo "$container|5|porthomclPairsBestHit|Start|$(date)" 
echo "$container|5|porthomclPairsBestHit|Start|$(date)" >> $LOGFILE

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

echo "$container|5|porthomclPairsBestHit|End|$(date)" 
echo "$container|5|porthomclPairsBestHit|End|$(date)" >> $LOGFILE

# #####################################
# ### Step 6: Finding Orthologs
read -p "Press any key to continue... " -n1 -s


echo "$container|6|Finding orthologs (porthomclPairsOrthologs)|Start|$(date)"
echo "$container|6|Finding orthologs (porthomclPairsOrthologs)|Start|$(date)" >> $LOGFILE

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

echo "$container|6|Finding orthologs (porthomclPairsOrthologs)|End|$(date)" 
echo "$container|6|Finding orthologs (porthomclPairsOrthologs)|End|$(date)" >> $LOGFILE

# #####################################
# ### Step 7:  Finding Paralogs
read -p "Press any key to continue... " -n1 -s



echo "$container|7.1|Finding Paralogs (awk)|Start|$(date)" 
echo "$container|7.1|Finding Paralogs (awk)|Start|$(date)" >> $LOGFILE

cd $container
mkdir 7.ogenes

# # genes in the second column

loop_start=1
loop_end=$((loop_start + NUM_THREADS - 1))
loop_end=$(($loop_end>$taxon_count?$taxon_count:$loop_end))

while [ $loop_start -lt $loop_max ]
do
	for i in $(seq $loop_start $loop_end)
	do
		{ TAXON_FILE="`sed -n $i\p taxon_list`"; awk -F'[|\t]' '{print $4 >> ("7.ogenes/"$3".og.tsv")}' 6.orthologs/$TAXON_FILE.ort.tsv } &
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

while [ $loop_start -lt $loop_max ]
do
	for i in $(seq $loop_start $loop_end)
	do
		{ TAXON_FILE="`sed -n $i\p taxon_list`";awk -F'[|\t]' '{print $2 >> ("7.ogenes/"$1".og.tsv")}' 6.orthologs/$TAXON_FILE.ort.tsv } &
	done
	wait
	loop_start=$((loop_end + 1))
	loop_end=$((loop_start + NUM_THREADS - 1))
	loop_end=($(($loop_end>$taxon_count?$taxon_count:$loop_end)))
	
done

echo "$container|7.1|Finding Paralogs (awk)|End|$(date)" 
echo "$container|7.1|Finding Paralogs (awk)|End|$(date)" >> $LOGFILE
###########
read -p "Press any key to continue... " -n1 -s

echo "$container|7.2|Finding Paralogs (porthomclPairsOrthologs)|Start|$(date)" 
echo "$container|7.2|Finding Paralogs (porthomclPairsOrthologs)|Start|$(date)" >> $LOGFILE

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


echo "$container|7.2|Finding Paralogs (porthomclPairsOrthologs)|Start|$(date)" 
echo "$container|7.2|Finding Paralogs (porthomclPairsOrthologs)|Start|$(date)" >> $LOGFILE



# #####################################
# ### Step 8:  Running MCL
read -p "Press any key to continue... " -n1 -s

echo "$container|8|Running MCL|Start|$(date)" 
echo "$container|8|Running MCL|Start|$(date)" >> $LOGFILE


cat $containe/6.orthologs/*.tsv >> $containe/8.all.ort.tsv

mcl $containe/8.all.ort.tsv  --abc -I 1.5 -t $NUM_THREADS -o $containe/8.all.ort.group

cat $containe/7.paralogs/*.tsv >> $containe/8.all.par.tsv

mcl $containe/8.all.par.tsv  --abc -I 1.5 -t $NUM_THREADS -o $containe/8.all.par.group


echo "$container|8|Running MCL|End|$(date)" 
echo "$container|8|Running MCL|End|$(date)" >> $LOGFILE