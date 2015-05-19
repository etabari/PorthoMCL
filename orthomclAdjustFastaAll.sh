#!/bin/bash

if [ -z "$1" ] 
then
	echo "------------------- orthomclAdjustFastaAll  help -------------------"
	echo
    echo "orthomclAdjustFastaAll will look into all the subfolders of supplied genome folder and adjust faa files to fasta using the strand as the genome identifier and 4 column"
    echo
    echo "if you don't understand the above statement please read orthomclAdjustFasta help"
    echo
    echo
    echo "orthomclAdjustFastaAll <genomesfolder>"
    echo
    echo "		Please provide the folder you have your genomes stored."
    echo
    echo "------------------- orthomclAdjustFasta  help -------------------"
    ./orthomclAdjustFasta
    exit 1
fi

GENOMES="$1"

for genome in $GENOMES/*
do
	for faa in $genome/*.faa
	do
		echo 'processing' $faa

		strand=$(basename "$faa")
		strand="${strand%.*}" # remove .faa

		./orthomclAdjustFasta $strand $faa 4 || exit 1
	done

done
