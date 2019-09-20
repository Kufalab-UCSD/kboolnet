#!/bin/bash
######################################################
# Adrian C
# Simple script to extract modules from rxncon file
# and make regulatory/species reaction graphs for each
# module
######################################################

##### Configuration #####
inFile="" # Google Drive URL of rxncon file
outDir= # Directory to which output will be written
modules=() # Space-separated list of modules
rxnconDir= # Directory where rxncon scripts are located
kboolnetDir= # Root directory of kboolnet repository

#### Setup ####
date=$(date '+%y%m%d')
rxnconExcel=$outDir/${date}_master.xlsx

mkdir -p $outDir
mkdir -p $outDir/excel
mkdir -p $outDir/regulatory
mkdir -p $outDir/sr
rm -f $outDir/excel/* $outDir/regulatory/* $outDir/sr/*

curl $inFile --output $rxnconExcel

for module in "${modules[@]}"
do
    echo "Extracting module $module"
    python3 $kboolnetDir/Python/extract_modules.py --file $rxnconExcel --modules $module --output $outDir/excel/$module.xlsx
    echo "Creating regulatory graph for module $module"
    python3 $rxnconDir/rxncon2regulatorygraph.py $outDir/excel/$module.xlsx --output $outDir/regulatory/${module}
    echo "Creating species-reaction graph for module $module"
    python3 $rxnconDir/rxncon2srgraph.py --output $outDir/sr/${module} $outDir/excel/$module.xlsx
done
