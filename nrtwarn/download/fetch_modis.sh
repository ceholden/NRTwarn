#!/bin/bash
#$ -V
#$ -l h_rt=24:00:00
#$ -j y

if [ -z $1 ]; then
    echo "Error - must specify a tile"
    exit 1
fi
tile=$1

if [ -z $2 ]; then
    echo "No years specified - collecting all years"
    years=$(/usr/bin/seq 2000 $(date +%Y))
else
    years=$2
fi

mkdir -p $1/download
cd $1/download

# Usage
# -v --> verbose
# -q --> check for existing files
# -p --> MODIS product
# -s --> sensor platform
# -y --> year
# -t --> tile
for year in $years; do
    get_modis.py -v -q -p MOD09GQ.005 -s MOLT -y $year -t $tile -o $year
    get_modis.py -v -q -p MYD09GQ.005 -s MOLA -y $year -t $tile -o $year
    get_modis.py -v -q -p MOD09GA.005 -s MOLT -y $year -t $tile -o $year
    get_modis.py -v -q -p MYD09GA.005 -s MOLA -y $year -t $tile -o $year
done

echo "Done!"
