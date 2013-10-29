#!/bin/bash
#USAGE: ./specify_contaminant.sh bast_db blast_threads contaminant_reads.fa

db=$1
threads=$2
reads=$3

dir=$(dirname $0)

blastn -num_threads $threads -db $db -query $reads -outfmt 6 -dust yes -out blast_out.txt
python $dir/parse_blast_out.py blast_out.txt
