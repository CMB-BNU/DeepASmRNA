#!/usr/bin/sh

#make blast database
makeblastdb -in $1 -dbtype nucl

#sequence alignment using blastn
blastn -query $1 -db $1 -strand plus -evalue 1E-10 -outfmt 5 -ungapped -num_threads 20 -out $1\out.xml 

conda create -n CoAST python=3.7
conda activate CoAST
conda install tensorflow=2.1

#predict AS transcript pair

python3 predictAS.py $1\out.xml $1\as.txt >$1\as.seq

#classify AS transcript pair

python classify.py $1\as.seq -m $2 -o result.txt
