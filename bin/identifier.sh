#!/usr/bin/sh

help() {
 
    #echo 'def help is running'
    #echo "终端输入脚本名称,参数数量,参数: $0,$#,$@"
    echo "please do as follow:"
    # echo "注意 -r -u -f -c 参数请勿同时使用"
    #echo "    -h          : help information"
    echo "    Userage : identifier.sh transcript.fasta model ### model, choosing from [arabidopsis, rice, human], arabidopsis or rice for plant, human for animal"
    echo "    All the output files located in the path of input fasta file which prefix is the name of the fasta file"
    exit 1
}
 
if [[ $# == 0 || "$1" == "-h"  || "$1" == "-help"  || "$1" == "--help" ]]; then
    help
    exit 1
fi
if [[ $# != 2 ]]; then
    echo "We need two parameters,1st for transcript.fasta,2th for predicting model"
    exit 1
fi
filename=$1
if [[ $# == 2 ]]; then
    if [[ ${filename:0-2} != "fa" && ${filename:0-5} != "fasta" ]]; then
        echo "We need two parameters,1st for transcript.fasta"
        exit 1
    fi
fi
if [[ $# == 2 && $2 != "arabidopsis" && $2 != "rice" && $2 != "human" ]]; then
    echo "We need two parameters,2th for predicting model, choosing from [arabidopsis, rice, human], arabidopsis or rice for plant, human for animal"
    exit 1
fi

echo "Step1: make blast database"
date
makeblastdb -in $1 -dbtype nucl

echo "Step2: sequence alignment using blastn"
date
blastn -query $1 -db $1 -strand plus -evalue 1E-10 -outfmt 5 -ungapped -num_threads 20 -out $1\out.xml 

#conda create -n CoAST python=3.7
#conda activate CoAST
#conda install tensorflow=2.1
#conda install biopython

echo "Step3: predict AS transcript pair"
date

python predictAS.py $1\out.xml $1 $1\as.txt >$1\as.seq

echo "Step4: classify AS transcript pair"
date

python classifyAS.py $1\as.seq -m $2 -o $1\as_type.txt

echo "done,thank u"
date
