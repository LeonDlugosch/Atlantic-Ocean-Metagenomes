#!/bin/bash
SPAdes=/user/kila0393/SPAdes_3.15.2/bin/spades.py
m=500
threads=16
k=21,33,55

while : ; do
    case $1 in
    	-i)
            shift
            inDir=$1
            shift
            ;;
        -m)
            shift
            m=$1
            shift
            ;;
        -t)
            shift
            t=$1
            shift
            ;;
        -n)
            shift
            n=$1
            shift
            ;;
        -k)
            shift
            k=$1
            shift
            ;;
        -h|--help)
            echo "PLACEHOLDER"
            exit
            ;;
        *)  if [ -z "$1" ]; then break; fi
            inDir=$1
            shift
            ;;
    esac
done

main_dir=.
out_dir=${main_dir}/Contigs
tmp_dir=${main_dir}/tmp_${n}

mkdir -p ${out_dir}
mkdir -p ${tmp_dir}

for s in $(cat Files_${n}.txt); do
          
  if [ -d ${tmp_dir} ]; then
    echo Deleting old contig files...
    rm -rf ${tmp_dir}
    mkdir -p ${tmp_dir}
  fi
	
python $SPAdes -1 ${inDir}/${s}_R1.fastq -2 ${inDir}/${s}_R2.fastq --meta --memory $m -k $k -o ${tmp_dir} -t $t
        
mv ${tmp_dir}/contigs.fasta ${out_dir}/${s}_contigs.fasta
mv ${tmp_dir}/scaffolds.fasta ${out_dir}/${s}_scaffolds.fasta
mv ${tmp_dir}/contigs.paths ${out_dir}/${s}_contigs.paths
mv ${tmp_dir}/scaffolds.paths ${out_dir}/${s}_scaffolds.paths
       
done
rm Files_${n}.txt
