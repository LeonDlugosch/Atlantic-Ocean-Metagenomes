#!/bin/sh
#
# Defaults
rDir=.
threads=30
cDir=0
idx=1
oDir=.
depth=0
mode=db
readmode=0
sortme=0
db=0
a_rm=no
# Colors
red='\033[0;31m'
green='\033[0;32m'
nc='\033[0m'
# Programms 
SumBAM=/bioinf/home/leon.dlugosch/Resources/metabat/jgi_summarize_bam_contig_depths

while : ; do
    case $1 in
        -r)
            shift
            rDir=$1
            shift
            ;;
        -o)
            shift
            oDir=$1
            shift
            ;;
        -depth)
            shift
            depth=1
            shift
            ;;
        -idx)
            shift
            idx=1
            shift
            ;;
        -db)
            shift
            db=$1
            shift
            ;;
        -se)
            shift
            readmode=1
            shift
            ;;
        -pe)
            shift
            readmode=0
            shift
            ;;
        -adaper_rm)
            shift
            a_rm=yes
            shift
            ;;
        -deplete_rRNA)
            shift
            sortme=yes
            shift
            ;;
        -c)
            shift
            mode=co
            cDir=$1
            shift
            ;;
        -t)
            shift
            threads=$1
            shift
            ;; 
        -h|--help)
            echo "This script uses Bowtie2 (very-sensitive mode) and samtools to map paired Illumina reads to assembled contigs and outputs .sam, .bam and, if the -idx option is set read abuandance files."
            echo "USAGE:"
            echo "-r or -rDir [PATH]: define path to a directory containing quality quality checked reads from Illumina sequencing (mandatory)"
            echo "-c or -cDir [PATH]: define path to a directory containing corresponding assembled contigs from the reads used for mapping (mandatory)"
            echo "NOTE: file names should be: XZY_R1.fastq XZY_R2.fastq for reads and XZY_contigs.fasta"
            echo "-o [PATH]: Directory in which results will be saved. Default: ./"
            echo "-read.mode [OPTION]: Should single-end [se] or paired-end [pe] reads be mapped to the database?"
            echo "-t [INT]: number of threads allocated for the process. Default: 4"
            echo "-idx: output mapped read abundance as .txt file"
            exit
            ;; 
        *)  
            if [ -z "$1" ]; then break; fi
            ERR=1
            shift
            ;;
    esac
done

echo "idx: " $idx
echo "rdir: " $rDir
echo "odir: "$oDir
echo "depth: "$depth
echo "db: "$db

############################# Checking & editing variables #############################
mkdir -p $oDir/temp
echo "$oDir"
if [[ "${rDir}" == 0 ]]; then
    echo "ERROR!"
    echo -e "${red}No read-directory defined. Script aborted.${nc}"
    echo -e "${red}Use -rDir /PATH/TO/READS to tell the script where they are.${nc}"
    exit
fi

if [[ ! -d $oDir ]]; then 
    echo "ERROR"
    echo -e "${red}Output directory does not exist. Set with -o or -oDIR /PATH/TO/OUTPUT/DIR${nc}"
    exit
fi

if [[ "${rDir: -1}" == "/" ]]; then
    rDir=${rDir::-1}
fi

if [[ "${oDir: -1}" == "/" ]]; then
    oDir=${oDir::-1}
fi
####################################### Adapter removal ########################################
if [[ "$a_rm" == "yes" ]]; then
mkdir -p $oDir/trimmo/Paired
mkdir -p $oDir/trimmo/Unpaired

for Sample in $(cat $Dir/Temp/Files.txt)
    do
    java -jar $trimmo PE \
                     -threads $threads \
                     $inDir/${Sample}_*1.fastq \
                     $inDir/${Sample}_*2.fastq \
                     $oDir/trimmo/Paired/${Sample}_R1.fastq \
                     $oDir/trimmo/Unpaired/${Sample}_SE_R1.fastq \
                     $oDir/trimmo/Paired/${Sample}_R2.fastq \
                     $oDir/trimmo/Unpaired/${Sample}_SE_R2.fastq \
                     SLIDINGWINDOW:4:15 \
                     MINLEN:100 \
                     ILLUMINACLIP:/bioinf/home/leon.dlugosch/Resources/Adapter/Contaminants_2.fas:2:30:10

done
rDir=$oDir/trimmo/Paired
fi
####################################### rRNA depletion ########################################

if [[ "${sortme}" != 0 ]]; then
    echo $sortme
    echo "rRNA depletion start..."
    ( cd $rDir && ls *.fastq ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d > $oDir/temp/files.txt
    for s in $(cat $oDir/temp/files.txt); do
    
    mkdir $oDir/SortmeRNA/Out -p
    mkdir $oDir/SortmeRNA/Unaligned -p

        if [[ "${readmode}" == 0 ]]; then
sortmerna   --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/bac/silva-bac-16s-id90.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/bac/silva-bac-23s-id98.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/euk/silva-euk-18s-id95.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/euk/silva-euk-28s-id98.fasta \
          	--ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/5s/rfam-5.8s-database-id98.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/5s/rfam-5s-database-id98.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/arc/silva-arc-16s-id95.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/arc/silva-arc-23s-id98.fasta \
            --reads $rDir/${s}_R1.fastq \
            --reads $rDir/${s}_R2.fastq \
            --workdir $oDir/SortmeRNA/Out/ \
            --other $oDir/SortmeRNA/Unaligned/ \
            --threads 1:1:$threads \
            --paired_out \
            --fastx \
            -e 0.00001 \
            -v
        fi

            
        if [[ "${readmode}" == 1 ]]; then
sortmerna --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/bac/silva-bac-16s-id90.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/bac/silva-bac-23s-id98.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/euk/silva-euk-18s-id95.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/euk/silva-euk-28s-id98.fasta \
          	--ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/5s/rfam-5.8s-database-id98.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/5s/rfam-5s-database-id98.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/arc/silva-arc-16s-id95.fasta \
            --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/arc/silva-arc-23s-id98.fasta \
            --reads $rDir/${s}_R1.fastq \
            --reads $rDir/${s}_R2.fastq \
            --workdir $oDir/SortmeRNA/Out/ \
            --other $oDir/SortmeRNA/Unaligned/ \
            --threads 1:1:$threads \
            --paired_out \
            --fastx \
            -e 0.00001 \
            -v 
        fi
done

rDir=$oDir/SortmeRNA/Unaligned

fi
####################################### Contig Mapping ########################################

if [[ "$mode" == "co" ]]; then

if [[ "${cDir}" == 0 ]]; then
    echo "ERROR!"
    echo -e "${red}No contig-directory defined. Script aborted.${nc}"
    echo -e "${red}Use -cDir /PATH/TO/CONTIGS to tell the script where they are.${nc}"
    exit
fi

if [[ "${cDir: -1}" == "/" ]]; then
    cDir=${cDir::-1}
fi

contigs=$( cd $cDir && ls *.fa* ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | wc -l
reads=$( cd $rDir && ls *.fa* ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d | wc -l
if [[ "$contigs" == "$reads" ]]; then
    echo -e "${green}Number of Read files matches number of contig files!${nc}"
else
    echo -e "${red}Number of Read files does not match number of contig files!${nc}"
    echo "Paired read libraries: "$reads
    echo "Contig files: "$contigs
    exit
fi

mkdir -p $oDir/bam
mkdir -p $oDir/sam
ls $rDir

( cd $rDir && ls *.fastq ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d > $oDir/temp/files.txt
for s in $(cat $oDir/temp/files.txt)
    do
echo -e "Bulding bowtie2 database from ${green}"${s}"_contigs.fasta${nc} ..."
bowtie2-build $cDir/${s}_contigs.fasta $oDir/temp/Bowtie2.db
DB=$oDir/temp/Bowtie2.db

echo
echo -e "Mapping reads from ${green}"$s"_R1.fastq${nc} and ${green}"$s"_R2.fastq${nc} to database..." 
bowtie2 --very-sensitive-local \
-x $DB \
-1 $rDir/${s}_R1.fastq \
-2 $rDir/${s}_R2.fastq \
-p $threads \
-S $oDir/temp/${s}.sam

echo -e "Reformating ${green}"$s".sam${nc} to ${green}"$s".bam${nc} ..."
samtools view -b -S $oDir/temp/${s}.sam > $oDir/bam/${s}.bam
echo -e "Deleting ${green}"$s".sam${nc}..."
rm $oDir/temp/${s}.sam
echo -e "Sorting ${green}"$s".bam${nc}..."
samtools sort $oDir/bam/${s}.bam > $oDir/bam/${s}.sorted.bam 
rm $oDir/bam/${s}.bam

if [[ "$idx" == 1 ]]; then
echo -e "Creating, sorting and indexing ${green}"$s".bam${nc} file..."
mkdir -p $oDir/map 
samtools index $oDir/bam/${s}.sorted.bam 
samtools idxstats $oDir/bam/${s}.sorted.bam  > $oDir/map/${s}.mapped.txt
fi

if [[ "$depth" == 1 ]]; then
echo -e "Generating depth-file from ${green}" ${s}".sorted.bam${nc}..."
mkdir -p $oDir/depth
$SumBAM --outputDepth $oDir/depth/${s}_depth.txt $oDir/bam/${s}.sorted.bam 
fi
echo ""
echo -e "Finished mapping of ${green}"$s"${nc}!"
done 
fi

####################################### Database Mapping ########################################

if [[ "${db}" != "0" ]]; then
mode=db
fi
if [[ "${mode}" == "db" ]]; then

# if [[ "${db}" == 0 ]]; then
#     echo "ERROR!"
#     echo -e "${red}No database file defined. Script aborted.${nc}"
#     echo -e "${red}Use -db /PATH/TO/FILE to tell the script where the database is.${nc}"
#     exit
# fi

mkdir -p $oDir/bam
mkdir -p $oDir/sam
mkdir -p $oDir/map

filetype=$(echo "${db}" | rev | cut -d '.' -f 1 | rev)
echo "$db"
echo $filetype
if [[ "$filetype" == "fasta" ]] || [[ "$filetype" == "fa" ]] || [[ "$filetype" == "fna" ]] || [[ "$filetype" == "faa" ]] || [[ "$filetype" == "fastq" ]]; then
mkdir -p $oDir/db
echo "Selected database is formated in fasta-format. Bowtie2 database will be created and can be found in the db-directory in the output folder."
bowtie2-build ${db} $oDir/db/Bowtie2.db
DB=$oDir/db/Bowtie2.db
fi

if [[ "$filetype" == "db" ]]; then
echo "Using provided Bowtie2 Database."
DB=${db}
fi
echo "$filetype"
if [[ "$filetype" != "fasta" ]] && [[ "$filetype" != "fa" ]] && [[ "$filetype" != "fna" ]] && [[ "$filetype" != "faa" ]] && [[ "$filetype" != "db" ]]; then
echo "ERROR!"
echo -e "${red}No suitable database provided. Script aborted.${nc}"
exit
fi
echo "Readmode is:" ${readmode}
####################################### PE read mapping ########################################

if [[ "${readmode}" == "0" ]]; then
( cd $rDir && ls *.fastq ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d > $oDir/temp/files.txt
for s in $(cat $oDir/temp/files.txt); do
echo -e "Mapping reads from ${green}"$s"_R1.fastq${nc} and ${green}"$s"_R2.fastq${nc} to database..." 
bowtie2 --very-sensitive-local \
-x $DB \
-1 $rDir/${s}_R1.fastq \
-2 $rDir/${s}_R2.fastq \
-p $threads \
-S $oDir/temp/${s}.sam

echo -e "Reformating ${green}"$s".sam${nc} to ${green}"$s".bam${nc} ..."
samtools view -b -S $oDir/temp/${s}.sam > $oDir/bam/${s}.bam
echo -e "Deleting ${green}"$s".sam${nc}..."
rm $oDir/temp/${s}.sam
echo -e "Sorting ${green}"$s".bam${nc}..."
samtools sort $oDir/bam/${s}.bam > $oDir/bam/${s}.sorted.bam 
rm $oDir/bam/${s}.bam
echo -e "Creating, sorting and indexing ${green}"$s".bam${nc} file..."
samtools index $oDir/bam/${s}.sorted.bam 
samtools idxstats $oDir/bam/${s}.sorted.bam  > $oDir/map/${s}.mapped.txt
echo -e "Deleting ${green}"$s".bam files${nc}..."
rm $oDir/bam/*.bam
rm $oDir/bam/*.bai
done 
fi

####################################### SE read mapping ########################################
echo "SE mapping" 
if [[ "${readmode}" == "1" ]]; then
echo "readmode se"
ls $rDir
( cd $rDir && ls *.fastq ) |  awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d> $oDir/temp/files.txt
echo "oDir: "$oDir
echo "DB: "$DB
echo "U: "$rDir/${s}
echo "p: "$threads 
echo "S: "$oDir/temp/${s}.sam
cat $oDir/temp/files.txt

for s in $(cat $oDir/temp/files.txt); do
echo -e "Mapping reads from ${green}"$s"${nc}"
bowtie2 --very-sensitive-local \
-x $DB \
-U $rDir/${s} \
-p $threads \
-S $oDir/temp/${s}.sam

echo -e "Reformating ${green}"$s".sam${nc} to ${green}"$s".bam${nc} ..."
samtools view -b -S $oDir/temp/${s}.sam > $oDir/bam/${s}.bam
echo -e "Deleting ${green}"$s".sam${nc}..."
rm $oDir/temp/${s}.sam
echo -e "Sorting ${green}"$s".bam${nc}..."
samtools sort $oDir/bam/${s}.bam > $oDir/bam/${s}.sorted.bam 
rm $oDir/bam/${s}.bam
echo -e "Creating, sorting and indexing ${green}"$s".bam${nc} file..."
samtools index $oDir/bam/${s}.sorted.bam 
samtools idxstats $oDir/bam/${s}.sorted.bam  > $oDir/map/${s}.mapped.txt
echo -e "Deleting ${green}"$s".bam files${nc}..."
rm $oDir/bam/*.bam
rm $oDir/bam/*.bai
done
fi
fi
##################################### Removing Temp-files #######################################
echo ""
echo "Finishing touches..."
#rm -rf $oDir/temp
echo "Done! :-)"