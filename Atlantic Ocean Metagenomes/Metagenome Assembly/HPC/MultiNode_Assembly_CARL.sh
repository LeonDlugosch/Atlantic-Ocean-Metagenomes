#!/bin/bash
inDir=0
njobs=2
p=mpcb.p
rt_d=10
rt_h=0
k=21,33,55
u=0
while : ; do
    case $1 in
    	-i)
            shift
            inDir=$1
            shift
            ;;
        -njobs)
            shift
            njobs=$1
            shift
            ;;
        -p)
            shift
            p=$1
            shift
            ;;
        -mail)
            shift
            mail=$1
            shift
            ;;
        -rt_d)
            shift
            rt_d=$1
            shift
            ;;
        -rt_h)
            shift
            rt_h=$1
            shift
            ;;
        -k)
            shift
            k=$1
            shift
            ;;
        -h|--help)
            echo "Splitting list of Metagenome Samples for SPAdes 3.15 assembly and submitting multiple Jobs to CARL"
            echo ""
            echo "USAGE:"
            echo "-i [PATH]: Directory containing forward and reverse reads. Filenames: *_R1.fastq, *_R2.fastq"
            echo "-njobs [INT]: Number of jobs to submit to CARL. Samples will be equally distibuted among jobs."
            echo "-p [OPTION]: CARL partition. default: mpcb.p (nodes: 128, max. threads: 16, max memory: 495GB), other options: mpcs.p (nodes: 158, max. threads: 24, max memory: 243GB), mpcp.p (nodes: 2, max. threads: 40, max memory: 1975G)"
            echo "-mail [EMAIL ADDRESS]: E-mail address for notification in case of FAIL or FINISH of jobs"
            echo "-rt_d [INT]: max. job runtime in Days. Jobs > 21d will not start! default: 10"
            echo "-rt_d [INT]: max. job runtime in hours - will be added to -rt_d. Jobs > 21d will not start! default: 0"
            echo "-k [INT,INT,INT...]: k-mer sizes for assembly. default: 21,33,55"
            exit
            ;;
        *)  if [ -z "$1" ]; then break; fi
            inDir=$1
            shift
            ;;
    esac
done

if [[ "$inDir" == 0 ]]; then
    echo "Splitting list of Metagenome Samples for SPAdes 3.15 assembly and submitting multiple Jobs to CARL"
    echo ""
    echo "USAGE:"
    echo "-i [PATH]: Directory containing forward and reverse reads. Filenames: *_R1.fastq, *_R2.fastq"
    echo "-njobs [INT]: Number of jobs to submit to CARL. Samples will be equally distibuted among jobs."
    echo "-p [OPTION]: CARL partition. default: mpcb.p (nodes: 128, max. threads: 16, max memory: 495GB), other options: mpcs.p (nodes: 158, max. threads: 24, max memory: 243GB), mpcp.p (nodes: 2, max. threads: 40, max memory: 1975G)"
    echo "-mail [EMAIL ADDRESS]: E-mail address for notification in case of FAIL or FINISH of jobs"
    echo "-rt_d [INT]: max. job runtime in Days. Jobs > 21d will not start! default: 10"
    echo "-rt_d [INT]: max. job runtime in hours - will be added to -rt_d. Jobs > 21d will not start! default: 0"
    echo "-k [INT,INT,INT...]: k-mer sizes for assembly. default: 21,33,55"

    exit
fi

if [[ "$p" == "mpcb.p" ]]; then
    m=495
    t=16
fi
if [[ "$p" == "mpcs.p" ]]; then
    m=243
    t=24
fi
if [[ "$p" == "mpcp.p" ]]; then
    m=1975
    t=40
fi

( cd $inDir && ls *.fastq ) | cut -f 1 -d . | cut -f 1,2 -d _ | uniq -d > ./Files.txt

nFiles=$(cat ./Files.txt | wc -l)
nSamples=$(echo "scale=1;($nFiles/$njobs)" | bc | awk '{print ($0-int($0)<0.0001)?int($0):int($0)+1}')

c=0
for ((i=1; i<=njobs; i++)); do
    if [[ "$i" == 1 ]]; then
        from=$(echo $i)
    else
        from=$(echo $c)
    fi

    if [[ "$i" < "$njobs" ]]; then
        to=$(echo "$(($i * $nSamples))")
    else
        to=$(echo "$nFiles")
    fi
    
    c=$(echo "$(($to + 1))")
    echo "FROM: $from"
    echo "TO: $to"
    awk "NR==$from, NR==$to" Files.txt
    awk "NR==$from, NR==$to" Files.txt > Files_${i}.txt
    echo "#!/bin/bash" > JOB_$i.slurm
    echo "#SBATCH --ntasks=1" >> JOB_$i.slurm
    echo "#SBATCH --cpus-per-task=$t" >> JOB_$i.slurm
    echo "#SBATCH --mem=${m}G" >> JOB_$i.slurm
    echo "#SBATCH --time=${rt_d}-${rt_h}:00:00" >> JOB_$i.slurm
    echo "#SBATCH --output=LOG.%A.out" >> JOB_$i.slurm
    echo "#SBATCH --error=ERRORS.%A.err" >> JOB_$i.slurm
    echo "#SBATCH --mail-type=START,END,FAIL" >> JOB_$i.slurm
    echo "#SBATCH --mail-user=$mail" >> JOB_$i.slurm
    echo "#SBATCH --partition=$p" >> JOB_$i.slurm
    echo "##### SBATCH --array=1-10%3" >> JOB_$i.slurm
    echo "bash SPAdes315_MG.sh -i $inDir -m $m -t $t -n $i -k $k" >> JOB_$i.slurm

    sbatch JOB_$i.slurm
    echo ""
    echo ""
done
rm Files.txt