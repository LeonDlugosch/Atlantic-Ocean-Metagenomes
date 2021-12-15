# Atlantic Ocean Metagenomes Reference Gene Cathalogue
## Assembly and analysis pipeline (under construction)

The assemblie pipeline is written for server-environment of the Institute for [Chemistry and Biology of the Marine Environment (ICBM)](https://uol.de/en/icbm) of the Carl-von-Ossietzky University Oldenburg, Germany as it uses Institute servers as well as the [HPC](https://uol.de/en/school5/sc/high-perfomance-computing/hpc-facilities) facility of the University of Oldenburg. Third party software may required to run the assembly.

Downstream analysis such as statistics etc. are included in the analysis directory.

This Page is under construction and a proper documentation will be added in the near future (until the end of 2021).

## About the dataset
Fractionated samples were taken on the RV Polarstern cruise ANT28-4 and ANT28-5 from March to May 2012 and metagenomes from 22 stations were generated using Illumina HiSeq2500 (click here for details).
<p align="center">
  <img src="images/GitMap_web.png">
</p>

In short, Illumina sequences were quality trimmed and residual adapter seqeunces were removed using [timmomatic](http://www.usadellab.org/cms/?page=trimmomatic):
```
java -jar timmomatic.jar PE \
          -phred33 \
          -threads 20 \
          ANT28_179_R1.fastq \
          ANT28_179_R1.fastq \
          ./01_QC/Paired/ANT28_179_R1.fastq \
          ./01_QC/Unpaired/ANT28_179_R1.fastq \
          ./01_QC/Paired/ANT28_179_R2.fastq \
          ./01_QC/Unpaired/ANT28_179_R2.fastq \
          ILLUMINACLIP:adapter.fa:2:30:10:2:true \
          SLIDINGWINDOW:4:15 \
          LEADING:3 \
          MINLEN:100 
```
and subsequently assembled using [metaSPAdes](https://cab.spbu.ru/software/meta-spades/) (HPC).
```
python spades.py -1 ./01_QC/Paired/ANT28_179_R1.fastq \
                 -2 ./01_QC/Paired/ANT28_179_R2.fastq \
                 --meta \
                 --memory 243 \
                 -k 11,33,55 \
                 -o ./tmp \
                 -t 16
```

## Data availability
Illumina seqeuncing data, assemblies, Atlantic Ocean Metagenome gene cathalogue (AOM-RGC), gene abundance tables and environmental data are available here:  

[Illumina sequencing data](https://www.ebi.ac.uk/ena/browser/view/PRJEB34453) </br>
[Environmental data](https://doi.pangaea.de/10.1594/PANGAEA.906247) </br>
[AOM-RGC, gene abundance tables, assemblies and predicted genes]() </br>

