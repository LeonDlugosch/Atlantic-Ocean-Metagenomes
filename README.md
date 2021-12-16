# Atlantic Ocean Metagenomes Reference Gene Cathalogue
## Assembly and analysis pipeline

The assembly pipeline is written for server-environment of the Institute for [Chemistry and Biology of the Marine Environment (ICBM)](https://uol.de/en/icbm) of the Carl-von-Ossietzky University Oldenburg, Germany as it uses Institute servers as well as the [HPC](https://uol.de/en/school5/sc/high-perfomance-computing/hpc-facilities) facility of the University of Oldenburg. Third party software may required to run the assembly.

Downstream analysis such as statistics etc. are included in the analysis directory.

This Page is under construction and a proper documentation will be added in the near future (until the end of 2021).

## About the dataset
Fractionated CTD samples from 20m depth were taken on the RV Polarstern cruise ANT28-4 and ANT28-5 from March to May 2012 and metagenomes from 22 stations were generated using Illumina HiSeq2500 ([Dlugosch et al.]() for details).
<p align="center">
  
  <img src="images/GitMap_web.png">
  
</p>
>Station from RV Polarstein cruises ANT28-4 and ANT28-5. Colour indicates [average chlorophyll a concentration from 2012]().
</br>
</br>
In short, Illumina sequences were quality trimmed and residual adapter seqeunces were removed using [timmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
and subsequently assembled using [metaSPAdes](https://cab.spbu.ru/software/meta-spades/) (HPC). Genes from contigs ware predicted using [Prodigal](https://github.com/hyattpd/Prodigal), dereplicated and clustered at 95% sequence identity [using usearch](https://www.drive5.com/usearch/manual/cmd_cluster_fast.html). Resulting seqeunces were classified taxonomically using [Kaiju](https://kaiju.binf.ku.dk/) with the ProGenomes and Refseq databases (seqeunce taxonomy is integrated during [dataset generation](https://github.com/LeonDlugosch/Atlantic-Ocean-Metagenomes/blob/main/Analysis/00_AOM_Taxonomy_Function_CountTable.R)). Sequences were functionally classified using [GhostKOALA](https://www.kegg.jp/ghostkoala/) and the [CAZyme database](http://www.cazy.org/) ([diamond](https://github.com/bbuchfink/diamond) blastx --more-sensitive). High quality reads were mappted to the AOM-RGC to using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). 

See [here]() for a more detailed manual of the assembly-pipeline. 

## Data availability
Illumina seqeuncing data, assemblies, Atlantic Ocean Metagenome gene cathalogue (AOM-RGC), gene abundance tables and environmental data are available here:  

[Illumina sequencing data](https://www.ebi.ac.uk/ena/browser/view/PRJEB34453) </br>
[Environmental data](https://doi.pangaea.de/10.1594/PANGAEA.906247) </br>
[AOM-RGC, gene abundance tables, assemblies and predicted genes]() </br>

