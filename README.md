![alt tag](ressources/CoCoLogo.PNG)
# **CoCo**: Count Corrector for embedded and multi-mapped genes.

CoCo is a pipeline designed to improve the evaluated abundance of embedded and multi-mapped genes from RNA-seq alignment data. CoCo is divided in two main parts which can be used together or seperately, depending on the user's preference. The first part is a correction of the gene annotation used by read assignment softwares such as featureCounts or HTSeq in order to correct the evaluated read counts for embedded genes such as snoRNA, that overlap features of their host gene's transcripts such as retained introns and exons. The second part of the correction distributes multi mapped reads in relation to the evaluated read counts obtained from single mapped reads.

We've shown that the combination of these two corrections gives the best correlation of the evaluated gene abundance from RNAseq data with abundance evaluated from other methods such as qPCR.



## **Getting Started**

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

CoCo is supported on Linux (Ubuntu, Fedora, CentOS).

### **Prerequisites**

Here are the tools that must be installed for CoCo to work:

* [Subread] (v1.5.2 or higher)(http://subread.sourceforge.net/), for the *featureCounts* function, which is used by correct_count to produce the read counts per genes. [Click here for installation guidelines](http://bioinf.wehi.edu.au/subread-package/)
* [BEDtools] (v.2.25.0 or higher)(http://bedtools.readthedocs.io/en/latest/), for the *intersect* function which is used by correct_annotation and for *genomecov* which is used by correct_bedgraph. [Click here for installation guidelines](http://bedtools.readthedocs.io/en/latest/content/installation.html)
* [samtools] (v1.3.1 or higher)(http://samtools.sourceforge.net/), for the *view* function which is used to verify bam integrity and to get max read length. [Click here for installation guidelines](http://www.htslib.org/download/)
* [pairedBamToBed12](https://github.com/Population-Transcriptomics/pairedBamToBed12), only required for the correct_bedgraph mode.


CoCo scripts are mostly python3, so be sure you have python3 installed and working. (tested for python 3.5)
CoCo also uses some python3 packages:
* [pandas] (v0.18.0 or higher)(http://pandas.pydata.org/), [Click here for installation guidelines](http://pandas.pydata.org/pandas-docs/stable/install.html)



### **Installing**

Installing CoCo is done in a few easy steps.

First, you must clone the git repository to your favorite directory.

```
cd /path/to/clone/
git clone http://gitlabscottgroup.med.usherbrooke.ca/scott-group/coco.git
```

Next, you have to add the /bin/ path of the coco repository to your PATH environment variable.
this is easily done by opening your .bashrc file in your $HOME repository and add this line at the bottom:

```
export PATH=/path/to/clone/coco/bin:$PATH
```

And then you're pretty much done! You should have access to the **coco** command from your terminal.

 ```
coco --help
>>	CoCo: Count Corrector for embedded and multi-mapped genes.
	Version: 0.2.7
	Usage: CoCo <Run mode> <Run mode specific arguments>

	Run modes:
		correct_annotation, CA, ca:	Produce modified annotation for embedded genes.
		correct_count, CC, cb:		Produce gene expression values, taking multi-mapped reads into account.
		correct_bedgraph, CB, cb:	Produce paired-end bedgraphs.

	-h, --help      show this help message and exit
	-v, --version   show CoCo version

	For more info, please refer to the README file.
```



## **Usage**

Once the repository is downloaded and its /bin/ path is added to your PATH environment variable, CoCo can be called directly from your terminal.
CoCo is divided into two main run modes: **correct_annotation** and **correct_count**, and one accessory mode: **correct_bedgraph**, each with their own parameters, and can be used as such:

```
coco [run mode] [args]
```

You can also run it directly from the bin:

```
cd path/to/coco/bin/
./coco [run mode] [args]
```

Here is how a user (you) can use the three CoCo modules to get to the desired output (whether it be the corrected gene expression values, the corrected paired-end bedgraph or both)

<img src="ressources/CoCoPipeline.PNG" alt="CoCo Pipeline" style="width: 700px;"/>

For detailed information about the usage of every run modes, please refer to the [MANUAL.md](MANUAL.md).



## **Authors**

* **Gabrielle Deschamps-Francoeur** - *Making of the multi-mapped read correction and collaborated on the correction for ambiguous reads.*
* **Vincent Boivin** - *Making of the annotation correction for ambiguous reads.* - [boiv2803](http://gitlabscottgroup.med.usherbrooke.ca/u/boiv2803)
* **Michelle Scott** - *Mastermind of the project and awesome PI* - [scottgroup](http://scottgroup.med.usherbrooke.ca/)

Questions should be directed to: _gabrielle.deschamps-francoeur@usherbrooke.ca_


## **License**

CoCo: Count Corrector for embedded and multi-mapped genes.
Copyright (C) 2017 Gabrielle Deschamps-Francoeur, Vincent Boivin & Michelle Scott

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

For more information about the GNU General Public License, see <http://www.gnu.org/licenses/>.



## **How to cite us**

Please cite as :

Deschamps-Francoeur G & Boivin V, Abou Elela S and Scott MS. CoCo: RNA-seq Read Assignment Correction for Nested Genes 
and Multimapped Reads. Bioinformatics. 2019. Available from : <https://doi.org/10.1093/bioinformatics/btz433>



## **Acknowledgments**

* Hat tip to Fabien Dupuis-Sandoval for starting this project from scratch.
* Many thanks to Jean-Michel Garant for his general wisdom.
* Thanks to Kamil Slowikowski [slowkow](https://gist.github.com/slowkow) for its GTF.py script that is used by correct_annotation.



## **References**
**CoCo**: Deschamps-Francoeur G & Boivin V, Abou Elela S and Scott MS. CoCo: RNA-seq Read Assignment Correction for Nested Genes 
and Multimapped Reads. Bioinformatics. 2019. Available from : <https://doi.org/10.1093/bioinformatics/btz433>

**featureCounts**: Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014;30(7):923-30.

**BEDtools**: Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841-2.

**samtools**: Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078-9.

