# **CoCo**: Count Corrector for embedded and multi-mapped genes.

CoCo is a pipeline designed to improve the evaluated abundances of embedded and multi-mapped genes from RNA-seq alignment data. CoCo is divided in two main parts which can be used together or seperately, depending on user's preferences. The first part being a correction of the gene annotation used by read assignment softwares such as featureCounts or HTSeq in order to correct the evaluated read counts for embedded genes such as snoRNA, that overlap features of their host gene's transcripts such as retained introns and exons. The second part of the correction distributes multi mapped reads in relation to the evaluated read counts obtained from single mapped reads.

We've shown that the combination of these two corrections gives the best correlation of the evaluated gene abundance from RNAseq data with abundances evaluated from other methods such as qPCR.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Here are the tools that must be installed for CoCo to work:
* [Subread](http://subread.sourceforge.net/), for the *featureCounts* function, which is used by CorrectCount to produce the read counts per genes. [Click here for installation guidelines](http://bioinf.wehi.edu.au/subread-package/)
* [BEDtools](http://bedtools.readthedocs.io/en/latest/), for the *intersect* function which is used by CorrectAnnotation and for *genomecov* which is used by CorrectBedgraph.[Click here for installation guidelines](http://bedtools.readthedocs.io/en/latest/content/installation.html)

CoCo scripts are mostly python3, so be sure you have python3 installed and working. (tested for python 3.5)
CoCo also uses some python3 packages:
* [pandas](http://pandas.pydata.org/), [Click here for installation guidelines](http://pandas.pydata.org/pandas-docs/stable/install.html)

### Installing

Installing CoCo is done in a few easy steps.

First, you must clone the git repositroy to your favourite directory.

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
	Version: 1.0.0
	Usage: CoCo <Run mode> <Run mode specific arguments>

	Run modes:
		CorrectAnnotation:	Produce modified annotation for embedded genes.
		CorrectCount:		Produce gene expression values, taking multi-mapped reads into account.
		CorrectBedgraph:	Produce paired-end bedgraphs.

	-h, --help	show this help message and exit

	For more info, please refer to the README file.
```

## Usage

Once the repository is downloaded and its /bin/ path is added to your PATH environment variable, CoCo can be called directly from your terminal.
CoCo is divided into two main run modes: **CorrectAnnotation** and **CorrectCount**, and one accessory mode: **CorrectBedgraph**.

### CorrectAnnotation
CorrectAnnotation should be used first in order to produce a modified annotation file from an input gene transfer format (.gtf) annotation file obtained from **Ensembl**.
Basic usage:
```
coco CorrectAnnotation path/to/your/annotation/hg38_annotation.gtf
```

The corrected annotation consists in a version of the original annotation where exon portions of host genes that overlap embedded genes are removed (see image below).

![alt tag](ressources/CorrectAnnotation.PNG)

The considered embedded gene biotypes are: snoRNA, scaRNA, snRNA, miRNA and tRNA. Therefore, **the original annotation file should have a "gene_biotype" entry for each genes**. This might not be included into the annotation files from other sources than Ensembl and so **we recommend that you get the annotation from Ensembl**.

### CorrectCount

CorrectCount is used to produce the corrected read counts per gene as well as their count per million (CPM) and approximated transcripts per million (TPM) values. This step should use the modified (.gtf) annotation file produced 


### CorrectBedgraph



## Running the tests

Explain how to run the automated tests for this system


### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```


## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Authors

* **Gabrielle Deschamps-Francoeur** - *Making of the multi-mapped read correction.*
* **Vincent Boivin** - *Making of the annotation correction for ambiguous reads.* - [boiv2803](http://gitlabscottgroup.med.usherbrooke.ca/u/boiv2803)
* **Michelle Scott** - *Mastermind of the project and awesome PI* - [scottgroup](http://scottgroup.med.usherbrooke.ca/)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to Fabien Dupuis-Sandoval for starting this project from scratch.
* Many thanks to Jean-Michel Garrant as well for his general wisdom.

## References
CoCo:
Subread package:
bedTools:

