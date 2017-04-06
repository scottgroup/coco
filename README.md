# CoCo: Count Corrector for embedded and multi-mapped genes.

CoCo is a pipeline designed to improve the evaluated abundances of genes from RNA-seq alignment data. CoCo is divided in two main parts which can be used together or seperately, depending on user's preferences. The first part being a correction of the gene annotation used by read assignment softwares such as featureCounts or HTSeq in order to correct the evaluated read counts for embedded genes such as snoRNA, that overlap features of their host gene's transcripts such as retained introns and exons. The second part of the correction distributes multi mapped reads in relation to the evaluated read counts obtained from single mapped reads.

We've shown that the combination of these two corrections gives the best correlation of the evaluated gene abundance from RNAseq data with abundances evaluated from other methods such as qPCR.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

CoCo has a few python and R dependencies.

```
Give examples
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

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

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Gabrielle Deschamps-Francoeur** - *Making of the multi-mapped read correction.
* **Vincent Boivin** - *Making of the annotation correction for ambiguous reads.* - [boiv2803](http://gitlabscottgroup.med.usherbrooke.ca/u/boiv2803)
* **Michelle Scott** - *Mastermind of the project and awesome PI* - [scottgroup](http://scottgroup.med.usherbrooke.ca/)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to Fabien Dupuis-Sandoval for starting this project from scratch
* Many thanks to Jean-Michel Garrant as well for his general wisdom.


