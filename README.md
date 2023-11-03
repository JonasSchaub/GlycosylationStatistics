[![DOI](https://zenodo.org/badge/284081201.svg)](https://zenodo.org/doi/10.5281/zenodo.7081511)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub issues](https://img.shields.io/github/issues/JonasSchaub/SugarRemoval.svg)](https://GitHub.com/JonasSchaub/SugarRemoval/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/JonasSchaub/SugarRemoval.svg)](https://GitHub.com/JonasSchaub/SugarRemoval/graphs/contributors/)
[![GitHub release](https://img.shields.io/github/release/JonasSchaub/SugarRemoval.svg)](https://github.com/JonasSchaub/SugarRemoval/releases/)

# Description and Analysis of Glycosidic Residues in the Largest Open Natural Products Database
##### Code for automated, systematic detection of sugar moieties in the COlleCtion of Open Natural prodUcTs (COCONUT) database

## Description
This repository contains Java source code for automatically detecting and analysing glyosidic moieties <i>in silico</i> in 
the largest open natural products database COCONUT, 
as described in [Schaub, J.; Zielesny, A.; Steinbeck, C.; Sorokina, M. Description and Analysis of Glycosidic Residues in the Largest 
Open Natural Products Database. Biomolecules 2021, 11, 486.](https://doi.org/10.3390/biom11040486),
using the [Sugar Removal Utility](https://doi.org/10.1186/s13321-020-00467-y).
<br>Additionally, similar analyses are done with datasets from the ZINC15 database, DrugBank, and ChEMBL.
<br>Python scripts and Jupyter notebooks for the curation of some used datasets and analysing and visualising the results 
are also supplied in this repository.
<p>
Please NOTE that the code in this repository is primarily supposed to show how the glycosylation statistics published
in the article linked above were generated and to allow reproduction and executing of the same analyses for other
datasets. It is not considered a software by itself. Hence, things like the publication of a Maven artifact for 
straightforward installation are not given here.
<br>

The [Sugar Removal Utility](https://github.com/JonasSchaub/SugarRemoval), however, can be installed as a Maven artifact in
a straightforward manner and used in your own scripts and workflows to analyse other datasets this way.

</p>

## Contents
### Source code for glycosylation statistics analysis
In the directory <i>/src/test/java/de/unijena/cheminf/deglycosylation/stats/</i> the class 
<i>GlycosylationStatisticsTest</i> can be found. It is a JUnit test class with multiple test methods that can be run in 
a script-like fashion to do the various analyses. Using an IDE like e.g. IntelliJ is recommended. Please note that some 
directories etc. will need to be adjusted and some datasets be put into the <i>/src/test/resources/</i> directory 
(see below) to run the tests yourself.<p>

The directory <i>/Python_scripts_and_notebooks/</i> contains a python script for picking a diverse subset of a larger 
datasets using the [RDKit MaxMin algorithm](http://www.rdkit.org/docs/GettingStartedInPython.html#picking-diverse-molecules-using-fingerprints).
For the reported analyses, it has been used to reduce in size the downloaded ZINC "in-vitro" subset while preserving 
diversity. Additionally, two Jupyter Notebooks can be found in this directory that have been used to analyse and visualise 
some of the test results. 

## Installation
This is a Maven project. In order to do the described analyses on your own, download or clone the repository and
open it in a Maven-supporting IDE (e.g. IntelliJ) as a Maven project and execute the pom.xml file. Maven will then take
care of installing all dependencies.
<br>To run the COCONUT-analysing tests, a MongoDB instance needs to be running on your platform and the COCONUT NP 
database imported to it. The respective MongoDB dump can be downloaded at 
[https://coconut.naturalproducts.net/download](https://coconut.naturalproducts.net/download).
<br>To run the Python scripts and Jupyter Notebooks, installing [Anaconda](https://www.anaconda.com) is recommended, to 
also ease the installation of required libraries, like the open-source cheminformatics software [RDKit](http://www.rdkit.org). 

## Required datasets
* **COCONUT**: To run the COCONUT-analysing tests, a MongoDB instance needs to be running on your platform and the COCONUT NP
  database imported to it. The respective MongoDB dump can be downloaded at
  [https://coconut.naturalproducts.net/download](https://coconut.naturalproducts.net/download). Please check and adjust 
  the credentials for the connection to MongoDB in the code and adjust them if needed. One test method also 
  analyses COCONUT in the form of an SDF. This file can also be obtained from the given webpage and needs to be placed 
  in the <i>/src/test/resources/</i> directory.
* **ZINC15**: A list of available ZINC15 subsets can be found [here](http://zinc15.docking.org/substances/subsets/). It is 
  recommended to use the program [wget](https://www.gnu.org/software/wget/) to download the subsets. All subsets were 
  downloaded as SMILES files.
    * ZINC "for-sale": A part of the ZINC "for-sale" subset was downloaded for the published analyses and further reduced
      in size using the <i>ZINC_for-sale_curation.py</i> script located in the <i>/Python_scripts_and_notebooks/</i> directory. 
      One test method curates the dataset further. After this is done, the curated datasets needs to be placed in the
      <i>/src/test/resources/</i> directory for it to be analysed by other test methods.
    * ZINC "in-vitro": One test method curates the dataset. After this is done, the curated datasets needs to be placed in the
      <i>/src/test/resources/</i> directory for it to be analysed by other test methods.
    * ZINC "biogenic": The ZINC "biogenic" dataset needs to be placed in the <i>/src/test/resources/</i> directory to be 
      used for the curation of the other datasets.
* **Manually curated review of bacterial natural products sugar moieties**: Two of the test methods do a substructure search 
  in COCONUT for sugar moieties reported in bacterial natural products, manually curated by 
  [Elshahawi et al.](https://doi.org/10.1039/C4CS00426D). This dataset is already supplied in this repository in the 
  <i>/src/test/resources/</i> directory. 
* **ChEMBL**: The ChEMBL 28 database is curated in one test method and analysed for glycosidic moieties in another. To run the
  curation test, the dataset has to be placed in the <i>/src/test/resources/</i> directory as an SDF. After curation, the 
  curated dataset has to be placed in the same directory.
* **DrugBank**: The DrugBank "all structures" dataset is curated in one test method and analysed for glycosidic moieties in another. To run the
  curation test, the dataset has to be placed in the <i>/src/test/resources/</i> directory as an SDF. After curation, the
  curated dataset has to be placed in the same directory.

## Dependencies
* Java Development Kit (JDK) version 17
    * [AdoptOpenJDK](https://adoptopenjdk.net) (as one possible source of the JDK)
* Chemistry Development Kit (CDK) version 2.8
    * [Chemistry Development Kit on GitHub](https://cdk.github.io/)
* Apache Maven version 4
    * [Apache Maven](http://maven.apache.org)
* JUnit version 4.13.2
    * [JUnit 4](https://junit.org/junit4/)
* Java MongoDB Driver version 3.11.1
    * [Java MongoDB Driver documentation](https://docs.mongodb.com/drivers/java/)
* Python version 4.7.1
    * [Anaconda](https://www.anaconda.com) (as a recommended source for the python interpreter and many more)
* RDKit: Open-Source Cheminformatics Software version 2020.03
    * [RDKit homepage](http://www.rdkit.org)

## References and useful links
**Glycosylation statistics of COCONUT publication**
* [Schaub, J., Zielesny, A., Steinbeck, C., Sorokina, M. Description and Analysis of Glycosidic Residues in the Largest Open Natural Products Database. Biomolecules 2021, 11, 486. https://doi.org/10.3390/biom11040486](https://doi.org/10.3390/biom11040486)

**Sugar Removal Utility**
* [Schaub, J., Zielesny, A., Steinbeck, C., Sorokina, M. Too sweet: cheminformatics for deglycosylation in natural products. J Cheminform 12, 67 (2020). https://doi.org/10.1186/s13321-020-00467-y](https://doi.org/10.1186/s13321-020-00467-y)
* [SRU Source code](https://github.com/JonasSchaub/SugarRemoval)
* [Sugar Removal Web Application](https://sugar.naturalproducts.net)
* [Source Code of Web Application](https://github.com/mSorok/SugarRemovalWeb)

**Chemistry Development Kit (CDK)**
* [Chemistry Development Kit on GitHub](https://cdk.github.io/)
* [Steinbeck C, Han Y, Kuhn S, Horlacher O, Luttmann E, Willighagen EL. The Chemistry Development Kit (CDK): An Open-Source Java Library for Chemo- and Bioinformatics. J Chem Inform Comput Sci. 2003;43(2):493-500.](https://dx.doi.org/10.1021%2Fci025584y)
* [Steinbeck C, Hoppe C, Kuhn S, Floris M, Guha R, Willighagen EL. Recent Developments of the Chemistry Development Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics. Curr Pharm Des. 2006; 12(17):2111-2120.](https://doi.org/10.2174/138161206777585274)
* [May JW and Steinbeck C. Efficient ring perception for the Chemistry Development Kit. J. Cheminform. 2014; 6:3.](https://dx.doi.org/10.1186%2F1758-2946-6-3)
* [Willighagen EL, Mayfield JW, Alvarsson J, Berg A, Carlsson L, Jeliazkova N, Kuhn S, Pluska T, Rojas-Chertó M, Spjuth O, Torrance G, Evelo CT, Guha R, Steinbeck C, The Chemistry Development Kit (CDK) v2.0: atom typing, depiction, molecular formulas, and substructure searching. J Cheminform. 2017; 9:33.](https://doi.org/10.1186/s13321-017-0220-4)
* [Groovy Cheminformatics with the Chemistry Development Kit](https://github.com/egonw/cdkbook)

**COlleCtion of Open NatUral producTs (COCONUT)**
* [COCONUT Online home page](https://coconut.naturalproducts.net)
* [Sorokina, M., Merseburger, P., Rajan, K. et al. COCONUT online: Collection of Open Natural Products database. J Cheminform 13, 2 (2021). https://doi.org/10.1186/s13321-020-00478-9](https://doi.org/10.1186/s13321-020-00478-9)
* [Sorokina, M., Steinbeck, C. Review on natural products databases: where to find data in 2020. J Cheminform 12, 20 (2020).](https://doi.org/10.1186/s13321-020-00424-9)

**ZINC** 
* [ZINC15 Homepage](http://zinc15.docking.org)
* [Sterling and Irwin, J. Chem. Inf. Model, 2015 http://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00559](https://doi.org/10.1021/acs.jcim.5b00559)

**MongoDB**
* [MongoDB homepage](https://www.mongodb.com)
* [Java MongoDB Driver documentation](https://docs.mongodb.com/drivers/java/)

**RDKit**
* [RDKit homepage](http://www.rdkit.org)
* [Getting started with the RDKit in Python](http://www.rdkit.org/docs/GettingStartedInPython.html)

**DrugBank**
* [DrugBank homepage](https://go.drugbank.com)
* [Wishart DS, Feunang YD, Guo AC, Lo EJ, Marcu A, Grant JR, Sajed T, Johnson D, Li C, Sayeeda Z, 
  Assempour N, Iynkkaran I, Liu Y, Maciejewski A, Gale N, Wilson A, Chin L, Cummings R, Le D, Pon A, 
  Knox C, Wilson M. DrugBank 5.0: a major update to the DrugBank database for 2018. Nucleic Acids Res. 
  2017 Nov 8. doi: 10.1093/nar/gkx1037.](https://doi.org/10.1093/nar/gkx1037)

**ChEMBL**
* [ChEMBL homepage](https://www.ebi.ac.uk/chembl/)
* [Gaulton A, Hersey A, Nowotka M, Bento AP, Chambers J, Mendez D, Mutowo P, Atkinson F, Bellis LJ, Cibrián-Uhalte E, 
  Davies M, Dedman N, Karlsson A, Magariños MP, Overington JP, Papadatos G, Smit I, Leach AR. (2017) 'The ChEMBL 
  database in 2017.' Nucleic Acids Res., 45(D1) D945-D954.](http://dx.doi.org/10.1093/nar/gkw1074)