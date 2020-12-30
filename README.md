# explicit

**Got a list of genes and you wonder what TFs regulate their expression? EXPLICIT is the right tool to try.**

The EXPLICIT approach has been developed to construct a gene expression predictor model for the plant species *Arabidopsis thaliana*. The predictor uses the expression of 1,678 transcription factor (TF) genes to predict the expression of 29,182 non-TF genes. It further enables downstream inference of TF regulators for genes and gene modules functioning in diverse plant pathways. Please check the [original paper](https://github.com/MaShisongLab/explicit#Reference) by Geng *et al.* for more details. The EXPLICIT package presented here enables users to 1. Infer TF regulators for their own gene modules; 2. Draw chord diagrams showing TF-target genes regulation for the modules; 3. Create custom gene expression predictors using their own gene expression data. (*Note: below is an example showing the analysis flow-chart for a gene module involved in vascular system development.*)

<a><img alt="The Analysis Work-flow" src="./data/working_flow.png" align="center" width="700" ></a>

## Table of Contents
- [Install](https://github.com/MaShisongLab/explicit#Install)
- [Usage](https://github.com/MaShisongLab/explicit#Usage)
   - Infer TF regulators for gene modules
   - Draw chord diagrams showing TF-target genes regulation for the modules
   - Create custom gene expression predictor
- [Additional Information](https://github.com/MaShisongLab/explicit#additional-information)
- [Reference](https://github.com/MaShisongLab/explicit#Reference)

## Install
This package requires [Perl](https://www.activestate.com/products/perl/downloads/), [R](https://www.r-project.org/), and the [circlize](https://www.rdocumentation.org/packages/circlize/) package in R. 

`circlize` can be installed within an R console via the command:

```R
install.packages("circlize")
```
[MATLAB](https://www.mathworks.com/products/matlab.html) is optional. Only required if you want to create your own predictor model using custom expression data. 

Once the required software is installed, just download or clone the whole package to a local computer and start using it from the package's home directory. 

## Usage

### 1. Infer TF regulators for gene modules

#### a. Prepare the module file
The file used to store gene modules information is "modules_to_analyze.txt". Edit the file with your own modules. The file has the following format, with the first column being gene ids and the second column being module names. The two columns are separated by a tab. For gene ids, only standard Arabidopsis AGI ids are currently supported. Multiple modules can be analyzed at the same time. Once finished editing, save the file without changing its name.
```
Gene_Name   ModuleID
AT1G25360   Module138
AT2G22340   Module138
AT5G75660   Moudle138
AT2G22130   Module139
AT4G12350   Module139
.........   ........
```
#### b. Conduct enrichment assay to identify TF regulators for the modules
The Perl script "getArabidopsisRegulatorTFs.pl" will do the job. It takes the modules from the file "modules_to_analyze.txt" to conduct an enrichment assay to identify potential TF regulators. Results are outputted to a file named "results.regulator.tfs.txt", which can be open in EXCEL.

The command to use:
```shell
perl getArabidopsisRegulatorTFs.pl
```
Here is an example of the output results:
![alt-text](./data/results_sample.png "Output result example")


### 2. Draw chord diagrams showing TF-target genes regulation for the modules
#### a. Obtain the chord-list file for a module of interest
The Perl script "getChordLists.pl" will extract the TF-target gene pairs from the "results.regulator.tfs.txt" for the module specified. By default, it will take the top 50 TFs and top 15 target genes. The results are outputted to a file named "chord.lists.txt", which will be used in the next step to draw a chord diagram. 
```shell
perl getChordLists.pl XXXXX
```
Replace `XXXXX` with the name of the module.
#### b. Draw the chord diagram according to the chord-list in R, using the `circlize` package
The `circlize` package in R will be used to draw the chord diagram showing the TF-target genes interaction, as specified in the "chord.lists.txt". Open an R console and navigate to the home directory of the explicit package, which contains the "chord.lists.txt" file. Within the R console, type the following commands:
```R
source("Rscripts.R")
library("circlize")
drawChordDiagram(chordfile = "chord.lists.txt", ratio = 1)
drawChordDiagram(chordfile = "chord.lists.txt", ratio = 0.6)
```
Ratio specifies the relative size of the target gene area occupies. You can repeat step a and b to draw diagrams for another module.

#### c. Use a single command in R to draw the diagram
One can also directly issue the following command within a R console to draw a Chord diagram for a module:
```R
source("Rscripts.R")
library("circlize")
directChordDiagram( module="XXXXX", ratio = 1, tfnum = 50, targetnum = 15)
```
`tfnum` and `targetnum` specify, respectively, the maximum number of TFs and target genes to be included within the chord Diagram.

### 3. Create custom gene expression predictor 
Currently, we have only the gene expression predictor model for *Arabidopsis thaliana*. We are working on predictor models for other species. At the same time, you can also create your won custom gene expression predictor. However, a large number of training samples are required for training the model. The number should be at least 5 - 10 times larger than the number of input TFs.

The MATLAB function <B>explicit</B>, as specified within the file "explicit.m", will be used to create the model. The file can be found within the data folder. MATLAB is required for this analysis.
```matlab
mdl = explicit( TF_expression, TG_expression, TF_name, TG_name)
```
`TF_expression`: the expression matrix for TF, with rows representing samples and columns representing genes <br>
`TG_expression`: the expression matrix for target genes, with rows representing samples and columns representing genes <br> 
`TF_name`: the names of the TF genes <br>
`TG_name`: the names of the target genes <br>

## Additional Information
Will update soon.

## Reference

Will update soon.
