# Tools_Expand_Genes_Network
 Tools used to analyze and expand gene networks of Vitis Vinifera and Human. All this tools used Python 3 and R.

## Prerequisites
  Python 3

  To install the required Python libraries and R packages:
  ```
   ./setup.sh
  ```

  In directory 'import_doc' is present the file 'vv_exprdata_2.csv'.  It is a large file uploaded with Git LFS that generate a pointer to the real file. To have the correct version of the file you need to:
   - download the repository with Git LFS, or,
   - get the file at this link (<!--insert link-->) and substitute to the file in 'import_doc' obtained with git clone.
<!-- Python 3
 Library of Python (install with pip3):
   * *datetime*
   ```
    pip3 install datetime
   ```
   * *pandas*
   ```
    pip3 install pandas
   ```
   * *rpy2*
   ```
    pip3 install rpy2
   ```
   * *matplotlib*
   ```
    pip3 install matplotlib
   ```
   * *matplotlib-venn*
   ```
    pip3 install matplotlib-venn
   ```
   * *networkx*
   ```
    pip3 install networkx
   ```
   * *numpy*
   ```
    pip3 install numpy
   ```
   * *scipy*
   ```
    pip3 install scipy
   ```
   * *rpack*
   ```
    pip3 install rectangle-packer
   ```


 For *biological_validation.py* is required R with [*topGO*](https://bioconductor.org/packages/release/bioc/html/topGO.html) and [DREME](http://meme-suite.org/doc/dreme.html) form MEME suite.
   Install *topGO* library:
   ```
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install("topGO")
   ```
   Install DREME: [DREME tool](http://meme-suite.org/doc/download.html)
-->

## TOOL 1: *managerList.py*
<!--Tool used to build the complete graph of interaction and find the LGN. In the human case study based on FANTOM we suggest to use all the isoforms of the genes for a complete analysis.
How to use:-->
Tool used to build the complete graph of interaction and find the LGN.
The output of this tool is the graphical representation (*.png*) of the complete graph discovered observing the file in input with corresponding legend (*.txt*). It is created the *.json* file compatible with Cytoscape. The tool generates a *.csv* file containing the textual representation of the graph, basically are written the edges of the graph. ManagerList returns the degree of the genes and the connected components of the graph.
The LGN, Local Gene Network of the previous graph, discovered using an adapted version of Charikar algorithm, corresponds to the *Core graph* files. The files for the LGN are the same of the complete graph previously described.
*ATTENTION*: In the human case study based on FANTOM we suggest to use all the isoforms of the genes for a complete analysis.
How to use:
```
Usage: python3 managerList.py PARAM [FILTERS]... -files [FILES]...
FILES can be a list of .csv or .zip. These are the expansion list from OneGenE
PARAM:
 '-vitis'  Lists of vitis
 '-human'  Lists of human. Need to be follow by:
  '-fantom'  Lists from Fantom DB
  '-TCGA'    Lists from TCGA DB
FILTERS:
 '-a'             Autosave image of graphs. If -a is present, it save automatically .png. USE IN MICROSOFT WINDOWS
 '-f [NUMBER]'    Ignored genes with frel<=NUMBER
 ONLY FOR VITIS GENES:
  '-t [PATTERN,...]'      Take genes that in 'functional annotation' or 'Network1' column there is at least one pattern
 ONLY FOR HUMAN GENES:
  '-i ['comp'/'notcomp']' Ignored edges between isoforms of same gene
```
Example with Cuticle genes:
```
python3 managerList.py -vitis -f 0.1 -a -files example_lists/Vitis_7genes_MYB-ERF/Example.zip
```
```
python3 managerList.py -human -fantom -f 0.5 -a -i comp -files example_lists/Human_22genes_cancer/Fantom/Fantom.zip
```
Help:
```
python3 managerList.py --help
```

## TOOL 2: *integrateCoupleGenes.py*
Tool used to expand the LGN adding to the graph the genes connected to the LGN with specific characteristic. The expansion can be done by:
  - apply a filter on the relative frequency;
  - apply a filter on the rank of the genes in the expansion list;
  - taking only the common genes;
  - apply a filter on the columns *Network1* and *Network2*
In output, the tool returns the graphical (*graphGenes.png*) and textual (*edges_graph.csv*) representation of the expanded graph  with the corresponding legend and the file *.json* to import in Cytoscape. Optionally, it can be obtained the histogram and the venn diagram. The histogram represents the genes selected for the expansion (blue columns) respect all the genes present in the expansion list of the starting gene (green columns). The venn diagram (graphical representation up to 5 genes in LGN, textual for any number of genes in LGN) represents how the genes in expansion list are shared between the LGN.
How to use:
```
Usage: python3 integrateCoupleGenes.py PARAM TYPEA [FILTERS]... -files [GENES] [FILES]... [ISOFORM]
PARAM:
 '-vitis'  Lists of vitis
 '-fantom' Lists of fantom DB
 '-TCGA'   Lists of TCGA DB
TYPEA:
 '-frel'        Build expansion network based on FREL. Required filter '-f'
 '-rank [INT]'  Build expansion network based on RANK. Take top genes
 '-shared'      Build expansion network based on SHARED GENES
FILTERS:
 '-a'           Autosave image of graphs. If -a is present, it save automatically .png. USE IN MICROSOFT WINDOWS
 '-c'           Add edges between associated genes
 '-e'           Print Venn Diagram and Histogram for complete analysis
 '-f [NUMBER]'  Ignored genes with frel<=NUMBER
GENES: file .csv with the genes to analyze. Example: ('CoupleGeneToIntegrate/coupleGene.csv')
FILES can be a list of .csv or .zip
ISOFORM: file .csv from the execution of ManagerList.py with the composition of edge gene-gene (tag '-i comp'). To use only with '-fantom'
```
Example with Cuticle genes:
```
python3 integrateCoupleGenes.py -vitis -shared -f 0.1 -e -a -files CoupleGeneToIntegrate/coupleGene0.csv example_lists/Vitis_7genes_MYB-ERF/Example.zip
```
```
python3 integrateCoupleGenes.py -fantom -shared -a -e -f 0.1 -files CoupleGeneToIntegrate/coupleGene.csv example_lists/Human_22genes_cancer/Fantom/Fantom.zip networkOutput/*dateExecutionTool1*/isoformInEdge.csv
```
Help:
```
python3 integrateCoupleGenes.py --help
```

## TOOL 3: *biological_validation.py*
Tool that execute biological validation. It perform two analysis:
  - Computes a Gene Ontology analysis on set of genes.
  - Prepare the fasta file for the execution of XSTREME tool.
Only for Vitis Vinifera.
How to use:
```
Usage: python3 biological_validation.py PARAM [FILTERS]... LIST_GENES COMPLETE_GENOME
PARAM:
 '-topGO' Execute GO validation
 '-dreme' Execute DREME analysis
FILTERS (ONLY FOR DREME):
 '-exe'  Execute dreme analysis in local. Without prepare fasta file for DREME website
LIST_GENES      List of genes in .csv file.
COMPLETE_GENOME COmplete file information for validation
 * For topGO: file map from vitis ID to GO ID. ('import_doc/V1_GOcomplete.txt')
 * For Dreme: complete list of genes in genome. ('import_doc/grape_1k_upstream.fasta')
```
Example topGO with Cuticle genes:
```
python3 biological_validation.py -topGO example_lists/Vitis_7genes_MYB-ERF/topGO_Vitis/MYB_ERF_topGO0.csv import_doc/V1_GOcomplete.txt
```
Example DREME with Cuticle genes:
```
python3 biological_validation.py -dreme example_lists/Vitis_7genes_MYB-ERF/topGO_Vitis/MYB_ERF_topGO0.csv import_doc/grape_1k_upstream.fasta
```
Help:
```
python3 biological_validation.py --help
```
