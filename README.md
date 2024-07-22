# ImmunQTL analysis
## -Map genetic effects on the cellular composition of tumor microenvironment 

[![N|Solid](https://cldup.com/dTxpPi9lDf.thumb.png)](https://nodesource.com/products/nsolid)

[![Build Status](https://travis-ci.org/joemccann/dillinger.svg?branch=master)](https://travis-ci.org/joemccann/dillinger)

This repository contains the code used to perform the analysis in the "Interpreting genetic control of cellular infiltration and its distinct role in human tumor microenvironment" paper.
The purpose of making the code available is for transparency and data analysis reproducibility.

## Overview

Tumor microenvironment (TME) is a complex assembly of tumor, immune, stromal, and extracellular components. Organization of these components at the cellular and tissue levels plays a crucial role in the effectiveness of antitumor immunity. Single-cell RNA sequencing (scRNA-seq) provides an unprecedented resolution to decipher the heterogeneous populations in TME, allowing identification of novel cell-types and discovery of unknown associations. With the extension of CIBERSORTx, we deconvolved TCGA bulk RNA-seq with the reference profiles generated from scRNA-seq data of the corresponding tumor type. Then, we performed a pan-cancer evaluation of the contribution of germline variation on tumor microenvironment composition, and identified immunity quantitative trait loci (immunQTLs) across 23 cancer types from TCGA. Furthermore, we identified associations between cell fractions and gene expression, and explored dynamic interactions between tumor and immune cells.


## Workflow

The analyses contained in our work:

- [Process the scRNA-seq data ] - The work flow to process the scRNA-seq, including quality control, doublet removal, cell clusering and cell annotation.
- [Run CIBERSORTx and assess their deconvolution results] - Apply two additional deconvolution methods, MuSiC and Bisque to compare the results.
- [Calculate immunQTLs] - Assess the genetic effect on cellular composition with the adjustment of covariates.
- [Infer intercellular communication using CellChat] - Predict cell-cell contact interactions,ECM-receptor interactions and signaling interactions.
- [Infer intercellular communication using NicheNet] - Predict the regulatory activity of ligands.
- [Explore gene-cell association] - Prioritize genes mediate immunQTL effect.
- [Exhibit the finding with plots] -Contain the code related to figure creation.



## Data source

All accession codes or web links for publicly available scRNA-seq datasets are described in the Supplementary Tables. RNA-seq profiles (level 3) and genotype data (level 2) of 7,707 samples were obtained from the TCGA data portal (https://portal.gdc.cancer.gov/). RNA-seq data conducted in tumor tissues and adjacent normal tissues from our 154 CRC samples have been deposited in the Genome Sequence Archive (GSA) under accession number HRA007986, which contains RAW sequencing data. Genotype data for our 154 CRC samples have been deposited in GSA under accession number GVM000801. GWAS summary statistics for 124 traits in UK Biobank were obtained from pan-cancer GWAS (https://github.com/Wittelab/pancancer_pleiotropy) and gwasATLAS (https://atlas.ctglab.nl/). Pan-cancer GWAS summary statistics for Asian was available through the Biobank Japan (http://jenger.riken.jp/en/). ChIP-seq peaks, TF-binding sites, and eCLIP-seq data among human cancer cell lines were downloaded from the ENCODE portal (https://www.encodeproject.org/data/annotations/). Expression profile data and drug sensitivity data of cancer cell lines were obtained from the GDSC (https://www.cancerrxgene.org/). 

## License

Miaolab

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [dill]: <https://github.com/joemccann/dillinger>
   [git-repo-url]: <https://github.com/joemccann/dillinger.git>
   [john gruber]: <http://daringfireball.net>
   [df1]: <http://daringfireball.net/projects/markdown/>
   [markdown-it]: <https://github.com/markdown-it/markdown-it>
   [Ace Editor]: <http://ace.ajax.org>
   [node.js]: <http://nodejs.org>
   [Twitter Bootstrap]: <http://twitter.github.com/bootstrap/>
   [jQuery]: <http://jquery.com>
   [@tjholowaychuk]: <http://twitter.com/tjholowaychuk>
   [express]: <http://expressjs.com>
   [AngularJS]: <http://angularjs.org>
   [Gulp]: <http://gulpjs.com>

   [PlDb]: <https://github.com/joemccann/dillinger/tree/master/plugins/dropbox/README.md>
   [PlGh]: <https://github.com/joemccann/dillinger/tree/master/plugins/github/README.md>
   [PlGd]: <https://github.com/joemccann/dillinger/tree/master/plugins/googledrive/README.md>
   [PlOd]: <https://github.com/joemccann/dillinger/tree/master/plugins/onedrive/README.md>
   [PlMe]: <https://github.com/joemccann/dillinger/tree/master/plugins/medium/README.md>
   [PlGa]: <https://github.com/RahulHP/dillinger/blob/master/plugins/googleanalytics/README.md>
