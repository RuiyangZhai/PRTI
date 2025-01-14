# __PRTI__
![PRTI](http://218.8.241.248:3838/PRTI/image/PRTI.png)

`PRTI` (Predicting patients Response to Taxane-based chemotherapy and I mmunotherapy) is an R package designed to predict tumor treatment response based on bulk and single-cell transcriptomics.We also provide an interactive online webserver with the same name: [PRTI](http://218.8.241.248:3838/PRTI/ "PRTI web").

## Installation

Please ensure that the dependent packages required by `PRTI` are installed correctly.
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("AUCell", "genefu", "GSVA", "pROC", "Seurat")
BiocManager::install(packages)

devtools::install_github("RuiyangZhai/PRTI")
```

## Example
The `PRTI` package provides the bulk RNA-seq and scRNA-seq data required for the example code. And the stRNA-seq data required for the example code can be downloaded from [here](https://www.alipan.com/s/gXX6keVexSg).

PredictICI and PredictTax are used in the same way. Just take PredictICI as an example.
### Bulk RNA-seq data
```r
# Load R package and internal data
library(PRTI)
data("gene_table",package = "PRTI")
data("pdata_table",package = "PRTI")


# Calculate the PredictICI score, where the Response parameter can be NULL
resultICI = PredictICI(expr = gene_table, Response = pdata_table$Response, verbose = T)

#画图代码
```
### scRNA-seq data
```r
# Load R package and internal data
library(PRTI)
load(system.file("extdata", "Seurat_obj.RData", package = "PRTI"))

# Calculate the PredictICI score
# The 'gsva', 'ssgsea' and 'AUCell' parameters can be selected for the 'method' parameter.
resultICI = scPredictICI(object = Seurat_obj, method = "AUCell", slot="data", asssy="RNA")

#画图代码
```
### stRNA-seq data
```r
# Load R package and data
library(PRTI)
Spatial_obj = readRDS("../Spatial.rds")

# Calculate the PredictICI score
# The 'gsva', 'ssgsea' and 'AUCell' parameters can be selected for the 'method' parameter.
resultICI = scPredictICI(object = Spatial_obj, method = "AUCell", slot="data", asssy="Spatial")

#画图代码
```
### Published signatures
`PRTI` also provides functions for calculating published signatures.
```r
# Load R package and internal data
library(PRTI)
data("gene_table",package = "PRTI")
data("pdata_table",package = "PRTI")

# Calculate a single signature
resultTLS = TLS.Sig(expr = gene_table, Response = pdata_table$Response, verbose = T)

# Calculate multiple signatures
Sig_list = c('PredictICI','PD1.Sig','GEP.Sig','CYT.Sig')
#Setting the 'sig_list' parameter to 'all' means that all signatures are calculated.
#Setting the parameter to 'Tax' or 'ICI' means that all signatures of the corresponding category are calculated.

results = calculateSig(expr = gene_table, Response = pdata_table$Response,Sig_list = Sig_list, verbose = T)

#画图代码
```

## Contact
Any technical question please contact Ruiyang Zhai (22b928040@stu.hit.edu.cn) or Te Ma (23b928040@stu.hit.edu.cn).
