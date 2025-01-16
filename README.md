# PRTI
![PRTI](https://github.com/RuiyangZhai/img/blob/main/PRTI.png?raw=true)

`PRTI` (Predicting patients Response to Taxane-based chemotherapy and I mmunotherapy) is an R package designed to predict tumor treatment response based on bulk and single-cell transcriptomics. We also provide an interactive online webserver with the same name: [PRTI](http://218.8.241.248:3838/PRTI/ "PRTI web").

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
The `PRTI` package provides the bulk RNA-seq and scRNA-seq data required for the example code. And the stRNA-seq data required for the example code can be downloaded from [here](https://github.com/RuiyangZhai/data/releases/download/PRTI/Spatial.rds).

PredictICI and PredictTax are used in the same way. Just take PredictICI as an example.
### Bulk RNA-seq data
```r
# Load R package and internal data
library(PRTI)
data("gene_table",package = "PRTI")
data("pdata_table",package = "PRTI")


# Calculate the PredictICI score, where the Response parameter can be NULL
resultICI = PredictICI(expr = gene_table, Response = pdata_table$Response, verbose = T)

# ROC curve
prediction = resultICI$prediction
prediction$Response = pdata_table$Response
roc.list = pROC::roc(Response ~ PredictICI + gene_up + gene_down, data = prediction,
                     levels=c("NR","R"),direction="<")

ggroc(roc.list,linewidth=0.8)+
  theme_bw()+
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=1, color = "black")+
  coord_equal()+
  guides(colour=guide_legend(title=NULL))+
  theme(panel.grid=element_blank())+
  scale_color_manual(values = c("#EC3232", "#0787C3", "#F6944B"))
```

<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/ROC.png?raw=true" width="500">
</div>

### scRNA-seq data
```r
# Load R package and internal data
library(PRTI)
load(system.file("extdata", "Seurat_obj.RData", package = "PRTI"))

# Calculate the PredictICI score
# The 'gsva', 'ssgsea' and 'AUCell' parameters can be selected for the 'method' parameter.
resultICI = scPredictICI(object = Seurat_obj, method = "AUCell", slot="data", asssy="RNA")

# Plots
Seurat::FeaturePlot(resultICI, features = c("PredictICI","gene_up","gene_down"),ncol =3)
Seurat::VlnPlot(resultICI, features = c("PredictICI","gene_up","gene_down"),pt.size=0,ncol =3)
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/FeaturePlot.png?raw=true" width="800">
  <img src="https://github.com/RuiyangZhai/img/blob/main/VlnPlot.png?raw=true" width="800">
</div>

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

# > resultTLS
# $prediction
# Pt1         Pt2         Pt3         Pt4         Pt5         Pt6         Pt7         Pt8         Pt9        Pt10        Pt11 
# 0.38388322  0.08368658 -0.01756699 -0.24335755  0.29058169 -0.09588680 -0.04006636 -0.03644579 -0.29494496  0.58040490 -0.35744601 
# Pt12        Pt13 
# -0.28338173  0.03053978 
# 
# $auc
# [1] 0.6904762


# Calculate multiple signatures
Sig_list = c("PredictICI","CYT.Sig","GEP.Sig","PDL1.Sig","PD1.Sig","CTLA4.Sig","TLS.Sig",
             "IFN.Y.10.Sig", "IFN.Y.28.Sig", "HOT.score.Sig","IFNy4gene.Sig","IIS.Sig",
            "Inflamed20genes.Sig","Roh.immune.Sig","GBP2.Sig")
#Setting the 'sig_list' parameter to 'all' means that all signatures are calculated.
#Setting the parameter to 'Tax' or 'ICI' means that all signatures of the corresponding category are calculated.
results = calculateSig(expr = gene_table, Response = pdata_table$Response,Sig_list = Sig_list, verbose = T)

# Cor table
prediction = results$prediction
cor_table = cor(prediction,use="complete.obs",method="spearman")

library(pheatmap)
pheatmap(cor_table, 
         border_color="white",
         treeheight_row = 0,
         treeheight_col = 0,
         angle_col = 90,
         color = colors <- c("#ECC042","#C9BD5C","#A6BA77","#83B892","#60B5AD","#4DA5B5","#4B89AB","#486DA2","#455198","#43358E"),
         breaks = seq(0, 1, by = 0.1),
         legend = TRUE)
```

<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/corHeatmap.png?raw=true" width="550">
</div>

## Contact
Any technical question please contact Ruiyang Zhai (22b928040@stu.hit.edu.cn) or Te Ma (23b928040@stu.hit.edu.cn).
