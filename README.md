# PRTI
![PRTI](https://github.com/RuiyangZhai/img/blob/main/PRTI.png?raw=true)

`PRTI` (Predicting patients Response to Taxane-based chemotherapy and Immunotherapy) is an R package designed to predict tumor treatment response based on bulk and single-cell transcriptomics. We also provide an interactive online webserver with the same name: [PRTI](http://218.8.241.248:3838/PRTI/ "PRTI web").

## Installation

Please ensure that the dependent packages required by `PRTI` are installed correctly.
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("AUCell", "genefu", "GSVA", "pROC", "Seurat")
BiocManager::install(packages)
devtools::install_github("Lothelab/CMScaller")

devtools::install_github("RuiyangZhai/PRTI")
```

## Example
The `PRTI` package provides the bulk RNA-seq data required for the example code. We also provided the download link of the [scRNA-seq](https://github.com/RuiyangZhai/data/releases/download/PRTI/Seurat_obj.RData) and [stRNA-seq](https://github.com/RuiyangZhai/data/releases/download/PRTI/Spatial.RData) data required for the example code.

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
library(pROC)
library(ggplot2)
prediction = resultICI$prediction
prediction$Response = pdata_table$Response
roc.list = roc(Response ~ PredictICI + gene_up + gene_down, data = prediction,
               levels=c("NR","R"),direction="<")
ci.auc.list = lapply(roc.list, function(x){round(ci.auc(x),3)})

labels = paste0(names(ci.auc.list)," AUC:",unlist(lapply(ci.auc.list,function(x){x[2]})),
            "(95%CI:",unlist(lapply(ci.auc.list,function(x){x[1]})),"-",
            unlist(lapply(ci.auc.list,function(x){x[3]})),")")

cols = c("#EC3232", "#0787C3", "#F6944B")
p = ggroc(roc.list,linewidth=0.8)+
    theme_bw()+
    geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=1, color = "black")+
    coord_equal()+
    guides(colour=guide_legend(title=NULL))+
    theme(panel.grid=element_blank())+
    scale_color_manual(values = cols)

for(j in 1:length(labels)) {
  p = p + annotate("text", x = 0.55, y = 0.15-(j*0.05), label = labels[j],hjust = 0, size = 4)
} 
p
```

<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/ROC.png?raw=true" width="600">
</div>

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
colors = colorRampPalette(c("#ECC042","#4FB4BB","#43358E"))(10)
pheatmap(cor_table, 
         border_color="white",
         treeheight_row = 0,
         treeheight_col = 0,
         angle_col = 90,
         color = colors,
         breaks = seq(0, 1, by = 0.1),
         legend = TRUE)
```

<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/corHeatmap.png?raw=true" width="550">
</div>

### scRNA-seq data
```r
# Load R package and internal data
library(PRTI)
load("Seurat_obj.RData")

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
load("Spatial.RData")

# Calculate the PredictICI score
# The 'gsva', 'ssgsea' and 'AUCell' parameters can be selected for the 'method' parameter.
resultICI = scPredictICI(object = Spatial_obj,method = "AUCell",slot="data",asssy="Spatial")

# Plots
library(ggplot2)

Seurat::SpatialPlot(Spatial_obj,alpha = 0,image.alpha =1)+theme(legend.position = "none")
Seurat::SpatialDimPlot(Spatial_obj,repel = F, label = F,label.size = 5,image.alpha =0.5)
Seurat::SpatialFeaturePlot(Spatial_obj, features = c("PredictICI"),image.alpha =0.5)+theme(legend.position = "right")
```
![Spatial](https://github.com/RuiyangZhai/img/blob/main/Spatial.png?raw=true)

### CIsubtype
Example data can be downloaded from [here](https://github.com/RuiyangZhai/data/releases/download/PRTI/CIsubtype.RData).
```r
# Load R package and data
library(PRTI)
load("CIsubtype.RData")

#Predicting CIsubtypes
CIsubtypes = CIsubtype(gene_table,threshold = 0.2,doPlot = TRUE)

#Plots
data_plot<-as.data.frame(table(CIsubtypes$prediction,pdata_table$Response))
colnames(data_plot)[1:2] = c("Group","Response")
data_plot$Percent = (data_plot$Freq/sum(data_plot$Freq))*100

chisq_data = reshape2::dcast(data_plot,Response ~ Group, value.var = "Freq")
rownames(chisq_data) = chisq_data[,1]
chisq_data = chisq_data[,-1]

library(ggplot2)
ggplot(dfdata, aes(x = Group, y = Percent, fill = Response)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(aes(label = paste0(Percent, "%")), 
            position = position_stack(vjust = 0.5), size = 5) +
  labs(x = "Group", y = "Percent",
       title = "CI subtype",
       subtitle =paste0("Chisq, p = ",chisq.test(chisq_data)[["p.value"]])) +
  scale_fill_manual(values = c("R" = "#56B4E9","NR" = "#E69F00")) +
  theme_classic() +
  theme(
    axis.line = element_line(size = 1),
    text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))
```
<div align=center>
  <img src="https://github.com/RuiyangZhai/img/blob/main/CIsubtype.png?raw=true" width="600">
</div>

## Contact
Any technical question please contact Ruiyang Zhai (zhairuiyang@foxmail.com) or Te Ma (mate.compbio@foxmail.com).
