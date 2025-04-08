#'Predicting the sensitivity of cancer patients to immune checkpoint inhibitors (ICI).
#'
#'This function is based on single sample gene set enrichment analysis (ssGSEA), which uses gene expression values to calculate the sensitivity score of samples to immune checkpoint inhibitors (ICI).
#' @importFrom pROC roc
#'
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step.
#' @param ... Other arguments passed on to gsva()'s params argument.
#' @export
#'
#' @examples
#' res = PredictICI(expr = gene_table, Response = pdata_table$Response)
#'
PredictICI <- function(expr,Response=NULL,verbose=FALSE,...) {
  gene_up=c('IFNG','LAG3','GBP1','CXCL9','TBX21','IGHV3-48','TRGC1','FASLG','APOL6',
            'PLA2G2D','IGHV3-11','JAKMIP1','IGLV3-21','CXCL11','PSME1','GBP4','CTLA4',
            'CD274','FBXO6','PSMB10','IGHV4-59','IGLV3-1','IL4I1','TAP1','IGKV3-11',
            'KLRC2','KLRK1','IGKV1-12','IL21','UBE2L6','IGHV4-34','H2AC21','CAPG',
            'P2RX5','PSMB9','TNIP3','CEP128','MYO7A','IL32','OASL','CLEC6A','KLRC4-KLRK1',
            'GBP2','AIM2','OAS3','FCRL6','HAPLN3','GPR25','H2BC10','FCRL3','GFI1','H1-3')
  gene_down=c('HSPG2','OGN','NINL','EFNA1','SMO','BCAM','EIF3E','ITGB8','PRKG1',
              'NFE2L2','NDST1','RPL23','TFF3','PARVA','EXTL2','CXADR','HSD11B2',
              'RPS13','SEMA6A','SMAD5','ID4','CD24','PBX1','COL4A5','PLXNB1','IGF1R',
              'PRKAB1','HEYL','JUP','CSRP1','SORBS1','TPT1','GNG12','WLS','NECTIN1',
              'MAP2','SOX6','NID1','EXT1','SGCD','SDC2','PLA2G12A','TSC22D1','PHLDB1',
              'DYNC2H1','NFIA','PTPRZ1','SORT1')

  predict_result = .predict_sig(expr = expr,gene_up=gene_up,gene_down = gene_down,
                                predict_name = "ICI",method = "ssGSEA",verbose=verbose,
                                nCores = 1,normAUC = TRUE,...)

  auc_value = NA
  if (!is.null(Response)) {
    auc_result = pROC::roc(Response,predict_result$PredictICI,levels=c("NR","R"),direction="<")
    auc_value = as.numeric(auc_result$auc)
    if (verbose) {
      print(auc_result$auc)
    }
  }
  out_put = list(prediction = predict_result,auc = auc_value)
  return(out_put)
}


#'Predicting the sensitivity of cancer patients to Taxanes.
#'
#'This function is based on gene set variation analysis (GSVA), which uses gene expression values to calculate the sensitivity score of samples to Taxanes.
#' @importFrom pROC roc
#'
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step.
#' @param ... Other arguments passed on to gsva()'s params argument.
#' @export
#'
#' @examples
#' res = PredictTax(expr = gene_table, Response = pdata_table$Response)
#'
PredictTax <- function(expr,Response=NULL,verbose=FALSE,...) {
  gene_up=c('GBP1','GTSE1','CTSC','FBXO5','GMNN','GBP5','NUP153','CXCL13','PDIA6',
            'LYZ','SNRPA','PLA2G7','CASP8AP2','TCF3','ZBED2','BATF2','TFDP1','IGKV4-1',
            'ELF5','BTN2A1','TFDP2','H3C3','GTF2E2','SMARCD2','SLC35B1','MGAT1','MED1',
            'CRLF1','MAP3K3')
  gene_down=c('SUSD4','S100G','UNC5B','SDC4','PPP2R2C','EIF4G3','TACC2','RPS4Y1',
              'LAMP5','CTTN','FAN1','HGSNAT','KIF3A','PLA2G12A','DOCK1','ULK2','ARHGEF12',
              'FLNB','CGA','ERLIN2','CELSR2','BAIAP3','HMGCL','PURA','FBXL16','SPATA18',
              'BHLHE40','GPRC5A','SLC7A2','MATN3')

  predict_result = .predict_sig(expr = expr,gene_up=gene_up,gene_down = gene_down,
                                predict_name = "Tax",method = "GSVA",verbose=verbose,
                                nCores = 1,normAUC = TRUE,...)

  auc_value = NA
  if (!is.null(Response)) {
    auc_result = pROC::roc(Response,predict_result$PredictTax,levels=c("NR","R"),direction="<")
    auc_value = as.numeric(auc_result$auc)
    if (verbose) {
      print(auc_result$auc)
    }
  }
  out_put = list(prediction = predict_result,auc = auc_value)
  return(out_put)
}


#' Function to compute the signatures.
#'
#' This function computes signatures scores from gene expression values.
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step.
#' @param Sig_list The signatures that need to be calculated.Setting the parameter to `all` means calculating all signatures.Setting the parameter to `Tax` or `ICI` means calculating all signatures for the corresponding category.
#' @export
#' @examples
#' Sig_list = c('GEP.Sig','DDIR.Sig','HOT.score.Sig','CYT.Sig')
#' res = calculateSig(expr = gene_table, Response = pdata_table$Response, Sig_list = Sig_list)
#'
calculateSig <- function(expr,Response=NULL,Sig_list,verbose=FALSE) {
  Tax_sig = c('PredictTax','TP53.Sig','Liu.2018.Sig','MammaPrint.Sig','GGI.Sig',
              'DDIR.Sig','RB.Sig','CIN70.Sig','Bai.2022.Sig','Zhu.2020.Sig',
              'CES.Sig','Birkbak.2018.Sig','RAD51.Sig','RB.loss.Sig','TYMS.Sig',
              'Tegafur_uracil.Sig','CD96.Sig','Juul.2010.Sig','Chen.2021.Sig','Weng.2022.Sig',
              'ILB.Sig','Oncotype.DX.Sig')

  ICI_sig = c('PredictICI','GEP.Sig','GBP2.Sig','Roh.immune.Sig','Traf3_KO.Sig','IFN.Y.28.Sig',
              'CTLA4.Sig','Atezolizumab.Sig','PDL1.Sig','CYT.Sig','CD8A.Sig','Hwang.2020.Sig',
              'IFN.Y.10.Sig','TSE.Sig','Inflamed20genes.Sig','NI.score.Sig','HOT.score.Sig',
              'SELECT1.Sig','CD96.Sig','SIA.Sig','IMPRES.Sig')
  Sig_table = data.frame(Signature = c(Tax_sig,ICI_sig),
                         Group = c(rep("Tax",length(Tax_sig)),
                                   rep("ICI",length(ICI_sig))))
  if ("all" %in% Sig_list) {
    Sig_temp = Sig_table
  }else if("Tax" %in% Sig_list){
    Sig_temp = Sig_table[Sig_table$Group=="Tax",]
  }else if("ICI" %in% Sig_list){
    Sig_temp = Sig_table[Sig_table$Group=="ICI",]
  }else{
    Sig_temp = Sig_table[Sig_table$Signature%in%Sig_list,]
  }
  Sig_temp = unique(Sig_temp$Signature)
  if (length(Sig_temp)<1) {
    print(Sig_table)
    stop("Please select the correct signature.")
  }

  predict_result = lapply(Sig_temp, function(x){
    if (verbose) {
      message(paste0("\n",x," is running\n"))
    }
    temp_fun = get(x)
    y = try(temp_fun(expr = expr,Response = Response,verbose = verbose),silent = !verbose)
    if ("try-error" %in% class(y)) {
      y=rep(NA,ncol(expr))
      names(y) = colnames(expr)
      y = list(prediction=y,auc=NA)
    }else{
      if (x=='PredictTax'|x=='PredictICI') {
        prediction = y[["prediction"]]
        prediction = prediction[,3]
        names(prediction) = rownames(y[["prediction"]])
        y[["prediction"]] = prediction
      }
    }
    return(y)
  })

  predict_table = as.data.frame(lapply(predict_result, function(x){x[1]}))
  colnames(predict_table) = Sig_temp

  auc_value = unlist(lapply(predict_result, function(x){x[2]}))
  names(auc_value) = Sig_temp

  out_put = list(prediction = predict_table,auc = auc_value)
  return(out_put)
}

#'PredictICI function for single cell and spatial transcriptome data
#'
#'Similar to the PredictICI function, but the input data is in Seurat format.
#' @import GSVA
#' @importFrom utils packageVersion
#' @import Seurat
#' @import SeuratObject
#' @param object An Seurat object.
#' @param method Scoring algorithm, including `GSVA`, `ssGSEA` and `AUCell`.
#' @param slot Specific assay data to get.
#' @param asssy Specific assay to get data from. Spatial transcriptome data may need to set to `Spatial`.
#' @param verbose Gives information about each calculation step.
#' @param nCores When the method is selected as' aucell ', this is one of the parameters of the AUCell_calcAUC function.
#' @param normAUC When the method is selected as' aucell ', this is one of the parameters of the AUCell_calcAUC function
#' @param ... Other arguments passed on to gsva()'s or AUCell_buildRankings()'s params argument.
#' @export
#'
#' @examples
#' # res = scPredictICI(object = Seurat_obj, method = "AUCell", slot="data", asssy="RNA")
#'
scPredictICI <- function(object,method="AUCell",slot="data",asssy=NULL,verbose=FALSE,
                           nCores = 1,normAUC = TRUE,...) {
  gene_up=c('IFNG','LAG3','GBP1','CXCL9','TBX21','IGHV3-48','TRGC1','FASLG','APOL6',
            'PLA2G2D','IGHV3-11','JAKMIP1','IGLV3-21','CXCL11','PSME1','GBP4','CTLA4',
            'CD274','FBXO6','PSMB10','IGHV4-59','IGLV3-1','IL4I1','TAP1','IGKV3-11',
            'KLRC2','KLRK1','IGKV1-12','IL21','UBE2L6','IGHV4-34','H2AC21','CAPG',
            'P2RX5','PSMB9','TNIP3','CEP128','MYO7A','IL32','OASL','CLEC6A','KLRC4-KLRK1',
            'GBP2','AIM2','OAS3','FCRL6','HAPLN3','GPR25','H2BC10','FCRL3','GFI1','H1-3')
  gene_down=c('HSPG2','OGN','NINL','EFNA1','SMO','BCAM','EIF3E','ITGB8','PRKG1',
              'NFE2L2','NDST1','RPL23','TFF3','PARVA','EXTL2','CXADR','HSD11B2',
              'RPS13','SEMA6A','SMAD5','ID4','CD24','PBX1','COL4A5','PLXNB1','IGF1R',
              'PRKAB1','HEYL','JUP','CSRP1','SORBS1','TPT1','GNG12','WLS','NECTIN1',
              'MAP2','SOX6','NID1','EXT1','SGCD','SDC2','PLA2G12A','TSC22D1','PHLDB1',
              'DYNC2H1','NFIA','PTPRZ1','SORT1')

  seurat_version <- utils::packageVersion("Seurat")
  if (seurat_version >= "5.0") {
    expr = Seurat::GetAssayData(object = object,assay = asssy,layer = slot)
  }else{
    expr = Seurat::GetAssayData(object = object,assay = asssy,slot = slot)
  }
  if (method == "GSVA") {
    message("Using the gsva algorithm for single-cell data may run for a long time")
  }
  predict_result = .predict_sig(expr = expr,gene_up = gene_up,gene_down = gene_down,
                                verbose = verbose,method = method,nCores = nCores,
                                normAUC = normAUC,predict_name="ICI",...)
  object <- Seurat::AddMetaData(object, metadata = predict_result)
  return(object)
}

#'PredictTax function for single cell and spatial transcriptome data
#'
#'Similar to the PredictTax function, but the input data is in Seurat format.
#' @import GSVA
#' @importFrom utils packageVersion
#' @import Seurat
#' @import SeuratObject
#' @param object An Seurat object.
#' @param method Scoring algorithm, including `GSVA`, `ssGSEA` and `AUCell`.
#' @param slot Specific assay data to get.
#' @param asssy Specific assay to get data from. Spatial transcriptome data may need to set to `Spatial`.
#' @param verbose Gives information about each calculation step.
#' @param nCores When the method is selected as' aucell ', this is one of the parameters of the AUCell_calcAUC function.
#' @param normAUC When the method is selected as' aucell ', this is one of the parameters of the AUCell_calcAUC function
#' @param ... Other arguments passed on to gsva()'s or AUCell_buildRankings()'s params argument.
#' @export
#'
#' @examples
#' # res = scPredictTax(object = Seurat_obj, method = "AUCell", slot="data", asssy="RNA")
#'
scPredictTax <- function(object,method="AUCell",slot="data",asssy=NULL,verbose=FALSE,
                           nCores = 1,normAUC = TRUE,...) {
  gene_up=c('GBP1','GTSE1','CTSC','FBXO5','GMNN','GBP5','NUP153','CXCL13','PDIA6',
            'LYZ','SNRPA','PLA2G7','CASP8AP2','TCF3','ZBED2','BATF2','TFDP1','IGKV4-1',
            'ELF5','BTN2A1','TFDP2','H3C3','GTF2E2','SMARCD2','SLC35B1','MGAT1','MED1',
            'CRLF1','MAP3K3')
  gene_down=c('SUSD4','S100G','UNC5B','SDC4','PPP2R2C','EIF4G3','TACC2','RPS4Y1',
              'LAMP5','CTTN','FAN1','HGSNAT','KIF3A','PLA2G12A','DOCK1','ULK2','ARHGEF12',
              'FLNB','CGA','ERLIN2','CELSR2','BAIAP3','HMGCL','PURA','FBXL16','SPATA18',
              'BHLHE40','GPRC5A','SLC7A2','MATN3')

  seurat_version <- utils::packageVersion("Seurat")
  if (seurat_version >= "5.0") {
    expr = Seurat::GetAssayData(object = object,assay = asssy,layer = slot)
  }else{
    expr = Seurat::GetAssayData(object = object,assay = asssy,slot = slot)
  }
  if (method == "GSVA") {
    message("Using the gsva algorithm for single-cell data may run for a long time")
  }
  predict_result = .predict_sig(expr = expr,gene_up = gene_up,gene_down = gene_down,
                                verbose = verbose,method = method,nCores = nCores,
                                normAUC = normAUC,predict_name="Tax",...)
  object <- Seurat::AddMetaData(object, metadata = predict_result)
  return(object)
}

#'Predicting chemoimmunotherapy response subtypes
#'
#'This function classifies samples into subtypes responding to chemoimmunotherapy based on the gene expression matrix.
#' @importFrom CMScaller ematAdjust ntp
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param threshold FDR thresholds used to determine CI-M subtypes.
#' @param nPerm an integer, number of permutations for pvalue estimation.
#' @param verbose Gives information about each calculation step.
#' @param doPlot logical, whether to produce prediction subHeatmap
#' @param ... Other arguments passed on to ematAdjust()'s or ntp()'s params argument.
#' @export
#'
#' @examples
#' res = CIsubtype(expr = gene_table,threshold = 0.2)
#'
CIsubtype <- function(expr,threshold=0.2,nPerm=1000,verbose=FALSE,doPlot=TRUE,...) {
  file_path = system.file("extdata", "CItype.csv", package = "PRTI")
  CItype = read.csv(file_path,header = T)
  emat = ematAdjust(expr, normMethod = "quantile",verbose = verbose,...)
  predictiCI = ntp(emat, CItype, doPlot=doPlot, nPerm=nPerm,verbose = verbose,...)
  predictiCI$prediction = as.character(predictiCI$prediction)
  predictiCI$prediction[predictiCI$FDR>=threshold] = "CI-M"
  return(predictiCI)
}
