#' Gene Expression Matrix
#'
#' A data.frame of expression values where rows correspond to genes and columns correspond to samples.
#'
#' @format A data frame with rows as genes and columns as samples.
#' @docType data
"gene_table"

#' Response Information
#'
#' A data frame containing the Response information for the samples.
#'
#' @format A data frame with the sample names and their response data.
#' @docType data
"pdata_table"


#' Internal function to calculate signature
#' @noRd
.checkGene <- function(expr,genes,verbose=FALSE) {
  gene_num = sum(genes%in%rownames(expr))
  gene_percent = round((gene_num/length(genes))*100,1)
  if (verbose) {
    cat(gene_num," genes are found ","(",gene_percent,"%).\n",sep ="")
  }
  if (gene_num<=1) {
    stop("Too few genes matched!")
  }
  if (gene_percent!=100) {
    warning("The prediction results may be inaccurate due to the lack of genes!")
  }
}

#' Internal function to calculate signature
#' @noRd
.calculateAUC <- function(y,Response,verbose=FALSE) {
  auc_value = NA
  if (!is.null(Response)) {
    auc_result = pROC::roc(Response,y,levels=c("NR","R"),direction="<")
    auc_value = as.numeric(auc_result$auc)
    if (verbose) {
      print(auc_result$auc)
    }
  }
  out_put = list(prediction = y,auc = auc_value)
  return(out_put)
}

#' Internal function to calculate signature
#' @noRd
.calculateSingle <- function(expr,gene,verbose=FALSE) {
  if (!gene%in%rownames(expr)) {
    stop("Too few genes matched!")
  }
  x = expr[rownames(expr)%in%gene,]
  y = as.numeric(x[1, ])
  names(y) = colnames(x)
  return(y)
}

#' Internal function to calculate signature
#' @noRd
.calculateMean <- function(expr,genes,verbose=FALSE) {
  .checkGene(expr = expr,genes = genes,verbose = verbose)
  x = expr[rownames(expr)%in%genes,]
  y=colMeans(x)
  return(y)
}

#' Internal function to calculate signature
#' @noRd
.calculateSum <- function(expr,genes,verbose=FALSE) {
  .checkGene(expr = expr,genes = genes,verbose = verbose)
  x = expr[rownames(expr)%in%genes,]
  y=colSums(x)
  return(y)
}

#' Internal function to calculate signature
#' @noRd
.calculateWeightedMean  <- function(expr,genes,coefs,verbose=FALSE,na.rm = FALSE) {
  .checkGene(expr = expr,genes = genes,verbose = verbose)

  index = which(genes %in% row.names(expr))
  genes = genes[index]
  coefs = coefs[index]
  x = expr[genes,]
  y = apply(x, 2, function(a){stats::weighted.mean(a,coefs,na.rm = na.rm)})
  return(y)
}

#' Internal function to calculate signature
#' @import GSVA
#' @importFrom utils packageVersion
#' @noRd
.calculateGSVA <- function(expr,genes,verbose=FALSE,...) {
  .checkGene(expr = expr,genes = genes,verbose = verbose)
  gsva_version <- utils::packageVersion("GSVA")
  if (gsva_version >= "1.50.0") {
    gsvapar = GSVA::gsvaParam(exprData = as.matrix(expr),geneSets = list(genes),kcdf='Gaussian',...)
    predict_result = GSVA::gsva(gsvapar,verbose=verbose)
  } else {
    predict_result = GSVA::gsva(expr = as.matrix(expr),gset.idx.list = list(genes),method='gsva',kcdf='Gaussian',verbose=verbose,...)
  }
  y <- as.numeric(predict_result)
  names(y)=colnames(predict_result)
  return(y)
}

#' Internal function to calculate signature
#' @import GSVA
#' @importFrom utils capture.output packageVersion
#' @noRd
.calculatessGSEA <- function(expr,genes,verbose=FALSE,...) {
  .checkGene(expr = expr,genes = genes,verbose = verbose)
  gsva_version <- utils::packageVersion("GSVA")
  if (gsva_version >= "1.50.0") {
    gsvapar = GSVA::ssgseaParam(exprData = as.matrix(expr),geneSets = list(genes),...)
    if (verbose) {
      predict_result <- GSVA::gsva(gsvapar,verbose=verbose)
    }else{
      capt = utils::capture.output(predict_result <- GSVA::gsva(gsvapar,verbose=verbose), file = NULL)
    }
  } else {
    if (verbose) {
      predict_result <- GSVA::gsva(expr = as.matrix(expr),gset.idx.list = list(genes),method="ssgsea",verbose=verbose,...)
    }else{
      capt = utils::capture.output(predict_result <- GSVA::gsva(expr = as.matrix(expr),gset.idx.list = list(genes),method="ssgsea",verbose=verbose,...), file = NULL)
    }
  }

  y <- as.numeric(predict_result)
  names(y)=colnames(predict_result)
  return(y)
}

#' Internal function to calculate signature
#' @noRd
.calculateREO <- function(expr,genes_a,genes_b,verbose=FALSE) {
  index = which(genes_a%in%rownames(expr)&genes_a%in%rownames(expr))

  gene_num = length(index)
  gene_percent = round((gene_num/length(genes_a))*100,1)
  if (verbose) {
    cat(gene_num," pairs of genes are found ","(",gene_percent,"%).\n",sep ="")
  }
  if (gene_num<=1) {
    stop("Too few genes matched!")
  }
  if (gene_percent!=100) {
    warning("The prediction results may be inaccurate due to the lack of genes!")
  }
  genes_a = genes_a[index]
  genes_b = genes_b[index]

  x_a = expr[rownames(expr)%in%genes_a,]
  x_b = expr[rownames(expr)%in%genes_a,]
  y = ifelse(x_a>x_b,1,0)
  y = ifelse(colSums(y)>=nrow(y)*0.5,1,0)
  return(y)
}

#' Geometric mean
#' @noRd
.gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

#' Predictive function
#' @import GSVA
#' @importFrom utils packageVersion capture.output
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC getAUC
#' @noRd
.predict_sig <- function(expr,gene_up,gene_down,predict_name,method=c('GSVA','ssGSEA','AUCell'),
                         verbose=FALSE,nCores = 1,normAUC = TRUE,...) {
  method <- match.arg(method)
  genes_list = list(gene_up=gene_up,gene_down=gene_down)

  gene_up_num = sum(genes_list$gene_up%in%rownames(expr))
  gene_down_num = sum(genes_list$gene_down%in%rownames(expr))
  gene_up_percent = round((gene_up_num/length(genes_list$gene_up))*100,1)
  gene_down_percent = round((gene_down_num/length(genes_list$gene_down))*100,1)

  if (verbose) {
    cat(gene_up_num," up-regulated genes in the ",predict_name,"_gene_list are found ","(",gene_up_percent,"%).\n",sep ="")
    cat(gene_down_num," down-regulated genes in the ",predict_name,"_gene_list are found ","(",gene_down_percent,"%).\n\n",sep ="")
  }
  if (gene_up_num<=1 | gene_down_num<=1) {
    stop("Too few genes matched!")
  }
  if (gene_up_percent!=100 | gene_down_percent!=100) {
    warning("The prediction results may be inaccurate due to the lack of genes!")
  }

  if (method=="AUCell") {
    cells_rankings <- AUCell::AUCell_buildRankings(exprMat = as.matrix(expr),splitByBlocks=TRUE,
                                                   verbose = verbose,...)
    cells_AUC <- AUCell::AUCell_calcAUC(geneSets = genes_list, rankings = cells_rankings,
                                        verbose = verbose,nCores = nCores,normAUC = normAUC)
    predict_result = data.frame(gene_up = as.numeric(AUCell::getAUC(cells_AUC)["gene_up", ]),
                                gene_down = as.numeric(AUCell::getAUC(cells_AUC)["gene_down", ]))
    predict_result[,paste0("Predict",predict_name)] = (predict_result$gene_up)-(predict_result$gene_down)
    rownames(predict_result) = cells_AUC@colData@rownames
  }else if(method=="GSVA"){
    gsva_version <- utils::packageVersion("GSVA")
    if (gsva_version >= "1.50.0") {
      gsvapar = GSVA::gsvaParam(exprData = as.matrix(expr),geneSets = genes_list,kcdf='Gaussian',...)
      predict_result = GSVA::gsva(gsvapar,verbose=verbose)
    }else{
      predict_result = GSVA::gsva(expr = as.matrix(expr),gset.idx.list = genes_list,method='gsva',kcdf='Gaussian',verbose=verbose,...)
    }
    predict_result = as.data.frame(t(predict_result))
    predict_result[,paste0("Predict",predict_name)] = as.numeric(predict_result$gene_up) - as.numeric(predict_result$gene_down)
  }else if(method=="ssGSEA"){

    gsva_version <- utils::packageVersion("GSVA")
    if (gsva_version >= "1.50.0") {
      gsvapar = GSVA::ssgseaParam(exprData = as.matrix(expr),geneSets = genes_list,...)
      if (verbose) {
        predict_result <- GSVA::gsva(gsvapar,verbose=verbose)
      }else{
        capt = utils::capture.output(predict_result <- GSVA::gsva(gsvapar,verbose=verbose), file = NULL)
      }
    }else{
      if (verbose) {
        predict_result <- GSVA::gsva(expr = as.matrix(expr),gset.idx.list = genes_list,method="ssgsea",verbose=verbose,...)
      }else{
        capt = utils::capture.output(predict_result <- GSVA::gsva(expr = as.matrix(expr),gset.idx.list = genes_list,method="ssgsea",verbose=verbose,...), file = NULL)
      }
    }

    predict_result = as.data.frame(t(predict_result))
    predict_result[,paste0("Predict",predict_name)] = as.numeric(predict_result$gene_up) - as.numeric(predict_result$gene_down)
  }
  return(predict_result)
}
