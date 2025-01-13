#' Function to compute the signature as published by Ahmed et al. in 2007.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Ahmed et al. in 2007.
#' @references Ahmed, Ahmed Ashour et al. “The extracellular matrix protein TGFBI induces microtubule stabilization and sensitizes ovarian cancers to paclitaxel.” Cancer cell vol. 12,6 (2007): 514-27. doi:10.1016/j.ccr.2007.11.014
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = TGFBI.Sig(expr = gene_table,Response = pdata_table$Response)
#'
TGFBI.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('TGFBI')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' A function for calculating the average value of target genes of docetaxel in DrugBank.
#'
#' A list of target genes with direct pharmacological action with docetaxel is found from DrugBank, and the average value of their expression is calculated.
#' @references Knox, Craig et al. “DrugBank 6.0: the DrugBank Knowledgebase for 2024.” Nucleic acids research vol. 52,D1 (2024): D1265-D1275. doi:10.1093/nar/gkad976
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = Docetaxel.Sig(expr = gene_table,Response = pdata_table$Response)
#'
Docetaxel.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('TUBB1','MAP2','MAP4','MAPT')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Malla et al. in 2020.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Malla et al. in 2020.
#' @references Malla, Sudhir B et al. “In-depth Clinical and Biological Exploration of DNA Damage Immune Response as a Biomarker for Oxaliplatin Use in Colorectal Cancer.” Clinical cancer research : an official journal of the American Association for Cancer Research vol. 27,1 (2021): 288-300. doi:10.1158/1078-0432.CCR-20-3237
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = DDIR.Sig(expr = gene_table,Response = pdata_table$Response)
#'
DDIR.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CXCL10','MX1','IDO1','IFI44L','CD2','GBP5','PRAME','ITGAL','LRP4',
             'APOL3','CDR1','FYB1','TSPAN7','RAC2','KLHDC7B','GRB14','MT-RNR1',
             'KIF26A','CD274','CD109','ETV7','MFAP5','OLFM4','PI15','FOSB','TAFA5',
             'NLRC5','PRICKLE1','EGR1','CLDN10','ADAMTS4','SP140L','ANXA1','RSAD2',
             'ESR1','IKZF3','OR2I1P','EGFR','NAT1','LATS2','CYP2B6','PTPRC','PPP1R1A',
             'AL137218.1')
  coefs <- c(0.023,0.0226,0.0221,0.0191,0.019,0.0181,0.0177,0.0176,-0.0159,0.0151,
             -0.0149,0.0149,-0.0148,0.0148,0.014,0.0137,-0.0136,-0.0136,0.0133,-0.0129,
             0.0124,-0.0121,-0.0117,-0.0115,-0.0111,-0.011,0.0101,-0.0089,-0.0086,-0.0086,
             -0.0085,0.0084,-0.0082,0.0081,0.0079,0.0073,0.007,-0.0066,0.0065,-0.0063,
             0.0061,0.0051,-0.0041,-0.0017)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}


#' Function to compute the signature as published by Foy et al. in 2022.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Foy et al. in 2022.
#' @references Foy, Jean-Philippe et al. “Immunologically active phenotype by gene expression profiling is associated with clinical benefit from PD-1/PD-L1 inhibitors in real-world head and neck and lung cancer patients.” European journal of cancer (Oxford, England : 1990) vol. 174 (2022): 287-298. doi:10.1016/j.ejca.2022.06.034
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @param ... Other arguments passed on to gsva()'s params argument.
#' @export
#' @examples
#' a = HOT.score.Sig(expr = gene_table,Response = pdata_table$Response)
#'
HOT.score.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('CCL19','CCR2','CCR4','CCR5','CD27','CD40LG','CD8A','CXCL10','CXCL11',
             'CXCL13','CXCL9','CXCR3','CXCR6','FASLG','FGL2','GZMA','GZMH','IDO1',
             'IFNG','IRF8','LAG3','LYZ','MS4A1','PDCD1','TBX21','TLR7','TLR8')

  y = .calculateGSVA(expr = expr,genes = genes,verbose = verbose,...)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}


#' Function to compute the signature as published by Şenbabaoğlu et al. in 2016.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Şenbabaoğlu et al. in 2016.
#' @references Şenbabaoğlu, Yasin et al. “Tumor immune microenvironment characterization in clear cell renal cell carcinoma identifies prognostic and immunotherapeutically relevant messenger RNA signatures.” Genome biology vol. 17,1 231. 17 Nov. 2016, doi:10.1186/s13059-016-1092-z
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @param ... Other arguments passed on to gsva()'s params argument.
#' @export
#' @examples
#' a = APM.Sig(expr = gene_table,Response = pdata_table$Response)
#'
APM.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2','TAPBP')

  y = .calculatessGSEA(expr = expr,genes = genes,verbose = verbose,...)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}


#' Function to compute the signature as published by Qi et al. in 2016.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Qi et al. in 2016.
#' @references Qi, Lishuang et al. “An individualised signature for predicting response with concordant survival benefit for lung adenocarcinoma patients receiving platinum-based chemotherapy.” British journal of cancer vol. 115,12 (2016): 1513-1519. doi:10.1038/bjc.2016.370
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = GPS_3.Sig(expr = gene_table,Response = pdata_table$Response)
#'
GPS_3.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes_a <- c('REXO5','TREM1','ERN1')
  genes_b <- c('ZMYND10','ANPEP','FA2H')

  y = .calculateREO(expr = expr,genes_a = genes_a,genes_b = genes_b,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by van't Veer et al. in 2002.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by van't Veer et al. in 2002.
#' @references van 't Veer, Laura J et al. “Gene expression profiling predicts clinical outcome of breast cancer.” Nature vol. 415,6871 (2002): 530-6. doi:10.1038/415530a
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @param ... Other arguments passed on to gene70()'s params argument.
#' @importFrom pROC roc
#' @importFrom genefu gene70
#' @importFrom utils data
#' @export
#' @examples
#' a = MammaPrint.Sig(expr = gene_table,Response = pdata_table$Response)
#'
MammaPrint.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  if (!exists("sig.gene70"))
    utils::data("sig.gene70",package ="genefu")
  genes = c('ALDH4A1','FGF18','BBC3','SCUBE2','RUNDC1','CCN4','GSTM3','ZNF385B',
            'RTN4RL1','ECI2','TGFB3','STK32B','MS4A7','AP2B1','DHX58','TMEM74B',
            'CCNE2','CENPA','PRC1','NMU','IGFBP5','PITRM1','PLAAT1','MSANTD3',
            'MCM6','CDCA7','RFC4','ORC6','SLC2A3','ADGRG6','DCK','DTL','COL4A2',
            'MELK','MTDH','UCHL5','RAB6B','GPR180','LPCAT1','CDC42BPA','NDC80',
            'GMPS','ECT2','MMP9','OXCT1','GNAZ','EXT1','CMC2','DIAPH3','QSOX2',
            'NUSAP1','TSPYL5')

  Entrez_ID = c(8659,8817,27113,57758,146923,8840,2947,151126,146760,10455,7043,
                55351,58475,163,79132,55321,9134,1058,9055,10874,3488,10531,57110,
                91283,4175,83879,5984,23594,6515,57211,1633,51514,1284,9833,92140,
                51377,51560,160897,79888,8476,10403,8833,1894,4318,5019,2781,2131,
                56942,81624,169714,51203,85453)

  .checkGene(expr = expr,genes = genes,verbose = verbose)
  x = t(expr[rownames(expr)%in%genes,])

  annot = data.frame(row.names = genes,probe = genes,EntrezGene.ID = Entrez_ID)
  annot = annot[annot$probe %in% colnames(x),]

  y <- genefu::gene70(data=x, annot=annot, do.mapping=TRUE,verbose = verbose,...)
  y = y$score
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Sotiriou et al. in 2006.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Sotiriou et al. in 2006.
#' @references Sotiriou, Christos et al. “Gene expression profiling in breast cancer: understanding the molecular basis of histologic grade to improve prognosis.” Journal of the National Cancer Institute vol. 98,4 (2006): 262-72. doi:10.1093/jnci/djj052
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @param ... Other arguments passed on to ggi()'s params argument.
#' @importFrom pROC roc
#' @importFrom genefu ggi
#' @importFrom utils data
#' @export
#' @examples
#' a = GGI.Sig(expr = gene_table,Response = pdata_table$Response)
#'
GGI.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  if (!exists("sig.ggi"))
    utils::data("sig.ggi",package ="genefu")
  genes = c('KPNA2','SLC7A5','MARS1','UBE2N','DDX39A','MYBL2','BIRC5','MCM2','PLK1',
            'FOXM1','CCNB2','UBE2S','CDC20','UBE2C','RNASEH2A','TIMELESS','SPAG5',
            'CDK1','LMNB1','MAD2L1','CCNA2','PTTG1','HMGB3','BUB1B','DLGAP5','ZWINT',
            'TRIP13','AURKA','NDC80','CDK2','GTSE1','KIF11','EXO1','TROAP','CDC25A',
            'NUDT1','FEN1','ESPL1','TTK','MELK','CENPA','CCNE2','CENPE','BLM','KIF14',
            'HMMR','CENPF','CCT5','TUBA1C','KIF2C','AURKB','BUB1','KIFC1','CDKN3',
            'RRM2','TPX2','TUBA1B','MKI67','MCM4','NCAPH','OIP5','H2AZ1','SHMT2',
            'GMPS','CCNB1','RAB5IF','PRC1','NUSAP1','KIF4A','CMC2','CEP55','ORMDL2',
            'NCAPG','HJURP','KIF20A','CENPU','DSCC1','KIF15','POLQ','CENPN','ASPM',
            'PARPBP','MCM10','CDCA3','CDCA8','PIMREG','DONSON','RACGAP1','FRY',
            'IFT88','CX3CR1','STARD13','LAMB2','TPT1','CYBRD1','SESN1','BBS1','IFT46',
            'PIGV','CFAP69','JHY','WDR19','SIRT3','MPHOSPH8')

  Entrez_ID = c(3838,8140,4141,7334,10212,4605,332,4171,5347,2305,9133,27338,991,
                11065,10535,8914,10615,983,4001,4085,890,9232,3149,701,9787,11130,
                9319,6790,10403,1017,51512,3832,9156,10024,993,4521,2237,9700,7272,
                9833,1058,9134,1062,641,9928,3161,1063,22948,84790,11004,9212,699,
                3833,1033,6241,22974,10376,4288,4173,23397,11339,3015,6472,8833,
                891,55969,9055,51203,24137,56942,55165,29095,64151,55355,10112,
                79682,79075,56992,10721,55839,259266,55010,55388,83461,55143,54478,
                29980,29127,10129,8100,1524,90627,3913,7178,79901,27244,582,56912,
                55650,79846,79864,57728,23410,54737)

  .checkGene(expr = expr,genes = genes,verbose = verbose)
  x = t(expr[rownames(expr)%in%genes,])

  annot = data.frame(row.names = genes,probe = genes,EntrezGene.ID = Entrez_ID)
  annot = annot[annot$probe %in% colnames(x),]

  y <- genefu::ggi(data=x, annot=annot, do.mapping=TRUE,verbose = verbose,...)
  y = y$score
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Paik et al. in 2004.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Paik et al. in 2004.
#' @references Paik, Soonmyung et al. “A multigene assay to predict recurrence of tamoxifen-treated, node-negative breast cancer.” The New England journal of medicine vol. 351,27 (2004): 2817-26. doi:10.1056/NEJMoa041588
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @param ... Other arguments passed on to oncotypedx()'s params argument.
#' @importFrom pROC roc
#' @importFrom genefu oncotypedx
#' @importFrom utils data
#' @export
#' @examples
#' a = Oncotype.DX.Sig(expr = gene_table,Response = pdata_table$Response)
#'
Oncotype.DX.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  if (!exists("sig.oncotypedx"))
    utils::data("sig.oncotypedx",package ="genefu")
  genes = c('MKI67','AURKA','BIRC5','CCNB1','MYBL2','MMP11','CTSV','GRB7','ERBB2',
            'ESR1','PGR','BCL2','SCUBE2','GSTM1','CD68','BAG1','ACTB','GAPDH','RPLP0',
            'GUSB','TFRC')

  Entrez_ID = c(4288,6790,332,891,4605,4320,1515,2886,2064,2099,5241,596,57758,
                2944,968,573,60,2597,6175,2990,7037)

  .checkGene(expr = expr,genes = genes,verbose = verbose)
  x = t(expr[rownames(expr)%in%genes,])

  annot = data.frame(row.names = genes,probe = genes,EntrezGene.ID = Entrez_ID)
  annot = annot[annot$probe %in% colnames(x),]

  y <- genefu::oncotypedx(data=x, annot=annot, do.mapping=TRUE,verbose = verbose,...)
  y = y$score
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Birkbak et al. in 2018.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Birkbak et al. in 2018.
#' @references Birkbak, N J et al. “Overexpression of BLM promotes DNA damage and increased sensitivity to platinum salts in triple-negative breast and serous ovarian cancers.” Annals of oncology : official journal of the European Society for Medical Oncology vol. 29,4 (2018): 903-909. doi:10.1093/annonc/mdy049
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = Birkbak.2018.Sig(expr = gene_table,Response = pdata_table$Response)
#'
Birkbak.2018.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('BLM','FANCI','BRCA1')
  if (!all(genes%in%rownames(expr))) {
    stop("Too few genes matched!")
  }
  x <- t(expr[genes,])
  y <- (x[,1]+x[,2])/x[,3]

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Zhang et al. in 2016.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Zhang et al. in 2016.
#' @references Zhang, Weiguo et al. “Centromere and kinetochore gene misexpression predicts cancer patient survival and response to radiotherapy and chemotherapy.” Nature communications vol. 7 12619. 31 Aug. 2016, doi:10.1038/ncomms12619
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = CES.Sig(expr = gene_table,Response = pdata_table$Response)
#'
CES.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CENPA','CENPN','CENPM','CENPK','CENPL','CENPW','CENPU','HJURP',
             'OIP5','ZWINT','NDC80','SPC24','SPC25','NUF2')
  y = .calculateSum(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Rooney et al. in 2015.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Rooney et al. in 2015.
#' @references Rooney, Michael S et al. “Molecular and genetic properties of tumors associated with local immune cytolytic activity.” Cell vol. 160,1-2 (2015): 48-61. doi:10.1016/j.cell.2014.12.033
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = CYT.Sig (expr = gene_table,Response = pdata_table$Response)
#'
CYT.Sig  <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('GZMA','PRF1')
  if (!all(genes%in%rownames(expr))) {
    stop("Too few genes matched!")
  }
  x <- expr[genes,]
  y <- apply(x,2,.gm_mean)

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' A function for calculating the average value of target genes of PD-L1 inhibitor in DrugBank.
#'
#' A list of target genes with direct pharmacological action with PD-L1 inhibitor is found from DrugBank, and the average value of their expression is calculated.
#' @references Knox, Craig et al. “DrugBank 6.0: the DrugBank Knowledgebase for 2024.” Nucleic acids research vol. 52,D1 (2024): D1265-D1275. doi:10.1093/nar/gkad976
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = PDL1.Sig(expr = gene_table,Response = pdata_table$Response)
#'
PDL1.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('CD274')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Liu et al. in 2012.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Liu et al. in 2012.
#' @references Liu, Jeff C et al. “Seventeen-gene signature from enriched Her2/Neu mammary tumor-initiating cells predicts clinical outcome for human HER2+:ERα- breast cancer.” Proceedings of the National Academy of Sciences of the United States of America vol. 109,15 (2012): 5832-7. doi:10.1073/pnas.1201105109
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = HTICS.Sig(expr = gene_table,Response = pdata_table$Response)
#'
HTICS.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('AURKB','CLDN8','NPY','ATP7B','CHAF1B','SCRN1','CCNA2','CCNB1','NRP1',
             'CD74','C1QB','CD72','VCAM1','ITGB2','CD180','CCR2','ST8SIA4')
  coefs <- c(1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Ayers et al. in 2017.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Ayers et al. in 2017.
#' @references Ayers, Mark et al. “IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade.” The Journal of clinical investigation vol. 127,8 (2017): 2930-2940. doi:10.1172/JCI91190
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = IFN.Y.10.Sig(expr = gene_table,Response = pdata_table$Response)
#'
IFN.Y.10.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('IFNG','STAT1','CCR5','CXCL9','CXCL10','CXCL11','IDO1','PRF1','GZMA','HLA-DRA')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Higgs et al. in 2018.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Higgs et al. in 2018.
#' @references Higgs, Brandon W et al. “Interferon Gamma Messenger RNA Signature in Tumor Biopsies Predicts Outcomes in Patients with Non-Small Cell Lung Carcinoma or Urothelial Cancer Treated with Durvalumab.” Clinical cancer research : an official journal of the American Association for Cancer Research vol. 24,16 (2018): 3857-3866. doi:10.1158/1078-0432.CCR-17-3451
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = IFNy4gene.Sig(expr = gene_table,Response = pdata_table$Response)
#'
IFNy4gene.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('IFNG','CD274','LAG3','CXCL9')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Ayers et al. in 2017.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Ayers et al. in 2017.
#' @references Ayers, Mark et al. “IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade.” The Journal of clinical investigation vol. 127,8 (2017): 2930-2940. doi:10.1172/JCI91190
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = GEP.Sig(expr = gene_table,Response = pdata_table$Response)
#'
GEP.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CCL5','CD27','CD274','CD276','CD8A','CMKLR1','CXCL9','CXCR6','HLA-DQA1',
             'HLA-DRB1','HLA-E','IDO1','LAG3','NKG7','PDCD1LG2','PSMB10','STAT1','TIGIT')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Ayers et al. in 2017.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Ayers et al. in 2017.
#' @references Ayers, Mark et al. “IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade.” The Journal of clinical investigation vol. 127,8 (2017): 2930-2940. doi:10.1172/JCI91190
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = IFN.Y.28.Sig(expr = gene_table,Response = pdata_table$Response)
#'
IFN.Y.28.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('IL2RG','CXCR6','CD3D','CD2','ITGAL','TAGAP','CIITA','HLA-DRA','PTPRC',
             'CXCL9','CCL5','NKG7','GZMA','PRF1','CCR5','CD3E','GZMK','IFNG','HLA-E',
             'GZMB','PDCD1','SLAMF6','CXCL13','CXCL10','IDO1','LAG3','STAT1','CXCL11')

  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Şenbabaoğlu et al. in 2016.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Foy et al. in 2016.
#' @references Şenbabaoğlu, Yasin et al. “Tumor immune microenvironment characterization in clear cell renal cell carcinoma identifies prognostic and immunotherapeutically relevant messenger RNA signatures.” Genome biology vol. 17,1 231. 17 Nov. 2016, doi:10.1186/s13059-016-1092-z
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @param ... Other arguments passed on to gsva()'s params argument.
#' @export
#' @examples
#' a = IIS.Sig(expr = gene_table,Response = pdata_table$Response)
#'
IIS.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes_list = list(
    CD8_T_cells_genes=c('ABT1','TLE5','APBA2','ARHGAP8','MAPKAPK5-AS1','TMEM259','TRMO','CAMLG','CD8A','CD8B','CDKN2AIP','DNAJB1','FLT3LG','GADD45A','GZMM','KLF9','LEPROTL1','LIME1','MYST3','PF4','PPP1R2','PRF1','PRR5','RBM3','SF1','SFRS7','SLC16A7','TBCC','THUMPD1','TMC6','TSC22D3','VAMP2','ZEB1','ZFP36L2','ZNF22','ZNF609','ZNF91'),
    T_cells_genes=c('BCL11B','CD2','CD28','CD3D','CD3E','CD3G','CD6','CD96','GIMAP5','ITM2A','LCK','NCALD','PRKCQ','SH2D1A','SKAP1','TRA','TRAC','TRAT1','TRBC1'),
    T_helper_cells_genes=c('ANP32B','ASF1A','ATF2','BATF','BORA','CD28','DDX50','FAM111A','FRYL','GOLGA8A','ICOS','ITM2A','LRBA','NAP1L4','NUP107','PHF10','PPP2R5C','RPA1','SEC24C','SLC25A12','SRSF10','TRA','UBE2L3','YME1L1'),
    Tcm_cells_genes=c('AQP3','ATF7IP','ATM','CASP8','CDC14A','CEP68','CLUAP1','CREBZF','CYLD','DOCK9','FAM153B','FOXP1','FYB1','HNRNPH1','INPP4B','KLF12','LOC441155','MAP3K1','KMT2A','N4BP2L2-IT2','NEFL','NFATC3','PCM1','PCNX1','PDXDC2P-NPIPB14P','PHC3','POLR2J2','PSPC1','REPS1','RPP38','SLC7A6','SNRPN','ST3GAL1','STX16','TIMM8A','TRAF3IP3','TXK','TXLNGY','USP9Y'),
    Tem_cells_genes=c('AKT3','SND1-IT1','CCR2','DDX17','EWSR1','FLI1','GDPD5','LTK','MEFV','NFATC4','PRKY','TBC1D5','TBCD','TRA','EZR'),
    Th1_cells_genes=c('APBB2','APOD','ATP9A','BST2','BTG3','CCL4','CD38','CD70','CMAHP','CSF2','CTLA4','DGKI','DOK5','DPP4','DUSP5','EGFL6','GGT1','HBEGF','IFNG','IL12RB2','IL22','LRP8','LRRN3','LTA','SGCB','SYNGR3','ZBTB32'),
    Th17_cells_genes=c('IL17A','IL17RA','RORC'),
    Th2_cells_genes=c('ADCY1','AHI1','ANK1','BIRC5','CDC25C','CDC7','CENPF','CXCR6','DHFR','EVI5','GATA3','GSTA4','HELLS','IL26','LAIR2','LIMA1','MB','MICAL2','NEIL3','PHEX','PMCH','PTGIS','SLC39A14','SMAD2','SNRPD1','WDHD1'),
    Treg_cells_genes=c('FOXP3'),
    aDC_genes=c('CCL1','EBI3','IDO1','LAMP3','OAS3'),
    B_cells_genes=c('ABCB4','BACH2','BCL11A','BLK','BLNK','CCR9','CD19','CD72','COCH','CR2','DTNB','FCRL2','GLDC','GNG7','HLA-DOB','HLA-DQA1','IGHA1','IGHG1','IGHM','IGKC','IGL','FAM30A','MEF2C','MICAL3','MS4A1','OSBPL10','PNOC','QRSL1','SCN3A','SLC15A2','SPIB','TCL1A','TNFRSF17'),
    Cytotoxic_cells_genes=c('APBA2','APOL3','CTSW','DUSP2','GNLY','GZMA','GZMH','KLRB1','KLRD1','KLRF1','KLRK1','NKG7','RORA','RUNX3','SIGIRR','WHAMMP3','ZBTB16'),
    DC_genes=c('CCL13','CCL17','CCL22','CD209','HSD11B1','NPR1','PPFIBP2'),
    Eosinophils_genes=c('ABHD2','ACACB','C9orf156','CAT','CCR3','CLC','CYSLTR2','EMR1','EPN2','GALC','GPR44','HES1','HIST1H1C','HRH4','CD101','IL5RA','KBTBD11','KCNH2','LRP5L','MYO15B','RCOR3','RNASE2','RRP12','SIAH1','SMPD3','SYNJ1','TGIF1','THBS1','THBS4','TIPARP','TKTL1'),
    iDC_genes=c('ABCG2','BLVRB','CARD9','CD1A','CD1B','CD1C','CD1E','CH25H','CLEC10A','CSF1R','CTNS','F13A1','FABP4','FZD2','GSTT1','GUCA1A','HS3ST2','LMAN2L','MMP12','MS4A6A','NUDT9','PDXK','PPARG','PREP','RAP1GAP','SLC26A6','SLC7A8','SYT17','TACSTD2','DCSTAMP','VASH1'),
    Macrophages_genes=c('APOE','ATG7','BCAT1','CCL7','CD163','CD68','CD84','CHI3L1','CHIT1','CLEC5A','COL8A2','COLEC12','CTSK','CXCL5','CYBB','DNASE2B','EMP1','FDX1','FN1','GM2A','GPC4','ANOS1','MARCO','ME1','MS4A4A','MSR1','PCOLCE2','PTGDS','RAI14','SCARB2','SCG5','SGMS1','SULT1C2'),
    Mast_cells_genes=c('ABCC4','ADCYAP1','CALB2','CEACAM8','CMA1','CPA3','CTSG','ELANE','GATA2','HDC','HPGD','HPGDS','KIT','LINC01140','MAOB','MLPH','MPO','MS4A2','NR0B1','PPM1H','PRG2','PTGS1','SCG2','SIGLEC6','SLC18A2','SLC24A3','TAL1','TPSAB1','TPSB2','VWA5A'),
    Neutrophils_genes=c('ALPL','BST1','CD93','CEACAM3','CREB5','CRISPLD2','CSF3R','CYP4F3','DYSF','FCAR','FCGR3B','CPPED1','FPR1','FPR2','G0S2','H2BC4','HPSE','CXCR1','CXCR2','KCNJ15','LILRB2','MGAM','MME','PDE4B','S100A12','SIGLEC5','SLC22A4','SLC25A37','TECPR2','TNFRSF10C','VNN3'),
    NK_CD56bright_cells_genes=c('DUSP4','FOXJ1','LPCAT4','MADD','MARCHF6','MPPED1','MUC3B','TRAPPC9','PLA2G6','RRAD','XCL1'),
    NK_CD56dim_cells_genes=c('S1PR5','TTC38','GTF3C1','GZMB','IL21R','KIR2DL3','KIR2DS1','KIR2DS2','KIR2DS5','KIR3DL1','KIR3DL2','KIR3DL3','KIR3DS1','SPON2','PMEPA1'),
    NK_cells_genes=c('ADARB1','AF107846','ALDH1B1','APBB2','ATL2','BCL2','CDC5L','FGF18','FUT5','FZR1','GAGE2A','IGFBP5','KANK2','LDB3','MAPRE3','MCM3AP','MRC2','NCR1','PDLIM4','PRX','PSMD4','RP5-886K2.1','SGMS1','SLC30A5','PPP4R3A','SPN','TBXA2R','TCTN2','TINAGL1','TRPV6','XCL1','XCL2','ZNF205','ZNF528','ZNF747'),
    pDC_genes=c('IL3RA')
  )

  .checkGene(expr = expr,genes = unlist(genes_list),verbose = verbose)

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

  y = colSums(predict_result,na.rm = T)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Şenbabaoğlu et al. in 2016.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Foy et al. in 2016.
#' @references Şenbabaoğlu, Yasin et al. “Tumor immune microenvironment characterization in clear cell renal cell carcinoma identifies prognostic and immunotherapeutically relevant messenger RNA signatures.” Genome biology vol. 17,1 231. 17 Nov. 2016, doi:10.1186/s13059-016-1092-z
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @param ... Other arguments passed on to gsva()'s params argument.
#' @export
#' @examples
#' a = TIS.Sig(expr = gene_table,Response = pdata_table$Response)
#'
TIS.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes_list = list(
    CD8_T_cells_sig=c('ABT1','TLE5','APBA2','ARHGAP8','MAPKAPK5-AS1','TMEM259','TRMO','CAMLG','CD8A','CD8B','CDKN2AIP','DNAJB1','FLT3LG','GADD45A','GZMM','KLF9','LEPROTL1','LIME1','MYST3','PF4','PPP1R2','PRF1','PRR5','RBM3','SF1','SFRS7','SLC16A7','TBCC','THUMPD1','TMC6','TSC22D3','VAMP2','ZEB1','ZFP36L2','ZNF22','ZNF609','ZNF91'),
    T_cells_sig=c('BCL11B','CD2','CD28','CD3D','CD3E','CD3G','CD6','CD96','GIMAP5','ITM2A','LCK','NCALD','PRKCQ','SH2D1A','SKAP1','TRA','TRAC','TRAT1','TRBC1'),
    T_helper_cells_sig=c('ANP32B','ASF1A','ATF2','BATF','BORA','CD28','DDX50','FAM111A','FRYL','GOLGA8A','ICOS','ITM2A','LRBA','NAP1L4','NUP107','PHF10','PPP2R5C','RPA1','SEC24C','SLC25A12','SRSF10','TRA','UBE2L3','YME1L1'),
    Tcm_cells_sig=c('AQP3','ATF7IP','ATM','CASP8','CDC14A','CEP68','CLUAP1','CREBZF','CYLD','DOCK9','FAM153B','FOXP1','FYB1','HNRNPH1','INPP4B','KLF12','LOC441155','MAP3K1','KMT2A','N4BP2L2-IT2','NEFL','NFATC3','PCM1','PCNX1','PDXDC2P-NPIPB14P','PHC3','POLR2J2','PSPC1','REPS1','RPP38','SLC7A6','SNRPN','ST3GAL1','STX16','TIMM8A','TRAF3IP3','TXK','TXLNGY','USP9Y'),
    Tem_cells_sig=c('AKT3','SND1-IT1','CCR2','DDX17','EWSR1','FLI1','GDPD5','LTK','MEFV','NFATC4','PRKY','TBC1D5','TBCD','TRA','EZR'),
    Th1_cells_sig=c('APBB2','APOD','ATP9A','BST2','BTG3','CCL4','CD38','CD70','CMAHP','CSF2','CTLA4','DGKI','DOK5','DPP4','DUSP5','EGFL6','GGT1','HBEGF','IFNG','IL12RB2','IL22','LRP8','LRRN3','LTA','SGCB','SYNGR3','ZBTB32'),
    Th17_cells_sig=c('IL17A','IL17RA','RORC'),
    Th2_cells_sig=c('ADCY1','AHI1','ANK1','BIRC5','CDC25C','CDC7','CENPF','CXCR6','DHFR','EVI5','GATA3','GSTA4','HELLS','IL26','LAIR2','LIMA1','MB','MICAL2','NEIL3','PHEX','PMCH','PTGIS','SLC39A14','SMAD2','SNRPD1','WDHD1'),
    Treg_cells_sig=c('FOXP3')
  )

  .checkGene(expr = expr,genes = unlist(genes_list),verbose = verbose)

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

  y = colSums(predict_result,na.rm = T)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Lv et al. in 2023.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Lv et al. in 2023.
#' @references Lv, Jiawei et al. “The tumor immune microenvironment of nasopharyngeal carcinoma after gemcitabine plus cisplatin treatment.” Nature medicine vol. 29,6 (2023): 1424-1436. doi:10.1038/s41591-023-02369-6
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @param ... Other arguments passed on to gsva()'s params argument.
#' @export
#' @examples
#' a = ILB.Sig(expr = gene_table,Response = pdata_table$Response)
#'
ILB.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('MS4A1','ICOSLG','CD27','IGHD','IGHM','CD69','ZFP36L2','CREM','CXCR4',
             'ATP1B3','NR4A2','RGS2','SAT1','EZR','ZNF331','PPP1R14B','GPR183',
             'LY9','YPEL5','JUND','SKIL')

  y = .calculatessGSEA(expr = expr,genes = genes,verbose = verbose,...)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Montironi et al. in 2023.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Montironi et al. in 2023.
#' @references Montironi, Carla et al. “Inflamed and non-inflamed classes of HCC: a revised immunogenomic classification.” Gut vol. 72,1 (2023): 129-140. doi:10.1136/gutjnl-2021-325918
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @param ... Other arguments passed on to gsva()'s params argument.
#' @export
#' @examples
#' a = Inflamed20genes.Sig(expr = gene_table,Response = pdata_table$Response)
#'
Inflamed20genes.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('CCL5','CD2','CD3D','CD48','CD52','CD53','CXCL9','CXCR4','FYB1','GZMA',
             'GZMB','GZMK','IGHG1','IGHG3','LAPTM5','LCP2','PTPRC','SLA','TRAC','TRBC2')

  y = .calculatessGSEA(expr = expr,genes = genes,verbose = verbose,...)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' A function for calculating the average value of target genes of CTLA4 inhibitor in DrugBank.
#'
#' A list of target genes with direct pharmacological action with CTLA4 inhibitor is found from DrugBank, and the average value of their expression is calculated.
#' @references Knox, Craig et al. “DrugBank 6.0: the DrugBank Knowledgebase for 2024.” Nucleic acids research vol. 52,D1 (2024): D1265-D1275. doi:10.1093/nar/gkad976
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = CTLA4.Sig(expr = gene_table,Response = pdata_table$Response)
#'
CTLA4.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('CTLA4')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Hugo et al. in 2016.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Hugo et al. in 2016.
#' @references Hugo, Willy et al. “Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma.” Cell vol. 165,1 (2016): 35-44. doi:10.1016/j.cell.2016.02.065
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = IPRES.Sig(expr = gene_table,Response = pdata_table$Response)
#'
IPRES.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('ANGPT2','AXL','CCL13','CCL2','CCL7','CDH1','FAP','FLT1','IL10','LOXL2','RORA','RORB','RORC',
             'TAGLN','TWIST2','VEGFA','VEGFC','WNT5A')
  x = expr[rownames(expr)%in%genes,]
  x = apply(expr,1,scale)
  x = t(x)
  colnames(x) = colnames(expr)
  y = .calculateMean(expr = x,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Rooney et al. in 2015.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Rooney et al. in 2015.
#' @references Rooney, Michael S et al. “Molecular and genetic properties of tumors associated with local immune cytolytic activity.” Cell vol. 160,1-2 (2015): 48-61. doi:10.1016/j.cell.2014.12.033
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = Juul.2010.Sig (expr = gene_table,Response = pdata_table$Response)
#'
Juul.2010.Sig  <- function(expr,Response=NULL,verbose=FALSE) {
  genes_a = c('BUB1B','CDK1','AURKB','TTK')
  genes_b = c('UGCG','CERT1')

  .checkGene(expr = expr,genes = c(genes_a,genes_b),verbose = verbose)

  if (sum(genes_a%in%rownames(expr))<=1 | sum(genes_b%in%rownames(expr))<=1) {
    stop("Too few genes matched!")
  }

  x_a <- expr[rownames(expr)%in%genes_a,]
  x_b <- expr[rownames(expr)%in%genes_b,]

  y=apply(x_a,2,.gm_mean) - apply(x_b,2,.gm_mean)

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}


#' Function to compute the signature as published by Zhang et al. in 2021.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Zhang et al. in 2021.
#' @references Zhang, Zhihui et al. “m6A regulators as predictive biomarkers for chemotherapy benefit and potential therapeutic targets for overcoming chemotherapy resistance in small-cell lung cancer.” Journal of hematology & oncology vol. 14,1 190. 10 Nov. 2021, doi:10.1186/s13045-021-01173-4
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = m6A.score.Sig(expr = gene_table,Response = pdata_table$Response)
#'
m6A.score.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c("ZCCHC4","IGF2BP3","ALKBH5","YTHDF3","METTL5","G3BP1","RBMX")
  coefs <- c(0.7942,-0.2645,-0.4484,-0.6853,0.4749,0.246,0.0911)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}


#' A function for calculating the average value of target genes of PD1 inhibitor in DrugBank.
#'
#' A list of target genes with direct pharmacological action with PD1 inhibitor is found from DrugBank, and the average value of their expression is calculated.
#' @references Knox, Craig et al. “DrugBank 6.0: the DrugBank Knowledgebase for 2024.” Nucleic acids research vol. 52,D1 (2024): D1265-D1275. doi:10.1093/nar/gkad976
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = PD1.Sig(expr = gene_table,Response = pdata_table$Response)
#'
PD1.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('PDCD1')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' A function for calculating the average value of target genes of Paclitaxel in DrugBank.
#'
#' A list of target genes with direct pharmacological action with Paclitaxel is found from DrugBank, and the average value of their expression is calculated.
#' @references Knox, Craig et al. “DrugBank 6.0: the DrugBank Knowledgebase for 2024.” Nucleic acids research vol. 52,D1 (2024): D1265-D1275. doi:10.1093/nar/gkad976
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = Paclitaxel.Sig(expr = gene_table,Response = pdata_table$Response)
#'
Paclitaxel.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('TUBB1','BCL2','MAP4','MAP2','MAPT')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' A function for calculating the average value of target genes of atezolizumab in DrugBank.
#'
#' A list of target genes with direct pharmacological action with atezolizumab is found from DrugBank, and the average value of their expression is calculated.
#' @references Knox, Craig et al. “DrugBank 6.0: the DrugBank Knowledgebase for 2024.” Nucleic acids research vol. 52,D1 (2024): D1265-D1275. doi:10.1093/nar/gkad976
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = Atezolizumab.Sig(expr = gene_table,Response = pdata_table$Response)
#'
Atezolizumab.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CD274','PDCD1')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Roh et al. in 2017.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Roh et al. in 2017.
#' @references Roh, Whijae et al. “Integrated molecular analysis of tumor biopsies on sequential CTLA-4 and PD-1 blockade reveals markers of response and resistance.” Science translational medicine vol. 9,379 (2017): eaah3560. doi:10.1126/scitranslmed.aah3560
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = Roh.immune.Sig (expr = gene_table,Response = pdata_table$Response)
#'
Roh.immune.Sig  <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('GZMA','GZMB','PRF1','GNLY','HLA-A','HLA-B','HLA-C','HLA-E','HLA-F','HLA-G',
             'HLA-H','HLA-DMA','HLA-DMB','HLA-DOA','HLA-DOB','HLA-DPA1','HLA-DPB1',
             'HLA-DQA1','HLA-DQA2','HLA-DQB1','HLA-DRA','HLA-DRB1','IFNG','IFNGR1',
             'IFNGR2','IRF1','STAT1','PSMB9','CCR5','CCL3','CCL4','CCL5','CXCL9',
             'CXCL10','CXCL11','ICAM1','ICAM2','ICAM3','ICAM4','ICAM5','VCAM1')

  .checkGene(expr = expr,genes = genes,verbose = verbose)

  x <- expr[rownames(expr)%in% genes,]
  y <- apply(x,2,.gm_mean)

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Mezheyeuski et al. in 2023.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Mezheyeuski et al. in 2023.
#' @references Mezheyeuski, Artur et al. “An immune score reflecting pro- and anti-tumoural balance of tumour microenvironment has major prognostic impact and predicts immunotherapy response in solid cancers.” EBioMedicine vol. 88 (2023): 104452. doi:10.1016/j.ebiom.2023.104452
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = SIA.Sig(expr = gene_table,Response = pdata_table$Response)
#'
SIA.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CD8A','C1QA','C1QB','C1QC')
  if (!all(genes%in%rownames(expr))) {
    stop("Too few genes matched!")
  }
  x <- t(expr[genes,])
  y <- x[,1]/rowMeans(x[,-1])

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}


#' Function to compute the signature as published by Cabrita et al. in 2020.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Cabrita et al. in 2020.
#' @references Cabrita, Rita et al. “Tertiary lymphoid structures improve immunotherapy and survival in melanoma.” Nature vol. 577,7791 (2020): 561-565. doi:10.1038/s41586-019-1914-8
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = TLS.Sig(expr = gene_table,Response = pdata_table$Response)
#'
TLS.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CD79B','CD1D','CCR6','LAT','SKAP1','CETP','EIF1AY','RBP5','PTGDS')
  if (sum(genes%in%rownames(expr))<=1) {
    stop("Too few genes matched!")
  }
  x = expr[rownames(expr)%in%genes,]
  x = apply(expr,1,scale)
  x = t(x)
  colnames(x) = colnames(expr)
  y = .calculateMean(expr = x,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Wang et al. in 2007.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Wang et al. in 2007.
#' @references Wang, Yunfei et al. “GBP2 is a prognostic biomarker and associated with immunotherapeutic responses in gastric cancer.” BMC cancer vol. 23,1 925. 2 Oct. 2023, doi:10.1186/s12885-023-11308-0
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = GBP2.Sig(expr = gene_table,Response = pdata_table$Response)
#'
GBP2.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('GBP2')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Niu et al. in 2022.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Niu et al. in 2022.
#' @references Niu, Decao et al. “The epiphany derived from T-cell-inflamed profiles: Pan-cancer characterization of CD8A as a biomarker spanning clinical relevance, cancer prognosis, immunosuppressive environment, and treatment responses.” Frontiers in genetics vol. 13 974416. 11 Aug. 2022, doi:10.3389/fgene.2022.974416
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = CD8A.Sig(expr = gene_table,Response = pdata_table$Response)
#'
CD8A.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('CD8A')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Carter et al. in 2006.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Carter et al. in 2006.
#' @references Carter, Scott L et al. “A signature of chromosomal instability inferred from gene expression profiles predicts clinical outcome in multiple human cancers.” Nature genetics vol. 38,9 (2006): 1043-8. doi:10.1038/ng1861
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = CIN70.Sig(expr = gene_table,Response = pdata_table$Response)
#'
CIN70.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('TPX2','PRC1','FOXM1','CDK1','TGIF2','MCM2','H2AZ1','TOP2A','PCNA','UBE2C',
             'MELK','TRIP13','NCAPD2','MCM7','RNASEH2A','RAD51AP1','KIF20A','CDC45','MAD2L1',
             'ESPL1','CCNB2','FEN1','TTK','CCT5','RFC4','ATAD2','CKAP5','NUP205','CDC20',
             'CKS2','RRM2','ELAVL1','CCNB1','RRM1','AURKB','MSH6','EZH2','CTPS1','DKC1',
             'OIP5','CDCA8','PTTG1','CEP55','H2AX','CMAS','NCAPH','MCM10','LSM4',
             'MTB','ASF1B','ZWINT','PBK','ZWILCH','CDCA3','ECT2','CDC6','UNG',
             'MTCH2','RAD21','ACTL6A','SRSF2','HDGF','NXT1','NEK2',
             'DHCR7','AURKA','NDUFAB1','NEMP1','KIF4A')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Ertel et al. in 2010.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Ertel et al. in 2010.
#' @references Ertel, Adam et al. “RB-pathway disruption in breast cancer: differential association with disease subtypes, disease-specific prognosis and therapeutic response.” Cell cycle (Georgetown, Tex.) vol. 9,20 (2010): 4153-63. doi:10.4161/cc.9.20.13454
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = RB.loss.Sig(expr = gene_table,Response = pdata_table$Response)
#'
RB.loss.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('AEBP2','ANAPC11','ANAPC5','ANGPTL2','ANLN','ANP32B','ARHGAP21','ASF1B',
             'ATM','BIRC5','BRCA1','BRCA2','NCAPH','BUB1','CASP8AP2','CBX2','CBX5',
             'CCNA2','CCNB1','CCNB2','CCNF','CD34','CDC20','CDC25C','CDC45','CDC6',
             'CDCA3','CDCA5','CDCA7','CDCA8','CDK2','CDKN1C','CDKN3','CENPA','CHAF1B',
             'CHEK1','CKAP2','CKS2','DCK','DDIT4','DEK','DLGAP5','DNAJC9','DNMT1','DOK1',
             'DTYMK','E2F1','ECT2','EGR1','EI24','EIF4A2','ERCC5','ETV4','EZH2','FBLN1',
             'FEN1','FEZ1','FOXM1','GMNN','GTSE1','H2AZ1','HAT1','HELB','HIP1R','HMGA2',
             'HMGB2','HMGB3','HMGN1','HMGN2','HNRNPC','HNRNPD','HNRNPR','HNRNPU','INCENP',
             'KIF11','KIF1C','KIF20A','KIF22','KIF23','KIF2C','KLF4','KPNA2','LBR',
             'LIG1','MAD2L1','MCM2','MCM3','MCM4','MCM5','MCM7','MDM2','MKI67','MRE11',
             'MSH2','MTCP1','NAP1L1','NASP','NEK2','NUSAP1','ORC6','PCNA','PDCD6IP',
             'PDGFRA','PERP','PHC1','PHF13','PLK1','PLK4','PLTP','PML','POLD1','PRC1',
             'PRDX4','PRIM1','PSIP1','RACGAP1','RAD21','RAD51','RAD51AP1','RBBP4','RBL1',
             'REV3L','RFC2','RFC5','ARHGEF28','RPA1','RRM1','RRM2','SIVA1','SLBP','SLK',
             'SMC2','SMC4','STMN1','TACC3','TARDBP','TCF19','TFDP1','TK1','TMPO',
             'TOP2A','TOPBP1','TRIP13','TTK','TUBGCP3','TYMS','UHRF1','UPF3B','USP1',
             'WAC','WDHD1','ZMAT3','XRCC1','AURKA','AURKB','CCNE1','CCNE2','CDT1','DBF4')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Lee et al. in 2021.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Rooney et al. in 2021.
#' @references Lee, Joo Sang et al. “Synthetic lethality-mediated precision oncology via the tumor transcriptome.” Cell vol. 184,9 (2021): 2487-2502.e13. doi:10.1016/j.cell.2021.03.030
#' @importFrom pROC roc
#' @importFrom stats quantile
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = SELECT1.Sig (expr = gene_table,Response = pdata_table$Response)
#'
SELECT1.Sig  <- function(expr,Response=NULL,verbose=FALSE) {
  genes_a = c('CXCL16','IL15RA','CD27','TNFRSF13C','TNFRSF13B','ICAM4','CD4','CD8A',
              'LTBR','IFITM2')
  genes_b = c('PDCD1','CD274')

  .checkGene(expr = expr,genes = c(genes_a,genes_b),verbose = verbose)

  if (sum(genes_a%in%rownames(expr))<=1 | sum(genes_b%in%rownames(expr))<1) {
    stop("Too few genes matched!")
  }

  x_a <- expr[rownames(expr)%in%genes_a,]
  x_b <- expr[rownames(expr)%in%genes_b,]

  t1 = apply(x_a, 1, function(g){stats::quantile(g,1/3)})
  x_a<-ifelse(x_a<t1,1,0)
  t2 = apply(x_b, 1, function(g){stats::quantile(g,0.3)})
  x_b[x_b<t2]=0

  y=colMeans(x_b)*colMeans(x_a)

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Lee et al. in 2021.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Rooney et al. in 2021.
#' @references Lee, Joo Sang et al. “Synthetic lethality-mediated precision oncology via the tumor transcriptome.” Cell vol. 184,9 (2021): 2487-2502.e13. doi:10.1016/j.cell.2021.03.030
#' @importFrom pROC roc
#' @importFrom stats quantile
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = SELECT3.Sig (expr = gene_table,Response = pdata_table$Response)
#'
SELECT3.Sig  <- function(expr,Response=NULL,verbose=FALSE) {
  genes_a = c('NFKB1','ZNF667','NMNAT1','GPR83','DAPK1','ELOVL7','TINAGL1','SLC5A8',
              'PCDHGA7','ZNF470','PIWIL4','SNX25','ZNF471','PCDHGA10','ZNF300','FBXO42',
              'GALNT15','CPM','ECHDC3','RBM7','POU2F2','CA13','SOX10','ARHGAP6','H6PD')
  genes_b = c('TOP2A')

  .checkGene(expr = expr,genes = c(genes_a,genes_b),verbose = verbose)

  if (sum(genes_a%in%rownames(expr))<=1 | sum(genes_b%in%rownames(expr))<1) {
    stop("Too few genes matched!")
  }

  x_a <- expr[rownames(expr)%in%genes_a,]
  x_b <- expr[rownames(expr)%in%genes_b,]

  t1 = apply(x_a, 1, function(g){stats::quantile(g,1/3)})
  x_a<-ifelse(x_a<t1,1,0)
  t2 = apply(x_b, 1, function(g){stats::quantile(g,0.3)})
  x_b[x_b<t2]=0

  y=colMeans(x_b)*colMeans(x_a)

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Lee et al. in 2021.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Rooney et al. in 2021.
#' @references Lee, Joo Sang et al. “Synthetic lethality-mediated precision oncology via the tumor transcriptome.” Cell vol. 184,9 (2021): 2487-2502.e13. doi:10.1016/j.cell.2021.03.030
#' @importFrom pROC roc
#' @importFrom stats quantile
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = SELECT4.Sig (expr = gene_table,Response = pdata_table$Response)
#'
SELECT4.Sig  <- function(expr,Response=NULL,verbose=FALSE) {
  genes_a = c('HNRNPR','FAM83G','DNAJC12','RALBP1','TXK','CSNK1G3','RPL4','NADK',
              'MYBPC2','AGTR2','SLC7A6','ACOX2','SOD2','STXBP1','LIMK2','KCTD10',
              'ARHGDIG','DIRAS1','POLA2','SLC2A5','SLC9A9','PMP22','SOD2','ECHDC3',
              'KIRREL1')
  genes_b = c('FOLR3','FOLR1','TOP1','TYMS','TOP1MT')

  .checkGene(expr = expr,genes = c(genes_a,genes_b),verbose = verbose)

  if (sum(genes_a%in%rownames(expr))<=1 | sum(genes_b%in%rownames(expr))<1) {
    stop("Too few genes matched!")
  }

  x_a <- expr[rownames(expr)%in%genes_a,]
  x_b <- expr[rownames(expr)%in%genes_b,]

  t1 = apply(x_a, 1, function(g){stats::quantile(g,1/3)})
  x_a<-ifelse(x_a<t1,1,0)
  t2 = apply(x_b, 1, function(g){stats::quantile(g,0.3)})
  x_b[x_b<t2]=0

  y=colMeans(x_b)*colMeans(x_a)

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Lee et al. in 2021.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Rooney et al. in 2021.
#' @references Lee, Joo Sang et al. “Synthetic lethality-mediated precision oncology via the tumor transcriptome.” Cell vol. 184,9 (2021): 2487-2502.e13. doi:10.1016/j.cell.2021.03.030
#' @importFrom pROC roc
#' @importFrom stats quantile
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = SELECT6.Sig (expr = gene_table,Response = pdata_table$Response)
#'
SELECT6.Sig  <- function(expr,Response=NULL,verbose=FALSE) {
  genes_a = c('NFKB1','ZNF667','NMNAT1','GPR83','DAPK1','ELOVL7','SMAP2','TINAGL1',
              'SLC5A8','PCDHGA7','ZNF470','PIWIL4','ZNF471','ERI1','KCNH4','AADACL2',
              'RPL11','RFX2','CPM','ECHDC3','RBM7','POU2F2','STXBP1','SOX10','AGTPBP1')

  genes_b = c('TOP2A','RRM1')

  .checkGene(expr = expr,genes = c(genes_a,genes_b),verbose = verbose)

  if (sum(genes_a%in%rownames(expr))<=1 | sum(genes_b%in%rownames(expr))<1) {
    stop("Too few genes matched!")
  }

  x_a <- expr[rownames(expr)%in%genes_a,]
  x_b <- expr[rownames(expr)%in%genes_b,]

  t1 = apply(x_a, 1, function(g){stats::quantile(g,1/3)})
  x_a<-ifelse(x_a<t1,1,0)
  t2 = apply(x_b, 1, function(g){stats::quantile(g,0.3)})
  x_b[x_b<t2]=0

  y=colMeans(x_b)*colMeans(x_a)

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Mariathasan et al. in 2018.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Mariathasan et al. in 2018.
#' @references Mariathasan, Sanjeev et al. “TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells.” Nature vol. 554,7693 (2018): 544-548. doi:10.1038/nature25501
#' @importFrom pROC roc
#' @importFrom stats prcomp
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = Pan_F_TBRS.Sig(expr = gene_table,Response = pdata_table$Response)
#'
Pan_F_TBRS.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('ACTA2','ACTG2','ADAM12','ADAM19','CNN1','COL4A1','CCN2','CTPS1','RFLNB',
             'FSTL3','HSPB1','IGFBP3','PXDC1','SEMA7A','SH3PXD2A','TAGLN','TGFBI','TNS1','TPM1')

  .checkGene(expr = expr,genes = genes,verbose = verbose)

  x = t(expr[rownames(expr)%in%genes,])
  x = scale(x)
  pca1 <- stats::prcomp(x)
  y <- pca1$x[,1]

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}


#' Function to compute the signature as published by Mariathasan et al. in 2018.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Mariathasan et al. in 2018.
#' @references Mariathasan, Sanjeev et al. “TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells.” Nature vol. 554,7693 (2018): 544-548. doi:10.1038/nature25501
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @param ... Other arguments passed on to gsva()'s params argument.
#' @export
#' @examples
#' a = T.Persistence.Sig(expr = gene_table,Response = pdata_table$Response)
#'
T.Persistence.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  gene_up=c('GBP1','GTSE1','CTSC','FBXO5','GMNN','GBP5','CXCL13',
            'PDIA6','SNRPA','PLA2G7','CASP8AP2','TCF3','ZBED2',
            'BATF2','TFDP1','IGKV4_1','ELF5','BTN2A1','H3C3',
            'GTF2E2','SMARCD2','SLC35B1','CRLF1','MAP3K3')
  gene_down=c('SUSD4','SDC4','EIF4G3','LAMP5','CTTN','FAN1',
              'HGSNAT','KIF3A','DOCK1','ULK2','ARHGEF12','FLNB',
              'CGA','ERLIN2','CELSR2','BAIAP3','SKP1','HMGCL',
              'PURA','FBXL16','SPATA18','BHLHE40','GPRC5A',
              'SLC7A2','MATN3')

  predict_result = .predict_sig(expr = expr,gene_up=gene_up,gene_down = gene_down,
                                predict_name = "T",method = "ssGSEA",verbose=FALSE,
                                nCores = 1,normAUC = TRUE,...)

  y = predict_result$PredictT
  names(y) = rownames(predict_result)

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}


#' Function to compute the signature as published by Risi et al. in 2019.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Risi et al. in 2019.
#' @references Risi, Emanuela et al. “An RB-1 loss of function gene signature as a tool to predict response to neoadjuvant chemotherapy plus anti-HER2 agents: a substudy of the NeoALTTO trial (BIG 1-06).” Therapeutic advances in medical oncology vol. 11 1758835919891608. 11 Dec. 2019, doi:10.1177/1758835919891608
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = RB.Sig(expr = gene_table,Response = pdata_table$Response)
#'
RB.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CDC45','TPX2','CDCA5','CCNA2','MYBL2','UBE2C','CDCA3','MCM10','FOXM1',
             'BIRC5','CENPA','AURKB','KIF2C','CDCA8','TICRR','ORC1','PLK1','EXO1',
             'RAD54L','CEP55','CDC20','CENPI','AURKA','TROAP','POLQ','KIF4A','CLSPN',
             'BUB1','BLM','RAD51','CCNB2','CENPO','CDT1','DEPDC1B','NCAPH','BUB1B',
             'CENPE','KIF4B','PKMYT1','DLGAP5','MELK','KIF20A','PTTG1','TRIP13','GTSE1',
             'PTTG3P','SPC25','CDC25A','CENPM','PTTG2','FAM83D','CENPN','ORC6','CHEK1',
             'CDKN3','KIF11','MTFR2','KIF23','PIF1','NDC80','ASPM','DEPDC1','NEIL3',
             'SGO1','KIFC1','MCM7','MKI67','TTK','KIF14','OIP5','NUSAP1','ASF1B',
             'PIMREG','MND1','STIL','RRM2','PRC1','ANLN','FANCI','SKA1','SKA3','MCM4',
             'ARHGAP11A','KIF15','AUNIP','CENPW','NUF2')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Zhu et al. in 2020.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Zhu et al. in 2020.
#' @references Zhu, Jieqiang et al. “Cancer genomics predicts disease relapse and therapeutic response to neoadjuvant chemotherapy of hormone sensitive breast cancers.” Scientific reports vol. 10,1 8188. 18 May. 2020, doi:10.1038/s41598-020-65055-4
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = Zhu.2020.Sig(expr = gene_table,Response = pdata_table$Response)
#'
Zhu.2020.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('TOP2A','AURKA','CKS2','CCNB2','CDK1','SLC19A1','E2F8','E2F1','PRC1','KIF11')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Oshi et al. in 2022.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Oshi et al. in 2022.
#' @references Oshi, Masanori et al. “Development of a novel BRCAness score that predicts response to PARP inhibitors.” Biomarker research vol. 10,1 80. 12 Nov. 2022, doi:10.1186/s40364-022-00427-8
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @param ... Other arguments passed on to gsva()'s params argument.
#' @export
#' @examples
#' a = BRCAness.Sig(expr = gene_table,Response = pdata_table$Response)
#'
BRCAness.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('ADSL','AP2M1','ATP1A4','CBX2','CCDC125','CDC45','CDCA2','CENPO','CREB3L2',
             'DNAJB11','ECE2','FAM136A','GORASP2','JOSD1','KCTD3','KRT222','LRRC28',
             'MAOA','MAPT','MCM5','OTX1','POLG','POLR2H','PPP1R14B','PRC1','PRSS53',
             'PSAT1','PSMD2','PTPRT','SNRPA1','ST3GAL4','TMEM150B','TRIP13','WARS1')

  y = .calculateGSVA(expr = expr,genes = genes,verbose = verbose,...)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' A function for calculating the average value of target genes of etoposide in DrugBank.
#'
#' A list of target genes with direct pharmacological action with etoposide is found from DrugBank, and the average value of their expression is calculated.
#' @references Knox, Craig et al. “DrugBank 6.0: the DrugBank Knowledgebase for 2024.” Nucleic acids research vol. 52,D1 (2024): D1265-D1275. doi:10.1093/nar/gkad976
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = TOP2A.Sig(expr = gene_table,Response = pdata_table$Response)
#'
TOP2A.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('TOP2A')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' Function to compute the signature as published by Liu et al. in 2018.
#'
#' This function computes signature scores from gene expression values following the algorithm used for the signature as published by Liu et al. in 2018.
#' @references Liu, Gang et al. “Seven Genes Based Novel Signature Predicts Clinical Outcome and Platinum Sensitivity of High Grade IIIc Serous Ovarian Carcinoma.” International journal of biological sciences vol. 14,14 2012-2022. 3 Nov. 2018, doi:10.7150/ijbs.28249
#' @importFrom pROC roc
#' @param expr A matrix or data.frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param Response A vector containing "R" and "NR" that matches the samples, used to calculate AUC. Default to NULL, and AUC will not be calculated.
#' @param verbose Gives information about each calculation step
#' @export
#' @examples
#' a = Liu.2018.Sig(expr = gene_table,Response = pdata_table$Response)
#'
Liu.2018.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CD3D','RAB36','GADD45G','CXCL13','TSPAN13','CXCL9','IL2RG')
  coefs <- c(-0.1074,-0.2262,-0.1740,-0.1017,-0.1197,-0.1933,-0.0706)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}
