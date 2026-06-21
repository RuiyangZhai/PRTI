#' @noRd
GEP.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CCL5','CD27','CD274','CD276','CD8A','CMKLR1','CXCL9','CXCR6','HLA-DQA1',
             'HLA-DRB1','HLA-E','IDO1','LAG3','NKG7','PDCD1LG2','PSMB10','STAT1','TIGIT')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
GBP2.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('GBP2')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
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

#' @noRd
IFN.Y.28.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('IL2RG','CXCR6','CD3D','CD2','ITGAL','TAGAP','CIITA','HLA-DRA','PTPRC',
             'CXCL9','CCL5','NKG7','GZMA','PRF1','CCR5','CD3E','GZMK','IFNG','HLA-E',
             'GZMB','PDCD1','SLAMF6','CXCL13','CXCL10','IDO1','LAG3','STAT1','CXCL11')

  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
CTLA4.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('CTLA4')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
Atezolizumab.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CD274','PDCD1')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
APM.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2','TAPBP')

  y = .calculatessGSEA(expr = expr,genes = genes,verbose = verbose,...)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
IFNy4gene.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('IFNG','CD274','LAG3','CXCL9')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
PDL1.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('CD274')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
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

#' @noRd
CD8A.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('CD8A')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
IFN.Y.10.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('IFNG','STAT1','CCR5','CXCL9','CXCL10','CXCL11','IDO1','PRF1','GZMA','HLA-DRA')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
TSE.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene_up <- c('IDO1','CXCL10','CXCL9','HLA-DRA','STAT1','IFNG','CD3D','IDO1','CIITA',
               'CD3E','CCL5','GZMK','CD2','HLA-DRA','CXCL13','IL2RG','NKG7','HLA-E','CXCR6',
               'LAG3','TAGAP','CXCL10','STAT1','GZMB','TIGIT','CD27','CD8A','PDCD1LG2','LAG3',
               'CD274','CXCR6','CMKLR1','NKG7','CCL5','PSMB10','IDO1','CXCL9','HLA-DQA1',
               'CD276','STAT1','HLA-DRB1','HLA-E','CD8A','CCL2','CCL3','CCL4','CXCL9','CXCL10',
               'ICOS','GZMK','IRF1','HLA-DMA','HLA-DMB','HLA-DOA','HLA-DOB','IFNG','CXCL9',
               'CD8A','GZMA','GZMB','CXCL10','PRF1','TBX21','HSPA1A','HSPA1B','HSPB1','DNAJB1',
               'GLRX','ACTB','HSP90AA1','SOD1','HSP90AB1','KLRD1','PHLDA1','ENTPD1','TPI1',
               'TUBA4A','HSPD1','CD69','DNAJB4','JUN','ADGRE5','HOPX','UBB','HSPA6','NDUFC2',
               'SPOCK2','CBX3','GIMAP7','SELENOH','LITAF','CALM3','PRDX1','GLUL','MRPL18','ATP5MC1',
               'ATP5F1C','LDHB','HSPH1','PPP2R5C','BAG3','NFKBIA','RHOB','ATP5F1A','ECH1','HNRNPA2B1',
               'ZFAND2A','CCL20','PSMA7','PRF1','RAC2','SMS','CLK1','GSTO1','STAT3','YWHAB','EID1',
               'HNRNPC','PTGER2','ATP5PO','TBCB','SUMO2','SLTM','SEC61B','LSP1','TMCO1','IL2RB',
               'LAMTOR2','S100A11','FAM177A1','EPSTI1','CNIH1','CCT3','ATP5PF','DNAJA4','TNFRSF1B',
               'JAK1','TXN2','G0S2','STING1','POLR2I','ITGAE','RNF181','GEM','HPRT1','ZNF706','ASPM',
               'TMEM71','CTSA','MAST4','THEMIS2','CDKN2A','BUB1','TRIB1','MICOS10','NDUFA6','CD82',
               'SAR1B','GADD45B','SGO2','UQCRC1','CHURC1','SOCS1','KIFC1','FBXO43','CENPO','TANK',
               'TMBIM6','TRAF1','ZNF544','PAK2','MDFIC','F2R','WDR1','ADPGK','RNF145','SLA','SRRT',
               'FAM174C','MRPS26','TBC1D20','STOM','PRMT1','COX16','ZNF440','SON','TRAF3IP3','TMED9',
               'HSD17B10','PRKCB','DUSP1','NDUFB5','LINC01480','COX6A1','CD27-AS1','ATP5PD','KARS1',
               'HLA-C','CD53','IMMP1L','ATP5MC2','INPP5F','ANKRD10','IFIT2','SSBP1','HLA-DRB1','NAMPT',
               'NFAT5','PPP1R12A','CRACDL','CCDC167','OXR1','ISOC1','ARL6IP1','GLO1','FIP1L1','SYNCRIP',
               'N4BP2','METTL23','KMT2E-AS1','POLR2G','PTPN22','LAX1','PRDM2','VDAC1','LYRM4','ISCU',
               'CYC1','LBH','DSTN','CCNG2','SH3KBP1','GNA15','MRPS21','SAMD9','SMAP2','ZC3H7A','TUBB2A',
               'CCT6A','C12orf75','TROAP','VAPA','GRB2','RSL1D1','MRPL4','BAZ1B','CCDC32','USP51',
               'AC026304.1','GPR171','SLC7A5','CSNK1A1','TIMP1','PRRG4','ZCCHC17','SEC11A','PIP4K2A',
               'CSTF2','PYURF','RASAL3','STIP1','CCDC74A','IL18R1','ELL2','AURKAIP1','HDAC1','CXCL13',
               'NR4A1','POLDIP2','HNRNPA3','ZC3HAV1L','HEY1','RTF1','ADAMTS17','MRPL42','LAMP1','DGKE',
               'CRYBG1','POLR2E','GGCT','PHACTR2','MAP3K8','NME3','ARPC5','GET3','KNL1','PARP14','PDAP1',
               'KDSR','ETV7','RWDD1','TRIM69','ID1','MT2A','NEDD9','AP2S1','GIMAP4','MPHOSPH8','PDCD6',
               'CGAS','EEF1G','OTULIN','ADPRM','DNAJA1','SRA1','ANKRD28','S1PR4','RNF114','BTLA','TAF7',
               'AL451085.2','MRPL37','PARP1','ADGRG1','DENND10','CAPZA2','PTGER4','INCENP','HIBCH','KIF4A',
               'CSF1','VIM','TYMS','TGIF1','C10orf67','KIF20B','RFK','MFNG','PRKAR1B','TNFSF10','ABCB1',
               'ARSB','AQP3','ORC1','MRPL54','OAT','BTBD10','AKT3','FAM72D','NEK11','DOK5','TRIM16','SERPINH1',
               'CCR7','CD3D','CD3E','CD3G','CD4','SELL','BCL6','CD3D','CD3E','CD3G','CD4','CD40LG','CD84',
               'CXCR5','ICOS','IL6R','PDCD1','SLAMF1','STAT3','TNFSF4','CCR1','CCR5','CD4','CSF2','CXCR3',
               'DPP4','HAVCR2','IFNA1','IFNGR1','IL2','KLRD1','TNF','TNFSF11','CCR3','CCR4','CCR7','CCR8',
               'CD4','CSF2','CXCR4','GATA3','HAVCR1','ICOS','IL10','IL13','IL1R1','IL4','IL5','IL6','PTGDR2',
               'CD3D','CD3E','CD3G','CD4','GATA3','IRF4','STAT6','CCR4','CCR6','CD38','CD3D','CD3E','CD3G',
               'CD4','IL17A','IL17F','IL1R1','IL21','IL22','KLRB1','LINC-ROR','STAT3','AHR','CCR10','CCR4',
               'CCR6','CD3D','CD3E','CD3G','CCR4','CD4','CNGB1','CTLA4','ENTPD1','FOXP3','IKZF2','IL2RA','ISG20',
               'ITGAE','LAG3','LRRC32','NT5E','SELL','TNFRSF18','TNFRSF4','CCL3','CCL4','CD2','CD3D','CD3E',
               'CD3G','CD8A','CD8B','CST7','GZMA','GZMB','IFNG','NKG7','PRF1','CCR7','LEF1','SELL','TCF7','ANXA2',
               'ANXA5','CD63','CD81','EFNB1','EFNB2','EFNB3','ICAM1','ICAM2','LGALS1','LGALS3','PECAM1','PF4',
               'TSPAN10','VCAM1','LTB4R','CCL1','CCL2','CCL4','CCL5','CCL7','CXCL10','CXCL11','CXCL13','CXCL16',
               'CXCL3','CXCL9','CXCR3','CXCR5','CXCR6','PF4','BTLA','VSIR','CD160','CD274','CTLA4','HAVCR2','LAG3',
               'LAIR1','LGALS9','PDCD1','PDCD1LG2','PVR','PVRIG','TIGIT','TNFRSF14','CD226','CD27','CD28','CD40',
               'CD40LG','CD80','CD86','ICOS','ICOSLG','NECTIN2','TNFRSF25','TNFRSF4','TNFRSF9','TNFSF4','TNFSF9',
               'CD14','CD302','CD68','CD74','GZMA','GZMB','GZMH','GZMK','GZMM','ITGAX')
  gene_down <- c('ACTA2','ACTG2','ADAM12','ADAM19','CNN1','COL4A1','CCN2','CTPS1','RFLNB','FSTL3','HSPB1','IGFBP3',
                 'PXDC1','SEMA7A','SH3PXD2A','TAGLN','TGFBI','TNS1','TPM1','COL1A1','COL1A2','COL6A1','COL6A2',
                 'COL6A3','DCN','FAP','THY1','DCN','PAPPA','SFRP4','THBS2','LY86','CXCL14','FOXF1','COL10A1','ACTG2',
                 'APBB1IP','SH2D1A','SULF1','MSR1','C3AR1','FAP','PTGIS','ITGBL1','BGN','CXCL12','ECM2','FCGR2A',
                 'MS4A4A','CCN4','COL1A2','MS4A6A','EDNRA','VCAM1','ADGRA2','SCUBE2','AIF1','HEPH','LUM','PTGER3',
                 'RUNX1T1','CDH5','PIK3R5','RAMP3','LDB2','COX7A1','EDIL3','DDR2','FCGR2B','PLPPR4','COL15A1','AOC3',
                 'ITIH3','FMO1','PRKG1','PLXDC1','VSIG4','COL6A3','SGCD','COL3A1','F13A1','OLFML1','IGSF6','COMP',
                 'HGF','GIMAP5','ABCA6','ITGAM','MAF','ITM2A','CLEC7A','ASPN','LRRC15','ERG','CD86','TRAT1','COL8A2',
                 'TCF21','CD93','CD163','GREM1','LMOD1','TLR2','ZEB2','C1QB','KCNJ8','KDR','CD33','RASGRP3','TNFSF4',
                 'CCR1','CSF1R','BTK','MFAP5','MXRA5','ISLR','ARHGAP28','ZFPM2','TLR7','ADAM12','OLFML2B','ENPP2',
                 'CILP','SIGLEC1','SPON2','PLXNC1','ADAMTS5','SAMSN1','CH25H','COL14A1','EMCN','RGS4','PCDH12',
                 'RARRES2','CD248','PDGFRB','C1QA','COL5A3','IGF1','SP140','TFEC','TNN','ATP8B4','ZNF423','FRZB',
                 'SERPING1','ENPEP','CD14','DIO2','FPR1','IL18R1','HDC','NME8','PDE2A','RSAD2','ITIH5','FASLG',
                 'MMP3','NOX4','WNT2','LRRC32','CXCL9','TENM4','FBLN2','EGFL6','IL1B','SPON1','CD200','FLNA','LOX',
                 'CTHRC1','COL6A3','MFAP5','TNC','FOXC2','COL6A2','CALD1','FBN1','LRRC7','IGF1','FN1','SERPINE2',
                 'ECM1','LAMA2','VIM','CALU','PTGER3','EMP3','ACTA2','EDNRA','EDNRB','FASLG','FGFR1','FGFR2','FGFR3',
                 'FGFR4','FLT1','FLT4','KDR','TNFSF10','TNFSF12','COL10A1','COL11A1','COL1A1','COL3A1','COL4A4',
                 'COL5A2','FN1','ITGA7','LAMB1','LAMB2','LAMB3','NOS1','NOS2','NOS3','PDGFRA','PDPN','SQSTM1',
                 'TNC','VIM','SKP1P2','ACVR1','ACVR1C','ACVR2A','ACVR2B','ACVRL1','AMH','AMHR2','BMP2','BMP4',
                 'BMP5','BMP6','BMP7','BMP8A','BMP8B','BMPR1A','BMPR1B','BMPR2','CDKN2B','CHRD','COMP','CREBBP',
                 'CUL1','DCN','E2F4','E2F5','EP300','FST','GDF5','GDF6','GDF7','ID1','ID2','ID3','ID4','IFNG',
                 'INHBA','INHBB','INHBC','INHBE','LEFTY1','LEFTY2','LTBP1','MAPK1','MAPK3','MYC','NODAL','NOG',
                 'PITX2','PPP2CA','PPP2CB','PPP2R1A','PPP2R1B','RBL1','RBL2','RBX1','RHOA','ROCK1','ROCK2','RPS6KB1',
                 'RPS6KB2','SKP1','SMAD1','SMAD2','SMAD3','SMAD4','SMAD5','SMAD6','SMAD7','SMAD9','SMURF1','SMURF2',
                 'SP1','TFDP1','TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','THBS1','THBS2','THBS3','THBS4','TNF',
                 'ZFYVE16','ZFYVE9')

  gene_up_num = sum(gene_up%in%rownames(expr))
  gene_down_num = sum(gene_down%in%rownames(expr))
  gene_up_percent = round((gene_up_num/length(gene_up))*100,1)
  gene_down_percent = round((gene_down_num/length(gene_down))*100,1)

  if (verbose) {
    cat(gene_up_num," up-regulated genes in the ",'TSE.Sig',"_gene_list are found ","(",gene_up_percent,"%).\n",sep ="")
    cat(gene_down_num," down-regulated genes in the ",'TSE.Sig',"_gene_list are found ","(",gene_down_percent,"%).\n\n",sep ="")
  }
  if (gene_up_num<=1 | gene_down_num<=1) {
    stop("Too few genes matched!")
  }
  if (gene_up_percent!=100 | gene_down_percent!=100) {
    warning("The prediction results may be inaccurate due to the lack of genes!")
  }

  expr_up = expr[rownames(expr)%in%gene_up,]
  expr_down = expr[rownames(expr)%in%gene_down,]

  y=colMeans(expr_up,na.rm = T)-colMeans(expr_down,na.rm = T)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
Hwang.2020.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('CBLB','CCR7','CD27','CD48','FOXO1','FYB1','HLA-B','HLA-G','IFIH1','IKZF4','LAMP3','NFKBIA','SAMHD1')

  y = .calculatessGSEA(expr = expr,genes = genes,verbose = verbose,...)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
Inflamed20genes.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('CCL5','CD2','CD3D','CD48','CD52','CD53','CXCL9','CXCR4','FYB1','GZMA',
             'GZMB','GZMK','IGHG1','IGHG3','LAPTM5','LCP2','PTPRC','SLA','TRAC','TRBC2')

  y = .calculatessGSEA(expr = expr,genes = genes,verbose = verbose,...)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
NI.score.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('SLC16A3','PLA2G2A','SLC2A1')
  coefs <- c(0.0009901,0.00005335,0.006663)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
PD1.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('PDCD1')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
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

#' @noRd
HOT.score.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('CCL19','CCR2','CCR4','CCR5','CD27','CD40LG','CD8A','CXCL10','CXCL11',
             'CXCL13','CXCL9','CXCR3','CXCR6','FASLG','FGL2','GZMA','GZMH','IDO1',
             'IFNG','IRF8','LAG3','LYZ','MS4A1','PDCD1','TBX21','TLR7','TLR8')

  y = .calculateGSVA(expr = expr,genes = genes,verbose = verbose,...)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
CD96.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('CD96')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
IMPRES.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes_a <- c('PDCD1','CD27','CTLA4','CD40','CD86','CD28','CD80','CD274','CD86','CD40',
               'CD86','CD40','CD28','CD40','TNFRSF14')
  genes_b <- c('TNFSF4','PDCD1','TNFSF4','CD28','TNFSF4','CD86','TNFRSF9','VSIR','HAVCR2',
               'PDCD1','CD200','CD80','CD276','CD274','CD86')

  y = .calculateREO(expr = expr,genes_a = genes_a,genes_b = genes_b,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

