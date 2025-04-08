#' @noRd
TP53.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('TP53')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
Liu.2018.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CD3D','RAB36','GADD45G','CXCL13','TSPAN13','CXCL9','IL2RG')
  coefs <- c(-0.1074,-0.2262,-0.1740,-0.1017,-0.1197,-0.1933,-0.0706)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
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

#' @noRd
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

#' @noRd
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

#' @noRd
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

#' @noRd
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
#' @noRd
Bai.2022.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CD48','GPR65','C3AR1','CD2','CD3E','ARHGAP9')
  coefs <- c(0.14771103,-0.19725147,-0.2050475,-0.15374589,-0.13619434,-0.1600497)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
Zhu.2020.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('TOP2A','AURKA','CKS2','CCNB2','CDK1','SLC19A1','E2F8','E2F1','PRC1','KIF11')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
CES.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CENPA','CENPN','CENPM','CENPK','CENPL','CENPW','CENPU','HJURP',
             'OIP5','ZWINT','NDC80','SPC24','SPC25','NUF2')
  y = .calculateSum(expr = expr,genes = genes,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
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

#' @noRd
RAD51.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('RAD51')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
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

#' @noRd
TYMS.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('TYMS')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
Tegafur_uracil.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('DPYD','TYMS')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
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
#' @noRd
Chen.2021.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('POLR2B','GAS6','CRY1','BCL2L1','ARG1','ORAI3','TRAF3','ZSWIM4','IRF1',
             'LEMD1','ACTB')
  coefs <- c(-0.172,0.057,-0.027,0.028,0.075,0.01,-0.0004,0.106,-0.080,0.07,0.044)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
Weng.2022.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('PTGES3','PAICS','GNPNAT1','PGM3','MTHFD2','DCK','MTAP','SLC25A36','GBE1','RRM2','KCTD3',
             'ACADSB','ABCD3','BCKDHB','PHOSPHO2','FUT4','EDEM3','NEU4','SLC16A1','ELOVL7','SLC6A8')
  coefs <- c(-0.1053,-0.1874,-0.0133,0.02893,0.08862,-0.0043,0.06428,0.22468,0.00308,0.0679,0.0029,
             -0.0436,-0.0339,-0.0137,-0.0525,-0.1185,0.00089,-0.0684,0.6165,-0.0162,0.04312)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
ILB.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('MS4A1','ICOSLG','CD27','IGHD','IGHM','CD69','ZFP36L2','CREM','CXCR4',
             'ATP1B3','NR4A2','RGS2','SAT1','EZR','ZNF331','PPP1R14B','GPR183',
             'LY9','YPEL5','JUND','SKIL')

  y = .calculatessGSEA(expr = expr,genes = genes,verbose = verbose,...)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
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
