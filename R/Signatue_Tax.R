#' @noRd
Oshi.2021.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CDCA8','MCM2','MCM6','MELK','DEK')
  coefs <- c(1.522,0.780,0.708,1.089,0.694)

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
Zhu.2020.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('TOP2A','AURKA','CKS2','CCNB2','CDK1','SLC19A1','E2F8','E2F1','PRC1','KIF11')
  y = .calculateMean(expr = expr,genes = genes,verbose = verbose)
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
TYMS.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('TYMS')
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
DTG_S.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('ADAM17','AMD1','LSS','NDC80','PLS3','YARS1','ACADVL','CANT1','CDK7','IGF1R',
             'IVD','PDK2')
  coefs <- c(1,1,1,1,1,1,-1,-1,-1,-1,-1,-1)

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
NEPAL.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  gene_up=c("ESPN", "NR0B2", "CLSPN", "TMEM61", "MLLT11", "NUF2", "ASPM", "PROX1", "RRM2", "NRXN1",
            "FAM161A", "FAM178B", "DAPL1", "CRYBA2", "SCG2", "DNER", "KIF1A", "CADPS", "TAGLN3", "PLS1",
            "ECT2", "PEX5L", "UCHL1", "CENPE", "CENPU", "CENPK", "CCNB1", "FOXD1", "PCSK1", "TUBB2B",
            "SCGN", "HOXA11", "HEPACAM2", "VGF", "NRCAM", "SYP", "BEX1", "ST18", "CA8", "CCNE2", "BAALC",
            "RIMS2", "CEL", "LCN15", "NAV2", "FAM111B", "ZWINT", "CDK1", "MKI67", "TROAP", "ESPL1", "ASCL1",
            "SLITRK6", "ZIC2", "SFTA3", "NKX2-1", "TRIM9", "CDKN3", "CHGA", "SCG3", "CCNB2", "PRC1", "ORC6",
            "TOX3", "GINS2", "TOP2A", "BIRC5", "NPTX1", "TYMS", "ADCYAP1", "CDH2", "NOL4", "CELF4", "SYT4",
            "ONECUT2", "CHGB", "INSM1", "FOXA2", "TPX2", "E2F1", "MYBL2", "UBE2C", "UHRF1", "ELAVL3",
            "CACNA1A", "DLL3", "TNNT1", "GTSE1", "NPPB", "INA")
  gene_down=c("SNAP23", "ABCC4", "ANPEP", "IQGAP2", "LTF", "SLC30A4", "TBC1D4", "RDH11", "GLUD1", "PPM1K",
              "MSMB", "KLK4", "SYNGR2", "MIPEP", "SOCS2", "RAB27A", "ADRB2", "STEAP2", "CCNG2", "SORD",
              "ITM2A", "NPY", "CCDC69", "RNF149", "ERGIC1", "ZFP36L2", "OR51E2", "MLPH", "KLK2", "SPRY1",
              "ACKR1", "TGM2", "TPSAB1", "MT1E", "TACC1", "NEDD4L", "SLC39A6", "SMARCA2", "ALDH1A3", "KLK3")

  predict_result = .predict_sig(expr = expr,gene_up=gene_up,gene_down = gene_down,
                                predict_name = "T",method = "ssGSEA",verbose=FALSE,
                                nCores = 1,normAUC = TRUE,...)

  y = predict_result$PredictT
  names(y) = rownames(predict_result)

  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
TOP2A.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('TOP2A')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
BRCAness.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('ADSL','AP2M1','ATP1A4','CBX2','CCDC125','CDC45','CDCA2','CENPO','CREB3L2',
             'DNAJB11','ECE2','FAM136A','GORASP2','JOSD1','KCTD3','KRT222','LRRC28',
             'MAOA','MAPT','MCM5','OTX1','POLG','POLR2H','PPP1R14B','PRC1','PRSS53',
             'PSAT1','PSMD2','PTPRT','SNRPA1','ST3GAL4','TMEM150B','TRIP13','WARS1')

  y = .calculateGSVA(expr = expr,genes = genes,verbose = verbose,...)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
Fu.2021.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('ADAMDEC1','CCL18','CD79A','CD96','CXCL13','DIRAS3','ERBB4','EVL','GAMT',
             'GBP1','GFRA1','GZMB','HSPB8','IGHM','IRS1','ITK','LOC102723479','MAPT',
             'PADI2','RLN2','SEL1L3','SERPINA5','STC2','STK32B','SYBU')
  coefs <- c(0.00321747620626765,0.0457079749167309,8.61152358256599,6.22205851428899,
             -0.585126092824241,-6.08198493202845,1.72908010036751,-1.70368931131805,
             -8.84896004120253,-0.764626193845283,-0.115908259316488,-0.0752619689246736,
             -1.28866942797256,-1.37319937849059,0.250096649476748,-2.30297033083433,
             0.385454564188641,0.286187494306212,0.783128470665541,-1.56204367828805,
             -2.98426861278556,0.25651424658033,0.430345120497431,-1.28399430856461,
             -0.706271090221699)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
DDIR.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('CXCL10','MX1','IDO1','IFI44L','CD2','GBP5','PRAME','ITGAL','LRP4',
             'APOL3','CDR1','FYB1','TSPAN7','RAC2','KLHDC7B','GRB14','MT-RNR1',
             'KIF26A','CD274','CD109','ETV7','MFAP5','OLFM4','PI15','FOSB','TAFA5',
             'NLRC5','PRICKLE1','EGR1','CLDN10','ADAMTS4','SP140L','ANXA1','RSAD2',
             'ESR1','IKZF3','OR2I1P','EGFR','NAT1','LATS2','CYP2B6','PTPRC',
             'PPP1R1A','AL137218.1')
  coefs <- c(0.023,0.0226,0.0221,0.0191,0.019,0.0181,0.0177,0.0176,-0.0159,0.0151,
             -0.0149,0.0149,-0.0148,0.0148,0.014,0.0137,-0.0136,-0.0136,0.0133,
             -0.0129,0.0124,-0.0121,-0.0117,-0.0115,-0.0111,-0.011,0.0101,-0.0089,
             -0.0086,-0.0086,-0.0085,0.0084,-0.0082,0.0081,0.0079,0.0073,0.007,
             -0.0066,0.0065,-0.0063,0.0061,0.0051,-0.0041,-0.0017)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
  y = ifelse(y>=0.3681,1,0)
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
