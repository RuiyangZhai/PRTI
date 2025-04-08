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
Traf3_KO.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  genes <- c('HLA-A','STUM','MAP3K14','NFKB2','NFKBIE','IRF8','KIAA1522','CXCL10','CYBA','LGALS3','CASP4',
           'COQ8A','NFKBIA','SLC2A6','TTN','BTN2A2','PLSCR1','TUBB3','SERPINA3','CD9','FRMPD1','CTSH',
           'GLIPR2','NDNF','RPPH1','ZBP1','CASP12','ANXA1','TNIP1','HPSE','IFNGR2','RELB','ANXA4','NUDT7',
           'RGS16','LSM5','SLC11A2','ECSCR','DNM1','IFI16','TNFRSF9','CCDC122','PSMB9','PSMB10','CSTB','JUNB',
           'TMLHE','TEX30','VAMP5','ZNHIT3','KLF6','PSTK','MED31','VAMP7','RPA3','DNAJB9','ZNF25','NQO1',
           'IFRD1','PLEKHA2','CAPG','CCL15','C1R','ARHGEF40','MGST1','CMTR2','SERPINB9','PYCARD','RAP1B',
           'KBTBD8','MRPL53','CDKN1B','MAD2L1','SPRYD7','ZSWIM7','TMOD1','JAGN1','CCDC58','HSCB','MRPL14',
           'VBP1','CYB561D2','ETFA','MOCOS','FBXL5','ID3','MRPS16','CIB1','CCT4','CNIH4','TMEM42','IL12RB1',
           'PEX7','MTRF1','GSTA5','CCNE2','ATP6V1F','ERGIC2','NAE1','ARHGDIB','PHB2','ISCA2','MRPS24',
           'GADD45GIP1','AMN1','HSPA13','FKBP11','TRNAU1AP','CEP44','MRPL35','NT5C3A','TSEN2','RBM7','ATF3',
           'UBE2V2','SEC13','PIGP','ZNF442','BIRC2','FASTKD2','PCGF1','BOLA3','SAT1','CHAC2','ZNF616','RPP30',
           'DGUOK','HSPE1','BIRC3','TAPBPL','POLR1C','ALG6','COMMD2','NUDCD2','USE1','TMED2','CHCHD6','RTRAF',
           'UFM1','SMS','RFC3','TRNT1','RPAIN','TANK','BET1','TRAF2','MED30','ALKBH7','SAR1B','TAF1D','SP140',
           'CHADL','MRPL22','KIT','PTPN18','PSMA5','DCAF13','TRMT10C','SEC61B','ECHS1','RNF7','NAMPT','PDCD10',
           'FAM104A','ATOX1','PPP1R35','ARL1','TAPBP','C4orf3','CAV2','CNEP1R1','B9D1','SDF2L1','COPS3','MANF',
           'MRPS10','TES','HINT1','EPS8','NFU1','GBP6','PEX13','TMEM219','APH1B','YIPF5','SSR4','HLA-E','LSM7',
           'PGRMC1','SELENOS','GNG10','C20orf144','B9D2','FGFR1OP2','NT5C','IDI1','KLHDC2','CNBP','RPL37','SERP1',
           'ITGA4','WWC1','GPATCH8','ARHGEF17','GAS7','EIF4G1','CYP26B1','GTF2I','MEF2D','FARP1','PKNOX2','NOLC1',
           'ARHGAP32','SETD1B','ENG','NOL6','PDXDC1','MTR','IK','CHSY1','SLC35F1','PHLDB2','TRAK1','DGKZ','MACF1',
           'SORT1','PLD1','ANKRD50','ARID5B','ATRNL1','DLC1','IRF2BPL','SARDH','CEP250','PRUNE2','RASGRP3','FLNA',
           'MARVELD1','PLXNB2','PTPN14','NEO1','CRIM1','CHN1','SCD','SAP130','PPP1R37','TRIM56','ZMAT3','LTBP3',
           'PABPC4','SSBP4','MSN','DGKD','ANP32A','ADAM23','ADGRG1','ZNF704','CEP170B','ANKRD34A','CNTFR','MAP4K4',
           'DBN1','PRRC2B','PKD1','PHYHIPL','FADS2','BAG6','NACC2','DAAM2','DCDC2','MEGF8','SORBS3','CDC42EP1',
           'AMOTL1','AFAP1L1','BAHCC1','SNX25','SALL2','DAB2IP','ANP32E','LPGAT1','EIF4EBP2','ARMCX4','CDON','BCOR',
           'DENND1A','RNF44','LAMC1','POU3F2','BRD1','HAP1','ANKRD52','EEF2K','RGS3','FLYWCH1','SMARCC2','PRSS12',
           'MTCH1','MEX3A','MCC','TNS2','PYGO1','SCAMP2','IGFBP4','RUSC2','KCTD17','CA14','PIGT','APLP2','CD276',
           'SORBS1','TNS3','TARBP1','PRDM16','ANGPTL2','SOGA1','PITPNM2','MARCKSL1','HUNK','TULP4','TLE3','C2CD2L',
           'ADAMTS4','TTYH3','ATF5','TNS1','SELENON','MAML3','SDC3','TENM4','ITPKB','NOTCH1','SLC25A23','DAG1',
           'TBC1D16','GJA3','PKDCC','PRSS23','DDAH1','NIBAN2','CEMIP2','TMEM132A','CLMN','CFAP54','SVIL','CORO2B',
           'CKB','R3HDM2','RNF144A','PTCH1','MASP1','FILIP1','CITED1','SDK1','RAB11FIP4','DPYSL3','IGF2R','CSNK1E',
           'SMC1A','PHF2','SIDT2','RASSF2','TUBB2B','BCAN','TCF7L1','COL5A3','PARVB','ABCA12','COL4A2','LMO1',
           'CTDSP2','CLSTN1','LIMCH1','ARFGEF3','SHISA2','CDKN1A','CNTN1','MAN2A2','ALDH2','ATG9A','ADAM12',
           'SLC44A3','MEST','PBXIP1','AKAP12','SULF2','CELSR2','SEMA6A','ITGB5','ENAH','CYFIP2','FRZB','NES',
           'KANK4','CLVS1','KIF26B','COL4A1','MBP','FSCN1','CARMIL2')

  coefs = c(1,0.943375093,0.804272405,0.720120934,0.640698128,0.630424328,0.608157684,0.544272792,0.532159393,
            0.531215795,0.520840774,0.460552502,0.420953427,0.408002948,0.395908256,0.395469859,0.390948623,
            0.390526316,0.387847044,0.379459542,0.376918931,0.374155768,0.367214776,0.365470419,0.361198574,
            0.356761919,0.353745121,0.350766403,0.350524303,0.347488941,0.347132736,0.34434446,0.34141393,
            0.340339804,0.339804148,0.329315815,0.325985135,0.325898119,0.325336397,0.324499612,0.316259904,
            0.314936487,0.309842169,0.309009936,0.308268615,0.307862923,0.304522758,0.302805739,0.301325143,
            0.299064087,0.297791898,0.29392998,0.290829882,0.290258195,0.289316377,0.288238541,0.287360232,
            0.286432405,0.285957059,0.285701457,0.285203091,0.283807482,0.281587645,0.28120642,0.280701615,
            0.280509214,0.28022405,0.279813375,0.279744829,0.279369213,0.278336657,0.277946689,0.277193237,
            0.276350703,0.27624081,0.27613027,0.275596238,0.275544823,0.275413987,0.274487972,0.273682648,
            0.273542408,0.273035924,0.272646625,0.27204669,0.271888312,0.271879474,0.271106003,0.27047889,
            0.270242089,0.269929499,0.268911899,0.268391971,0.267584121,0.267381477,0.266385903,0.266310737,
            0.266188416,0.265640755,0.265008779,0.264913372,0.264752001,0.264340724,0.264299562,0.263589789,
            0.262709062,0.262384345,0.261760328,0.261563015,0.261166676,0.260807396,0.260759435,0.260179497,
            0.259033493,0.25880954,0.258791993,0.258783589,0.258562853,0.25830344,0.258073625,0.257427004,
            0.257324587,0.256786195,0.256779704,0.256392331,0.25635175,0.255731119,0.255686229,0.255669221,
            0.255610234,0.255485249,0.255344119,0.255287236,0.255265598,0.254890794,0.254651687,0.254447155,
            0.254433181,0.254428945,0.25436729,0.253828355,0.253638236,0.253165735,0.253140874,0.252605901,
            0.252542665,0.252423434,0.251870424,0.251081404,0.250926542,0.250891298,0.250503105,0.250500316,
            0.250464886,0.250460472,0.249974794,0.24986619,0.249671982,0.249391234,0.249336993,0.249327792,
            0.24930927,0.249251125,0.249192983,0.248918493,0.248511083,0.248500651,0.247804488,0.247712128,
            0.247505216,0.247101039,0.246670954,0.246252623,0.245977306,0.245895986,0.245794537,0.245474105,
            0.245394328,0.245219106,0.244917122,0.244843738,0.244781808,0.244448993,0.244333521,0.244132583,
            0.244051423,0.243718775,0.243364015,0.242502926,0.242402941,0.241309475,0.240993824,0.24098583,
            0.240860632,0.240438074,0.240425575,0.24028504,0.240080919,0.24002373,0.239874418,-0.188965866,
            -0.189061606,-0.189255496,-0.18941669,-0.189445872,-0.189508736,-0.189531439,-0.18954102,-0.189542685,
            -0.189701888,-0.189962255,-0.190143624,-0.190270726,-0.190448867,-0.190638381,-0.190984312,-0.191157664,
            -0.191678877,-0.191806897,-0.191838283,-0.192295219,-0.192569938,-0.193183908,-0.193302091,-0.193385674,
            -0.193464429,-0.193565451,-0.193741849,-0.194257485,-0.194548542,-0.195263624,-0.19548664,-0.195673022,
            -0.196051138,-0.196165182,-0.196211047,-0.197015369,-0.19703565,-0.197461013,-0.197532185,-0.197856393,
            -0.197857328,-0.198097431,-0.198560116,-0.198577933,-0.198846033,-0.198951262,-0.199290755,-0.200545528,
            -0.20083059,-0.2008471,-0.20091939,-0.200940577,-0.201172658,-0.201208005,-0.201638812,-0.20179751,
            -0.201874701,-0.202209131,-0.2025856,-0.202609771,-0.203590045,-0.20366417,-0.203814593,-0.203884482,
            -0.204188703,-0.204363634,-0.204462095,-0.205169962,-0.20554628,-0.205546771,-0.2055818,-0.206241024,
            -0.206601188,-0.206610622,-0.206740165,-0.207237021,-0.207311413,-0.207692202,-0.207891078,-0.208187191,
            -0.20886227,-0.209273924,-0.209425024,-0.209679368,-0.209699326,-0.209949049,-0.210388527,-0.210390984,
            -0.210692077,-0.211482266,-0.211836637,-0.213041125,-0.213143453,-0.213767017,-0.214245493,-0.214678726,
            -0.214796769,-0.215339437,-0.215410962,-0.216445447,-0.216636226,-0.217039937,-0.217696204,-0.217718074,
            -0.217979148,-0.219295247,-0.219857985,-0.219958871,-0.220507027,-0.220787748,-0.221179354,-0.221675221,
            -0.221877128,-0.222276312,-0.222849923,-0.223471158,-0.223553474,-0.224039041,-0.224505567,-0.225775707,
            -0.22583029,-0.226641115,-0.228078109,-0.230130619,-0.230526481,-0.231093483,-0.231303636,-0.232810115,
            -0.232831782,-0.232895887,-0.233123887,-0.233256061,-0.233721665,-0.234521522,-0.235390503,-0.23567749,
            -0.23596224,-0.236299594,-0.236642363,-0.237254468,-0.237552736,-0.237779275,-0.238526887,-0.238701829,
            -0.241239117,-0.241662212,-0.241851499,-0.24447621,-0.245315105,-0.246055267,-0.24615747,-0.24679295,
            -0.247393833,-0.248018168,-0.248769114,-0.249643361,-0.250016206,-0.251863468,-0.253545529,-0.254824585,
            -0.257269769,-0.258047353,-0.258235492,-0.260981529,-0.261758959,-0.261842775,-0.266858623,-0.267791312,
            -0.269341965,-0.269359499,-0.269914956,-0.276485868,-0.277065175,-0.278285746,-0.279862155,-0.28046863,
            -0.286096727,-0.290928124,-0.293397581,-0.293520169,-0.304598862,-0.305677859,-0.306663425,-0.308581978,
            -0.310296019,-0.314626367,-0.318127354,-0.319658407,-0.320747952,-0.321010065,-0.335054674,-0.338127165,
            -0.352128315,-0.363885683,-0.36704499,-0.438361051,-0.500138526,-0.548936104,-0.567430393)

  y = .calculateWeightedMean(expr = expr,genes = genes,coefs = coefs,verbose = verbose)
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
Hwang.2020.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('CBLB','CCR7','CD27','CD48','FOXO1','FYB1','HLA-B','HLA-G','IFIH1','IKZF4','LAMP3','NFKBIA','SAMHD1')

  y = .calculatessGSEA(expr = expr,genes = genes,verbose = verbose,...)
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
HOT.score.Sig <- function(expr,Response=NULL,verbose=FALSE,...) {
  genes <- c('CCL19','CCR2','CCR4','CCR5','CD27','CD40LG','CD8A','CXCL10','CXCL11',
             'CXCL13','CXCL9','CXCR3','CXCR6','FASLG','FGL2','GZMA','GZMH','IDO1',
             'IFNG','IRF8','LAG3','LYZ','MS4A1','PDCD1','TBX21','TLR7','TLR8')

  y = .calculateGSVA(expr = expr,genes = genes,verbose = verbose,...)
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
CD96.Sig <- function(expr,Response=NULL,verbose=FALSE) {
  gene <- c('CD96')
  y = .calculateSingle(expr = expr,gene = gene,verbose = verbose)
  out_put = .calculateAUC(y = y,Response = Response,verbose = verbose)
  return(out_put)
}

#' @noRd
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

