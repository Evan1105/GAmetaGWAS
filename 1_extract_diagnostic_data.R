#!/usr/bin/Rscript

# TITLE: maternal GWAS meta-analysis of GestAge, 2018-07; JB.
# GOAL: generate initial report for all cohorts about delivered data files
# also, export the data in a RData format, which is faster to laod.
# script name: "1_extract_diagnostic_data.R" (used to be called "1_digest_SNPREPORT.R")

# run from: /home/jonas/meta/
# run like:  ./1_extract_diagnostic_data.R 2>&1 | tee -a met/log_file.txt

library(dplyr)
library(hexbin)

raw_dir = "/home/jonas/meta/raw/"  # here the raw gz files are stored
met_dir = "/home/jonas/meta/met/" # here the SNP data (meta data) will be stored
dia_dir = "/home/jonas/meta/dia/" # here the digested diagnostic reports are stored
#raw_dir = "~/Biostuff/mount_ws2/metaß/raw/" # when local
#met_dir = "~/Biostuff/mount_ws2/meta/met/" # when local
#dia_dir = "~/Biostuff/mount_ws2/meta/dia/" # when local

file_list = list.files(path = raw_dir,pattern =".*RESULTS.*txt.gz",
                       recursive = F,include.dirs = F,full.names = F)

### EXCEPRIONS TO BE OMITED FOR NOW!!! ***
file_list = file_list[-grep("STORK|EGCUT|23andMe",file_list)]
file_list = file_list[which(file_list!="NFBC1966-allPTD-RESULTS-18042018.txt.gz")]
file_list = file_list[which(file_list!="NFBC1966-earlyPTD-RESULTS-18042018.txt.gz")]
file_list = file_list[which(file_list!="NFBC1966-GAnrm-RESULTS-18042018.txt.gz")]
file_list = file_list[which(file_list!="NFBC1966-GAraw-RESULTS-18042018.txt.gz")]
file_list = file_list[which(file_list!="NFBC1966-postTerm-RESULTS-18042018.txt.gz")]

# .. also, pay atention to Viva files which are missing, and STORK files, that are not up to date now.
# .. also, 23andMe files need reformating.


cat("\n\n\n\t\t starting... \n\n")

for (file_name in file_list) { 
        message(file_name)
        
        message("\t loading file..")
        dat = read.table(gzfile(paste(raw_dir,file_name,sep="")), #nrows=1e3, #****
                         h=T,na.strings = c("NA",".","-"),stringsAsFactors = F) 
        
        #col_ix = which(colnames(dat) %in% c("SNPID","RSID","CHR","POS","EFF_ALLELE","REF_ALLELE"))
        #o = dat[,col_ix]
        
        message("\t saving in RData format..")
        RData_file_name = paste(unlist(strsplit(file_name,"\\.txt\\.gz")),".RData",sep="")
        save(list = c("dat"),file = paste(met_dir,RData_file_name,sep="")) # save in .RData data
        
        # 0) CHECK FOR MISSING COLUMNS AND IF NEEDED - CREATE DUMMY REPLACEMENTS
        message("\t step 0: missing columns..")
        if (!"SNPID" %in% colnames(dat)) {
                warning("'SNPID' column was not found. will be replaced with 'RSID'",immediate. = T)
                dat$SNPID = dat$RSID  # assuming RSID EXISTS
        }
        
        if (!"P_VAL_DOM" %in% colnames(dat)) {
                warning("'P_VAL_DOM' column was not found. will be replaced with dummy values",immediate. = T)
                dat$P_VAL_DOM = sample(c(0.2,0.5,0.8),nrow(dat),replace=T) # dummy values
        }
        
        if (!"P_VAL_REC" %in% colnames(dat)) {
                warning("'P_VAL_REC' column was not found. will be replaced with dummy values",immediate. = T)
                dat$P_VAL_REC = sample(c(0.2,0.5,0.8),nrow(dat),replace=T) # dummy values
        }
        
        if (!"INFO" %in% colnames(dat)) {
                warning("'INFO' column was not found. will be replaced with dummy values",immediate. = T)
                dat$INFO = sample(c(0.2,0.5,0.8),nrow(dat),replace=T) # dummy values
        }
        
        if (!"IMPUTED" %in% colnames(dat)) {
                warning("'IMPUTED' column was not found. will be replaced with dummy values",immediate. = T)
                dat$IMPUTED = sample(c(3,4),nrow(dat),replace = T) # dummy values
        }
        
        if (!"N" %in% colnames(dat)) {
                warning("'N' column was not found. will be replaced with dummy values",immediate. = T)
                dat$N = 10 # dummy values
        }
        
        # if columns are not missing, but velues are:
        if (all(is.na(dat$P_VAL_DOM))) dat$P_VAL_DOM = sample(c(0.2,0.5,0.8),nrow(dat),replace=T) # faka data, signals a problem!
        if (all(is.na(dat$P_VAL_REC))) dat$P_VAL_REC = sample(c(0.2,0.5,0.8),nrow(dat),replace=T) # fake data, signals a problem
        if (all(is.na(dat$INFO))) dat$INFO = sample(c(0.2,0.5,0.8),nrow(dat),replace=T) # fake data, signals a problem
        
        # 1) general info
        message("\t step 1: sample the file")
        set.seed(100)
        r1_1 = dat[sample(nrow(dat),5,replace=F),]
        r1_2 = nrow(dat)
        
        # 2) SNPID
        message("\t step 2: SNPID..")
        r2_1 = sum(is.na(dat$SNPID))
        r2_2 = table(nchar(dat$SNPID),useNA = "a")
        r2_3 = table(substr(dat$SNPID,1,2),useNA = "a")
        
        snps = gsub("[0-9]+","num",dat$SNPID) # replace numbers with "num"
        snps = gsub("[A,T,C,G]+","ATCG",snps) # abbreviate 1000G SNP nomenclature
        snps = gsub("<(.*?)>","<..>",snps) # abbreviate 1000G SNP nomenclature
        tbl = table(snps,useNA = "a")
        dtf = data.frame(SNPID_pattern=names(tbl),count=as.numeric(tbl),stringsAsFactors = F)
        dtf$frequency = paste(round(dtf$count / sum(dtf$count) * 100,5),"%",sep="")
        r2_4 = dtf[order(dtf$count,decreasing = T),] # note the "NA" and "<NA>"!
        rm(snps,tbl,dtf)
        
        # 3) RSID
        message("\t step 3: RSID..")
        r3_1 = sum(is.na(dat$RSID))
        r3_2 = table(nchar(dat$RSID),useNA = "a")
        r3_3 = dat[which(nchar(dat$RSID)<=3),][1:10,]
        r3_4 = table(substr(dat$RSID,1,2),useNA = "a") ## beware of character "NA" !
        
        snps = gsub("[0-9]+","num",dat$RSID) # replace numbers with "num"
        snps = gsub("[A,T,C,G]+","ATCG",snps) # abbreviate 1000G SNP nomenclature
        snps = gsub("<(.*?)>","<..>",snps) # abbreviate 1000G SNP nomenclature
        tbl = table(snps,useNA = "a")
        dtf = data.frame(RSID_pattern=names(tbl),count=as.numeric(tbl),stringsAsFactors = F)
        dtf$frequency = paste(round(dtf$count / sum(dtf$count) * 100,5),"%",sep="")
        r3_5 = dtf[order(dtf$count,decreasing = T),] # note the "NA" and "<NA>"!
        rm(snps,tbl,dtf)
        
        # 4) CHR
        message("\t step 4: CHR..")
        r4_1 = table(dat$CHR,useNA = "a")
        #barplot(sort(r4_1,decreasing = T),cex.names = 0.6)
        obs_chr = sort(unique(dat$CHR))
        exp_chr = as.character(1:23)
        r4_2 = obs_chr[which(!obs_chr %in% exp_chr)] # not expected but found chromosome name
        r4_3 = exp_chr[which(!exp_chr %in% obs_chr)] # expected but not found chromosome name
        rm(obs_chr,exp_chr)
        
        # 5) fourth column
        message("\t step 5: POS..")
        r5_1 = sum(is.na(dat$POS))
        r5_2 = table(nchar(dat$POS),useNA = "a")
        r5_3 = table(substr(dat$POS,1,2),useNA = "a")
        
        ### ADDITIONAL
        
        # 6) is SNPID determined correctly?
        message("\t step 6..")
        chrpos = paste(dat$CHR,dat$POS,sep=":")
        r6_1 = table(dat$SNPID==chrpos,useNA = "a")  # **
        rm(chrpos)
        
        # 7) detailed view of marker info PER CHROMOSOME
        message("\t step 7..")
        # CAVEAT!  unique is not exactly unique, and duplicated is not exactly duplicated
        r7_1 = group_by(dat,CHR) %>% summarise(n=n(),
                                             unqSNPID=sum(!duplicated(SNPID)),
                                             dupSNPID=sum(duplicated(SNPID)),
                                             misSNPID=sum(is.na(SNPID)),
                                             unqRSID=sum(!duplicated(RSID)),
                                             dupRSID=sum(duplicated(RSID)),
                                             misRSID=sum(is.na(RSID)),
                                             minPoz=min(POS),
                                             maxPoz=max(POS),
                                             misPoz=sum(is.na(POS))) %>% ungroup()
        #tbl_4 = tbl_4[order(tbl_4$n,decreasing = T),]
        #barplot(tbl_4$minPoz,names.arg = tbl_4$CHR)
        #barplot(tbl_4$maxPoz,names.arg = tbl_4$CHR)
        
        ##### 8) duplicated SNPIDs
        message("\t step 8..")
        snps = unique(dat$SNPID[which(duplicated(dat$SNPID))]) # names of duplicated SNP
        tmp = dat[which(dat$SNPID %in% snps),]
        tmp = tmp[order(tmp$SNPID),]
        r8_1 = head(tmp)
        r8_2 = table(table(tmp$SNPID)) # number of times they were duplicated
        r8_3 = length(snps)
        rm(snps,tmp)

        ##### 9) duplicated RSIDs
        message("\t step 9..")
        snps = unique(dat$RSID[which(duplicated(dat$RSID))])
        tmp = dat[which(dat$RSID %in% snps),]
        tmp = tmp[order(tmp$RSID),]
        r9_1 = head(tmp)
        r9_2 = table(table(tmp$RSID)) # number of times they were duplicated
        r9_3 = length(snps)
        rm(snps,tmp)

        ##### 10) duplicated CHR:POS
        message("\t step 10..")
        chrpos = paste(dat$CHR,dat$POS,sep=":")
        snps = unique(chrpos[which(duplicated(chrpos))])
        chrpos = chrpos[which(chrpos %in% snps)]
        r10_2 = table(table(chrpos)) # number of times (replications) they were duplicated
        r10_3 = length(snps) # number of unique duplicated-names
        rm(chrpos,snps)
        
        ####################################################### more detailed
        dat$P_VAL_ADD = as.numeric(dat$P_VAL_ADD)
        dat$P_VAL_DOM = as.numeric(dat$P_VAL_DOM)
        dat$P_VAL_REC = as.numeric(dat$P_VAL_REC)
        dat$INFO = as.numeric(dat$INFO)

        short_file_name = paste(unlist(strsplit(file_name,"-"))[1:2],collapse="-")
        
        ##### 11) BETA
        message("\t step 11: BETA..")
        r11_1_hist = hist(dat$BETA_ADD,breaks=100,col="grey",xlab="BETA_ADD",main=short_file_name)
        r11_2_numr = sum(is.na(dat$BETA_ADD)) # **
        r11_3_qntl = quantile(dat$BETA_ADD,probs = seq(0.01,0.99,0.01),na.rm = T)
        
        ##### 12) SERR
        message("\t step 12: SERR..")
        r12_1_hist = hist(dat$SE_ADD,breaks=100,col="grey",xlab="SE_ADD",main=short_file_name)
        r12_2_numr = sum(is.na(dat$SE_ADD))
        r12_3_qntl = quantile(dat$SE_ADD,probs = seq(0.01,0.99,0.01),na.rm = T)
        
        ##### 13) PVAL_ADD
        message("\t step 13: P_VAL_ADD..")
        r13_1_hist = hist(dat$P_VAL_ADD,breaks=100,col="grey",xlab="P_VAL_ADD",main=short_file_name)
        r13_2_numr = sum(is.na(dat$P_VAL_ADD))
        r13_3_qntl = quantile(dat$P_VAL_ADD,probs = seq(0.01,0.99,0.01),na.rm = T)
        
        ##### 14) PVAL_DOM
        message("\t step 14: P_VAL_DOM..")
        r14_1_hist = hist(dat$P_VAL_DOM,breaks=100,col="grey",xlab="P_VAL_DOM",main=short_file_name)
        r14_2_numr = sum(is.na(dat$P_VAL_DOM))
        r14_3_qntl = quantile(dat$P_VAL_DOM,probs = seq(0.01,0.99,0.01),na.rm = T)
        
        ##### 15) PVAL_REC
        message("\t step 15: P_VAL_REC..")
        r15_1_hist = hist(dat$P_VAL_REC,breaks=100,col="grey",xlab="P_VAL_REC",main=short_file_name)
        r15_2_numr = sum(is.na(dat$P_VAL_REC))
        r15_3_qntl = quantile(dat$P_VAL_REC,probs = seq(0.01,0.99,0.01),na.rm = T)
        
        ##### 16) EAF_CASES
        message("\t step 16: EAF_CASES..")
        if ("EAF_CASES" %in% colnames(dat)) {
        if(any(!is.na(dat$EAF_CASES))) {
        r16_1_hist = hist(dat$EAF_CASES,breaks=100,col="grey",xlab="EAF_CASES",main=short_file_name)
        r16_3_qntl = quantile(dat$EAF_CASES,probs = seq(0.01,0.99,0.01),na.rm = T)
        } else {
                r16_1_hist = NULL
                r16_3_qntl = NULL
        }
        r16_2_numr = sum(is.na(dat$EAF_CASES))
        } else {
                r16_1_hist = NULL
                r16_3_qntl = NULL
                r16_2_numr = NULL
        }
        
        
        ##### 17) EAF_CONTR
        message("\t step 17: EAF_CONTR..")
        if ("EAF_CONTR" %in% colnames(dat)) {
        r17_1_hist = hist(dat$EAF_CONTR,breaks=100,col="grey",xlab="EAF_CONTR",main=short_file_name)
        r17_2_numr = sum(is.na(dat$EAF_CONTR))
        r17_3_qntl = quantile(dat$EAF_CONTR,probs = seq(0.01,0.99,0.01),na.rm = T)
        } else {
                r17_1_hist = NULL
                r17_2_numr = NULL
                r17_3_qntl = NULL
        }
        
        ##### 18) INFO
        message("\t step 18: INFO..")
        r18_1_hist = hist(dat$INFO,breaks=100,col="grey",xlab="INFO",main=short_file_name)
        r18_2_numr = sum(is.na(dat$INFO))
        r18_3_qntl = quantile(dat$INFO,probs = seq(0.01,0.99,0.01),na.rm = T)
        
        ##### 19) COMPOSITE MEASURE FOR P CONSISTENCY
        message("\t step 19..")
        betase = dat$BETA_ADD/dat$SE_ADD
        r19_1_plot = hexbin(x=betase,y=dat$P_VAL_ADD,xbins = 100,
                            xlab = "BETA/SE",ylab="P_VAL_ADD")
        rm(betase)
        r19_2_plot = hexbin(x=dat$P_VAL_ADD,y=dat$P_VAL_DOM,xbins = 100,
                            xlab = "P_VAL_ADD",ylab="P_VAL_DOM")
        r19_3_plot = hexbin(x=dat$P_VAL_ADD,y=dat$P_VAL_REC,xbins = 100,
                            xlab = "P_VAL_ADD",ylab="P_VAL_REC")
        r19_4_plot = hexbin(x=dat$P_VAL_DOM,y=dat$P_VAL_REC,xbins = 100,
                            xlab = "P_VAL_DOM",ylab="P_VAL_REC")
        
        
        ##### 20) COMPOSITE MEASURE FOR EFF-ALLELE-FREQ CONSISTENCY
        message("\t step 20..")
        if (all(c("EAF_CASES","EAF_CONTR") %in% colnames(dat))) {
        if(any(!is.na(dat$EAF_CASES))) { ##  not perfect solution! ***
        r20_1_plot = hexbin(x=dat$EAF_CASES,y=dat$EAF_CONTR,xbins = 100,
                            xlab = "EAF_CASES",ylab="EAF_CONTR")
        } else {
                r20_1_plot = NULL
        }
        } else {
                r20_1_plot = NULL
        }
        
        
        ####################################################### 
        
        ##### 21) STRAND
        message("\t step 21: STRAND..")
        r21_1_tbl = table(dat$STRAND,useNA = "a")
        
        ##### 22) BUILD
        message("\t step 22: BUILD..")
        if ("BUILD" %in% colnames(dat)) {
        r22_1_tbl = table(dat$BUILD,useNA = "a")
        } else { r22_1_tbl = NULL }
        
        ##### 23) N SAMPLES
        message("\t step 23: N SAMPLES..")
        if ("N" %in% colnames(dat)) {
        r23_1_tbl = table(dat$N,useNA = "a")
        } else { r23_1_tbl = NULL }
        
        ##### 24) IMPUTED
        message("\t step 24: IMPUTED..")
        if ("IMPUTED" %in% colnames(dat)) {
        r24_1_tbl = table(dat$IMPUTED,useNA = "a")
        } else { r24_1_tbl = NULL }
        
        ##### 25) IMPUTED-INFO COMPOSITE
        message("\t step 25..")
        if ("IMPUTED" %in% colnames(dat)) {
        infomiss = factor(is.na(dat$INFO)) # ,levels = c("TRUE","FALSE")
        r25_1_tbl = table(IMPUTED=dat$IMPUTED,INFOMISS=infomiss,useNA = "a")
        rm(infomiss)
        } else {
                r25_1_tbl = NULL
        }
        
        ##### 26) BETA-SE COMPOSITE
        message("\t step 26..")
        betamis = factor(is.na(dat$BETA_ADD),levels = c("TRUE","FALSE"))
        serrmis = factor(is.na(dat$SE_ADD),levels = c("TRUE","FALSE"))
        r26_1_tbl = table(BETAMISS=betamis,SERRMISS=serrmis)
        rm(serrmis)
        
        ##### 27) BETA-PADD COMPOSITE
        message("\t step 27..")
        #betamis = factor(is.na(dat$BETA_ADD),levels = c("TRUE","FALSE"))
        paddmis = factor(is.na(dat$P_VAL_ADD),levels = c("TRUE","FALSE"))
        r27_1_tbl = table(BETAMISS=betamis,PADDMISS=paddmis)
        rm(betamis)
        
        ##### 28) PVAL-REC_DOM COMPOSITE
        message("\t step 28..")
        precmis = factor(is.na(dat$P_VAL_REC),levels = c("TRUE","FALSE"))
        pdommis = factor(is.na(dat$P_VAL_DOM),levels = c("TRUE","FALSE"))
        r28_1_tbl = table(PRECMISS=precmis,PDOMMISS=pdommis)
        
        ##### 29) PVAL-ADD-REC COMPOSITE
        message("\t step 29..")
        #paddmis = factor(is.na(dat$P_VAL_ADD),levels = c("TRUE","FALSE"))
        #precmis = factor(is.na(dat$P_VAL_REC),levels = c("TRUE","FALSE"))
        r29_1_tbl = table(PADDMISS=paddmis,PRECMISS=precmis)
        
        ##### 30) PVAL-ADD-DOM COMPOSITE
        message("\t step 30..")
        #paddmis = factor(is.na(dat$P_VAL_ADD),levels = c("TRUE","FALSE"))
        #pdommis = factor(is.na(dat$P_VAL_DOM),levels = c("TRUE","FALSE"))
        r30_1_tbl = table(PADDMISS=paddmis,PDOMMISS=pdommis)
        rm(paddmis,precmis,pdommis)
        
        ####################################################  most signif
        
        ##### 31) extreme p-vlaues
        message("\t step 31..")
        cix = which(colnames(dat) %in% c("P_VAL_ADD","P_VAL_DOM","P_VAL_REC"))
        minP = apply(dat[,cix],1,function(x) ifelse(all(is.na(x)),NA,min(x,na.rm=T)))
        rix = which(minP<1e-6)
        if (length(rix)>0) {
        r31_1_dtf = dat[rix,]
        } else {
                r31_1_dtf = NULL
        }
        rm(cix,minP,rix)
        
        ####################################################  save
        message("\t saving..")
        diagn_rprt_file_name = paste(unlist(strsplit(file_name,"\\.txt\\.gz")),"_diagnRprt.RData",sep="")
        save(list=c("r1_1","r1_2",
                    "r2_1","r2_2","r2_3","r2_4",
                    "r3_1","r3_2","r3_3","r3_4","r3_5",
                    "r4_1","r4_2","r4_3",
                    "r5_1","r5_2","r5_3",
                    "r6_1",
                    "r7_1",
                    "r8_1","r8_2","r8_3",
                    "r9_1","r9_2","r9_3",
                    "r10_2","r10_3",
                    "r11_1_hist","r11_2_numr","r11_3_qntl",
                    "r12_1_hist","r12_2_numr","r12_3_qntl",
                    "r13_1_hist","r13_2_numr","r13_3_qntl",
                    "r14_1_hist","r14_2_numr","r14_3_qntl",
                    "r15_1_hist","r15_2_numr","r15_3_qntl",
                    "r16_1_hist","r16_2_numr","r16_3_qntl",
                    "r17_1_hist","r17_2_numr","r17_3_qntl",
                    "r18_1_hist","r18_2_numr","r18_3_qntl",
                    "r19_1_plot","r19_2_plot","r19_3_plot","r19_4_plot",
                    "r20_1_plot","r21_1_tbl","r22_1_tbl","r23_1_tbl",
                    "r24_1_tbl","r25_1_tbl","r26_1_tbl","r27_1_tbl",
                    "r28_1_tbl","r29_1_tbl","r30_1_tbl","r31_1_dtf"),
             file=paste(dia_dir,diagn_rprt_file_name,sep="")) # save diagnostic reports
}


