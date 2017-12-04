
### script used to load and digest GWAS catalog
### 2017 Dec 04. Jonas B.
### located in: setwd("~/Dropbox/GIT/GWAS_CATALOG")

### basic idea: 
## 1) only relevant data is left
## 2) multiple SNPs (per input line) get expanded to multiple lines (output)
## 3) same is done for CHR_ID and CHR_POS
## 4) SNPxSNP interactions are eliminated
## 5) chrs and pozs are converted to numeric values


### load GWAS catalog
load("~/Dropbox/GIT/GWAS_CATALOG/gwas_catalog_v1.0-associations_e90_r2017-11-27.Rdata")

# leave only most valuable columns
use_col_names = c("PUBMEDID","DATE","DISEASE.TRAIT",
                  "REGION","CHR_ID","CHR_POS",
                  "REPORTED.GENE.S.","MAPPED_GENE",
                  "STRONGEST.SNP.RISK.ALLELE","SNPS",
                  "P.VALUE","OR.or.BETA")
a = a[,use_col_names] 

# get rid of GxE interaction studies (1)
bad_rows = grep("x",a$SNPS)
length(bad_rows)
a = a[-bad_rows,]; rm(bad_rows)

# get rid of GxE interaction studies (2)
bad_rows = grep("SNP x SNP interaction",a$DISEASE.TRAIT)
length(bad_rows)
a = a[-bad_rows,]; rm(bad_rows)

## recode CHR_ID
a$CHR_ID[which(a$CHR_ID=="")] = NA
table(a$CHR_ID)

## recode CHR_POS
a$CHR_POS[which(a$CHR_POS=="")] = NA
table(substr(a$CHR_POS,1,1))


#### reformat in case multiple SNPs were reported
table(nchar(a$SNPS))  # there are multiple SNPs
rix = grep(";|,",a$SNPS) # two scenarios of how SNPs are separated
tmp = a[rix,]  # save for digestion
a = a[-rix,]   # remove and leave for further merge with digested bit


ins = NULL  # cummulative insertion that will be added to "a" dataframe later
prb = NULL  # report of problematic instances
for (i in 1:nrow(tmp)) {
        snps = unlist(strsplit(tmp$SNPS[i],", |; "))
        snps = snps[which(substr(snps,1,2)=="rs")] # get rid of HLA_DRB1_0901 and similar ! ***
        pozs = unlist(strsplit(tmp$CHR_POS[i],",|;|, |; "))
        chrs = unlist(strsplit(tmp$CHR_ID[i],",|;|, |; "))
        
        if( (length(snps)!=length(pozs)) | (length(snps)!=length(chrs))) {
                pozs = rep(NA,length(snps)) # too much unclear
                chrs = rep(chrs[1],length(snps)) # chr names are always(?) unique per row of input
                prb = rbind(prb,tmp[i,c("PUBMEDID","DISEASE.TRAIT","CHR_ID","CHR_POS","SNPS")])
        }
        xtr = data.frame(tmp[i,],dum=rep(NA,length(snps)),stringsAsFactors = F)
        xtr$SNPS = snps
        xtr$CHR_POS = pozs
        xtr$CHR_ID = chrs
        ins = rbind(ins,xtr)
        rm(snps,chrs,pozs,xtr)
}  # ignore warnings

# preview encountered problems:
print(prb)

# doublecheck SNP ID types
table(substr(ins$SNPS,1,2)) # all "rs"
table(substr(a$SNPS,1,2))   # mostly "rs"

# remove dummy column
ins = ins[,-ncol(ins)]

# merge with the the rest
a = rbind(a,ins); rm(ins)

# SNP names final fixes
a$SNPS[which(a$SNPS=="NR")] = NA
table(nchar(a$SNPS),useNA = "a")
a[which(nchar(a$SNPS)>20),]

# CHR names final fixes
table(a$CHR_ID,useNA = "a")
a$CHR_ID[which(a$CHR_ID=="X")] = 23
a$CHR_ID = as.numeric(a$CHR_ID)

# CHR positions final fixes
table(substr(a$CHR_POS,1,1),useNA = "a")
table(nchar(a$CHR_POS),useNA = "a")
a$CHR_POS = as.numeric(a$CHR_POS)




