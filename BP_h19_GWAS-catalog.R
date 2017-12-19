# update GWAS Catalog SNPs by assigning them hg19 positions from UCSC
# Pol & Jonas. 2017 Dec 15-20

#load("~/Desktop/Work/GWAS_catalog/gwas_catalog_v1.0-associations_e90_r2017-11-27.Rdata") # @Pol
load("~/Dropbox/GIT/GWAS_CATALOG/gwas_catalog_v1.0-associations_e90_r2017-11-27.Rdata") # @Jon

# preview SNP names
table(substr(a$SNPS,1,2))

# extract all lines with informative SNP names
a = a[grep("rs",a$SNPS),]

# extract SNP names
snps = unique((unlist(a$SNPS)))

# doublecheck whether SNP names look OK
table(nchar(snps)) # -> there are many written in a string (with multiple SNPs)
snps[which(nchar(snps)>15)] # non-standard ones

# examples of common problems:
"rs6858430 x rs4800250"
"rs8453; rs926938" # snps[grep(";",snps)]
"HLA-DQB1*02:01, rs558702" # snps[grep(",",snps)]

### separate into batches:

# pure good-old-fashioned rs SNPs
ix = grep("^rs[0-9]+$",snps)
gr1 = snps[ix]
snps = snps[-ix]

# those separated by ";"
ix = grep("^rs[0-9]+;.*rs[0-9]+$",snps)
gr2 = snps[ix]
snps = snps[-ix]
gr2 = unique(unlist(lapply(gr2,function(x) unlist(strsplit(x,"; ")))))

# those separated by "x"
ix = grep("^rs[0-9]+ x *rs[0-9]+$",snps)
gr3 = snps[ix]
snps = snps[-ix]
gr3 = unique(unlist(lapply(gr3,function(x) unlist(strsplit(x," x ")))))

# those separated by ","
ix = grep(",",snps)
gr4 = snps[ix]
snps = snps[-ix]
gr4 = unique(unlist(lapply(gr4,function(x) unlist(strsplit(x,", ")))))

# remaining ones separated by ";"
ix = grep(";",snps)
gr5 = snps[ix]
snps = snps[-ix]
gr5 = unique(unlist(lapply(gr5,function(x) unlist(strsplit(x,"; ")))))

# there is one ambiguous/nonstandard SNP name remaining
snps  # we can safely get rid of this one..

# pre-pre-final list of SNP names
snps = unique(c(gr1,gr2,gr3,gr4,gr5))

# number of still problematic instances
length(snps) - length(grep("rs[0-9]+$",snps))

# pre-final list of SNP names
ix = grep("rs[0-9]+$",snps)
snps = snps[ix]

# still problematic instances
table(nchar(snps)) # "rs8"? wtf..
snps = snps[which(nchar(snps)>3)]

# final check
table(nchar(snps)) # looks good!

# final-final
snps = unique(snps)

# export to UCSC browser
#write.table(snps,"~/Desktop/Work/GWAS_catalog/temp_GWAScatalog_snpList.txt", row.names=F, quote=F) # @Pol
write.table(snps,"~/Downloads/temp_GWAScatalog_snpList.txt", row.names=F, quote=F) # @Jon


#=====================================================
...
... upload to UCSC and extract GRCh37/hg19 coordinates
...
#=====================================================

# load the export from UCSC
#s = read.table("~/Desktop/Work/GWAS_catalog/temp_GWAScatalog_snpPositions_hg19.txt",h=T,sep="\t",stringsAsFactors = F) # @Pol
s = read.table("~/Downloads/temp_GWAScatalog_snpPositions_hg19.txt",h=T,stringsAsFactors = F,sep="\t") # @Jon

# eliminate nonsense
s = s[which(nchar(s$chrom) %in% c(4,5)),] # no haplo blocks (non-chromosomes)

# recode
s$chrom = gsub('chr', '', s$chrom)
s$chrom[which(s$chrom %in% c("x","X"))] = 23
s$chrom = as.numeric(s$chrom)

# reformat
s = s[,c("chrom","name","chromEnd")]

# save
#save(list="s", file="~/Desktop/Work/GWAS_catalog/e90_r2017-11-27_hg19positions.Rdata") # @Pol
save(list="s", file="~/Dropbox/GIT/GWAS_CATALOG/e90_r2017-11-27_hg19positions.Rdata") # @Jon

