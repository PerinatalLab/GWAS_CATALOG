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
length(snps) - length(grep("rs[0-9]+$",snps)) # HLA region SNPs

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
# (due to the UCSC retardedness, must be done in batches)
write.table(snps[1:10000],"~/Downloads/temp_GWAScatalog_snpList_1.txt", row.names=F,col.names = F, quote=F) # @Jon
write.table(snps[10001:20000],"~/Downloads/temp_GWAScatalog_snpList_2.txt", row.names=F,col.names = F, quote=F) # @Jon
write.table(snps[20001:30000],"~/Downloads/temp_GWAScatalog_snpList_3.txt", row.names=F,col.names = F, quote=F) # @Jon
write.table(snps[30001:40000],"~/Downloads/temp_GWAScatalog_snpList_4.txt", row.names=F,col.names = F, quote=F) # @Jon
write.table(snps[40001:length(snps)],"~/Downloads/temp_GWAScatalog_snpList_5.txt", row.names=F,col.names = F, quote=F) # @Jon

#=====================================================
...
... upload to 5 batches to UCSC and extract GRCh37/hg19 coordinates
...

# batch1. 22 of the 10000 given identifiers have no match in table snp150, field name.
# batch2. 66 of the 10000 given identifiers have no match in table snp150, field name. 
# batch3. 123 of the 10000 given identifiers have no match in table snp150, field name.
# batch4. 228 of the 10000 given identifiers have no match in table snp150, field name. 
# batch5. 5 of the 867 given identifiers have no match in table snp150, field name. 

... unhash headers
#=====================================================

# load the export from UCSC
s1 = read.table("~/Downloads/part1_temp.txt",h=T,stringsAsFactors = F,sep="\t") # @Jon
s2 = read.table("~/Downloads/part2_temp.txt",h=T,stringsAsFactors = F,sep="\t") # @Jon
s3 = read.table("~/Downloads/part3_temp.txt",h=T,stringsAsFactors = F,sep="\t") # @Jon
s4 = read.table("~/Downloads/part4_temp.txt",h=T,stringsAsFactors = F,sep="\t") # @Jon
s5 = read.table("~/Downloads/part5_temp.txt",h=T,stringsAsFactors = F,sep="\t") # @Jon

# combine
s = rbind(s1,s2,s3,s4,s5)

# inspect
table(nchar(s$chrom))
s$chrom[nchar(s$chrom)==13][1:10]
s$chrom[nchar(s$chrom)==14][1:10]
s$chrom[nchar(s$chrom)==15][1:10]

# leave only normal chromosome names
s = s[which(nchar(s$chrom) %in% c(4,5)),] # no haplo blocks (non-chromosomes)

# preview
table(s$chrom)

# recode chromosome names
s$chrom = gsub('chr', '', s$chrom)
s$chrom[which(s$chrom %in% c("x","X"))] = 23
s$chrom[which(s$chrom %in% c("y","Y"))] = 24
s$chrom = as.numeric(s$chrom)

# reformat
s = s[order(s$chrom,s$chromEnd),] # thus Y (least important) will be in the end
head(s)

# check for duplications
table(table(s$name))
s[which(s$name %in% s$name[which(duplicated(s$name))]),]  # pseudoautosomal? rs306896

# remove duplications
s = s[which(!duplicated(s$name)),] # Y chr will be removed, OK for me

# meta
nrow(s)

### final checks

# inspect SNP names
table(table(s$name)) # ok
table(nchar(s$name)) # ok

# inspect certainty of position
dif = s$chromEnd-s$chromStart
table(dif) # non =1 values might have problems in merging with HRC imputed data

# inspect: same position, different SNP names
pos = paste(s$chrom,s$chromEnd,sep="_")
table(table(pos))
s[which(pos %in% pos[duplicated(pos)]),] # might need actual allele names to resolve this


# save
save(list="s", file="~/Dropbox/GIT/GWAS_CATALOG/e90_r2017-11-27_hg19positions.Rdata") # @Jon


# summary:
# 444 SNP names were not found (out of 40867)
# total number of unique remaining SNP names: 40422
# Jonas B. 2017 Dec 19