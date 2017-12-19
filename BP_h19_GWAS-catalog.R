load("~/Desktop/Work/GWAS_catalog/gwas_catalog_v1.0-associations_e90_r2017-11-27.Rdata")

table(substr(a$SNPS,1,2))
a$SNPS[which(substr(a$SNPS,1,2)!="rs")]

a= a[grep("rs",a$SNPS),]

d= unique((unlist(a$SNPS)))

#========================================================================================================
# Run from here

write.table(a,"~/Desktop/Work/GWAS_catalog/snp_list.txt", row.names=F, quote=F) # Use this to retreive bp 
# from UCSC GRCh37/hg19 

s= read.table("~/Desktop/Work/GWAS_catalog/SNP_coordinates.txt",h=T)
s$name= as.character(s$name)
s$chrom= as.character(s$chrom)
s= s[(grep("_", s$chrom,invert=T)),]
s= s[!duplicated(s$name),]
s= s[,c(1,3,4)]
colnames(s)= c("chr","bp_hg19","snp")

s$chr=gsub('chr', '', s$chr)

write.table(s,"~/Desktop/Work/GWAS_catalog/snplist_bp_hg19.txt", row.names=F, quote=F, sep="\t")

