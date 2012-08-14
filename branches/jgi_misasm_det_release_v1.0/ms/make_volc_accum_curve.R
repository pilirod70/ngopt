#!/usr/bin/env Rscript

aco = sort(read.table("dat/Volc/A5/dsm3757.a5.contigs.len")[,1],decreasing=TRUE)
#acb = sort(read.table("dat/Volc/A5/dsm3757.a5.contigs.broken.len")[,1],decreasing=TRUE)
#acr.ec = sort(read.table("dat/Volc/A5/dsm3757.a5.ec.crude.scaffolds.len")[,1],decreasing=TRUE)
#afi.ec = sort(read.table("dat/Volc/A5/dsm3757.a5.ec.final.scaffolds.len")[,1],decreasing=TRUE)
acr = sort(read.table("dat/Volc/A5/dsm3757.a5.crude.scaffolds.len")[,1],decreasing=TRUE)
afi = sort(read.table("dat/Volc/A5/dsm3757.a5.final.scaffolds.len")[,1],decreasing=TRUE)
sco = sort(read.table("dat/Volc/SOAP/volc_soap.contig.len")[,1],decreasing=TRUE)
ssc = sort(read.table("dat/Volc/SOAP/volc_soap.scafSeq.len")[,1],decreasing=TRUE)

file = "volc_accum_plot.pdf"
pdf(file,width=14)
paste(file)

col=rainbow(5)
col[2] = "darkorange"
tmp = col[2]
col[2] = col[4]
col[4] = tmp
plot(cumsum(aco),type='l',col=col[1],xlab="N largest sequences ", ylab="No. of bases (Mb)",yaxt='n')
#plot(cumsum(ssc),type='l',col=col[5],xlab="N largest sequences ", ylab="No. of bases (Mb)",yaxt='n')
lbl=c(0,1,2,3,4)
at=lbl*1000000
axis(2,at=at,labels=lbl)
points(cumsum(acr),type='l',col=col[2])
points(cumsum(afi),type='l',col=col[3])
points(cumsum(sco),type='l',col=col[4])
points(cumsum(ssc),type='l',col=col[5])
#points(cumsum(aco),type='l',col=col[1])

names=c("A5 - ctg","A5 - scaf","A5 - scaf-QC","SOAP - ctg","SOAP - scaf")

legend("bottomright", legend=names,col=col,lty=1)


dev.off();

