#!/usr/bin/env Rscript

aco = sort(read.table("dat/Volc/A5/dsm3757.a5.contigs.len")[,1],decreasing=TRUE)
acb = sort(read.table("dat/Volc/A5/dsm3757.a5.contigs.broken.len")[,1],decreasing=TRUE)
acr = sort(read.table("dat/Volc/A5/dsm3757.a5.crude.scaffolds.len")[,1],decreasing=TRUE)
afi = sort(read.table("dat/Volc/A5/dsm3757.a5.final.scaffolds.len")[,1],decreasing=TRUE)
sco = sort(read.table("dat/Volc/SOAP/volc_soap.contig.len")[,1],decreasing=TRUE)
ssc = sort(read.table("dat/Volc/SOAP/volc_soap.scafSeq.len")[,1],decreasing=TRUE)

pdf("volc_accum_plot.pdf")
paste("volc_accum_plot.pdf")

col=rainbow(6)
tmp = col[2]
col[2] = col[3]
col[3] = tmp
#main="Sequence length accumulation",
plot(cumsum(aco),type='l',col=col[1],xlab="No. sequences ", ylab="No. of bases (Mb)",yaxt='n')
lbl=c(0,1,2,3,4)
at=lbl*1000000
axis(2,at=at,labels=lbl)
points(cumsum(acb),type='l',col=col[2])
points(cumsum(acr),type='l',col=col[3])
points(cumsum(afi),type='l',col=col[4])
points(cumsum(sco),type='l',col=col[5])
points(cumsum(ssc),type='l',col=col[6])

names=c("vA5ctg","vA5ctgQC","vA5-noQC","vA5","vSOAPctg","vSOAP")

legend("bottomright", legend=names,col=col,lty=1)


dev.off();

