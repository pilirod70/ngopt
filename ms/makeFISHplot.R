#!/usr/bin/env Rscript
library('fields')
nseg=360
eX=c(1,6.5)
eY=c(1.5,3)
c=c(3,3.33,3.67,4)
X=c(c)
Y=c(c+1)
par(mar=(c(4,4.5,4,2)+0.1))
pdf("FISHplot.png")
plot(X,Y,xlim=c(0,8),ylim=c(0,7),pch=16, xlab="Contig X", ylab="Contig Y",xaxt='n',cex.lab=2.0,cex=1.5) 
axis(1,at=c(0:8),labels=c(0:8),cex=1.5)
xline(4.5,lty=2)
yline(3.5,lty=2)
x.cent=3.5
y.cent=4.5
r=2
xx <- x.cent + 0.5*cos( seq(0,2*pi, length.out=nseg) )
yy <- y.cent + 1*sin( seq(0,2*pi, length.out=nseg) )
angle=7*pi/4
xr = (xx-x.cent)*cos(angle)-(yy-y.cent)*sin(angle)
yr = (xx-x.cent)*sin(angle)+(yy-y.cent)*cos(angle)
lines(xr+x.cent,yr+y.cent,col='darkgray',cex=1.5)
points(eX,eY,col='darkorange',pch=16,cex=1.5)
legend(5.75,7,c("Real","Erroneous"),pch=16,pt.cex=1.5,col=c("black","darkorange"))
dev.off()
