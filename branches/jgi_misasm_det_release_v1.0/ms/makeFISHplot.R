#!/usr/bin/env Rscript
library('fields')
xrad=0.45
yrad=0.9
nseg=360
eX=c(1,6.5)
eY=c(1.5,3)
c=c(3,3.33,3.67,4)
X=c(c)
Y=c(c+1)
par(mar=(c(4,4.5,4,2)+0.1))
pdf("FISHplot_XY.pdf")
plot(X,Y,xlim=c(0,8),ylim=c(0,7),pch=16, xlab="Contig X", ylab="Contig Y",xaxt='n',cex.lab=2.0,cex=1.5) 
axis(1,at=c(0:8),labels=c(0:8),cex=1.5)
xline(4.5,lty=2)
yline(3.5,lty=2)
x.cent=mean(X)
y.cent=mean(Y)
xx <- x.cent + xrad*cos( seq(0,2*pi, length.out=nseg) )
yy <- y.cent + yrad*sin( seq(0,2*pi, length.out=nseg) )
angle=7*pi/4
xr = (xx-x.cent)*cos(angle)-(yy-y.cent)*sin(angle)
yr = (xx-x.cent)*sin(angle)+(yy-y.cent)*cos(angle)
lines(xr+x.cent,yr+y.cent,col='darkgray',cex=1.5)
points(eX,eY,col='darkorange',pch=16,cex=1.5)
legend(5.75,7,c("Real","Erroneous"),pch=16,pt.cex=1.5,col=c("black","darkorange"))

dev.off();
pdf("FISHplot_XX.pdf")
X1=c(5,5.33,5.67,6)
X2=c(X1+2)
plot(X1,X2,xlim=c(0,8),ylim=c(0,8),pch=16,xlab="Contig X", ylab="Contig X",xaxt='n',cex.lab=2.0,cex=1.5,yaxt='n') 
axis(1,at=c(0:8),labels=c(0:8),cex=1.5)
axis(2,at=c(0:8),labels=c(0:8),cex=1.5)
xline(4.5,lty=2)
yline(4.5,lty=2)
x.cent=mean(X1)
y.cent=mean(X2)
xx <- x.cent + xrad*cos( seq(0,2*pi, length.out=nseg) )
yy <- y.cent + yrad*sin( seq(0,2*pi, length.out=nseg) )
angle=7*pi/4
xr = (xx-x.cent)*cos(angle)-(yy-y.cent)*sin(angle)
yr = (xx-x.cent)*sin(angle)+(yy-y.cent)*cos(angle)
lines(xr+x.cent,yr+y.cent,col='darkgray',cex=1.5)

dev.off();
pdf("FISHplot_YY.pdf")
Y1=c(0,0.33,0.67,1)
Y2=c(Y1+2)
plot(Y1,Y2,xlim=c(0,7),ylim=c(0,7),pch=16,xlab="Contig Y", ylab="Contig Y",xaxt='n',cex.lab=2.0,cex=1.5) 
axis(1,at=c(0:7),labels=c(0:7),cex=1.5)
xline(3.5,lty=2)
yline(3.5,lty=2)
x.cent=mean(Y1)
y.cent=mean(Y2)
xx <- x.cent + xrad*cos( seq(0,2*pi, length.out=nseg) )
yy <- y.cent + yrad*sin( seq(0,2*pi, length.out=nseg) )
angle=7*pi/4
xr = (xx-x.cent)*cos(angle)-(yy-y.cent)*sin(angle)
yr = (xx-x.cent)*sin(angle)+(yy-y.cent)*cos(angle)
lines(xr+x.cent,yr+y.cent,col='darkgray',cex=1.5)


dev.off()

