#  Supplementaryscript2.r
############################
#Traditional method-DeltaCT#
############################
#first calculate the mean between both technical replicates

deltaCtData <- RTqPCRdata[seq(1,nrow(RTqPCRdata),2),]

deltaCtData$ct <- colMeans(matrix(RTqPCRdata$ct,nrow=2),na.rm=T)

#Opt - Wat

target.opt <-  deltaCtData$ct[(which(deltaCtData$treat == "opt" & deltaCtData$gene == TARGET))]

refa.opt <-  deltaCtData$ct[(which(deltaCtData$treat == "opt" & deltaCtData$gene == REFA))]

refb.opt <-  deltaCtData$ct[(which(deltaCtData$treat == "opt" & deltaCtData$gene == REFB))]

target.wat <- deltaCtData$ct[(which(deltaCtData$treat == "wat" & deltaCtData$gene == TARGET))]

refa.wat <- deltaCtData$ct[(which(deltaCtData$treat == "wat" & deltaCtData$gene == REFA))]

refb.wat <- deltaCtData$ct[(which(deltaCtData$treat == "wat" & deltaCtData$gene == REFB))]

opt <- target.opt - 0.5*(refa.opt + refb.opt)

wat <- target.wat - 0.5*(refa.wat + refb.wat)

diff.opt.wat <- opt - wat

t.test.opt.wat <- t.test(opt, wat)

#aca - cat

target.aca <-  deltaCtData$ct[(which(deltaCtData$cultivar == "aca" & deltaCtData$gene == TARGET))]

refa.aca <-  deltaCtData$ct[(which(deltaCtData$cultivar == "aca" & deltaCtData$gene == REFA))]

refb.aca <-  deltaCtData$ct[(which(deltaCtData$cultivar == "aca" & deltaCtData$gene == REFB))]

target.cat <- deltaCtData$ct[(which(deltaCtData$cultivar == "cat" & deltaCtData$gene == TARGET))]

refa.cat <- deltaCtData$ct[(which(deltaCtData$cultivar == "cat" & deltaCtData$gene == REFA))]

refb.cat <- deltaCtData$ct[(which(deltaCtData$cultivar == "cat" & deltaCtData$gene == REFB))]

aca <- target.aca - 0.5*(refa.aca + refb.aca)

cat <- target.cat -0.5*(refa.cat + refb.cat)

diff.aca.cat <- aca - cat

t.test.aca.cat <- t.test(aca, cat)

#aca.opt - aca.wat

target.aca.opt <-  deltaCtData$ct[(which(deltaCtData$cultivar == "aca" & deltaCtData$gene == TARGET & deltaCtData$treat == "opt"))]

refa.aca.opt <-  deltaCtData$ct[(which(deltaCtData$cultivar == "aca" & deltaCtData$gene == REFA & deltaCtData$treat == "opt"))]

refb.aca.opt <-  deltaCtData$ct[(which(deltaCtData$cultivar == "aca" & deltaCtData$gene == REFB & deltaCtData$treat == "opt"))]

target.aca.wat <- deltaCtData$ct[(which(deltaCtData$cultivar == "aca" & deltaCtData$gene == TARGET & deltaCtData$treat == "wat"))]

refa.aca.wat <- deltaCtData$ct[(which(deltaCtData$cultivar == "aca" & deltaCtData$gene == REFA & deltaCtData$treat == "wat"))]

refb.aca.wat <- deltaCtData$ct[(which(deltaCtData$cultivar == "aca" & deltaCtData$gene == REFB & deltaCtData$treat == "wat"))]

aca.opt <- target.aca.opt - 0.5*(refa.aca.opt + refb.aca.opt)

aca.wat <- target.aca.wat -0.5*(refa.aca.wat + refb.aca.wat)

diff.aca.opt.wat <- aca.opt - aca.wat

t.test.aca.opt.wat <- t.test(aca.opt, aca.wat)

#cat.opt - cat.wat

target.cat.opt <-  deltaCtData$ct[(which(deltaCtData$cultivar == "cat" & deltaCtData$gene == TARGET & deltaCtData$treat == "opt"))]

refa.cat.opt <-  deltaCtData$ct[(which(deltaCtData$cultivar == "cat" & deltaCtData$gene == REFA & deltaCtData$treat == "opt"))]

refb.cat.opt <-  deltaCtData$ct[(which(deltaCtData$cultivar == "cat" & deltaCtData$gene == REFB & deltaCtData$treat == "opt"))]

target.cat.wat <- deltaCtData$ct[(which(deltaCtData$cultivar == "cat" & deltaCtData$gene == TARGET & deltaCtData$treat == "wat"))]

refa.cat.wat <- deltaCtData$ct[(which(deltaCtData$cultivar == "cat" & deltaCtData$gene == REFA & deltaCtData$treat == "wat"))]

refb.cat.wat <- deltaCtData$ct[(which(deltaCtData$cultivar == "cat" & deltaCtData$gene == REFB & deltaCtData$treat == "wat"))]

cat.opt <- target.cat.opt - 0.5*(refa.cat.opt + refb.cat.opt)

cat.wat <- target.cat.wat -0.5*(refa.cat.wat + refb.cat.wat)

diff.cat.opt.wat <- cat.opt - cat.wat

t.test.cat.opt.wat <- t.test(cat.opt, cat.wat)

#create a vector with the mean of the -diff (or delta delta Ct)

diff.expression.means <- c(mean(-diff.opt.wat),
mean(-diff.aca.cat),
mean(-diff.aca.opt.wat),
mean(-diff.cat.opt.wat))

names(diff.expression.means) <- c("opt - wat", "aca - cat", "aca.opt - aca.wat", "cat.opt - cat.wat ")

lower.conf.int <- c(-t.test.opt.wat$conf.int[1],
-t.test.aca.cat$conf.int[1],
-t.test.aca.opt.wat$conf.int[1],
-t.test.cat.opt.wat$conf.int[1])

upper.conf.int <- c(-t.test.opt.wat$conf.int[2],
-t.test.aca.cat$conf.int[2],
-t.test.aca.opt.wat$conf.int[2],
-t.test.cat.opt.wat$conf.int[2])

pdf("deltaDeltaCtLog2FC.pdf", h=8,w=8)
barplot <- barplot(diff.expression.means,beside=T, ylim=c(-5,5),las=1,yaxt='n',cex.names=1.1)


title(ylab=expression("Log"[2]*"FC"),line=2)
axis(2, at =c(seq(-15,15,1)),las =2,cex=.5)

abline(h=0)
abline(h=1,lty=3,lwd=2,col="red")
abline(h=-1,lty=3,lwd=2,col="red")
arrows(x0 = barplot, y0 = lower.conf.int, x1 = barplot,y1=upper.conf.int,code=3,angle=90,length=0.05,col="#964841",lwd=2)
dev.off()

FC.diff.expression.means <- 2^diff.expression.means

FC.lower.conf.int <- 2^upper.conf.int

FC.upper.conf.int <- 2^lower.conf.int

pdf("deltaDeltaCtFC.pdf", h=8,w=8)
barplot <- barplot(FC.diff.expression.means,beside=T, ylim=c(0,5),las=1,yaxt='n',cex.names=1)

title(ylab=expression("FC"),line=2)
axis(2, at =c(seq(0,10,1)),las =2,cex=.5)

abline(h=0)
abline(h=1,lty=3,lwd=2,col="red")
abline(h=.5,lty=3,lwd=2,col="blue")
arrows(x0 = barplot, y0 = FC.lower.conf.int, x1 = barplot,y1=FC.upper.conf.int,code=3,angle=90,length=0.05,col="#964841",lwd=2)
dev.off()


############################################
#make a single plot with Delta Ct and model#
############################################

names <- c("opt - wat", "aca - cat", "aca.opt - aca.wat", "cat.opt - cat.wat ")

pdf("LMM.versus.DeltaCt.pdf",h=9,w=9)
par(mfrow = c(2,2))
#log2FC lMM model
barplot <- barplot(confidence.intervals[,1],beside=T, ylim=c(

if(min(confidence.intervals[,2])>=-5){-5}else{min(confidence.intervals[,2])*1.2},
if(max(confidence.intervals[,3])<=5){5}else{max(confidence.intervals[,3])*1.2}

),las=1,yaxt='n',names="",main="LMM model",cex.main=2)

text( x= barplot,y=if(min(confidence.intervals[,2])>=-5){-6.5}else{min(confidence.intervals[,2])*1.4},
labels = names,srt=30,xpd=T)

axis(2, at =c(seq(-100,100,1)),las =2,cex=.5)

title(ylab=expression("Log"[2]*"FC"),line=2,cex.lab=1.5)

abline(h=0)
abline(h=1,lty=3,lwd=2,col="red")
abline(h=-1,lty=3,lwd=2,col="red")
arrows(x0 = barplot, y0 = confidence.intervals[,2], x1 = barplot,y1=confidence.intervals[,3],code=3,angle=90,length=0.05,col="#964841",lwd=2)

#FC lMM model

barplot <- barplot(FCCI[,1],beside=T, ylim=c(0,
if(max(FCCI[,3]) <= 5){5}else
{max(FCCI[,3])*1.20
}),
las=1,yaxt='n',cex.names=1,main="LMM model",names="",cex.main=2)

title(ylab=expression("FC"),line=2,cex.lab=1.5)

axis(2, at =c(seq(0,100)),las =2,cex=.5)

text( x= barplot,y=if(max(FCCI[,3]) <= 5 ){-.7}else if(max(FCCI[,3]) <= 15){- 1.4}else{-3.5},
labels = names,srt=30,xpd=T)

abline(h=0)
abline(h=1,lty=3,lwd=2,col="red")
abline(h=.5,lty=3,lwd=2,col="blue")
arrows(x0 = barplot, y0 = FCCI[,2], x1 = barplot,y1=FCCI[,3],code=3,angle=90,length=0.05,col="#964841",lwd=2)

FC.strapolator=FALSE
# Delta Delta Ct Log2FC
barplot <- barplot(diff.expression.means,beside=T, ylim=c(
if(min(upper.conf.int)>=-5){-5}else{min(upper.conf.int)*1.2},
if(max(lower.conf.int)<=5){5}else{max(upper.conf.int)*1.2}
),las=1,yaxt='n',cex.names=1.1,main=expression(paste(Delta,bold("Ct method"))),names="",cex.main=2)

text( x= barplot,y=if(min(upper.conf.int)>=-5){-6.5}else{min(upper.conf.int)*1.4},
labels = names,srt=30,xpd=T)

title(ylab=expression("Log"[2]*"FC"),line=2,cex.lab=1.5)

axis(2, at =c(seq(-100,100,1)),las =2,cex=.5)

abline(h=0)
abline(h=1,lty=3,lwd=2,col="red")
abline(h=-1,lty=3,lwd=2,col="red")
arrows(x0 = barplot, y0 = lower.conf.int, x1 = barplot,y1=upper.conf.int,code=3,angle=90,length=0.05,col="#964841",lwd=2)

#Ct methdo FC

barplot <- barplot(FC.diff.expression.means,beside=T, ylim=c(0,
if(max(FC.upper.conf.int) <= 5){5}else{max(FC.upper.conf.int)*1.2})
,las=1,yaxt='n',cex.names=1,main=expression(paste(Delta,bold("Ct method"))),names="",cex.main=2)

text( x= barplot,y=if(max(FC.upper.conf.int) <= 5 ){-.7}else if(max(FC.upper.conf.int) <= 15){- 1.4}else{-3.5},
labels = names,srt=30,xpd=T)

title(ylab=expression("FC"),line=2,cex.lab=1.5)
axis(2, at =c(seq(0,100,1)),las =2,cex=.5)

abline(h=0)
abline(h=1,lty=3,lwd=2,col="red")
abline(h=.5,lty=3,lwd=2,col="blue")
arrows(x0 = barplot, y0 = FC.lower.conf.int, x1 = barplot,y1=FC.upper.conf.int,code=3,angle=90,length=0.05,col="#964841",lwd=2)
dev.off()


