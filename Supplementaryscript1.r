#  Supplementaryscript1.r
############################
#         LMM model        #
############################
RTqPCRdata <- read.table("rawData.tab",header=T)

TARGET="target"
REFA = "MDH"
REFB = "UBQ"

t=15
j=1000
var_Target <- c()

for (i in seq(1,j,1)){
var_Target <- c(var_Target,var(sample(RTqPCRdata$ct[RTqPCRdata$gene==TARGET & !is.na(RTqPCRdata$ct) ],t,replace=T)))
}

var_RefA <- c()

for (i in seq(1,j,1)){
var_RefA <- c(var_RefA,var(sample(RTqPCRdata$ct[RTqPCRdata$gene==REFA & !is.na(RTqPCRdata$ct)],t,replace=T)))
}

var_RefB <- c()

for (i in seq(1,j,1)){
var_RefB <- c(var_RefB,var(sample(RTqPCRdata$ct[RTqPCRdata$gene==REFB & !is.na(RTqPCRdata$ct)],t,replace=T)))
}

pdf("variancePlots.pdf",h=(15.87/2.5),w=(9.05/2.5))
par(mfrow=c(3,1))
plot(density(var_Target,na.rm=T),lwd=3,main="Target Bootstrapped Variance",xlim=c(0,max(var_Target,na.rm=T)))
abline(v=var(RTqPCRdata$ct[RTqPCRdata$gene==TARGET],na.rm=T),col="red",lwd=2)
axis(3,at=round(var(RTqPCRdata$ct[RTqPCRdata$gene==TARGET],na.rm=T),3),col="red",lwd=2,labels=F)

mtext(3,at=round(var(RTqPCRdata$ct[RTqPCRdata$gene==TARGET],na.rm=T),3),text=round(var(RTqPCRdata$ct[RTqPCRdata$gene==TARGET],na.rm=T),3),col="red",line=.46,cex=.7)

plot(density(var_RefA),lwd=3,main="Reference A Bootstrapped Variance",xlim=c(0,max(var_Target,na.rm=T)))
abline(v=var(RTqPCRdata$ct[RTqPCRdata$gene==REFA],na.rm=T),col="red",lwd=2)
axis(3,at=round(var(RTqPCRdata$ct[RTqPCRdata$gene==REFA],na.rm=T),3),col="red",lwd=2,labels=F)

mtext(3,at=round(var(RTqPCRdata$ct[RTqPCRdata$gene==REFA],na.rm=T),3),text=round(var(RTqPCRdata$ct[RTqPCRdata$gene==REFA],na.rm=T),3),col="red",line=.46,cex=.7)

plot(density(var_RefB,na.rm=T),lwd=3,main="Reference B Bootstrapped Variance",xlim=c(0,max(var_Target,na.rm=T)))

abline(v=var(RTqPCRdata$ct[RTqPCRdata$gene==REFB],na.rm=T),col="red",lwd=2)
axis(3,at=round(var(RTqPCRdata$ct[RTqPCRdata$gene==REFB],na.rm=T),3),col="red",lwd=2,labels=F)

mtext(3,at=round(var(RTqPCRdata$ct[RTqPCRdata$gene==REFB],na.rm=T),3),text=round(var(RTqPCRdata$ct[RTqPCRdata$gene==REFB],na.rm=T),3),col="red",line=.46,cex=.7)


dev.off()

#install.packages("lme4")
#install.packages("multcomp")

library(lme4)
library(multcomp)

RTqPCRdata$interaction <- paste(RTqPCRdata$cultivar,RTqPCRdata$gene,RTqPCRdata$treat, sep=".")

mod <- lmer(ct ~ 0 + interaction + (1|sample) , data=RTqPCRdata )

summary(mod)

shapiro.test(resid(mod))

residuals <- resid(mod)[seq(1,length(resid(mod)),2)]

pdf("normalityTests.pdf",h=5,w=5)
par(mfrow=c(2,2))
hist(residuals,main="")
plot(residuals ~ predict(mod)[seq(1,length(resid(mod)),2)],xlab="model predictions",main="residuals vs. predicted")
plot(density(residuals),main="",)

qqnorm(residuals,main=paste("p=",round(as.numeric((shapiro.test(residuals)[2])),3)))

dev.off()

contr<- as.matrix(read.delim("contrast.table.tab",row.names=1))

resp.contr <- glht(mod, linfct = -contr)

summary(resp.contr)

confidence.intervals <- confint(resp.contr)$confint

#plot the Log2fold change

pdf("Log2FC.pdf", h=8,w=8)
barplot <- barplot(confidence.intervals[,1],beside=T, ylim=c(-5,5),las=1,yaxt='n',cex.names=1.1)


title(ylab=expression("Log"[2]*"FC"),line=2)
axis(2, at =c(seq(-15,15,1)),las =2,cex=.5)

abline(h=0)
abline(h=1,lty=3,lwd=2,col="red")
abline(h=-1,lty=3,lwd=2,col="red")
arrows(x0 = barplot, y0 = confidence.intervals[,2], x1 = barplot,y1=confidence.intervals[,3],code=3,angle=90,length=0.05,col="#964841",lwd=2)
dev.off()

FCCI <- 2^confidence.intervals

pdf("FC.pdf", h=8,w=8)
barplot <- barplot(FCCI[,1],beside=T, ylim=c(0,20),las=1,yaxt='n',cex.names=1)


title(ylab=expression("FC"),line=2)
axis(2, at =c(seq(0,10,1)),las =2,cex=.5)

abline(h=0)
abline(h=1,lty=3,lwd=2,col="red")
abline(h=.5,lty=3,lwd=2,col="blue")
arrows(x0 = barplot, y0 = FCCI[,2], x1 = barplot,y1=FCCI[,3],code=3,angle=90,length=0.05,col="#964841",lwd=2)
dev.off()

multicompTestResut <- summary(resp.contr)

capture.output(multicompTestResut, file = "multicompTestResut.txt")
