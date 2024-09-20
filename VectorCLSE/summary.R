
load("~/Desktop/Changepoint/aop/aop_singleCP/CLS/notes/simulation/MultiInnov_rho0.5_iii1.RData")

ciT = 0.8615
thetaT = 0.43075
rhoT = 0.5
parT = c(ciT, thetaT, rhoT)

EstTable = matrix(NA, nrow=4, ncol=length(parT))
colnames(EstTable) = expression(c[2], beta[0], rho)
row.names(EstTable) = c("True", "Mean", "Rel.bias", "sd")

EstTable[1,] <- round(parT, digits=4)
EstTable[2,] <- round(colMeans(res[,7:9]), digits=4)
EstTable[3,] <- round((EstTable[1,] - EstTable[2,])/EstTable[1,], digits=4)
EstTable[4,] <- round(apply(res[,7:9], 2, sd), digits=4)
EstTable



tmp <- apply(res[,7:9], 2, range)
tmp

setEPS()
postscript("notes/simulation/FigureParamEst_rho0.5.eps", width = 12, height = 8)
mycexfont <- 2
myupperY  <- 100
# mylabelX  <- 0.15
mylabelY  <- 90
par(cex.lab=1.2, cex.axis=1.2, mfrow=c(2,3))


# (1)
mylabelX <- 0.1
mybreaks <- seq(from=floor(tmp[1,1]*100)/100,
                to=ceiling(tmp[2,1]*100)/100,
                by=0.01)
hist(res$ci2, 
     # xlim=c(0.6, 1.1),
     ylim=c(0,myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(ciT+0.15, mylabelY, labels = expression('c'[2]==0.8615), cex=mycexfont)
abline(v=EstTable[1,1], col="red", lwd=1.5)
box(which = "plot", lty = "solid")
# (2)
mybreaks <- seq(from=floor(tmp[1,2]*100)/100,
                to=ceiling(tmp[2,2]*100)/100,
                by=0.01)
hist(res$beta0, 
     # xlim=c(EstTable[2,2]-myXlen/2, EstTable[2,2]+myXlen/2),
     # xlim=c(0.2, 0.7),
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(thetaT+0.2, mylabelY, labels = expression(beta[0]==0.4308), cex=mycexfont)
abline(v=EstTable[1,2], col="red", lwd=1.5)
box(which = "plot", lty = "solid")
# (3)
mybreaks <- seq(from=floor(tmp[1,3]*100)/100,
                to=ceiling(tmp[2,3]*100)/100,
                by=0.01)
hist(res$rho1, 
     # xlim=c(EstTable[2,2]-myXlen/2, EstTable[2,2]+myXlen/2),
     # xlim=c(0.2, 0.7),
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(rhoT+0.1, mylabelY, labels = expression(beta[1]==0.5), cex=mycexfont)
abline(v=EstTable[1,3], col="red", lwd=1.5)
box(which = "plot", lty = "solid")
# (4)
qqnorm(res$ci2, pch = 1, main=expression('c'[2]), frame = TRUE, cex.main=mycexfont)
qqline(res$ci2, col = "red", lwd = 1.5)
# (5)
qqnorm(res$beta0, pch = 1, main=expression(alpha[0]), frame = TRUE, cex.main=mycexfont)
qqline(res$beta0, col = "red", lwd = 1.5)
# (6)
qqnorm(res$rho1, pch = 1, main=expression(rho), frame = TRUE, cex.main=mycexfont)
qqline(res$rho1, col = "red", lwd = 1.5)

dev.off()


load("~/Desktop/Changepoint/aop/aop_singleCP/CLS/notes/simulation/MultiInnov_rho-0.5_iii1.RData")

ciT = 0.8615
thetaT = 0.43075
rhoT = -0.5
parT = c(ciT, thetaT, rhoT)

EstTable = matrix(NA, nrow=4, ncol=length(parT))
colnames(EstTable) = expression(c[2], beta[0], rho)
row.names(EstTable) = c("True", "Mean", "Rel.bias", "sd")

EstTable[1,] <- round(parT, digits=4)
EstTable[2,] <- round(colMeans(res[,7:9]), digits=4)
EstTable[3,] <- round((EstTable[1,] - EstTable[2,])/EstTable[1,], digits=4)
EstTable[4,] <- round(apply(res[,7:9], 2, sd), digits=4)
EstTable


tmp <- apply(res[,7:9], 2, range)
tmp

setEPS()
postscript("notes/simulation/FigureParamEst_rho-0.5.eps", width = 12, height = 8)
mycexfont <- 2
myupperY  <- 100
# mylabelX  <- 0.15
mylabelY  <- 90
par(cex.lab=1.2, cex.axis=1.2, mfrow=c(2,3))


# (1)
mylabelX <- 0.1
mybreaks <- seq(from=floor(tmp[1,1]*100)/100,
                to=ceiling(tmp[2,1]*100)/100,
                by=0.01)
hist(res$ci2, 
     # xlim=c(0.6, 1.1),
     ylim=c(0,myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(ciT+0.1, mylabelY, labels = expression('c'[2]==0.8615), cex=mycexfont)
abline(v=EstTable[1,1], col="red", lwd=1.5)
box(which = "plot", lty = "solid")
# (2)
mybreaks <- seq(from=floor(tmp[1,2]*100)/100,
                to=ceiling(tmp[2,2]*100)/100,
                by=0.01)
hist(res$beta0, 
     # xlim=c(EstTable[2,2]-myXlen/2, EstTable[2,2]+myXlen/2),
     # xlim=c(0.2, 0.7),
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(thetaT+0.1, mylabelY, labels = expression(beta[0]==0.4308), cex=mycexfont)
abline(v=EstTable[1,2], col="red", lwd=1.5)
box(which = "plot", lty = "solid")
# (3)
mybreaks <- seq(from=floor(tmp[1,3]*100)/100,
                to=ceiling(tmp[2,3]*100)/100,
                by=0.01)
hist(res$rho1, 
     # xlim=c(EstTable[2,2]-myXlen/2, EstTable[2,2]+myXlen/2),
     # xlim=c(0.2, 0.7),
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(rhoT+0.1, mylabelY, labels = expression(beta[1]==-0.5), cex=mycexfont)
abline(v=EstTable[1,3], col="red", lwd=1.5)
box(which = "plot", lty = "solid")
# (4)
qqnorm(res$ci2, pch = 1, main=expression('c'[2]), frame = TRUE, cex.main=mycexfont)
qqline(res$ci2, col = "red", lwd = 1.5)
# (5)
qqnorm(res$beta0, pch = 1, main=expression(alpha[0]), frame = TRUE, cex.main=mycexfont)
qqline(res$beta0, col = "red", lwd = 1.5)
# (6)
qqnorm(res$rho1, pch = 1, main=expression(rho), frame = TRUE, cex.main=mycexfont)
qqline(res$rho1, col = "red", lwd = 1.5)

dev.off()



