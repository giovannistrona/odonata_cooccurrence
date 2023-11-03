library(RColorBrewer)
library(BayesLogit)
library(Hmsc)
library(tidyr)
library(reshape2)
library(plyr)
library(abind)
library(viridis)
library(ggplot2)

data<-read.csv('odonata_for_hmsc.csv',header=T)
#colnames(data)[1:30]
data<-data[,-9] #remove hfi
data<-data[is.finite(rowSums(data)),]
data<-data[data$interp_t/100<40,]
data<-data[data$interp_t/100>0,]
data<-data[,which(colSums(data)!=0)]

data<-data[which(rowSums(data[,c(14:ncol(data))])>1),]

XData<-data[,c(1:13)]
Y<-as.matrix(data[,c(14:ncol(data))])
Y<-Y[,which(colSums(Y)>9)]

xycoords<-XData[,c(1,2)]
xycoords[,1]<-xycoords[,1]+runif(nrow(xycoords))/100000

XData<-XData[,3:13]
row.names(XData)<-1:dim(XData)[1]
row.names(Y)<-1:dim(XData)[1]
row.names(xycoords)<-1:dim(XData)[1]

n<-dim(XData)[1]	#number of sites

studyDesign <- data.frame(Site = as.factor(row.names(XData)),
                          Year = as.factor(XData$year))


thin = 100
samples = 250
nChains = 4
transient = as.integer(.5*(samples*thin))

#Define random levels
rL1=HmscRandomLevel(sData=xycoords)
rL2=HmscRandomLevel(units=unique(studyDesign$year))

XFormula= ~ as.factor(lu) +
            elev +
            poly(interp_t, degree = 2, raw=TRUE) +
            poly(interp_p, degree = 2, raw=TRUE) +
            tree_canopy +
            grass +
            soil +
            water +
            artificial


m <- Hmsc(Y=Y, XData=XData, 
          XFormula=XFormula,
          studyDesign=studyDesign,
          ranLevels=list("Site"=rL1,"Year"=rL2),
          distr="probit")



a=Sys.time()
m <- sampleMcmc(m, samples = samples, thin = thin, transient = transient,
nChains = nChains, nParallel = nChains, updater=list(GammaEta=FALSE),verbose=1)
print (paste('analysis_time',Sys.time()-a))

####
post=convertToCodaObject(m)
predY=computePredictedValues(m)
MF = evaluateModelFit(hM = m, predY = predY)
omega_cor <- computeAssociations(m)

save(m,file='odonata_hmsc_m.Rdata')
save(MF,file='odonata_hmsc_MF.Rdata')
save(post,file='odonata_hmsc_post.Rdata')
save(omega_cor,file='odonata_hmsc_omega_cor.Rdata')


cooc_m<-omega_cor[[1]]$mean
cooc_s<-omega_cor[[1]]$support

not_pos<-which((cooc_s>0.95)*(cooc_m>0)!=TRUE)
not_neg<-which((cooc_s<0.05)*(cooc_m<0)!=TRUE)

pos_sig<-cooc_m
pos_sig[not_pos]<-0

neg_sig<-cooc_m
neg_sig[not_neg]<-0

to_plot<-pos_sig+neg_sig




library(corrplot)
library(viridis)
mean_cooc<-omega_cor[[1]]$mean
pdf('libellule_co-occurrence.pdf')
par(oma=c(0,0,1,0))
corrplot(mean_cooc,is.corr = F,
         #addCoef.col = 'black',
         mar = c(1,0,2,0),
         main = 'mean co-occurrence',
         addCoefasPercent = F,
         number.cex = 0.3,
         number.digits = 2,  
         col=colorRampPalette(c("darkblue","white","darkgreen"))(200),diag=F,method='color',
         addgrid.col=NA,order='hclust',hclust.method="complete",tl.col='black',tl.cex=0.2,
         cl.lim = c(-1,1))

supp_cooc<-omega_cor[[1]]$support

corrplot(supp_cooc,is.corr = F,
         mar = c(0,0,2,0),
         #addCoef.col = 'black',
         main = 'support',
         addCoefasPercent = F,
         number.cex = 0.3,
         number.digits = 2,
         col=colorRampPalette(c("darkblue","darkblue","darkblue","white","darkgreen"))(200),diag=F,method='color',
         addgrid.col=NA,order='hclust',hclust.method="complete",tl.col='black',tl.cex=0.2,
         cl.lim = c(0,1))


corrplot(to_plot,is.corr = F,
         #addCoef.col = 'black',
         mar = c(0,0,2,0),
         main = 'supported co-occurrence',
         addCoefasPercent = F,
         number.cex = 0.3,
         number.digits = 2,  
         col=colorRampPalette(c("darkblue","white","darkgreen"))(200),diag=F,method='color',
         addgrid.col=NA,order='hclust',hclust.method="complete",tl.col='black',tl.cex=0.2,
         cl.lim = c(-1,1))

dev.off()


write.table(mean_cooc,'mean_co-occurrence.csv',sep=',',quote = F)
write.table(supp_cooc,'support.csv',sep=',',quote = F)
write.table(to_plot,'supported_co-occurrence.csv',sep=',',quote = F)


#model fit
load('./hmsc_files/odonata_hmsc_MF.Rdata')
load('./hmsc_files/odonata_hmsc_post.Rdata')
load('./hmsc_files/odonata_hmsc_m.Rdata')


sink('fit.txt')
print(paste('meanTjurR2',mean(MF$TjurR2)))
print(paste('meanRMSE',mean(MF$RMSE)))
print(paste('meanAUC',mean(MF$AUC)))
sink()

pdf('FigS2.pdf',height=6,width=6)
par(mfrow=c(2,2))
hist(effectiveSize(post$Beta), xlab="ess(beta)",main='',cex.axis=1.2,cex.lab=1.2,las=1)
hist(gelman.diag(post$Beta, multivariate=FALSE)$psrf, xlab="psrf(beta)",main='',cex.axis=1.2,cex.lab=1.2,las=1)
hist(MF$AUC,12,xlab='AUC',main='',cex.axis=1.2,cex.lab=1.2,las=1)
hist(MF$TjurR2,12,xlab='TjurR2',main='',cex.axis=1.2,cex.lab=1.2,las=1)

dev.off()


# interaction_list<-melt(omega_cor[[1]]$mean)
# interaction_list <- data.frame(interaction_list)
# g <- graph_from_data_frame(interaction_list, directed = F)


pdf("betapost.pdf")
plot(post$Beta)
dev.off()

beta=getPostEstimate(m, "Beta")
pdf("Beta_0.95.pdf")
plotBeta(m, beta, supportLevel=.95,spNamesNumbers=c(F,F))
dev.off()

predY = apply(abind(predY,along=3),c(1,2),mean)
pdf("TjurR2.pdf")
plot(colSums(m$Y)/m$ny,MF$TjurR2,main=paste("Explanatory Tjur R2. Mean = ", round(mean(MF$TjurR2,na.rm = TRUE),2),".", sep=""), xlab = "Prevalence")
dev.off()


pdf("VP.pdf")
group=c(1,2,3,3,4,4,5,6,7,8,9)
groupnames = c('lu','elev','temp','prec','tree','grass','soil','water','artificial')
VP = computeVariancePartitioning(hM = m, group = group, groupnames = groupnames)
plotVariancePartitioning(m, VP,col=brewer.pal(n = 6, name = "Dark2"))
dev.off()





#########
cooc_res<-read.csv('supported_co-occurrence.csv',header=T,row.names=1)
ncomb<-nrow(cooc_res)*(ncol(cooc_res)-1)

sum(cooc_res>0.5)/2/ncomb
sum(cooc_res<(-0.5))/2/ncomb

hist(as.matrix(cooc_res))


###plot VP as boxplots
###sum the land use categories

vp_bp_data<-data.frame(t(VP[[1]]))
vp_bp_data$luc_frac<-vp_bp_data$tree+vp_bp_data$grass+vp_bp_data$soil+vp_bp_data$artificial+vp_bp_data$water
vp_bp_data<-data.frame('habitat' = vp_bp_data$lu,
                       'elevation' = vp_bp_data$elev,
                       'temperature' = vp_bp_data$temp,
                       'precipitation' = vp_bp_data$prec,
                       'land use' = vp_bp_data$luc_frac,
                       'random year' = vp_bp_data$Random..Year,
                       'random_site' = vp_bp_data$Random..Site)


tot_var<-vp_bp_data$habitat+
         vp_bp_data$elevation+
         vp_bp_data$temperature+
         vp_bp_data$precipitation+
         vp_bp_data$land.use


mean(tot_var)
mean(vp_bp_data$random.year+vp_bp_data$random_site)


vp_bp_data$total_variables<-tot_var
vp_bp_data$total_re<-vp_bp_data$random.year+vp_bp_data$random_site
pdf('VP_boxplots.pdf',height=6,width=9)
boxplot(100*vp_bp_data,horizontal = T,las = 1,xlab = 'explained variance (%)')
dev.off()

