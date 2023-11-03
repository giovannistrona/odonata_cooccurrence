###random forest model 1
library(viridis)
library(randomForest)
library(akima)
library(fields)



###model regressor
data_f<-read.csv('complete_dataset_rf.csv',header=T)

res_continuous<-c()
oobs_continuous<-c()
for (rep in 1:100){
  rand_ids<-sample(1:nrow(data_f),nrow(data_f))
  data_train<-data_f[rand_ids[1:round(nrow(data_f)*0.8)],]
  data_test<-data_f[rand_ids[(round(nrow(data_f)*0.8)+1):nrow(data_f)],]
  rf <-randomForest(cooc~.,data=data_train, mtry=48, importance=TRUE,ntree=500)
  res_continuous<-rbind(res_continuous,cbind(predict(rf,data_test),data_test$cooc))
  oobs_continuous<-c(oobs_continuous,mean(rf$rsq*100)) 
  print (rep)
  }


cor.test(res_continuous[,1],res_continuous[,2])[[4]]**2
plot(res_continuous[,1],res_continuous[,2],pch='.')
dim(res_continuous)


###pos + neg
data_f<-read.csv('complete_dataset_rf.csv',header=T)

tre<-0.8
res_both<-c()

for (tre in seq(0.4,0.9,0.05)){
  for (rep in 1:100){
    data<-data_f
    data$cooc<-1*(abs(data$cooc)>tre)
    n<-min(table(data$cooc))
    data_bal<-data[c(sample(which(data$cooc==0),n),sample(which(data$cooc==1),n)),]
    rf <-randomForest(as.factor(cooc)~.,data=data_bal, mtry=48, importance=TRUE,ntree=500)
    res_both<-rbind(res_both,c(tre,rf[[5]][,3]))
    print (c(rep,'both',res_both[nrow(res_both),]))}
    }



###positive
data_f<-read.csv('complete_dataset_rf.csv',header=T)
res_pos<-c()
tre<-0.9
for (tre in seq(0.4,0.9,0.05)){
  for (rep in 1:100){
    data<-data_f
    data$cooc<-1*(data$cooc>tre)
    n<-min(table(data$cooc))
    data_bal<-data[c(sample(which(data$cooc==0),n),sample(which(data$cooc==1),n)),]
    rf <-randomForest(as.factor(cooc)~.,data=data_bal, mtry=48, importance=TRUE,ntree=500)
    res_pos<-rbind(res_pos,c(tre,rf[[5]][,3]))
    print (c(rep,'pos',res_pos[nrow(res_pos),]))}
  }



###negative
data_f<-read.csv('complete_dataset_rf.csv',header=T)
res_neg<-c()
for (tre in seq(0.4,0.9,0.05)){
  for (rep in 1:100){
    data<-data_f
    data$cooc<-1*(data$cooc<(-tre))
    n<-min(table(data$cooc))
    data_bal<-data[c(sample(which(data$cooc==0),n),sample(which(data$cooc==1),n)),]
    rf <-randomForest(as.factor(cooc)~.,data=data_bal, mtry=48, importance=TRUE,ntree=500)
    res_neg<-rbind(res_neg,c(tre,rf[[5]][,3]))
    print (c(rep,'neg',res_neg[nrow(res_neg),]))}
  }


###save results
colnames(res_continuous)<-c('predicted','observed')
colnames(res_both)<-c('treshold','I','II')
colnames(res_pos)<-c('treshold','I','II')
colnames(res_neg)<-c('treshold','I','II')

write.table(res_continuous,'res_continuous.csv',
            row.names=F,col.names=T,sep=',',quote=F)

write.table(res_both,'res_unsigned.csv',
            row.names=F,col.names=T,sep=',',quote=F)

write.table(res_pos,'res_positive.csv',
            row.names=F,col.names=T,sep=',',quote=F)

write.table(res_neg,'res_negative.csv',
            row.names=F,col.names=T,sep=',',quote=F)


###variable importance
###continuous
rf_cont <-randomForest(cooc~.,data=data_f, mtry=48, importance=TRUE,ntree=500)

tre<-0.5
##unsigned
data<-data_f
data$cooc<-1*(abs(data$cooc)>tre)
n<-min(table(data$cooc))
data_bal<-data[c(sample(which(data$cooc==0),n),sample(which(data$cooc==1),n)),]
rf_un <-randomForest(as.factor(cooc)~.,data=data_bal, mtry=48, importance=TRUE,ntree=500)

##positive
data<-data_f
data$cooc<-1*(data$cooc>tre)
n<-min(table(data$cooc))
data_bal<-data[c(sample(which(data$cooc==0),n),sample(which(data$cooc==1),n)),]
rf_pos <-randomForest(as.factor(cooc)~.,data=data_bal, mtry=48, importance=TRUE,ntree=500)

##negative
data<-data_f
data$cooc<-1*(data$cooc<(-tre))
n<-min(table(data$cooc))
data_bal<-data[c(sample(which(data$cooc==0),n),sample(which(data$cooc==1),n)),]
rf_neg <-randomForest(as.factor(cooc)~.,data=data_bal, mtry=48, importance=TRUE,ntree=500)


cont_vi<-rf_cont$importance
un_vi<-rf_un$importance
pos_vi<-rf_pos$importance
neg_vi<-rf_neg$importance


write.table(cont_vi,'var_imp_continuous.csv',
            row.names=T,col.names=T,sep=',',quote=F)

write.table(un_vi,'var_imp_unsigned.csv',
            row.names=T,col.names=T,sep=',',quote=F)

write.table(pos_vi,'var_imp_positive.csv',
            row.names=T,col.names=T,sep=',',quote=F)

write.table(neg_vi,'var_imp_negative.csv',
            row.names=T,col.names=T,sep=',',quote=F)



vi_data<-list(cont_vi,pos_vi,neg_vi)
title<-c('continuous cooccurrence','positive cooccurrence','negative cooccurrence')
pdf('var_imp.pdf',width=9,height=12)
par(mfrow=c(3,2))
for (n in 1:3){
  vi_mse<-(vi_data[[n]][1:24,1]+vi_data[[n]][25:48,1])/2
  barplot(rev(rev(sort(vi_mse))[1:5]),horiz=T,las=1,cex.names=1,xlab='MSE increment',main=title[n])
  
  vi_np<-(vi_data[[n]][1:24,2]+vi_data[[n]][25:48,2])/2
  barplot(rev(rev(sort(vi_np))[1:5]),horiz=T,las=1,cex.names=1,xlab='node purity increment',main=title[n])
}
dev.off()


####
pdf('continuous_model.pdf',width = 4.5,height=4.5)
plot(res_continuous[,1],res_continuous[,2],pch=16,col=plasma(1,alpha=0.05)[1],xlim=c(-1,1),
     ylab='observed',xlab='predicted',las=1,cex.lab=1.3,cex.axis=1.3,cex=0.5)

abline(0,1,col='red',lwd=2)
cor.test(res_continuous[,1],res_continuous[,2])[[4]]**2
dev.off()

###
x<-seq(0.4,0.9,0.05)
res<-list(res_pos,res_neg)
title<-c('positive cooccurrence','negative cooccurrence')
pdf('models.pdf',height=3,width=9)
par(mfrow=c(1,3))
###type I: non-cooc detected; type II: cooc not detected
bp_pos<-apply(1-res[[1]][res[[1]][,1]==0.5,c(2,3)],1,mean)
bp_neg<-apply(1-res[[2]][res[[2]][,1]==0.5,c(2,3)],1,mean)

print (c('positive',round(mean(bp_pos),2),round(sd(bp_pos),3)))
print (c('negative',round(mean(bp_neg),2),round(sd(bp_neg),3)))


boxplot(bp_pos,bp_neg,names=c('positive','negative'),
        ylab='',xlab='accuracy',
        las=1,cex.axis=1.2,cex.lab=1.2,outline=F,horizontal = T)



for (n in 1:2){
  y1<-aggregate(res[[n]][,2]~res[[n]][,1], FUN='mean')[,2]
  y2<-aggregate(res[[n]][,3]~res[[n]][,1], FUN = 'mean')[,2]
  ci1<-aggregate(res[[n]][,2]~res[[n]][,1], FUN='sd')[,2]
  ci2<-aggregate(res[[n]][,3]~res[[n]][,1], FUN = 'sd')[,2]

  plot(x,y1,ylim=c(min((y1-ci1),(y2-ci2)),max((y1+ci1),(y2+ci2))),
       type='n',xlab='threshold',ylab='errors',las=1,cex.axis=1.2,cex.lab=1.2,
       main = title[n])

  polygon(c(rev(x), x), c(rev(y1-ci1),(y1+ci1)),col = plasma(2,alpha=0.5)[1], border = NA)
  lines(x,y1,col = plasma(3)[1],lwd=2)

  polygon(c(rev(x), x), c(rev(y2-ci2),(y2+ci2)),col = plasma(3,alpha=0.5)[2], border = NA)
  lines(x,y2,col = plasma(3)[2],lwd=2)
  legend(0.4,max(y1+ci1,y2+ci2),c('type I','type II'),col=plasma(3)[1:2],pch=15,pt.cex=2)

}        


dev.off()



###confusion matrices
library(ggplot2)
library(gridExtra)

predicted <- factor(c(0, 1, 0, 1))
observed <- factor(c(0, 0, 1, 1))

conf_mat_plots<-list()
for (n in 1:2){
  Y<-c(mean(1-res[[1]][res[[n]][,1]==0.5,2]),
                mean(res[[n]][res[[n]][,1]==0.5,2]),
                mean(res[[n]][res[[n]][,1]==0.5,3]),
                mean(1-res[[n]][res[[n]][,1]==0.5,3]))



print(title[n])
print (c(mean(res[[n]][res[[n]][,1]==0.5,2]),
  sd(res[[n]][res[[n]][,1]==0.5,2]),
  mean(res[[n]][res[[n]][,1]==0.5,3]),
  sd(res[[n]][res[[n]][,1]==0.5,3])))

  
  # ###obs/pred
# abs/abs; abs/pres
# pres/abs; pres/pres 

  df <- data.frame(predicted, observed, Y)
  conf_mat_plots[[n]]<-ggplot(data =  df, mapping = aes(x = predicted, y = observed)) +
    geom_tile(aes(fill = Y), colour = "white") +
    geom_text(aes(label = round(Y,2)), vjust = 1,colour = "white") +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_bw() + theme(legend.position = "none")+ggtitle(title[n])+
    theme(plot.title = element_text(hjust = 0.5))
}

pdf('confusion_matrices.pdf',width=6,height=3)
grid.arrange(conf_mat_plots[[1]],conf_mat_plots[[2]],nrow=1)
dev.off()




###
a<-read.csv('supported_co-occurrence.csv',header=T,row.names=1)
dim(a)
sc<-0
pos<-0
neg<-0
for (i in 1:nrow(a)){
  for (j in 1:ncol(a)){
    if (i<j){
      sc<-sc+1
      if (a[i,j]>0.5){pos<-pos+1}
      if (a[i,j]<(-0.5)){neg<-neg+1}
        }
      }
    }


100*pos/sc
100*neg/sc


#####
a<-read.csv('cooc_vs_dist.csv',header=T)
pdf('figS3.pdf',width=5,height=5)
plot(a$cooc,a$dist,pch='.',las=1,cex.axis=1.2,cex.lab=1.2,xlab='co-occurrence',ylab='functional distance')
dev.off()

cor.test(a$cooc,a$dist)[[4]]**2

a$dist<-(a$dist-min(a$dist,na.rm=T))/(max(a$dist,na.rm=T)-min(a$dist,na.rm=T))
plot(a$cooc,a$dist,pch='.',las=1,cex.axis=1.2,cex.lab=1.2,xlab='co-occurrence',ylab='functional distance')

abline(lm(a$dist~a$cooc))

pdf('diagram_bi_vs_ef.pdf',width=4.5,height=4.5)
plot(0,0,type='n',xlim=c(-1,1),ylim=c(0,1),ylab='functional distance',xlab='co-occurrence',
     las=1,cex.axis=1.5,cex.lab=1.5)

lines(c(-1,1),c(0,1),col='blue',lwd=3)
lines(c(-1,1),c(1,0),col='forestgreen',lwd=3)
dev.off()

library(hexbin)
library(RColorBrewer)
library(viridis)
# Create data
# Make the plot
bin<-hexbin(a$cooc, a$dist, xbins=30)
my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
pdf('hexplot.pdf',height=4.5,width=6)
plot(bin, main="" , colramp=colorRampPalette(hcl.colors(12)),
     ylab='functional distance',xlab='co-occurrence'
     ) #, legend=F ) 
     
dev.off()
     



cors<-c()
for (i in 5:ncol(a)){
  cors<-c(cors,lm(a$cooc~a[,i])$coefficients[2])
}

xmax<-max(abs(cors),na.rm=T)
xmin<-(-xmax)

pdf('slopes.pdf',height=4.5,width=4.5)
hist(cors,8,las=1,cex.axis=1.2,cex.lab=1.2,xlab='slope',ylab='frequency',main='',
     xlim=c(xmin,xmax))
dev.off()
