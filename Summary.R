####### Summarize Salmon Data############################################
#######Collapse sites and years into per Region survival stats#########
S.sum<-numeric(12)
R2.sum<-numeric(12)
for(i in 1:12){
  colname<-paste("X",i,sep="")
  indices<-which(salmon.data[,colname]==1)
  S.sum[i]<-sum(salmon.data$S[indices])
  R2.sum[i]<-sum(salmon.data$R2[indices])
}
salmon.sum<-cbind(S.sum,R2.sum)
