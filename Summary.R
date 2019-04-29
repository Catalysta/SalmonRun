####### Summarize Salmon Data############################################
#######Collapse sites and years into per Region survival stats#########
#However must account for overdispersion from different sites with
#potentially different viability
S.sum<-numeric()
R2.sum<-numeric()
reg<-numeric()
year<-numeric()
k=1
for(r in levels(data$Region)){
  years<-unique(data[which(data$Region==r),]$BY)
  for (y in years){
    S.sum[k]<-sum(data[which(data$Region==r & data$BY==y),"S"])
    R2.sum[k]<-sum(data[which(data$Region==r & data$BY==y),"R2"])
    reg[k]<-r
    year[k]<-y
    k<-(k+1)
  }
}
result<-cbind(reg,year,S.sum,R2.sum)
View(result)


#Try collapsing years as well
S.sum<-numeric()
R2.sum<-numeric()
reg<-numeric()
k=1
for(r in levels(data$Region)){
    S.sum[k]<-sum(data[which(data$Region==r),"S"])
    R2.sum[k]<-sum(data[which(data$Region==r& !is.na(data$R2)),"R2"])
    reg[k]<-r
    k<-(k+1)
  }
result<-as.data.frame(cbind(reg,S.sum,R2.sum))
for (s in levels(result$reg)){
  region<- ifelse(result$reg==s,1,0)
  result<-cbind(result,region)
}
colnames(result)[2:15]<-c("S","R2","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12")
levels(result$reg)<-levels(data$Region)
#write.table(data,"salmon.txt",sep="\t",row.names=FALSE)
write.table(result,"salmon.txt",sep="\t",row.names=FALSE)


