# read in .dat map coordinates
map.coords<- read.table('map_coords.dat', header = FALSE)

#plot
#resize.win(20,15)

#postscript("coastpic.eps")#, width=16, height=5)

plot(map.coords[,1],map.coords[,2], type='l', col=gray(0.3), xlab="Degrees of longitude (W)", ylab="Degrees of latitude (N)")
leg.text<-c("Pink salmon spawning sites")#,"Chum","Sockeye")
legend(190,50,leg.text, pch=c(20,21,3),col='black')
# read in .csv data
stock.info<-read.csv('stockinfo.csv')
stock.info$Longitude<-360-stock.info$Longitude
dimnames(stock.info)[[2]]

# plot actual stock locations
lines(stock.info$Longitude[which(stock.info$Species=='pink')],stock.info$Latitude[which(stock.info$Species=='pink')],type='p',pch=20,cex=1)
#lines(stock.info$Longitude[which(stock.info$Species=='chum')],stock.info$Latitude[which(stock.info$Species=='chum')],type='p',pch=21,cex=1.4,lwd=1)
#lines(stock.info$Longitude[which(stock.info$Species=='sockeye')],stock.info$Latitude[which(stock.info$Species=='sockeye')],type='p',pch=3,cex=1.4,lwd=1)

#dev.off()

