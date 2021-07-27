############################################################################################################################################################################
#Activate necessary packages
############################################################################################################################################################################

library("gplots")
library("RColorBrewer")
library("ggplot2")


############################################################################################################################################################################
#File input
############################################################################################################################################################################
setwd("~/MAMBE/Prep module for MT/Microscopy/20210310_woMMC/WT/T4h")
fold<-"Membranes*"
folder<-paste("../",Sys.glob(fold),sep="")

setwd("~/MAMBE/Prep module for MT/Microscopy/20210310_woMMC/WT/T4h/Membranes01")  

PlotType<-"CellNorm"  #choose between CellNorm and PopNorm

############################################################################################################################################################################

chan<-"Blue*"
maxim<-0
mayim<-0
minyim<-17000

for (f in folder){
  setwd(f)
  RED<-Sys.glob(chan)
  for(k in 1:length(RED)){
    filer<-read.table(RED[k], sep="\t",dec=".", stringsAsFactors=FALSE,header=TRUE)
    if (max(filer$x)>maxim){
      maxim<-max(filer$x)
    }
    if (PlotType=="CellNorm"){
      if (max(filer$y)>mayim){                 #activate if normalize for maximum overall fluorescence 
        mayim<-max(filer$y)                    #
        }                                        #
      if (min(filer$y)<minyim){                #activate if normalize for maximum overall fluorescence 
        minyim<-min(filer$y)                   #
    }                                        #
    }
  }}

a<-seq(0,maxim, by=0.1)                      #define binwidth (0.02 um for PALM, 0.2 for EPI)
REDTOT<-data.frame()

counter<-1
for (f in folder){
  setwd(f)
  RED<-Sys.glob(chan)
  for(k in 1:length(RED)){
    file<-read.table(RED[k], sep="\t",dec=".", stringsAsFactors=FALSE,header=TRUE) # sep="," for ratio
    file<-as.data.frame(cbind(file$x,file$y))
    colnames(file)[1]<-"x"
    colnames(file)[2]<-"y"
    file[3]<-length(file$x)+(0.0001*k) #risk that two cells of the same length from different folders come from a file with same name...think of solution
    file[4]<-counter
    counter<-counter+1
    REDTOT<-rbind(REDTOT,file)
  }
  REDTOT1<-REDTOT[order(REDTOT$V3),]
}

df2<-data.frame()
for(k in 1:length(unique(REDTOT1$V4))){
  print(k)
  file<-REDTOT1[REDTOT1$V4==unique(REDTOT1$V4)[k],]
  if (PlotType=="CellNorm"){
    mayim<-max(file$y)            ###if active it normalizes on max cell intensity
    minyim<-min(file$y)           ###if active it normalizes on max cell intensity
  }
  mayim_c<-(mayim-minyim)
  print(mayim)
  for (i in 1:(length(a)-1)){
    if (a[i]<max(file$x)){
      m<-mean(file[file$x>=a[i] & file$x<a[i+1],]$y)
      if (is.nan(m)) {
        m<-mean(file[file$x>=a[i-1] & file$x<a[i+1],]$y)
      }
      m_c<-(m-minyim) #correct by minimum
      b<-replicate(1+(m_c*100/mayim_c),k) #for epi: m*100/mayim OR 1+(m_c*100/mayim_c)  #for palm: 1+((m/mayim)*1000))  
      c<-replicate(1+(m_c*100/mayim_c),a[i])  #for epi: m*100/mayim OR 1+(m_c*100/mayim_c)  #for palm: 1+((m/mayim)*1000))
      df1<-as.data.frame(c)
      df1[2]<-b
      df2<-rbind(df2,df1)
    }
  }
}

write.table(REDTOT1, file="../WT_CellNormMinMax_profiles_ordered.txt",sep=",",row.names = FALSE, quote=FALSE)   ###save file containing the fluorescence profiles ordereed by cell length
write.table(df2, file="../WT_CellNormMinMax_profiles_ordered_matrix.txt",sep=",",row.names = FALSE, quote=FALSE)  ###save matrix used to represent the fluorescence profiles as demographs via hist2d()
write.table(maxim, file="../WT_CellNormMinMax_maxim.txt",sep=",",row.names = FALSE, quote=FALSE)  ###save matrix used to represent the fluorescence profiles as demographs via hist2d()

#REDTOT1<-read.table(file="../Red_profiles_ordered.txt",sep=",", header=TRUE)   ###save file containing the fluorescence profiles ordereed by cell length
#df2<-read.table(file="../Red_profiles_ordered_matrix.txt",sep=",", header=TRUE)  ###save matrix used to represent the fluorescence profiles as demographs via hist2d()
#maxim<-read.table(file="../Red_maxim.txt",sep=",", header=TRUE)  ###save matrix used to represent the fluorescence profiles as demographs via hist2d()
#a<-seq(0,maxim[1,1], by=0.1)

rf1<-colorRampPalette(c("black","blue","red")) ####Define color scale
r <- rf1(256)

rf1<-colorRampPalette(c("black","white")) ####Define color scale
r <- rf1(256)

png(file="../Demograph.png",height=3000,width=2500,res=600)   ####Demograph
h2<-hist2d(df2, nbins=c(length(a)-1,length(unique(REDTOT1$V4))), col=r)
dev.off()


######################Peaks analysis
#########Establish FindPeak Function

findpeaks <- function(y, span = NULL)
{
  if (is.null(span)) span <- round(.2 * length(y))
  
  z <- embed(y, span)
  s <- span %/% 2
  v <- max.col(z, ties.method = "first") == 1 + s
  
  which(c(rep(FALSE, s), v, rep(FALSE, s)))
}
#############################################
setwd("~/MAMBE/Prep module for MT/Microscopy/20210310_woMMC/WT/T4h")  
fold<-"Membranes*"
folder<-paste("../",Sys.glob(fold),sep="")
setwd("~/MAMBE/Prep module for MT/Microscopy/20210310_woMMC/WT/T4h/Membranes01")  
 

chan<-"Red*"
chan1<-"Blue*"

cellsTOT<-data.frame()
#f<-folder[1]

for (f in folder){
  setwd(f)
  RED<-Sys.glob(chan)
  BLUE<-Sys.glob(chan1)
  cells<-as.data.frame(seq(from=1, to=length(RED),by=1))
  colnames(cells)[1]<-"Index"
  cells[2]<-0
  colnames(cells)[2]<-"Length"
  cells[3]<-0
  colnames(cells)[3]<-"RedPeak"
  cells[4]<-0
  colnames(cells)[4]<-"BluePeak"
  for(k in 1:length(RED)){
    filer<-read.table(RED[k], sep="\t",dec=".", stringsAsFactors=FALSE,header=TRUE)
    fileb<-read.table(BLUE[k], sep="\t",dec=".", stringsAsFactors=FALSE,header=TRUE)
    
    png(file=paste("Cell",k,".png",sep=""),height=2000,width=3000,res=600)  ####Demograph
    print(ggplot(data=fileb, aes(x=x, y=y))+
            geom_line(col="blue")+
            geom_point(data=fileb[findpeaks(fileb$y, span=10),][fileb[findpeaks(fileb$y, span=10),]$y-min(fileb$y)>(0.2*(max(fileb$y)-min(fileb$y))),], aes(x=x, y=y), col="blue")+
            geom_line(data=filer, aes(x=x, y=y), col="red")+
            geom_point(data=filer[findpeaks(filer$y, span=10),][filer[findpeaks(filer$y, span=10),]$y-min(filer$y)>(0.2*(max(filer$y)-min(filer$y))),], aes(x=x, y=y), col="red")+
            theme_bw())
    dev.off()
    
    cells[k,2]<-max(filer$x)
    TempPeakRed<-findpeaks(filer$y, span=10)
    cells[k,3]<-length(filer[TempPeakRed,][filer[TempPeakRed,]$y-min(filer$y)>(0.2*(max(filer$y)-min(filer$y))),]$x)
    TempPeakBlue<-findpeaks(fileb$y, span=10)
    cells[k,4]<-length(fileb[TempPeakBlue,][fileb[TempPeakBlue,]$y-min(fileb$y)>(0.2*(max(fileb$y)-min(fileb$y))),]$x)
    
  }
  cellsTOT<-rbind(cellsTOT,cells)
  cellsTOT$Index<-seq(1,length(cellsTOT$Index), by=1)
}
write.table(cellsTOT, file="../CellLength_and_Peaks_M_mc.txt",sep=",",row.names = FALSE, quote=FALSE)   ###save file containing the fluorescence profiles ordereed by cell length

####plot single profiles and peaks
#for (f in folder){
#  setwd(f)
#  RED<-Sys.glob(chan)
#  BLUE<-Sys.glob(chan1)
#for(w in 1:length(RED)){
#fileb<-read.table(BLUE[w], sep="\t",dec=".", stringsAsFactors=FALSE,header=TRUE)
#filer<-read.table(RED[w], sep="\t",dec=".", stringsAsFactors=FALSE,header=TRUE)

#4h_20210209_new_mito
RESK<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/20210209/ResK/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB1<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/20210209/DipB1/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB2<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/20210209/DipB2/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB3<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/20210209/DipB3/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)


#2h_20210209_new_mito
RESK<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/20210209_newMMC/ResK/T2h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB1<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/20210209_newMMC/DipB1/T2h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB2<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/20210209_newMMC/DipB2/T2h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB3<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/20210209_newMMC/DipB3/T2h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)

#4h_20210202_nomito
RESK<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/Analysed/20210202/RES/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB1<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/Analysed/20210202/dipB1/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB2<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/Analysed/20210202/DipB2/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB3<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/Analysed/20210202/DipB3/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)

#4h_20210204_mito
RESK<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/Analysed/20210204/ResK/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB1<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/Analysed/20210204/DipB1/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB2<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/Analysed/20210204/DipB2/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)
dipB3<-read.table("X:/Giacomo Giacomelli/DIPS_project/Students/Bente/Analysed/20210204/DipB3/T4h/CellLength_and_Peaks_M_mc.txt", sep=",", header=TRUE)

ggplot(data=RESK, aes(x=Length, y=..density..))+
  geom_histogram(fill="royalblue4", alpha=1, binwidth = 0.25)+
  geom_histogram(data=dipB3, aes(x=Length, y=..density..), fill="orange", alpha=0.5, binwidth=0.25)+
  geom_histogram(data=dipB2, aes(x=Length, y=..density..), fill="yellow", alpha=0.5, binwidth=0.25)+
  geom_histogram(data=dipB1, aes(x=Length, y=..density..), fill="red", alpha=0.5, binwidth=0.5)+
  theme_bw()

shapiro.test(RESK$Length)
kruskal.test(RESK$Length,dipB1$Length)
wilcox.test(RESK$Length,dipB1$Length)

k



ggplot(data=RESK, aes(x=RedPeak, y=..density..))+
  geom_histogram(fill="royalblue4", alpha=1, binwidth = 1)+
  #geom_histogram(data=dipB3, aes(x=RedPeak, y=..density..), fill="orange", alpha=0.5, binwidth=1)+
  geom_histogram(data=dipB2, aes(x=RedPeak, y=..density..), fill="yellow", alpha=0.5, binwidth=1)+
  #geom_histogram(data=dipB1, aes(x=RedPeak, y=..density..), fill="red", alpha=0.5, binwidth=1)+
  theme_bw()

ggplot(data=RESK, aes(x=BluePeak, y=..density..))+
  geom_histogram(fill="royalblue4", alpha=1, binwidth = 1)+
  #geom_histogram(data=dipB3, aes(x=BluePeak, y=..density..), fill="orange", alpha=0.5, binwidth=1)+
  geom_histogram(data=dipB2, aes(x=BluePeak, y=..density..), fill="yellow", alpha=0.5, binwidth=1)+
  #geom_histogram(data=dipB1, aes(x=BluePeak, y=..density..), fill="red", alpha=0.5, binwidth=1)+
  theme_bw()

ggplot(data=RESK, aes(x=Length, y=BluePeak))+
  geom_point(col="royalblue4", alpha=1)+
  geom_point(data=dipB3, aes(x=Length, y=BluePeak+0.05), col="orange", alpha=0.5)+
  geom_point(data=dipB2, aes(x=Length, y=BluePeak+0.10), col="yellow", alpha=0.5)+
  geom_point(data=dipB1, aes(x=Length, y=BluePeak+0.15), col="red", alpha=0.5)+
  theme_bw()

ggplot(data=RESK, aes(x=Length, y=RedPeak))+
  geom_point(col="royalblue4", alpha=1)+
  geom_point(data=dipB3, aes(x=Length, y=RedPeak+0.05), col="orange", alpha=0.5)+
  #geom_point(data=dipB2, aes(x=Length, y=RedPeak+0.10), col="yellow", alpha=0.5)+
  geom_point(data=dipB1, aes(x=Length, y=RedPeak+0.15), col="red", alpha=0.5)+
  theme_bw()

ggplot(data=RESK, aes(x=BluePeak, y=RedPeak))+
  geom_point(col="royalblue4", alpha=1)+
  geom_point(data=dipB3, aes(x=BluePeak+0.05, y=RedPeak+0.05), col="orange", alpha=0.5)+
  geom_point(data=dipB2, aes(x=BluePeak+0.10, y=RedPeak+0.10), col="yellow", alpha=0.5)+
  geom_point(data=dipB1, aes(x=BluePeak+0.15, y=RedPeak+0.15), col="red", alpha=0.5)+
  theme_bw()



shapiro.test(RESK$Length)
shapiro.test(dipB3$Length)

RESK[5]<-"ResK"
colnames(RESK)[5]<-"Strain"
dipB3[5]<-"dipB3"
colnames(dipB3)[5]<-"Strain"
dipB2[5]<-"dipB2"
colnames(dipB2)[5]<-"Strain"
dipB1[5]<-"dipB1"
colnames(dipB1)[5]<-"Strain"

TOT<-rbind(RESK, dipB1, dipB2, dipB3)

wilcox.test(Length ~ Strain, data=TOT) 
wilcox.test(RedPeak ~ Strain, data=TOT) 
wilcox.test(BluePeak ~ Strain, data=TOT) 
library("pgirmess")
kruskalmc(Length ~ Strain, data=TOT)



#png(file=paste("Cell",w,".png",sep=""),height=2000,width=3000,res=600)  ####Demograph
#print(ggplot(data=fileb, aes(x=x, y=y))+
#  geom_line(col="blue")+
#  geom_point(data=fileb[findpeaks(fileb$y, span=10),][fileb[findpeaks(fileb$y, span=10),]$y-min(fileb$y)>(0.4*(max(fileb$y)-min(fileb$y))),], aes(x=x, y=y), col="blue")+
#  geom_line(data=filer, aes(x=x, y=y), col="red")+
#  geom_point(data=filer[findpeaks(filer$y, span=10),][filer[findpeaks(filer$y, span=10),]$y-min(filer$y)>(0.4*(max(filer$y)-min(filer$y))),], aes(x=x, y=y), col="red")+
#  theme_bw())
#dev.off()
#}
#}

png(file=("../Length_blue_peaks.png"),height=2000,width=3000,res=600) 
ggplot(data=cellsTOT[cellsTOT$BluePeak==1,], aes(x=Length))+
  geom_histogram(binwidth=0.2, fill="royalblue4")+
  geom_histogram(data=cellsTOT[cellsTOT$BluePeak==2,],binwidth=0.2, alpha=0.5, fill="orange")
dev.off()

png(file=("../Length_red_peaks.png"),height=2000,width=3000,res=600) 
ggplot(data=cells[cells$RedPeak==0,], aes(x=Length))+
  geom_histogram(binwidth=0.2, fill="royalblue4")+
  geom_histogram(data=cells[cells$RedPeak==1,],binwidth=0.2, alpha=0.5, fill="orange")+
  geom_histogram(data=cells[cells$RedPeak==2,],binwidth=0.2, alpha=0.5, fill="red")
dev.off()

length(cells[cells$RedPeak==0,]$Index)
length(cells[cells$BluePeak==1,]$Index)


ggplot(data=filer, aes(x=x, y=y))+
  geom_line()+
  geom_point(data=filer[findpeaks(filer$y, span=10),][filer[findpeaks(filer$y, span=10),]$y-min(filer$y)>(0.4*(max(filer$y)-min(filer$y))),], aes(x=x, y=y), col="red")


ggplot(data=cellsTOT, aes(x=BluePeak, y=..density..))+
  geom_histogram(binwidth=1, fill="royalblue4")+
  geom_histogram(data=cellsTOT, aes(x=RedPeak, y=..density..), binwidth=1, alpha=0.5, fill="orange")



ggplot(data=cells, aes(x=RedPeak))+
  geom_histogram(binwidth = 1)+
  xlim(-1,5)

ggplot(data=cells, aes(x=BluePeak))+
  geom_histogram(binwidth = 1)+
  xlim(-1,5)

plot(filer$x, filer$y, type = "l", col = "gray")
pk.pos <-findpeaks(filer$y , span=5)
abline(v = filer$x[pk.pos], col = 4)
pks <- fitpeaks(filer$y, pk.pos)
apply(pks, 1,
      function(pkmodel) {
        lines(filer$x,
              dnorm(1:length(filer$x), pkmodel["rt"], pkmodel["sd"]) *
                pkmodel["area"],
              col = 2)
        invisible()
      })

getAllPeaks()