############################################################################################################################################
#The aim of this script is to automatically curate manually drawn cell profiles.
#This script define the beginning and end of a cell based on the membrane staining fluorescence profile.
#The outer limits of a cells are here defined as the maximum fluorescence intensity within the first and last 20% of the profile data.
############################################################################################################################################

#Set the working directory within one of the "Profiles" folder containing the "Blue*.txt", "Gray*.txt" and "Red*.txt" files (See output from "ProfilingCells_EPI.ijm"

fold<-"Profiles*"  #Input folders
folder<-paste("../",Sys.glob(fold), sep="")
fold1<-"Membranes*"  #Output folders
outfolder<-paste("../",Sys.glob(fold1),"/", sep="")

chan<-"Red*"
chan2<-"Blue*"

for (f in 1:length(folder)){
  setwd(folder[f])
  RED<-Sys.glob(chan)
  BLUE<-Sys.glob(chan2)
  
  for(k in 1:length(RED)){
    filer<-read.table(RED[k], sep="\t",dec=".", stringsAsFactors=FALSE,header=TRUE)
    fileb<-read.table(BLUE[k], sep="\t",dec=".", stringsAsFactors=FALSE,header=TRUE)
    start<-filer[filer$x<=(0.2*max(filer$x)) & filer$y==max(filer[filer$x<=(0.2*max(filer$x)),]$y),]$x
    end<-filer[filer$x>=(0.8*max(filer$x)) & filer$y==max(filer[filer$x>=(0.8*max(filer$x)),]$y),]$x
    filer_new<-filer[filer$x>=start & filer$x<=end,]
    filer_new$x<-filer_new$x-start
    fileb_new<-fileb[fileb$x>=start & fileb$x<=end,]
    fileb_new$x<-fileb_new$x-start
    write.table(filer_new, file=paste(outfolder[f],RED[k], sep=""),row.names = FALSE, quote=FALSE, sep="\t")
    write.table(fileb_new, file=paste(outfolder[f],BLUE[k], sep=""),row.names = FALSE, quote=FALSE, sep="\t")
    }}

