
# read in data folder holding Agilent text files and the annotation file;
# also set results folder

results.name<-"basic"
folder='test/'

if (exists("folder")) {

 if (folder!="") { 

  temp=readline(paste("Data folder ('return' for ",folder," ): ",sep=""))

  if (temp!="") { folder=temp }

 } else {  folder=readline("Data folder: ") }

} else { folder=readline("Data folder: ")}



if (exists("results.name")) {

 if (results.name!="") {

  temp=readline(paste("Data norm results name ('return' for ",results.name," ): ",sep=""))

  if (temp!="") { results.name=temp }

 } else { results.name=readline("Data norm results name: ") }

} else { results.name=readline("Data norm results name: ") }



# v3: uses AgilentNorm_v3 which attempts to find euploid regions for lowess
#     and Segment (which now accounts for replicate probes)
# v4: puts gc-norm first and uses euploid regions only for R~I loess;
#     prints out pseudoimage; also better normalization plots

resultFolder<- paste(folder,results.name,"/",sep="")
if (file.access(resultFolder)!=0) { print("unable to find result directoty"); stop() }

if (file.exists(paste(resultFolder,"norm.out.D",sep=""))) {
	temp=readline("Redo analysis (R) or append new samples (A):")
	if (temp=="R" | temp=="r") { append=NA }
	if (temp=="A" | temp=="a") { append=resultFolder}
	if (!exists("append")) { stop() }
} else { append=NA }


source("agilentNorm_v4.R")

norm.out<-agilentNorm(folder,append=append,write.gpr=FALSE)


qa<-norm.out$qa
uml<-norm.out$pre.gcnorm


write.table(norm.out$results,file=paste(resultFolder,"normalized.xls",sep=""),sep="\t",row.names=FALSE,na="")
write.table(uml,file=paste(resultFolder,"pre.gcnormalized.xls",sep=""),sep="\t",row.names=FALSE,na="")


temp=t(cbind(rownames(qa),qa))
write.table(temp,file=paste(resultFolder,"QA.preSeg.xls",sep=""),sep="\t",col.names=FALSE,na="")

save(norm.out,file=paste(resultFolder,"norm.out.D",sep=""))


source("segment.R")

segment.out<-segment.normalized(norm.out,resultFolder,append=append)
save(segment.out,file=paste(resultFolder,"segment.out.D",sep=""))

temp=segmentQA(segment.out$segment.results$output)
qa.s=rbind(qa,temp$data)

temp=t(cbind(rownames(qa.s),qa.s))
write.table(temp,file=paste(resultFolder,"QA.postSeg.xls",sep=""),sep="\t",col.names=FALSE,na="")

map=norm.out$map; chr=map$Ch; pos=map$Pos

rml=segment.out$rml
mml=segment.out$mml
sml=segment.out$sml
smlgf=segment.out$smlgf

profiles.html(resultFolder)
chr.html(resultFolder)



