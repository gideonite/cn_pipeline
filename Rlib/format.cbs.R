#Set up input files
#Input is
#cbs.output - the standard output from running CBS
#segments.p.output - the output from running segments.p
#outputfile - name of output file from running format.cbs
#header - whether the first line in the output should be the column names

format.cbs <- function(cbs.output,segments.p.output,outputfile,header=TRUE)
  {
    unique.names <- unique(segments.p.output[,1])
    n <- length(unique.names)

    chroms.vector <- cbs.output$data[,1]
    positions.vector <- cbs.output$data[,2]
    
    for(i in 1:n)  
      {
        new.output <- segments.p.output[segments.p.output[,1]==unique.names[i],]
        p <- nrow(new.output)

        new.data <- cbs.output$data[,i+2]
        previous.chrom <- new.output[1,2]
        previous.start <- new.output[1,3]
        previous.end <- new.output[1,4]
        previous.markers <- new.output[1,5]
        which.na <- which(is.na(new.data))
        chroms.na <- chroms.vector[which.na]
        positions.na <- positions.vector[which.na]
        count.probes <- rep(0,p)
        sum.missing <- sum(chroms.na==previous.chrom & positions.na>=previous.start & positions.na<=previous.end)
        count.probes[1] <- previous.markers+sum.missing
#Add count of informative markers
        for(j in 2:p)
          {
            new.chrom <- new.output[j,2]
            new.start <- new.output[j,3]
            new.end <- new.output[j,4]
            new.markers <- new.output[j,5]        
            sum.missing <- sum(chroms.na==new.chrom & positions.na>=new.start & positions.na<=new.end)
            count.probes[j] <- new.markers+sum.missing
          }
        new.output <-  cbind(new.output[,1:4],count.probes,new.output[,5:6],new.output[,8:10])
#Look in gaps between chromosomes
        new.matrix <- NULL
        previous.chrom <- new.output[1,2]
        previous.end <- new.output[1,3]
        previous.end <- new.output[1,4]
        which.missing <- which(chroms.na==previous.chrom & positions.na<=previous.start)
        if(length(which.missing)>0)
          {
            new.row <- c(unique.names[i],previous.chrom,min(positions.na[which.missing]),max(positions.na[which.missing]),length(which.missing),0,rep(NA,4))
            if(length(new.matrix)==0) new.matrix <- new.row
            else new.matrix <- rbind(new.matrix,new.row)
          }    
        for(j in 2:p)      
          {
            new.chrom <- new.output[j,2]
            new.start <- new.output[j,3]
            new.end <- new.output[j,4]
            if (new.chrom==previous.chrom)
              {
                which.missing <- which(chroms.na==new.chrom & positions.na<=new.start & positions.na>=previous.end)
                if(length(which.missing)>0)
                  {
                    new.row <- c(unique.names[i],new.chrom,min(positions.na[which.missing]),max(positions.na[which.missing]),length(which.missing),0,rep(NA,4))
                    if(length(new.matrix)==0) new.matrix <- new.row
                    else new.matrix <- rbind(new.matrix,new.row)
                  }
                previous.chrom <- new.chrom
                previous.end <- new.end
              }
            else if(new.chrom!=previous.chrom)
              {
                which.missing <- which(chroms.na==previous.chrom & positions.na>=previous.end)
                if(length(which.missing)>0)
                  {
                    new.row <- c(unique.names[i],previous.chrom,min(positions.na[which.missing]),max(positions.na[which.missing]),length(which.missing),0,rep(NA,4))
                    if(length(new.matrix)==0) new.matrix <- new.row
                    else new.matrix <- rbind(new.matrix,new.row)
                  }
                which.missing <- which(chroms.na==new.chrom & positions.na<=new.start)
                if(length(which.missing)>0)
                  {
                    new.row <- c(unique.names[i],new.chrom,min(positions.na[which.missing]),max(positions.na[which.missing]),length(which.missing),0,rep(NA,4))
                    if(length(new.matrix)==0) new.matrix <- new.row
                    else new.matrix <- rbind(new.matrix,new.row)
                  }
                previous.chrom <- new.chrom
                previous.end <- new.end                        
              } 
          }
        which.missing <- which(chroms.na==new.chrom & positions.na>=new.end)
        if(length(which.missing)>0)
          {
            new.row <- c(unique.names[i],new.chrom,min(positions.na[which.missing]),max(positions.na[which.missing]),length(which.missing),0,rep(NA,4))
            if(length(new.matrix)==0) new.matrix <- new.row
            else new.matrix <- rbind(new.matrix,new.row)
          }
# Now merge the two
        if(length(new.matrix)>0)
          {
            new.matrix <- matrix(new.matrix,ncol=10)
            new.matrix.chroms <- new.matrix[,2]
#            new.matrix.chroms <- as.numeric(new.matrix[,2])            
            new.matrix.ends <- as.numeric(new.matrix[,4])
            for(j in 1:nrow(new.matrix))
              {
                new.chrom <- new.matrix.chroms[j]
                new.end <- new.matrix.ends[j]
                if(new.chrom==new.output[1,2] & new.end<=as.numeric(new.output[1,3]))
#                if(new.chrom==as.numeric(new.output[1,2]) & new.end<=as.numeric(new.output[1,3]))                  
                  {
                    new.output <- rbind(new.matrix[j,],new.output)
                  }
                else if (new.chrom==new.output[nrow(new.output),2] & new.end>=as.numeric(new.output[nrow(new.output),3]))
#                else if (new.chrom==as.numeric(new.output[nrow(new.output),2]) & new.end>=as.numeric(new.output[nrow(new.output),3]))                  
                  {
                    new.output <- rbind(new.output,new.matrix[j,])
                  }
                else
                  {
                    which.chrom <- which(new.output[,2]==new.chrom)
                    no.match <- TRUE
                    k <- 1
                    while(no.match & k<=length(which.chrom))
                      {
                        new.position <- which.chrom[k]
                        matrix.start <- as.numeric(new.output[new.position,3])
                        if(new.end<=matrix.start)
                          {
                            new.output <- rbind(new.output[1:(new.position-1),],new.matrix[j,],new.output[new.position:nrow(new.output),])
                            no.match <- FALSE
                          }
                        k <- k+1
                      }
                    if(no.match)
                      {
                        new.output <- rbind(new.output[1:new.position,],new.matrix[j,],new.output[(new.position+1):nrow(new.output),])
                      }
                  }
              }
          }
        left.threecolumns <- c(NA,NA,NA)
        for(j in 2:nrow(new.output))
          {
            if(new.output[j-1,2]==new.output[j,2])
              {
                left.threecolumns <- rbind(left.threecolumns,as.numeric(new.output[j-1,8:10]))
              }
            else
              {
                left.threecolumns <- rbind(left.threecolumns,c(NA,NA,NA))
              }
          }
        new.output <- cbind(new.output[,1:7],left.threecolumns,new.output[,8:10])
        if(i==1) output <- new.output
        else output <- rbind(output,new.output)
      }
    names(output) <- c("sample","chrom","loc.start","loc.end","num.mark","num.informative","seg.mean","pval","l.lcl","l.ucl","r.pval","r.lcl","r.ucl")

    write.table(output,outputfile,sep="\t",row.names=FALSE,col.names=header,quote=FALSE)

  }

