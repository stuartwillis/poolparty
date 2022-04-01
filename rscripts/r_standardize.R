#!bin/R

#Dependencies
#install.packages("tidyr")
#install.packages("matrixStats")

  suppressMessages(require(data.table))
	suppressMessages(require(matrixStats))
	suppressMessages(require(tidyr))
	

	#getwd()

#Get info for bash script
	args <- commandArgs()
	name <- args[6]
		outdir <- args[7]
	outname<- gsub( ".sync", "_norm.sync", name )
		outname2 <- gsub( ".sync", "_poly_sites.txt", name )
		alertname <- gsub( ".sync", "", name )
		alertname <- sub('.*/', '', alertname )
  #Reading file
  gpops <- fread(file=name,stringsAsFactors=FALSE, showProgress=FALSE)
  #Break into genotypes and heading
  counts= gpops[, 4:(ncol(gpops))]
  gpops = gpops[, 1:3]  
  #counts <- as.matrix(counts)

  #Get right dir
#  activedir = getwd()
  tfile= paste0(outdir,"/","tempRf.tem")
  
  #Write temp file of allele counts
  write.table(counts, file= tfile, sep=':',
                          row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  #Read counts back in separate columns
  counts <- fread(tfile,stringsAsFactors=FALSE, sep=':', showProgress=FALSE)

  #Remove temporary file
  invisible(file.remove(tfile))
 
  #Create rownames for freq file
#  gpops3 <- unite(gpops, V1, c(V1,V2), sep="_", remove=FALSE)
  
  #Determine number of individuals 
  npops = ncol(counts)/6
  
  #summary info from file
  totcol = ncol(counts)
  thetseq = seq(6,totcol,6)
  thebseq = seq(1,totcol,6)
#  L = nrow(counts)

   #Break up Genotypes for each individual into separate lists
#    Var <- list()
#      for (i in (1:npops)) {
#        Var[[i]] <-counts[, thebseq[i]:thetseq[i]]
#      }
  CalcZ <- list()
  for (i in (1:npops)) {
    CalcZ[[i]] <-counts[, thebseq[i]:thetseq[i]]
  }

    rm(counts)
     invisible(gc())
    alert1= paste0("R ALERT: Determining .sync stats for ", alertname)
    print (alert1)
    #Make Matrix for quick calculations, and get stats 
    MVar <- list()
    MCON <- list()# read proportion of allele with fewer reads, homozygous sites are thus "1", with unscored sites recorded as "2"
    MSUM <- list()# summed read counts
    MCNT <- list()# number of alleles scored
    MMIN <- list()# read count of the allele with fewer reads

#   for (i in (1:npops)) {
#      MVar[[i]] <- as.matrix(Var[[i]])
#      is.na( MVar[[i]])<- MVar[[i]]==0
#      #MVar[[i]][3,]
#      MCON[[i]] <- rowMins(MVar[[i]], na.rm =TRUE) / rowSums(MVar[[i]], na.rm=TRUE)
#      MCON[[i]][!is.finite(MCON[[i]])] <- 2
#      MSUM[[i]] <- rowSums(MVar[[i]], na.rm =TRUE)
#      MCNT[[i]] <- rowSums(MVar[[i]] > 0, na.rm =TRUE)
#      MMIN[[i]] <- rowMins(MVar[[i]], na.rm =TRUE)
#      gc()
#    }
#   rm(MVar)
    for (i in (1:npops)) {
      MVar[[i]] <- as.matrix(CalcZ[[i]])
      is.na( MVar[[i]])<- MVar[[i]]==0
      #MVar[[i]][3,]
      MCON[[i]] <- rowMins(MVar[[i]], na.rm =TRUE) / rowSums(MVar[[i]], na.rm=TRUE)
      MCON[[i]][!is.finite(MCON[[i]])] <- 2
      MSUM[[i]] <- rowSums(MVar[[i]], na.rm =TRUE)
      MCNT[[i]] <- rowSums(MVar[[i]] > 0, na.rm =TRUE)
      MMIN[[i]] <- rowMins(MVar[[i]], na.rm =TRUE)
  #    gc()
    }
    rm(MVar)
    invisible(gc())
    #Apply Base Column Names
    #Var<-lapply(Var, function(x) {
      CalcZ<-lapply(CalcZ, function(x) {
       colnames(x) <- paste(c('A','T','C','G','N', 'del') )
      return(x) })

print("R ALERT: Performing corrections on .sync file")
  #Replace based on criteria
    #1 If an individual has more than 3 alleles at a position, set all alleles to 0 weight (expected error)
    #2 If single allele contribution (MCON) < 10.1 % make the minor allle 0 (expected error)
    #3 If number of unique alleles (MCNT) =1, then make 2 weight for single allele 
    #4 If MCNT =2, make 1 weight for each of the heterozygous alleles
    #5 If coverage (MSUM) =1, Make allele weight 1 
   
  #Percent variation allowed in alleles
  perV=0.101 #probably should be 0.101 to be consistent with above
#  MCON # read proportion of allele with fewer reads, homozygous sites are thus "1", with unscored sites recorded as "2"
#  MSUM # summed read counts
#  MCNT # number of alleles scored
#  MMIN # read count of the allele with fewer reads
  CalcZ <-lapply(seq_along(CalcZ),function(i) {rapply(CalcZ[[i]],function(x) ifelse(x>0 & MCNT[[i]] > 2,0,x), how = "replace")})#replace counts for 3+ allele sites with 0
  CalcZ <-lapply(seq_along(CalcZ),function(i) {rapply(CalcZ[[i]],function(x) ifelse( x == MMIN[[i]] & MCON[[i]] < perV,0,x), how = "replace")})#if allele balance is under threshold (perV), give allele with fewer reads count of 0
  CalcZ <-lapply(seq_along(CalcZ),function(i) {rapply(CalcZ[[i]],function(x) ifelse(x>0 & MCNT[[i]] == 1 ,2,x), how = "replace")})#where only a single allele is scored, record count of 2
  CalcZ <-lapply(seq_along(CalcZ),function(i) {rapply(CalcZ[[i]],function(x) ifelse(x>0 & MCNT[[i]] == 2 & MCON[[i]] < perV,2,x), how = "replace")})#where two alleles are scored and BELOW allele balance threshold, give allele with more reads count of 2
  #CalcZ <-lapply(seq_along(CalcZ),function(i) {rapply(CalcZ[[i]],function(x) ifelse(x>0 & MCNT[[i]] == 2 & MCON[[i]] < perV,1,x), how = "replace")})#where two alleles are scored and BELOW allele balance threshold, give allele with more reads count of 1; SCW says "YES!"
  CalcZ <-lapply(seq_along(CalcZ),function(i) {rapply(CalcZ[[i]],function(x) ifelse(x>0 & MCNT[[i]] == 2 & MCON[[i]] > perV,1,x), how = "replace")}) #where two alleles are scored and ABOVE allele balance threshold, give each allele count of 1
  CalcZ <-lapply(seq_along(CalcZ),function(i) {rapply(CalcZ[[i]],function(x) ifelse( x >0  & MSUM[[i]] == 1,1,x), how = "replace")}) #if only a single read is recorded, give that allele a count of 1
  CalcZ <- lapply(seq_along(CalcZ),function(i) {as.data.frame(CalcZ[[i]])})

 #PolY <-  Reduce("+", MCNT) #add a list of vectors or matrices together, row/column-wise
 ##problem is that if most individuals are homozygotes, i.e. scored at less than 2 alleles, even if 1+ individuals are scored at 3 alleles, the sum may still be less that 2*no. individuals, so not be identified as over that threshold for paralogous
 ##solution: add a step that converts all 1 counts to 2
 MCNTboost<-lapply(MCNT,function(x) ifelse(x==1,2,x))
 PolY <-  Reduce("+", MCNTboost)
  PolY <- as.data.frame(cbind(gpops, PolY))
  rm(gpops)
  rm(MCNTboost)
  rm(MCNT)
  rm(MCON)
  rm(MMIN)
  rm(MSUM)
  invisible(gc())
    FaiL <- npops * 2
    #sum(PolY > 4)
      #PolY<- PolY[ PolY > 4] #reduced data fram of sites with individuals scored at 3+ alleles
      PolY<- PolY[PolY$PolY>FaiL,] #reduced data fram of sites with individuals scored at 3+ alleles; SCW fix
      #nrow(PolY)
      PolY <- PolY[,1:2]
        alert2<- paste0("R ALERT: ", nrow(PolY), " SNPs in ", alertname, " had 3 or more alleles. These will be blacklisted.")
        print (alert2)
  
  #Combine all of the standardized individual .syncs
  CalcZ <- as.data.frame(Reduce("+", CalcZ))  #add a list of vectors or matrices together, row/column-wise

#  CalcZ <- as.data.frame(CalcZ)
    colnames(CalcZ) <- c('A','B','C','D','E','F')
    CalcZ <- unite(CalcZ, V1, c(A,B,C,D,E,F), sep=":", remove=TRUE)#TRUE?
 
       #write out CalcZ$V1
  print("R ALERT: R is writing files to disk")
  write.table(CalcZ$V1, file=outname, sep='\t',row.names= FALSE, quote=FALSE, col.names= FALSE)
  write.table(PolY, file=outname2, sep='\t',row.names= FALSE, quote=FALSE, col.names= FALSE)
  rm(CalcZ)
  rm(PolY)

