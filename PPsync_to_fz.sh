#!/bin/bash

#get name of sync and MAF

if [ -z "$2" ]
then
echo "Not enough arguments."
echo "Correct usage: SCRIPT.sh FILE.sync MAF"
echo "Output will be called FILE.fz."
exit 1
#else
fi

WD=( ` pwd `)

echo $1
echo $2
echo $WD


echo "#!bin/R" > fix_frequency.R
echo "suppressMessages(require(tidyr))" >> fix_frequency.R
echo "suppressMessages(require(data.table))" >> fix_frequency.R
echo "suppressMessages(require(stringr))" >> fix_frequency.R
echo "suppressMessages(require(matrixStats))" >> fix_frequency.R
echo "args <- commandArgs()" >> fix_frequency.R
echo "name <- args[6]" >> fix_frequency.R
echo "outdir <- args[7]" >> fix_frequency.R
echo "outdirF <- args[8]" >> fix_frequency.R
echo "MAF<- as.numeric(args[9])" >> fix_frequency.R
echo "alertname <- sub('.*/', '', name )" >> fix_frequency.R
echo "alertname2 <- gsub( \".sync\", \"\", alertname)" >> fix_frequency.R
echo "alert1<-paste0(\"R ALERT: Reading file \", alertname, \" into R\")" >> fix_frequency.R
echo "print(alert1)" >> fix_frequency.R
echo "gpops <- fread(file=name,stringsAsFactors=FALSE, showProgress=FALSE)" >> fix_frequency.R
echo "counts= gpops[, 4:(ncol(gpops))]" >> fix_frequency.R
echo "gpops = gpops[, 1:3] " >> fix_frequency.R
echo "colnames(gpops) <- c(\"Chr\", \"Pos\", \"Ref\")" >> fix_frequency.R
echo "counts <- as.matrix(counts)" >> fix_frequency.R
echo "activedir = getwd()" >> fix_frequency.R
echo "tfile= paste0(outdir,\"/\",\"tempRf.tem\")" >> fix_frequency.R
echo "print(\"R ALERT: Formatting file and calculating summary stats\")" >> fix_frequency.R
echo "npops = ncol(counts)" >> fix_frequency.R
echo "L = nrow(counts)" >> fix_frequency.R
echo "Var<-apply(counts,2,function(x)data.frame(matrix(unlist(lapply(x,function(y)as.numeric(unlist(strsplit(y,\":\"))))),ncol=6,byrow=T),stringsAsFactors = FALSE))" >> fix_frequency.R
echo "names(Var)<-NULL" >> fix_frequency.R
echo "rm(counts) ; invisible(gc())" >> fix_frequency.R
echo " if (is.na(args[10])) {" >> fix_frequency.R
echo "      print(\"R ALERT: Calculating depth of coverage at each position\")" >> fix_frequency.R
echo "      Dcov <- data.frame(matrix(, nrow = L, ncol = npops))" >> fix_frequency.R
echo "      for (i in (1:length(Var))) {" >> fix_frequency.R
echo "         Dcov[[i]] <- rowSums(Var[[i]])" >> fix_frequency.R
echo "         gc()" >> fix_frequency.R
echo "      }" >> fix_frequency.R
echo "      indel.N = list()" >> fix_frequency.R
echo "      for (i in (1: npops)) {" >> fix_frequency.R
echo "           indel.N[[i]] <- as.matrix(rowSums(Var[[i]][,5:6]))" >> fix_frequency.R
echo "      }" >> fix_frequency.R
echo "      indel.N <- as.data.frame(indel.N)" >> fix_frequency.R
echo "      colnames(indel.N) <- c(1:npops)" >> fix_frequency.R
echo "      colnames(Dcov) <- c(1:npops)" >> fix_frequency.R
echo "      Dcov\$Total <- rowSums(Dcov)" >> fix_frequency.R
echo "      Dcov <- cbind(gpops[,1:2], Dcov, indel.N)" >> fix_frequency.R
echo "      outname7 <- paste0(outdirF, alertname2, \"_coverage.txt\")" >> fix_frequency.R
echo "      write.table(Dcov, file=outname7," >> fix_frequency.R
echo "             row.names=FALSE, col.names=TRUE, quote=FALSE)  " >> fix_frequency.R
echo "      rm(Dcov) ; invisible(gc())" >> fix_frequency.R
echo "      rm(indel.N) ; invisible(gc())" >> fix_frequency.R
echo "	}" >> fix_frequency.R
echo "alertname <- sub('.*/', '', alertname )" >> fix_frequency.R
echo "Var<-lapply(Var, function(x) {" >> fix_frequency.R
echo "   colnames(x) <- paste(c('A','T','C','G', 'N', 'del') )" >> fix_frequency.R
echo "return(x) })" >> fix_frequency.R
echo " VarN<-lapply(Var, function(x) {" >> fix_frequency.R
echo "  colnames(x) <- paste(c('1','2','3','4', '5', '6') )" >> fix_frequency.R
echo "  return(x) })" >> fix_frequency.R
echo "print(\"R ALERT: Noting potential paralogs (>3 alleles per position)\")" >> fix_frequency.R
echo "  MCNT <- list()" >> fix_frequency.R
echo "  MVar <- list()" >> fix_frequency.R
echo "  for (i in (1: npops)) {" >> fix_frequency.R
echo "   MVar[[i]] <- as.matrix(Var[[i]][,1:4])" >> fix_frequency.R
echo "    MCNT[[i]] <- rowSums(MVar[[i]] >0 )" >> fix_frequency.R
echo "    MCNT[[i]] <- as.data.frame(MCNT[[i]])" >> fix_frequency.R
echo "    MCNT[[i]] <-rapply(MCNT[[i]], function(x) ifelse(x > 2,1,0), how = \"replace\")" >> fix_frequency.R
echo "    gc()" >> fix_frequency.R
echo "  }" >> fix_frequency.R
echo "PolY <-  Reduce(\"+\", MCNT)" >> fix_frequency.R
echo "rm(MCNT) ; invisible(gc())" >> fix_frequency.R
echo "rm(MVar) ; invisible(gc())" >> fix_frequency.R
echo "PolY <- cbind(gpops[,1:2], PolY)" >> fix_frequency.R
echo "    PolY1<- PolY[,1:2][ PolY[[3]] > 0]" >> fix_frequency.R
echo "    PolY2<- PolY[,1:2][PolY[[3]] >= npops*0.5]" >> fix_frequency.R
echo "    PolY3<- PolY[,1:2][PolY[[3]] >= npops]" >> fix_frequency.R
echo "outname4 <- paste0(outdirF, alertname2, \"_poly_one.txt\")" >> fix_frequency.R
echo "outname5 <- paste0(outdirF, alertname2, \"_poly_half.txt\")" >> fix_frequency.R
echo "outname6 <- paste0(outdirF, alertname2, \"_poly_all.txt\")" >> fix_frequency.R
echo "  write.table(PolY1, file=outname4," >> fix_frequency.R
echo "              row.names=FALSE, col.names=FALSE, quote=FALSE)" >> fix_frequency.R
echo "  write.table(PolY2, file=outname5," >> fix_frequency.R
echo "              row.names=FALSE, col.names=FALSE, quote=FALSE) " >> fix_frequency.R
echo "  write.table(PolY3, file=outname6," >> fix_frequency.R
echo "              row.names=FALSE, col.names=FALSE, quote=FALSE)" >> fix_frequency.R
echo "  rm(PolY) ; invisible(gc())" >> fix_frequency.R
echo "  rm(PolY1) ; invisible(gc())" >> fix_frequency.R
echo "  rm(PolY2) ; invisible(gc())" >> fix_frequency.R
echo "  rm(PolY3) ; invisible(gc())" >> fix_frequency.R
echo "print(\"R ALERT: Finished writing potential paralogs (>3 alleles per position)\")" >> fix_frequency.R
echo "MaxC <- list()" >> fix_frequency.R
echo "  for (i in (1: npops)) {" >> fix_frequency.R
echo "     MaxC[[i]] <- colnames(Var[[i]])[max.col(Var[[i]],ties.method=\"first\")]" >> fix_frequency.R
echo "     gc()" >> fix_frequency.R
echo "  }" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 1 Complete\")" >> fix_frequency.R
echo "MaxN <- list()" >> fix_frequency.R
echo "  for (i in (1: npops)) {" >> fix_frequency.R
echo "     MaxN[[i]] <- colnames(VarN[[i]])[max.col(VarN[[i]],ties.method=\"first\")]" >> fix_frequency.R
echo "     gc()" >> fix_frequency.R
echo "  }" >> fix_frequency.R
echo "rm(VarN) ; invisible(gc())" >> fix_frequency.R
echo "  print(\"R ALERT: Tabulating alleles: Step 2 Complete\")" >> fix_frequency.R
echo "  MaxC <- as.data.frame(MaxC)" >> fix_frequency.R
echo "  colnames(MaxC) <- c(1:npops)" >> fix_frequency.R
echo "  MaxC <- as.matrix(MaxC)" >> fix_frequency.R
echo "  print(\"R ALERT: Tabulating alleles: Step 3a Complete\")" >> fix_frequency.R
echo "  MaxN <- as.data.frame(MaxN)" >> fix_frequency.R
echo "  colnames(MaxN) <- c(1:npops)" >> fix_frequency.R
echo "  MaxN <- as.matrix(MaxN)" >> fix_frequency.R
echo "  MaxN <- data.frame(apply(MaxN, 2, function(x) as.numeric(as.character(x))))" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 3b Complete\")" >> fix_frequency.R
echo " Sumz<- lapply(Var, function(x) rowSums(x))" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 4 Complete\")" >> fix_frequency.R
echo "Var<-lapply(Var,function(x)as.matrix(x))" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 5 Complete\")" >> fix_frequency.R
echo "print(\"R ALERT: Calculating Allele Frequencies\")" >> fix_frequency.R
echo "Maxz<- lapply(Var, function(x) rowMaxs(x))" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 6 Complete\")" >> fix_frequency.R
echo "Freqz <- mapply(\"/\",Maxz,Sumz,SIMPLIFY = FALSE)  " >> fix_frequency.R
echo "rm(Maxz) ; invisible(gc())" >> fix_frequency.R
echo "rm(Sumz) ; invisible(gc())" >> fix_frequency.R
echo "rm(Var) ; invisible(gc())" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 7 Complete\")" >> fix_frequency.R
echo "Freqz <- as.data.frame(Freqz)" >> fix_frequency.R
echo "colnames(Freqz) <- c(1:npops)" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 8 Complete\")" >> fix_frequency.R
echo "Freqz <- as.matrix(Freqz)" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 9 Complete\")" >> fix_frequency.R
echo "modefunc <- function(x){" >> fix_frequency.R
echo "  tabresult <- tabulate(x)" >> fix_frequency.R
echo "  themode <- which(tabresult == max(tabresult))" >> fix_frequency.R
echo "  if(sum(tabresult == max(tabresult))>1) themode <- themode[[1]]" >> fix_frequency.R
echo "  return(themode)" >> fix_frequency.R
echo "}" >> fix_frequency.R
echo "MaxN <- apply(MaxN, 1, modefunc)" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 11 Complete\")" >> fix_frequency.R
echo "MaxN <- as.data.frame(chartr(\"123456\", \"ATCGND\", MaxN))" >> fix_frequency.R
echo "ComPos <- as.matrix(MaxN)" >> fix_frequency.R
echo "rm(MaxN) ; invisible(gc())" >> fix_frequency.R
echo "colnames(ComPos) <- \"A1\"" >> fix_frequency.R
echo "idx <- (!apply(MaxC, 2, function(x) x == ComPos))" >> fix_frequency.R
echo "rm(MaxC) ; invisible(gc())" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 12 Complete\")" >> fix_frequency.R
echo "Freqz[idx] <- 1- Freqz[idx]" >> fix_frequency.R
echo "rm(idx) ; invisible(gc())" >> fix_frequency.R
echo "Freqz<- round(Freqz, 3)" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 13 Complete\")    " >> fix_frequency.R
echo "print(\"R ALERT: Determining Positions that fail MAF\")" >> fix_frequency.R
echo "mAx <- rowMaxs(Freqz, na.rm = TRUE)" >> fix_frequency.R
echo "mIn<- rowMins(Freqz, na.rm = TRUE)" >> fix_frequency.R
echo "rAnge <- as.vector(mAx - mIn)" >> fix_frequency.R
echo "rAnge <- as.matrix(rAnge)" >> fix_frequency.R
echo "mIn <- as.matrix(mIn)" >> fix_frequency.R
echo "mAx <- as.matrix(mAx)" >> fix_frequency.R
echo "rownames(rAnge) <- c()" >> fix_frequency.R
echo "colnames(rAnge) <- \"MAFT\"" >> fix_frequency.R
echo "rownames(mAx) <- c()" >> fix_frequency.R
echo "colnames(mAx) <- \"MAX\"" >> fix_frequency.R
echo "rownames(mIn) <- c()" >> fix_frequency.R
echo "colnames(mIn) <- \"MIN\"" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 14 Complete\")" >> fix_frequency.R
echo "MAFFAIL <- cbind (gpops[,1:2],Freqz, rAnge, mAx, mIn)" >> fix_frequency.R
echo "rm(mAx) ; invisible(gc())" >> fix_frequency.R
echo "rm(mIn) ; invisible(gc())" >> fix_frequency.R
echo "rm(rAnge) ; invisible(gc())" >> fix_frequency.R
echo "MAFFAIL <- as.data.frame(MAFFAIL)" >> fix_frequency.R
echo "MAFFAIL <- subset(MAFFAIL[,1:2], ((MAFFAIL\$MAX >= (1-MAF) | MAFFAIL\$MIN <= MAF) & (MAFFAIL\$MAFT <= MAF)) )" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 15 Complete\")" >> fix_frequency.R
echo "proP<- round((nrow(MAFFAIL) / nrow(Freqz) * 100),2)" >> fix_frequency.R
echo "paste0(\"R ALERT: \", nrow(MAFFAIL), \" SNPS (\", proP, \"%) do not pass additional population MAF threshold of \", MAF, \" and have been noted\")" >> fix_frequency.R
echo "outname3 <- paste0(outdirF, alertname2, \"_MAF_fail.txt\")" >> fix_frequency.R
echo "for(r in 1:nrow(MAFFAIL)){" >> fix_frequency.R
echo "  write.table(MAFFAIL[r,], file=outname3,row.names=FALSE, col.names=FALSE, quote=FALSE,append=TRUE)" >> fix_frequency.R
echo "}" >> fix_frequency.R
echo "rm(MAFFAIL) ; invisible(gc())" >> fix_frequency.R
echo "Freqz <- cbind(gpops[,1:3],ComPos,Freqz)" >> fix_frequency.R
echo "print(\"R ALERT: Tabulating alleles: Step 16 Complete\")" >> fix_frequency.R
echo "rm(ComPos) ; invisible(gc())" >> fix_frequency.R
echo "rm(gpops) ; invisible(gc())" >> fix_frequency.R
echo " outname1 <- paste0(outdir, alertname2, \".fz\")" >> fix_frequency.R
echo "print(\"R ALERT: Writing output files\")" >> fix_frequency.R
echo " write.table(Freqz[1,], file=outname1,row.names=FALSE, col.names=TRUE, quote=FALSE)" >> fix_frequency.R
echo "for(r in 2:nrow(Freqz)){" >> fix_frequency.R
echo "  write.table(Freqz[r,], file=outname1,row.names=FALSE, col.names=FALSE, quote=FALSE,append=TRUE)" >> fix_frequency.R
echo "}" >> fix_frequency.R
echo "rm(Freqz) ; invisible(gc())" >> fix_frequency.R
echo "print(\"R ALERT: Frequency.R completed correctly\")" >> fix_frequency.R


rin=$WD/$1
rout=$WD/
rout2=$WD/filters/
MAF=$2
echo $rin
echo $rout
echo $rout2
echo $MAF
Rscript $WD/fix_frequency.R $rin $rout $rout2 $MAF

rm $WD/fix_frequency.R &> /dev/null
