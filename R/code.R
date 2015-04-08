
############################################
# NUCLEOTIDE LEVEL                         #
# 1) Dominant (function: nDOM)             #
# 2) Recessive (function: nRec)            #
############################################
# function 1: nDOM (require: stringr!)
# requirements input files:
# Chromosome, Region and Allele, for family: also Zygosity
#' Variant filtering with dominant inheritance
#' @description nDom should be used to filter variants under a dominant mode of 
#'   inheritance.
#' @param x the case(s).
#' @param y optional, the control(s).
#' @param family optional, "P-F" or "Ps-F" mode, indicating a (healthy) 
#'   parent(s)- (affected) progeny relation.
#' @details nDom checks for common variants within cases. If variants from 
#'   control(s) are present, they can be used to filter variants from the cases.
#'   In addition, various detectance and penetrance parameters can be specified 
#'   during the processing.
#'   
#'   If family is specified, penetrance is assumed to be complete. In addition,
#'   the causal variant is assumed to be heterozygous in the case.
#' @author Bart Broeckx
#' @examples
#' output <- nDom("CLCfile1", "CLCfile2"  )
#' output
#' output <- nDom(c("CLCfile1","CLCfile2"))
#' output
#' 


nDom <- function(x,y,family)
{ # general prints
  Text <- paste("Number of cases: ", length(x), "\n", sep="")
  cat(Text)
  if (!missing(y)){
    Text <- paste("Number of controls: ", length(y), "\n", sep="")
  } else {
    Text <- paste("Number of controls: 0", "\n")
  }
  cat(Text)
  if (!missing(family)) { 
    if (family == "P-F" && length(x) == "1" && length(y) == "1" ){
      Text <- paste("P-F mode: 1 healthy parent, 1 affected progeny", "\n", sep="")
      cat(Text)
    } else if (family == "Ps-F" &&  length(x) == "1" && length(y) == "2" ) {
      Text <- paste("Ps-F mode: 2 healthy parents, 1 affected progeny", "\n", sep="")
      cat(Text)
    } else {
      stop("family not correct: check number of cases/controls or selected inheritance") 
    }
  }
  
  # first x-vector (cases!!)
  output=NULL
  for (i in 1:length(x)){
    nameValue<-get(x[[i]])
    nameValue$comb <-paste(nameValue$Chromosome,nameValue$Region, nameValue$Allele, sep=" ")
    output <-rbind(output, nameValue)
    Text <- paste("Finished case",i, "\n" )
    cat(Text)
  }
  Freq <-as.data.frame(table(output$comb), stringsAsFactors=FALSE)
  
  # followed by y vector (controls!!)
  if (!missing(y) ) {
    output=NULL
    for (i in 1:length(y)){
      nameValuey<-get(y[[i]])
      nameValuey$comb <-paste(nameValuey$Chromosome,nameValuey$Region, nameValuey$Allele, sep=" ")
      output <-rbind(output, nameValuey)
      Text <- paste("Finished control",i, "\n" )
      cat(Text)
    }
    Freq2 <-as.data.frame(table(output$comb), stringsAsFactors=FALSE)
    
    # next, which penetrance level is required?
    if (missing(family)) {
      Text <- paste("Which penetrance level is required?", "\n" )
      cat(Text)
      for (i in 0:length(y)) {
        Text <- paste("level (",i,"): ", round(length(x)/(length(x)+i)*100),
                      "%", " (",length(x),"/",length(x),"+",i,")" ,"\n", sep="")
        cat(Text)
      }
      level <- as.numeric(readline("Choose level:"))
      id <-which(Freq2$Freq >= (level+1)) 
    } else if (family == "P-F" | family == "Ps-F") { 
      id <- which(Freq2$Freq >= "1") 
    }
    # exlcusion of those variants dependent on the penetrance level
    if (length(id) > 0) { Freq2 <- Freq2[id,]
                          
                          # exclusion of those variants present in the controls
                          id <- which(Freq$Var1 %in% Freq2$Var1)
                          if (length(id) > 0) {Freq<- Freq[-id,]}
    }
  }
  # finally, which detectance level is required?
  if (missing(family)) {
    
    Text <- paste("Which detectance level is required?", "\n" )
    cat(Text)
    for (i in 1:length(x)) {
      Text <- paste("level (",i,"): ", i/length(x)*100,"%"," (",i,"/",length(x),")","\n", sep="")
      cat(Text)
      
    }
    level <- readline("Choose level:")
    id <-which(Freq$Freq >= level)
  } else if (family == "P-F" | family == "Ps-F") { 
    id <- which(nameValue$comb %in% Freq$Var1)
    Zygosity <-nameValue$Zygosity
    Zygosity <-Zygosity[id]
    Freq <-cbind(Freq, Zygosity)              
    id <-which(Freq$Zygosity== "Heterozygous")          
    Freq <-Freq[id,]
    id <- which(Freq$Freq >= "1")   
  }
  
  if (length(id) > 0){
    q <-as.data.frame(stringr::str_split_fixed(Freq[id,1]," ",3), stringsAsFactors=FALSE) 
    q <-cbind(q,Freq[id,2])
    colnames(q) <- c("Chromosome", "Region", "Allele","Number of Samples")
    
    q
  } else {
    Text <- "No variants retained"
    cat(Text)}
  
  
}

##############################################
# function 2: nREC (require stringr!)
# requirements input files:
# Chromosome, Region and Allele AND Zygosity
#' Variant filtering with recessive inheritance
#' @description nRec should be used to filter variants under a recessive mode of 
#'   inheritance.
#' @param x the case(s).
#' @param y optional, the control(s).
#' @param family optional, "P-F" or "Ps-F" mode, indicating a (healthy) 
#'   parent(s)- (affected) progeny relation.
#' @details nRec checks for common variants that are homozygous in the cases. If variants from  
#'   control(s) are present, they can be used to filter variants from the cases. However, only homozygous variants are used to filter.
#'   In addition, various detectance and penetrance parameters can be specified 
#'   during the processing.
#'   
#'   If family is specified, penetrance is assumed to be complete. In addition,
#'   the causal variant is assumed to be homozygous in the case and heterozygous in the parent(s).
#' @author Bart Broeckx
#' @examples
#' output <- nRec("CLCfile1", "CLCfile2")
#' output
#' output <- nRec(c("CLCfile1","CLCfile2"))
#' output


nRec <- function(x,y, family)
{
  Text <- paste("Number of cases: ", length(x), "\n", sep="")
  cat(Text)
  if (!missing(y)){
    Text <- paste("Number of controls: ", length(y), "\n", sep="")
  } else {
    Text <- paste("Number of controls: 0", "\n")
  }
  cat(Text)
  if (!missing(family)) {
    if (family == "P-F" && length(x) == "1" && length(y) == "1" ){
      Text <- paste("P-F mode: 1 healthy parent, 1 affected progeny", "\n", sep="")
      cat(Text)
    } else if (family == "Ps-F" &&  length(x) == "1" && length(y) == "2" ) {
      Text <- paste("Ps-F mode: 2 healthy parents, 1 affected progeny", "\n", sep="")
      cat(Text)
    } else {
      stop("family not correct: check number of cases/controls or selected inheritance") 
    }
  }
  # first x-vector (cases!!)
  output=NULL
  for (i in 1:length(x)){
    nameValue<-get(x[[i]])
    id <- which(nameValue$Zygosity == "Homozygous")
    nameValue <- nameValue[id,]
    nameValue$comb <-paste(nameValue$Chromosome,nameValue$Region, nameValue$Allele, sep=" ")
    output <-rbind(output, nameValue)
    Text <- paste("Finished case",i, "\n" )
    cat(Text)
  }
  
  Freq <-as.data.frame(table(output$comb), stringsAsFactors=FALSE)
  
  # followed by y vector (controls!!)
  if (!missing(y) ) {
    output=NULL
    for (i in 1:length(y)){
      nameValuey<-get(y[[i]])
      id <- which(nameValuey$Zygosity == "Homozygous")
      nameValuey <- nameValuey[id,]
      nameValuey$comb <-paste(nameValuey$Chromosome,nameValuey$Region, nameValuey$Allele, sep=" ")
      output <-rbind(output, nameValuey)
      Text <- paste("Finished control",i, "\n" )
      cat(Text)
    }
    Freq2 <-as.data.frame(table(output$comb), stringsAsFactors=FALSE)
    
    # next, which penetrance level is required?
    if (missing(family)) {
      Text <- paste("Which penetrance level is required?", "\n" )
      cat(Text)
      for (i in 0:length(y)) {
        Text <- paste("level (",i,"): ", round(length(x)/(length(x)+i)*100),
                      "%", " (",length(x),"/",length(x),"+",i,")" ,"\n", sep="")
        cat(Text)
      }
      level <- as.numeric(readline("Choose level:"))
      id <-which(Freq2$Freq >= (level+1)) 
      
    } else if (family == "P-F" | family == "Ps-F") { 
      id <- which(Freq2$Freq >= "1") 
    }
    
    # exclusion of those variants present in the controls
    if (length(id) > 0) { Freq2 <- Freq2[id,]
                          
                          # exclusion of those variants present in the controls
                          id <- which(Freq$Var1 %in% Freq2$Var1)
                          if (length(id) > 0) {Freq<- Freq[-id,]}
    }
  }
  
  # finally, which detectance level is required?
  if (missing(family)) {
    
    Text <- paste("Which detectance level is required?", "\n" )
    cat(Text)
    for (i in 1:length(x)) {
      Text <- paste("level (",i,"): ", i/length(x)*100,"%"," (",i,"/",length(x),")","\n", sep="")
      cat(Text)
      
    }
    level <- readline("Choose level:")
    id <-which(Freq$Freq >= level)
  } else if (family == "P-F" ) { 
    output=NULL
    for (i in 1:length(y)){
      nameValuez<-get(y[[i]])
      id <- which(nameValuez$Zygosity == "Heterozygous")
      nameValuez <- nameValuez[id,]
      nameValuez$comb <-paste(nameValuez$Chromosome,nameValuez$Region, nameValuez$Allele, sep=" ")
      output <-rbind(output, nameValuez)
    }  
    Freq2 <-as.data.frame(table(output$comb), stringsAsFactors=FALSE)
    id <- which(Freq$Var1 %in% Freq2$Var1)
    
  } else if (family == "Ps-F"){
    output=NULL
    for (i in 1:length(y)){
      nameValuez<-get(y[[i]])
      id <- which(nameValuez$Zygosity == "Heterozygous")
      nameValuez <- nameValuez[id,]
      nameValuez$comb <-paste(nameValuez$Chromosome,nameValuez$Region, nameValuez$Allele, sep=" ")
      output <-rbind(output, nameValuez)
    }  
    Freq2 <-as.data.frame(table(output$comb), stringsAsFactors=FALSE)
    id <- which(Freq2$Freq == "2")
    Freq2 <- Freq2[id,]
    id <- which(Freq$Var1 %in% Freq2$Var1)
  }
  
  if (length(id) > 0){
    q <-as.data.frame(stringr::str_split_fixed(Freq[id,1]," ",3), stringsAsFactors=FALSE) 
    q <-cbind(q,Freq[id,2])
    colnames(q) <- c("Chromosome", "Region", "Allele","Number of Samples")
    q
  } else {
    Text <- "No variants retained"
    cat(Text)}
  
}


############################################
# GENE LEVEL (or other grouping var)       #
# 1) CLCfile                               #
# 2) Dominant (function: gDOM)             #
# 3) Recessive (function: gRec)            #
# 4) vcffile                               #
############################################
# Function 1: CLCfile (require stringr!)
# preparation of CLC files for gDom and gRec
# should be done SAMPLE PER SAMPLE!
#' Prepare CLC output for filtering
#' 
#' @description CLCfile should be used to process CLC output prior downstream
#'   analysis with gDom and gRec.
#' @param x the data frame to be processed.
#' @param y the column name containing the grouping variable.
#' @param multiple does one line contain >1 group? If true, specify split.
#' @param split the split that separates the group in one line.
#' @details This function is specifically designed to process the output from 
#'   CLC genomics workbench. It should be used before the variant filters gDom 
#'   and gRec are used to emphasize the grouping (e.g. gene, exon) to be used in
#'   downstream analysis. As often is the case, a variant might belong to >1 
#'   group (e.g. 2 genes). If so, the file is automatically reformatted to 
#'   assure that each line contains only one variant per group, if multiple is 
#'   TRUE and the split is specified. If multiple is FALSE, only the first part 
#'   of the y-column is maintained.This function should not be used for nDom and
#'   nRec.
#' @author Bart Broeckx 
#' @examples
#' CLCfile1proc <-CLCfile(CLCfile1, "Coding.region.change", TRUE, "; ")
#' CLCfile2proc <-CLCfile(CLCfile2, "Coding.region.change", TRUE, "; ")

CLCfile <- function(x,y, multiple,split)
{ # splitting up lines  
  if (!missing(multiple)) {
    if (multiple == "TRUE") {
      if (!missing(split)) {
        group <-x[,y]  
        group <-unlist(strsplit(group, split))
        repeats  <- 1+stringr::str_count(x[,y], split)
        x <-x[rep(seq_len(nrow(x)),repeats),]
        group <-gsub(":.*", "", group)
        output <-cbind(x,group)
      }
    }
  }
  # without split 
  if (missing(multiple) | multiple == "FALSE") {
    group <-x[,y]
    group <-gsub(":.*", "", group)
    output<-cbind(x,group)
  }
  output
}

##############################################
# Function 2: gDOm (require stringr!)
# filtering of previously processed files
# can be done for groups of samples
#' Variant filtering with dominant inheritance considering functional unit.
#' @description gDom should be used to filter variants under a dominant mode of 
#'   inheritance when common functional units (a gene, an exon, ...) should be considered.
#' @param x the case(s).
#' @param y optional, the control(s).
#' @details gDom checks for common functional units within cases. If variants from 
#'   control(s) are present, they can be used to filter variants from the cases.
#'   In addition, various detectance and penetrance parameters can be specified 
#'   during the processing. 
#'   
#'   Before using this function, each file from each case should have been processed by CLCfile or VCFfile with annot.  
#' @author Bart Broeckx
#' @examples
#' CLCfile1proc <-CLCfile(CLCfile1, "Coding.region.change", TRUE, "; ")
#' CLCfile2proc <-CLCfile(CLCfile2, "Coding.region.change", TRUE, "; ")
#'  output <- gDom("CLCfile1proc", "CLCfile2proc")
#'  output
#'  output <- gDom(c("CLCfile1proc","CLCfile2proc"))
#'  output




gDom <- function(x,y)
{ # first some general remarks
  Text <- paste("Number of cases: ", length(x), "\n", sep="")
  cat(Text)
  if (!missing(y)){
    Text <- paste("Number of controls: ", length(y), "\n", sep="")
    cat(Text)
    # followed by y vector (controls!!)
    outputy=NULL
    for (i in 1:length(y)){
      nameValuey<-get(y[[i]])
      nameValuey$comb <-paste(nameValuey$Chromosome,nameValuey$Region, nameValuey$Allele, sep=" ")
      id <- which(duplicated(nameValuey$comb))
      if (length(id) >0){
        nameValuey <- nameValuey[-id,]
      }
      outputy <-rbind(outputy, nameValuey)
      Text <- paste("Finished control",i, "\n" )
      cat(Text)
    }
    Freq2 <-as.data.frame(table(outputy$comb), stringsAsFactors=FALSE)
    # next, which penetrance level is required?
    Text <- paste("Which penetrance level is required?", "\n" )
    cat(Text)
    for (i in 0:length(y)) {
      Text <- paste("level (",i,"): ", round(length(x)/(length(x)+i)*100),
                    "%", " (",length(x),"/",length(x),"+",i,")" ,"\n", sep="")
      cat(Text)
    }
    levelpenet <- as.numeric(readline("Choose level:"))
    id <-which(Freq2$Freq >= (levelpenet+1)) 
    
    
    # exclusion of those variants dependent on the penetrance level
    Freq2 <- Freq2[id,]
    
  }else {
    Text <- paste("Number of controls: 0", "\n")
  }
  
  
  
  # filtering:
  # first x-vector (cases!!)
  output=NULL
  variantoutput=NULL
  for (i in 1:length(x)){
    nameValue<-get(x[[i]])
    nameValue$comb <-paste(nameValue$Chromosome,nameValue$Region, nameValue$Allele, sep=" ")
    
    
    # exclusion of those variants present in the controls
    if (!missing(y) && nrow(Freq2) != "0") {
      
      id <- which(nameValue$comb %in% Freq2$Var1)
      if (length(id) > 0) {nameValue<- nameValue[-id,]}
    }
    # intermezzo: grouping
    Freq <-as.data.frame(table(nameValue[,"group"]), stringsAsFactors=FALSE)
    id <- which(Freq$Freq >=1)
    Freq <- Freq[id,]
    output <-c(output, Freq$Var1)
    # variantsretained
    nameValue$overall <-paste(nameValue$comb, nameValue[,"group"], sep=" / ")
    variantoutput<-c(variantoutput, nameValue$overall)
    
    Text <- paste("Finished case",i, "\n" )
    cat(Text)
    
  }
  Freq <- as.data.frame(table(output))
  variants <-as.data.frame(table(variantoutput))
  # which detectance level is required?
  
  Text <- paste("Which detectance level is required?", "\n" )
  cat(Text)
  for (i in 1:length(x)) {
    Text <- paste("level (",i,"): ", i/length(x)*100,"%"," (",i,"/",length(x),")","\n", sep="")
    cat(Text)
  }
  level <- readline("Choose level:")
  
  
  # followed by final filtering
  
  id <-which(Freq$Freq >= level)
  if (length(id) > 0){
    q <-Freq[id,"output"]
    r <-as.data.frame(stringr::str_split_fixed(variants[,1]," / ",2), stringsAsFactors=FALSE) 
    
    id <- which(r[,2] %in% q)
    r <- r[id,]
    r<-cbind(as.data.frame(stringr::str_split_fixed(r[,1]," ",3), stringsAsFactors=FALSE),r[,2]) 
    r <- cbind(r,variants[id,"Freq"])
    colnames(r) <- c("Chromosome", "Region", "Allele", "Group", "Number of samples")
    r
    Text <- paste("Number of variants retained: ", nrow(r), "\n", sep="")
    cat(Text)
    Text <- paste("Number of genes retained: ", nrow(as.data.frame(table(r$Group))), "\n", sep="")
    cat(Text)
    r
    } else {
    Text <- "None retained"
    cat(Text)}
  }

##############################################
# Function 3: gRec
# to be used after filtering of previously processed files
#'Variant filtering with recessive inheritance considering functional unit.
#'@description gRec should be used to filter variants under a recessive mode of 
#'  inheritance when common functional units (a gene, an exon, ...) should be 
#'  considered.
#'@param x the case(s).
#'@param y optional, the control(s).
#'@param list logical, specify whether a list should be returned (TRUE) or not
#'  (false).
#'@details gRec checks for common functional units within cases. For a unit to 
#'  be retained, it should have at least one homozygous non-reference variant or
#'  at least be compound heterozygous.  If variants from control(s) are present,
#'  they can be used to filter variants from the cases. Homozygous variants 
#'  present in the controls will be used to filter variants from the cases 
#'  immediately. Pairwise combinations of heterozygous variants present in the
#'  controls can also be used to filter variants in the cases. In addition,
#'  various detectance and penetrance parameters can be specified during the
#'  processing.
#'  
#'  Before using this function, each file from each case should have been 
#'  processed by CLCfile or VCFfile with annot.
#'  
#'  The parameter list allows you to specify whether the output should contain
#'  only the retained variants and genes (FALSE) only or together with the
#'  compound heterozygous variants per sample (TRUE). This parameter is only
#'  useful when you have both case(s) and control(s) and when compound
#'  heterozygous variants were present in the controls. For further processing
#'  with commonvar, list should be FALSE.
#'@author Bart Broeckx
#' @examples 
#' CLCfile1proc <-CLCfile(CLCfile1, "Coding.region.change", TRUE, "; ")
#' CLCfile2proc <-CLCfile(CLCfile2, "Coding.region.change", TRUE, "; ")
#'output <- gRec("CLCfile1proc", "CLCfile2proc", TRUE)
#'output
#'output <- gRec(c("CLCfile1proc","CLCfile2proc"),, TRUE)
#'output
gRec <- function(x,y,list)
{ # first some general remarks
  if (missing(list)) {
    stop("specify list")
  }
  Text <- paste("Number of cases: ", length(x), "\n", sep="")
  cat(Text)
  if (!missing(y)){
    Text <- paste("Number of controls: ", length(y), "\n", sep="")
    cat(Text)
    # followed by y vector (controls!!)
    homo=NULL
    hetero=NULL
    for (i in 1:length(y)){
      nameValuey<-get(y[[i]])
      nameValuey$comb <-paste(nameValuey$Chromosome,nameValuey$Region, nameValuey$Allele, sep=" ")
      id <- which(nameValuey$Zygosity =="Homozygous")
      homozygous <- unique(nameValuey[id, "comb"])
      heterozygous <- nameValuey[-id,]
      # following line has been added
      heterozygous$group <- factor(heterozygous$group)
      # previous line has been added
      pairwisecontrol = NULL
      for (z in 1:length(levels(heterozygous[, "group"]))) {
        Test <- subset(heterozygous[, "comb"], heterozygous[,"group"]==levels(heterozygous[,"group"])[z])
        if (length(Test)>1){
          a <-combn(Test,2, paste, collapse=" / ")
          pairwisecontrol <- c(pairwisecontrol,a)}
        
      }
      hetero <-c(hetero, unique(pairwisecontrol))
      homo <-c(homo, homozygous)
      
      Text <- paste("Finished control",i, "\n" )
      cat(Text)
    }
    Freqhomo <-as.data.frame(table(homo), stringsAsFactors=FALSE)
    Freqhetero <-as.data.frame(table(hetero), stringsAsFactors=FALSE)
    
    # next, which penetrance level is required?
    Text <- paste("Which penetrance level is required?", "\n" )
    cat(Text)
    for (i in 0:length(y)) {
      Text <- paste("level (",i,"): ", round(length(x)/(length(x)+i)*100),
                    "%", " (",length(x),"/",length(x),"+",i,")" ,"\n", sep="")
      cat(Text)
    }
    levelpenet <- as.numeric(readline("Choose level:"))
    id <-which(Freqhomo$Freq >= (levelpenet+1)) 
    id2 <-which(Freqhetero$Freq >= (levelpenet+1)) 
    
    # exclusion of those variants dependent on the penetrance level
    Freqhomo <- Freqhomo[id,]
    Freqhetero <- Freqhetero[id2,]
    
  }else {
    Text <- paste("Number of controls: 0", "\n")
  }
  
  # filtering:
  # first x-vector (cases!!)
  output=NULL
  variantoutput=NULL
  pairwisecombinations <- vector(mode="list", length=length(x))
  names(pairwisecombinations) <- x
  for (s in 1:length(x)){
    nameValue<-get(x[[s]])
    nameValue$comb <-paste(nameValue$Chromosome,nameValue$Region, nameValue$Allele, sep=" ")
    # two step filtering: 
    # filter out all the homozygous variants present in the controls:
    if (!missing(y) && nrow(Freqhomo) != "0") {
      
      id <- which(nameValue$comb %in% Freqhomo$homo)
      if (length(id) > 0) {nameValue<- nameValue[-id,]}
    }
    # give homozy value 2, hetero value 1 for further filtering
    id <- which(nameValue$Zygosity == "Homozygous")
    Zygo <- rep(2, length(id))
    homo <- cbind(nameValue[id,], Zygo)
    Zygo <- rep(1, nrow(nameValue)-nrow(homo))
    hetero <- cbind(nameValue[-id,], Zygo)
    
    
    
    # second filtering: filter out pairwise combinations
    # generate pairwise combn
    
    if (!missing(y) && nrow(Freqhetero) != "0") { 
      pairwisecase = NULL
      # the following line has been added:
      hetero$group<- factor(hetero$group)
      # the previous line has been added
      for (t in 1:length(levels(hetero[, "group"]))) {
        Test <- subset(hetero[, "comb"], hetero[,"group"]==levels(hetero[,"group"])[t])
        if (length(Test)>1){
          a <-combn(Test,2, paste, collapse=" / ")
          pairwisecase <- c(pairwisecase,a)}
        
      }
      # filterstep:
      
      id <- which(pairwisecase %in% Freqhetero$hetero)
      if (length(id) == 0) {
        heteroretained <- unique(unlist(strsplit(pairwisecase, " / ")))
        id <- which(hetero$comb %in% heteroretained)
        hetero <- hetero[id,]
      } else if(length(id) > 0) {
        if (length(unique(id)) != length(pairwisecase)) {
          # the next line comma has been removed
          pairwisecase<- pairwisecase[-id]
          # the previous line comma has been removed
          heteroretained <- unique(unlist(strsplit(pairwisecase, " / ")))
          id <- which(hetero$comb %in% heteroretained)
          hetero <- hetero[id,]
        } else if (length(unique(id)) == length(pairwisecase)){
          hetero=NULL
          pairwisecase <- "none"
        }
      }
    }
    endvariants <-rbind(homo, hetero)
    if (!missing(y) && nrow(Freqhetero) != "0") {
      pairwisecombinations[[s]] <- pairwisecase }
    
    
    # followed by grouping per "gene" (defined as group)
    grouping = NULL
    for (i in 1:length(levels(endvariants[,"group"]))) {
      Test <- subset(endvariants, endvariants[,"group"]==levels(endvariants[,"group"])[i])
      if (nrow(Test) > 0) {
        count <-sum(Test$Zygo)
        Gene <-levels(endvariants[,"group"])[i]
        c <- data.frame(count,Gene)
        grouping <- rbind(grouping,c)}
      
    }
    
    # followed by selection of those "genes" with count >=2
    id <- which(grouping$count >= 2)
    grouping <- grouping[id,2]
    
    # selection of those genes with a count of >=2 from the original dataframe:
    id <- which(endvariants[,"group"] %in% grouping)
    endvariants <-endvariants[id,]
    
    Freq <-as.data.frame(table(endvariants[,"group"]), stringsAsFactors=FALSE)
    id <- which(Freq$Freq > 0)
    Freq <- Freq[id,]
    output <-c(output, Freq$Var1)
    # variantsretained
    endvariants$overall <-paste(endvariants$comb, endvariants[,"group"], sep=" / ")
    variantoutput<-c(variantoutput, endvariants$overall)
    
    Text <- paste("Finished case",s, "\n" )
    cat(Text)
    
  } # end of for loop x
  Freq <- as.data.frame(table(output))
  variants <-as.data.frame(table(variantoutput))
  
  # which detectance level is required?
  
  Text <- paste("Which detectance level is required?", "\n" )
  cat(Text)
  for (i in 1:length(x)) {
    Text <- paste("level (",i,"): ", i/length(x)*100,"%"," (",i,"/",length(x),")","\n", sep="")
    cat(Text)
  }
  level <- readline("Choose level:")
  
  
  # followed by final filtering
  
  id <-which(Freq$Freq >= level)
  if (length(id) > 0){
    q <-Freq[id,"output"]
    r <-as.data.frame(stringr::str_split_fixed(variants[,1]," / ",2), stringsAsFactors=FALSE) 
    
    id <- which(r[,2] %in% q)
    r <- r[id,]
    filteringpairwisecombinations<- r[,1]
    r<-cbind(as.data.frame(stringr::str_split_fixed(r[,1]," ",3), stringsAsFactors=FALSE),r[,2]) 
    r <- cbind(r,variants[id,"Freq"])
    colnames(r) <- c("Chromosome", "Region", "Allele", "Group", "Number of samples")
    r
    Text <- paste("Number of variants retained: ", nrow(r), "\n", sep="")
    cat(Text)
    Text <- paste("Number of genes retained: ", nrow(as.data.frame(table(r$Group))), "\n", sep="")
    cat(Text)
    
    if (list == "TRUE" && !missing(y) && nrow(Freqhetero) != "0"){
      resultlist <- vector(mode="list", length=length(pairwisecombinations)+1)
      resultlist[[1]] <- r
      
      varend <- vector(mode="list", length=length(x))
      
      for (i in 1:length(pairwisecombinations)){
        var<-pairwisecombinations[[i]]
        vardef=NULL
        for (w in 1:length(var)){
          varret<- var[w]
          substr<-unlist(strsplit(varret, " / "))
          presentone<-substr[1] %in% filteringpairwisecombinations
          presenttwo <-substr[2] %in% filteringpairwisecombinations
          if (presentone == TRUE | presenttwo == TRUE) {
            vardef <- c(vardef,varret)
          }
        }
        if (length(vardef) > 0) {
          varend[[i]] <- vardef
        } else {
          varend[[i]] <- "none" 
        }
      }
      
      resultlist[2:length(resultlist)] <- varend
      allnames<- c("endresult", names(pairwisecombinations))
      names(resultlist)<- allnames
      resultlist
      
    } else if (list == "FALSE" | missing(y)) {r }
  } else {
    Text <- "None retained"
    cat(Text)}
  
}



##############################################
# Function 4: vcffile
# to be used after filtering of previously processed files

#' Prepare VCF file for filtering.
#' @description VCFfile should be used to process VCF files prior downstream 
#'   analysis with any other function.
#' @param x the VCF file to be processed.
#' @param sample specify the column name containing the variant info for the 
#'   sample of interest.
#' @param filter logical, should the variants have a certain quality value in the 7th 
#'   mandatory column of the VCF format? 
#' @param value optional, which variants should be retained? (e.g. PASS)
#' @details This function should be used to process VCF files to make them 
#'   compatible for further filtering. If wanted, a preliminary filtering can be
#'   performed by specifiying filter = TRUE and a value in "value". As VCF files
#'   can contain variant information from >1 sample, the column name of the
#'   sample of interest should be specified.
#' @author Bart Broeckx
#'   @examples 
#'   # test 1: filter = TRUE
#'   x <- test
#'   sample <- "V10"
#'   filter <- TRUE
#'   value <- "PASS"
#'   VCFfile(x,sample, filter, value)
#'   
#'   # test 2: filter = FALSE
#'   x <- test
#'   sample <- "V10"
#'   filter <- FALSE
#'   VCFfile(x,sample, filter)


VCFfile <-function(x,sample,filter,value){ 
  id <- which(is.na(x[, sample]))
  if (length(id) >=1) {
    x <- x[-id,]
  } 
  if(filter == TRUE){
    if(!missing(value)){
      id <- which(x[,7] == value)
      x <-x[id,]
    }
  }
  # standard fields
  y <- x[,c(1,2,5,9)]
  # add the field containing the sample
  x <- x[,sample]
  x <- cbind(y,x)
  # start reformatting
  output = NULL
  for(i in 1:nrow(x)){
    var <- x[i,]
    format <-unlist(strsplit(as.character(var[,4]), ":"))
    id <-which(format == "GT")
    variant <-unlist(strsplit(as.character(var[,"x"]), ":"))
    variant <- variant[id]
    splitvariants <-unlist(strsplit(variant, "/"))
    variant1 <- as.numeric(splitvariants[1])
    variant2 <- as.numeric(splitvariants[2])
    Zygosity<- ifelse(variant1 == variant2,"Homozygous", "Heterozygous")
    allele <-unlist(strsplit(as.character(var[,3]),","))
    allele1 <- allele[variant1]
    allele2 <- allele[variant2]
    
    final2  <- cbind(var[,c(1:2)], allele2, Zygosity)
    colnames(final2) <- c("Chromosome", "Region", "Allele", "Zygosity")
    if (Zygosity == "Homozygous") {
      end <- final2
      
    } else  if (Zygosity == "Heterozygous"){ 
      if (variant1 != "0") {
        final1  <-cbind(var[,c(1:2)], allele1, Zygosity)   
        colnames(final1) <- c("Chromosome", "Region", "Allele", "Zygosity")
        end <-rbind(final1,final2)
      } else if (variant1 == "0"){
        end <- final2
      }
    } 
    
    
    output <- rbind(output, end)
  }
  output
}  

############################################
# Annotation of a file and variant filter  #
# 1) Annot                                 #
# 2) Varfilter                             #
############################################
# werkende ANNOT functie.
#'Annotate variant data frame.
#'
#'@description annot can be used to annotate a data frame containing sequencing 
#'  variants.
#'@param x the data frame to be annotated.
#'@param y the data frame containing the annotation.
#'@param type specify whether y is a BED or GTF file.
#'@param nomatch should variants located outside any of the regions found in the
#'  annotation data frame be removed (FALSE) or withheld (TRUE).
#'@param CLC logical, origin of the data frame to be annotated is CLC genomics 
#'  workbench?
#'@details This function can be used to annotate data frames. For further 
#'  processing by gDom, gRec and commonvar, annotation is necessary. For nDom 
#'  and nRec, annotation is optional. Both variant files from CLC Genomics 
#'  Workbench and VCF files (after VCFfile) can be processed. If a variant can 
#'  be allocated to >1 group (a gene, exon or other), the variant information is
#'  spread over n rows (with n the number of groups).
#'  
#'  If the file originates from the CLC genomics workbench (directly or after
#'  processing with nDom or nRec), always specify CLC is TRUE. This is important
#'  as CLC specifies chromosomal location for MNVs in a different way than VCF.
#'@author Bart Broeckx
#'  @examples data(CLCfile1)
#'  data(genBED)
#'  # produces an annotated CLC file with all variants retained.
#'  AnnotCLCfile1 <- annot(CLCfile1, genBED, type="BED", nomatch=TRUE, CLC=TRUE)
#'  # produces an annotated CLC file with only the variants that were allocated 
#'  # to a gene being retained.
#'  AnnotCLCfile2 <- annot(CLCfile1, genBED, type="BED", nomatch=FALSE, CLC=TRUE)
annot <- function(x, y, type, nomatch, CLC){
  if (missing(x) | missing(y) | missing(type)){
    stop("specify x, y and type")
  } else {
    # preparatory phase 
    if  (CLC == TRUE){
      
      idregion <-grep("[:^:]", x$Region)
      if (length(idregion) >=1){
        originalregion <- x[idregion,"Region"]
        x$Region <- gsub("[:^:].*", "", x$Region)  
      }
      
      idregion2 <- grep("[:..:]", x$Region)
      if (length(idregion2) >=1){
        originalregion2 <- x[idregion2, "Region"]
        x$Region <- gsub("[:..:].*", "", x$Region)
        
      }
      
    }
    
    x[,1]<-as.factor(x[,1])
    x[,2]<- as.numeric(x[,2])
    
    
    if (type == "BED") {
      y[,2]<- y[,2] + 1 
    }
    output <- data.frame()
    for (i in 1:nrow(x)){
      linei <- x[i,]
      id <- which(linei[,1] == y[,1])
      subset <- y[id,]
      if (type == "BED") {
        id <- which(linei[,2] >= subset[,2] & linei[,2] <= subset[,3]) 
        if  (CLC == TRUE){
          if (length(idregion) >=1 | length(idregion2) >=1 ){
            for (f in 1:length(idregion)){
              if (i == idregion[f]){ 
                linei[,2]<-originalregion[f]
                
              }  
            }
            for (g in 1:length(idregion2)){
              if (i == idregion2[g]){ 
                linei[,2]<-originalregion2[g]
                
              }  
            }
          }
        }
        
        if (length(id)>0){
          annotation<-as.vector(unlist(subset[id, 4]))
          nrep <- length(annotation)  
          replicated<- linei[rep(seq_len(nrow(linei)),nrep),]
          result <-cbind(replicated, annotation)
          
          output <-rbind(output, result)
        } else {
          if (nomatch == "TRUE"){
            annotation  <- "NOT FOUND"
            result <- cbind(linei, annotation)
            output <-rbind(output, result)
          } 
        }
      } else if (type== "GTF"){
        id <- which(linei[,2] >= subset[,4] & linei[,2] <= subset[,5]) 
        if  (CLC == TRUE){
          if (length(idregion) >=1 | length(idregion2) >=1 ){
            for (f in 1:length(idregion)){
              if (i == idregion[f]){ 
                linei[,2]<-originalregion[f]
                
              }  
            }
            for (g in 1:length(idregion2)){
              if (i == idregion2[g]){ 
                linei[,2]<-originalregion2[g]
                
              }  
            }
          }        
        }
        if (length(id)>0){
          annotation<-as.vector(unlist(subset[id, 10]))
          nrep <- length(annotation)  
          replicated<- linei[rep(seq_len(nrow(linei)),nrep),]
          result <-cbind(replicated, annotation)
        } else {
          if (nomatch == "TRUE"){
            annotation  <- "NOT FOUND"
            result <- cbind(linei, annotation)
            output <-rbind(output, result)
          } 
        }
      }
      
    }
    
  }
  
  output
}

#########################################
#' Preparing variant database.
#'
#'@description prepvar is used to prepare a variant database for further
#'  filtering with varfilter.
#'@param y the data frame containing the variants.
#'@param MAF optional numeric, the MAF of the non-reference variants to be 
#'  retained for filtering.
#'@param reference the column containing the reference sequence in the y data 
#'  frame.
#'@details This function prepares the variants from a database like dbSNP to be
#'  used for filtering variants in cases. MAF is optional and specifies a
#'  minimum allele frequency for the database variants to be retained for
#'  variant filtering. The input of y should be the standard output of the UCSC
#'  table browser ("all fields from selected table"). For large files, running
#'  time can be quite long. An alternative supporting parallel computing is
#'  prepvarpar. For the actual filtering, use varfilter.
#'@author Bart Broeckx
#'@examples
#' y <- SNP
#' reference <- "refNCBI"
#' MAF <- 0.03
#' a <-prepvar(y, MAF, reference)
prepvar <- function(y, MAF, reference){
  if (missing(y) | missing(reference)) {
    stop("specify y and reference")}
  filtering <- data.frame()
  y <- y[, -c(1, 4:7, 10:22, 24, 26)]
  y[, 2] <- y[, 2] + 1
  for (i in 1:nrow(y)) {
    subset <- y[i, ]
    allele <- unlist(strsplit(as.character(subset[, "alleles"]), 
                              ","))
    id <- which(allele != subset[, reference])
    nonrefallele <- allele[id]
    allelefreq <- unlist(strsplit(as.character(subset[, "alleleFreqs"]), 
                                  ","))
    nonreffreq <- allelefreq[id]
    nonreffreq <- as.numeric(nonreffreq)
    nrep <- length(nonreffreq)
    subset <- subset[rep(seq_len(nrow(subset)), nrep), c(1:2)]
    subsetend <- cbind(subset, nonrefallele, nonreffreq)
    filtering <- rbind(filtering, subsetend)
    if (i ==round(nrow(y)/20) | i ==round(nrow(y)/10) | i ==round(nrow(y)/5) |
          i ==round(nrow(y)/3.33) |i ==round(nrow(y)/2.5) | i ==round(nrow(y)/2) | 
          i ==round(nrow(y)/1.66) | i ==round(nrow(y)/1.42) | i ==round(nrow(y)/1.25) | 
          i ==round(nrow(y)/1.11)) {
      Text <- paste(round(i/nrow(y)*100), "procent prepared", "\n", sep=" ")
      cat(Text)
    }
  }
  Text <- paste("Number of variants for filtering: ", nrow(filtering), 
                "\n", sep = "")
  cat(Text)
  if (!missing(MAF)) {
    id <- which(filtering$nonreffreq >= MAF)
    if (length(id) >= 1) {
      filtering <- filtering[id, ]
      Text <- paste("Number of variants with > MAF: ", 
                    nrow(filtering), "\n", sep = "")
      cat(Text)
    }
    else if (length(id) < 1) {
      Text <- paste("No variants for filtering retained!", 
                    "\n", sep = "")
      cat(Text)
    }
  }
  filtering
}

#####################################################
#' Filtering variants from a database.
#'
#'@description varfilter is used to filter case variants using a variant
#'  database (like dbSNP).
#'@param x the data frame to be filtered
#'@param y the data frame containing the variants.
#'@details After preprocessing the variant database with prepvar or
#'  prepvarpar, the actual filtering is done with varfilter.
#'@author Bart Broeckx
#'@examples
#' out <-VCFfile(test, "V10", TRUE, "PASS")
#' a <-prepvar(SNP, 0.03,"refNCBI")
#' filtered <-varfilter(out,a)
varfilter<- function(x,y){
if (missing(x) | missing(y)) {
  stop("specify x and y")}
seqvariants <- paste(x$Chromosome, x$Region, x$Allele, sep = " ")
filtering <- y
filtervariants <- paste(filtering$chrom, filtering$chromStart, 
                        filtering$nonrefallele, sep = " ")
id <- which(seqvariants %in% filtervariants)
orcount <- nrow(x)
Text <- paste("Number of variants to be filtered: ", orcount, 
              "\n", sep = "")
cat(Text)
if (length(id) >= 1) {
  x <- x[-id, ]
}
if (nrow(x) >= 1) {
  Text <- paste("Number of variants from x after filtering retained: ", 
                nrow(x), "\n", sep = "")
  cat(Text)
}
else {
  Text <- paste("No variants left from x after filtering", 
                "\n", sep = "")
  cat(Text)
}
x
}




############################################
# Final putting everything together:       #
# commonvar                                #
############################################
#'Combining the output of several filters.
#'
#'@description commonvar should be used to combine the output of nDom/nRec and
#'  gDom/gRec filters.
#'@param x the data frames to be filtered (> 1).
#'@param group the column names of the columns containing the functional units.
#'@details This function allows the combination of nDom/nRec and gDom/gRec
#'  output. It retains all functional units in common, whithout assumptions on
#'  zygosity (as this should already have been taken care of in the previous
#'  filtering steps). Within each functional unit, the variants present in each
#'  x object are retained. Detectance can be specified and will be used on the
#'  functional group instead of the individual variant.
#'@author Bart Broeckx
#'@examples 
#'output <- nRec("CLCfile1", "CLCfile2"  )
#'output
#'CLCfile1proc <-CLCfile(CLCfile1, "Coding.region.change", TRUE, "; ")
#'CLCfile2proc <-CLCfile(CLCfile2, "Coding.region.change", TRUE, "; ")
#'output2 <- gDom("CLCfile1proc", "CLCfile2proc") 
#'output2
#'Annotoutput2 <- annot(output,genBED, type="BED", nomatch=FALSE, CLC=TRUE)
#'x <- c("output2", "Annotoutput2")
#'group <- c("Group", "annotation")
#'a <-commonvar(x,group)

commonvar <- function(x, group){
  
  # QC
  if (missing(x) | missing(group)) {
    stop("specify x and group")
  }
  if (length(x) != length(group)) {
    stop("x and group not equal length")
  }
  if (length(x) < 2) {
    stop("length(x) should be > 1")
  }
  
  # actual function
  common <- vector()
  for (i in 1:length(x)){
    nameValue<-get(x[[i]])
    grouper <- group[i] 
    var <-as.character(unique(nameValue[,grouper]))
    common <- c(common,var)
  }
  
  commongroup <-as.data.frame(table(common))
  Text <- paste("Which detectance level is required?", "\n" )
  cat(Text)
  for (i in 1:length(x)) {
    Text <- paste("level (",i,"): ", i/length(x)*100,"%"," (",i,"/",length(x),")","\n", sep="")
    cat(Text)
    
  }
  level <- readline("Choose level:")
  id <-which(commongroup$Freq >= level)
  if (length(id) > 0){
    result <-as.character(commongroup[id,-2])
    Text <- paste(length(result), "groups are retained", "\n" )
    cat(Text)
    
    
    output <- vector()
    if (length(result) >=1 ) {
      for (i in 1:length(x)){
        nameValue<-get(x[[i]])
        grouper <- group[i] 
        id <- which(nameValue[,grouper] %in% result)
        selection <- nameValue[id, c("Chromosome", "Region", "Allele", grouper)]
        pasted<-paste(selection$Chromosome, selection$Region, selection$Allele, selection[,grouper], sep=" ")
        output <-c(output, pasted)
      }
    }
    res <-as.data.frame(table(output))
    res$output <- as.character(res$output)
    final<- cbind(as.data.frame(stringr::str_split_fixed(res[,1]," ",4)),res$Freq)
    colnames(final) <- c("Chromosome", "Region", "Allele", "Group", "Freq")
    final
  }else {
    Text <- paste("Nothing retained", "\n" )
    cat(Text)
  }
}

#' CLCfile1 documentation
#' A dataset containing 8 variants (output CLC genomics workbench)
#' 
"CLCfile1"
#' CLCfile2 documentation
#' A dataset containing 9 variants (output CLC genomics workbench)
#' 
"CLCfile2"
#' SNP documentation
#' A dataset containing 472 variants from the dbSNP database
#' 
#' @source \url{http://genome.ucsc.edu/cgi-bin/hgTables?command=start}
"SNP"
#'genBED documentation
#' A bed file containing genes used to annotate files
#' 
#' @source \url{http://genome.ucsc.edu/cgi-bin/hgTables?command=start}
"genBED"
#' test documentation
#' An example of a VCF file.
"test"


#' Preparing variant database (parallel).
#'
#'@description prepvarpar is used to prepare a large variant database for 
#'  further filtering with varfilter.
#'@param y the data frame containing the variants.
#'@param MAF optional numeric, the MAF of the non-reference variants to be 
#'  retained for filtering.
#'@param reference the column containing the reference sequence in the y data 
#'  frame.
#'@param nproc the number of cores or workers used for parallel computing.
#'@details This function prepares the variants from a large database like dbSNP
#'  to be used for filtering variants in cases. MAF is optional and specifies a 
#'  minimum allele frequency for the database variants to be retained for 
#'  variant filtering. The input of y should be the standard output of the UCSC 
#'  table browser ("all fields from selected table"). This function allows 
#'  parallel computing, leading to substantial time improvement when entire 
#'  databases are used. If parallel computing is not necessary or not possible, 
#'  use prepvar.Prior to running this function, call function registerDoParallel
#'  from the doParallel package to permit parallel computing.
#'  @examples
#' #library(doParallel)
#' #registerDoParallel()
#' # nproc <-getDoParWorkers()
#' # a <-prepvarpar(SNP, 0.03,"refNCBI", 1)
prepvarpar <- function(y, MAF, reference, nproc){
  "%dopar%" <- NA
  rm("%dopar%")
  if (missing(y) | missing(reference)) {
    stop("specify y and reference")}
  y <- y[, -c(1, 4:7, 10:22, 24, 26)]
  y[, 2] <- y[, 2] + 1
  indic <-rep(1:nproc, length.out=nrow(y))
  indic <- sort(indic)
  indic <- as.factor(indic)
  y <- cbind(y,indic)
  splits <- iterators::isplit(y,y$indic)
  subsety=NULL
  result<- foreach::foreach(subsety=splits, .combine='rbind') %dopar% {
    filtering <- data.frame()
    for (i in 1:nrow(subsety$value)) {
      subset <- subsety$value[i, ]
      allele <- unlist(strsplit(as.character(subset[, "alleles"]), 
                                ","))
      id <- which(allele != subset[, reference])
      nonrefallele <- allele[id]
      allelefreq <- unlist(strsplit(as.character(subset[, "alleleFreqs"]), 
                                    ","))
      nonreffreq <- allelefreq[id]
      nonreffreq <- as.numeric(nonreffreq)
      nrep <- length(nonreffreq)
      subset <- subset[rep(seq_len(nrow(subset)), nrep), c(1:2)]
      subsetend <- cbind(subset, nonrefallele, nonreffreq)
      filtering <- rbind(filtering, subsetend)
    }
    prelist <-list(value=filtering)
    presubset<- subsety[-1]
    result <- c(prelist,presubset)
  }
  empty=NULL
  for (i in 1:nproc) {
    prepdata<-as.data.frame(result[[i]])
    empty <- rbind(empty,prepdata)
    
  }
  filtering<- empty 
  Text <- paste("Number of variants for filtering: ", nrow(filtering), 
                "\n", sep = "")
  cat(Text)
  if (!missing(MAF)) {
    id <- which(filtering$nonreffreq >= MAF)
    if (length(id) >= 1) {
      filtering <- filtering[id, ]
      Text <- paste("Number of variants with > MAF: ", 
                    nrow(filtering), "\n", sep = "")
      cat(Text)
    }
    else if (length(id) < 1) {
      Text <- paste("No variants for filtering retained!", 
                    "\n", sep = "")
      cat(Text)
    }
  }
  filtering
}
