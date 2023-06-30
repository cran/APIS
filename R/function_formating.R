#' Format data from genotype to recommanded format for APIS
#'
#' @param data genotype data
#' @param SampleName names of each individuals in the data (ordered)
#' @param marker_name name of each snp (ordered)
#' @param ploidy ploidy level of the genotype data
#' @param marker_type marker type (snp or microsatellite)
#'
#' @return a list of 2 : a dataframe of the data formatted and a dataframe with some metrics about each markers
#'
#' @keywords internal
#' @noRd
#'

Run_formating=function(data,SampleName,marker_name,ploidy,marker_type){
  new_data = as.data.frame(matrix(nrow = length(SampleName),ncol = length(marker_name)))
  rownames(new_data)=SampleName
  colnames(new_data)=marker_name
  MAF=c()
  CR=c()
  for (k in 1:(ncol(data)/ploidy)){
    new_data[,k]=apply(X = data[,(ploidy*(k-1)+1):(ploidy*k)],MARGIN = 1,FUN = paste,collapse = "/")
    if (marker_type=="SNP"){
      x=as.vector(as.matrix(data[,(ploidy*(k-1)+1):(ploidy*k)]))
      g=unique(x)[!is.na(unique(x))]
      if (length(g)==2){ # nombre dallele == 2
        f1=length(x[x==g[1] & !is.na(x)])/length(x[!is.na(x)])
        f2=length(x[x==g[2] & !is.na(x)])/length(x[!is.na(x)])
        MAF=c(MAF,min(f1,f2))
      } else if (length(g)>2){
        warning("Max allele is 2 and you have more for your marker : MAF will not be calculated !")
        MAF=c(MAF,NA)
      } else {
        MAF=c(MAF,0)
      }
    } else { # marker_type == 'microsat'
      # Nothing equivalent except heterozygous rate ?
    }
    CR=c(CR,length(which(!is.na(x)))/length(x))
  }
  # Remove duplicate individual by selecting the best CR
  count_na = function(x){
    round(1-(length(which(x==paste0(rep("NA",ploidy),collapse = "/")))/length(x)),3)
  }
  CR_ind=data.frame(SampleName=SampleName,CR=apply(X = new_data,MARGIN = 1,FUN = count_na))
  db_ind = c()
  db_nam = c()
  for (pat in c("_2",".2","_BIS",".BIS")){ # add more if necessary
    indice = regexpr(pattern = pat,text = CR_ind$SampleName,fixed = TRUE) # fixed : so that . is not a regular expression replacing all character
    db_ind = c(db_ind,CR_ind$SampleName[which(indice!=-1)])
    db_nam = c(db_nam,CR_ind$SampleName[which(CR_ind$SampleName %in% substr(x = CR_ind$SampleName[which(indice!=-1)],
                                                                            start = 1,
                                                                            stop = nchar(CR_ind$SampleName[which(indice!=-1)])-nchar(pat,)))])
  }
  if (length(db_nam)>1){
    for (k in 1:length(db_ind)){
      cr1=CR_ind$CR[CR_ind$SampleName==db_ind[k]]
      cr2=CR_ind$CR[CR_ind$SampleName==db_nam[k]]
      if (cr1>cr2){
        new_data = new_data[-which(rownames(new_data)==db_nam[k]),] # suppress
        rownames(new_data)[which(rownames(new_data)==db_ind[k])] = db_nam[k] # rename because there is a _2 or else
      } else {
        new_data = new_data[-which(rownames(new_data)==db_ind[k]),] # just suppress, the good name is still here
      }
    }
  }
  if (marker_type=="SNP"){
    df_SNP = data.frame(MarkerName=marker_name,MAF=MAF,CR=CR,toKeep=TRUE) # toKeep TRUE -> can't know what filter user wants
  } else { # marker_type == 'microsat'
    df_SNP = data.frame(MarkerName=marker_name,CR=CR,toKeep=TRUE)
  }

  return(list(new_data,df_SNP))
}
