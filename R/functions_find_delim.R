#' Find delimiter
#'
#' @param string a string where it is needed to find delimiter (file import for example)
#'
#' @return the delimiter
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' string = c("What;is;my;delimiter?")
#' delim = find_delim(string)
#' delim
#'
find_delim <- function(string) {
  non_alnum <- gsub("[[:alnum:]]", "", string)
  non_alnum_count <- table(strsplit(non_alnum, "")[[1]])
  if (length(non_alnum_count)>0){
    df <- data.frame(Charac = names(non_alnum_count), Occu = as.numeric(non_alnum_count), stringsAsFactors = FALSE)
    poss_charac = c(" ",";",",","\t","  ")
    df=df[df$Charac %in% poss_charac,]
    if (nrow(df)>0){
      del=df$Charac[which(df$Occu==max(df$Occu))]
      return(del)
    } else {
      return("") #default for read.table
    }
  } else {
    return("") #default for read.table
  }
}
