#' Format log_file
#'
#' Gives names to columns and formats variables as numeric
#'
#' @param log_file the log_file
#'
#' @return log_file updated
#'
#' @keywords internal
#' @noRd

format_logfile = function(log_file){
  colnames(log_file) = c('offspring', 'sire_1', 'dam_1', 'probability_1', 'mismatch_1',
                         'sire_2', 'dam_2', 'probability_2', 'mismatch_2',
                         'sire_3', 'dam_3', 'probability_3', 'mismatch_3',
                         'delta_1_2', 'delta_2_3')
  log_file$probability_1  = as.numeric(log_file$probability_1)
  log_file$probability_2  = as.numeric(log_file$probability_2)
  log_file$probability_3  = as.numeric(log_file$probability_3)
  log_file$delta_1_2      = as.numeric(log_file$delta_1_2)
  log_file$delta_2_3      = as.numeric(log_file$delta_2_3)

  log_file$mismatch_1 = as.integer(log_file$mismatch_1)
  log_file$mismatch_2 = as.integer(log_file$mismatch_2)
  log_file$mismatch_3 = as.integer(log_file$mismatch_3)
  return(log_file)
}
