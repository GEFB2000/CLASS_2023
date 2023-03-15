import_peaks <- function(consensus_file_path = broadpeakfilepath) {
  peak_files <- list.files(consensus_file_path, full.names = T, pattern = ".broadPeak")
  dbp_name <- sapply(peak_files, function(x){
    y <-  str_extract(x, "([^\\/]+$)")
    # NOTE ALTERNATIVE Gsub
    # gsub("_peaks.broadPeak", "", y)
    paste(unlist(strsplit(y, "_"))[c(1,2)], collapse = "_") 
    
  })
  peak_list <- c()
  
  # the for loop !
  for(i in 1:length(peak_files)) {
    # Import peaks
    peaks <- rtracklayer::import(peak_files[i])
    # Append this GRanges object to the of the list.
    peak_list <- c(peak_list, peaks)
    # Name the list elements by their TF name.
    names(peak_list)[length(peak_list)] <- dbp_name[i]
  }
  return(peak_list)
}

consensus_from_reduced <- function(dbp, peak_list) {
  dbp_peaks <- peak_list[grepl(as.character(dbp), names(peak_list))]
  suppressWarnings(all_peaks <- GenomicRanges::reduce(unlist(as(dbp_peaks, "GRangesList"))))
  all_peaks <- all_peaks[grepl("chr", seqnames(all_peaks))]
  
  # peak_exists <- lapply(dbp_peaks, function(x) {
  #   as.numeric(countOverlaps(all_peaks, x) > 0))
  # }) %>%
  # bind_rows() OR bind_cols()
  peak_exists <- matrix(NA, nrow = length(all_peaks), ncol = length(dbp_peaks))
  for(i in 1:length(dbp_peaks)) {
    suppressWarnings(peak_exists[,i] <- as.numeric(countOverlaps(all_peaks, dbp_peaks[[i]]) > 0))
  }
  # filter to consensus requiring peaks to be in all replicates
  dbp_consensus <- all_peaks[rowSums(peak_exists) == ncol(peak_exists)]
  # Required only two replicates == dbp_consensus <- all_peaks[rowSums(peak_exists) > 1]
  return(dbp_consensus)
}