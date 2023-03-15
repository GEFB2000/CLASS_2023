# import peaks/consensus_from_reduced
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

# usc formatting
ucsc_formating <- function(consensusFilePath = consensusFilePath, export_path = export_path) {
  
  consensus_file_list <- list.files(consensusFilePath, full.names = T, pattern = ".bed")
  
  dbps <- sapply(consensus_file_list, function(x) {
    y <- str_extract(x, "([^\\/]+$)")
    unlist(strsplit(y, "_"))[1]})
  
  peaks <- lapply(consensus_file_list, read.table, col.names = c("chr", "start", "end", "name", "score", "strand"))
  names(peaks) <- dbps
  print(length(peaks))
  canonical_chr <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")
  peaks <- lapply(peaks, function(x) x %>% filter(chr %in% canonical_chr))
  
  headers <- paste0("track type=bedGraph name=", names(peaks))
  new_filenames <- paste0("analysis/00_consensus_peaks/test/", names(peaks), ".bed")
  
  for(i in 1:length(peaks)) {
    # Write the header line
    writeLines(headers[[i]], new_filenames[[i]])
    # Append the broadPeak table data
    
    write.table(peaks[[i]], new_filenames[[i]],
                sep = "\t", col.names = FALSE, row.names = FALSE,
                quote = FALSE, append = TRUE)
  }
  
  return(c("done?"))
}