#  3-run_and_extract_LDhat.R
# version 1.6
# 12 July 2024
# Zane Libke
# CMEE Masters 2023-24 Imperial College London

# script containing function to run ldhat
library(parallel)

#PARSE commandline args (POP_NAME)
args <- commandArgs(trailingOnly = TRUE)
POP_NAME <<- args[1]
#METHOD <- args[2] work on this later - would be nice to be able to specify method and use different methods?

#GLOBAL PARAMETERS
window_width <<- as.integer(args[2])
sites_per_window <<- as.integer(args[3])
chunk <<- as.integer(args[4])
windows_per_chunk <<- as.integer(args[5])
n <<- as.integer(args[6])
randsamp <<- as.logical(args[7])
numcores <<- detectCores() - 2
parameters <- list(window_width=window_width, sites_per_window=sites_per_window, chunk=chunk, windows_per_chunk=windows_per_chunk, n=n, randsamp=randsamp, POP_NAME=POP_NAME, numcores=numcores)

#establish names for rdata files produced by 2-convert 
chr_3R <-paste("../data/3R_", POP_NAME, ".Rdata", sep = "")
chr_3L <-paste("../data/3L_", POP_NAME, ".Rdata", sep = "")
chr_2R <-paste("../data/2R_", POP_NAME, ".Rdata", sep = "")
chr_2L <-paste("../data/2L_", POP_NAME, ".Rdata", sep = "")
chr_X <-paste("../data/X_", POP_NAME, ".Rdata", sep = "")

#run LDhat on all windows of specified size
LDhat_runall <- function(iterator) {
  print(paste("RUNNING CHUNK:", iterator))
  # THE UNIQUE window IDs RUN BY THIS NODE
  window<-(iterator-1)*chunk
  window<-window+(1:chunk)
  window <- window[window %in% above_spw_IDs]
  len_window <- as.integer(length(window))
  if (len_window<1){return()} #if no windows reach the wpc threshold, return
  if (len_window < windows_per_chunk){ #if num windows above threshold < wpc, wpc <- num windows
    wpc <- len_window
  } else {wpc <- windows_per_chunk} #run the same amnt as inputted to windows_per_chunk param
  
  # if random sample mode on, shuffle window, but only if above wpc (also handle edge case of windows_per_chunk = 1)
  if (randsamp == TRUE){
    if(len_window>wpc && wpc>1){
      window <- sample(window)
    }
  }
  
  #init iterator_results_df
  iterator_results_df <- data.frame(
    POP_NAME=rep(NA_character_, chunk),
    chr=rep(NA_character_, chunk), 
    window_start=rep(NA_real_, chunk),
    window_end=rep(NA_real_, chunk),
    rho=rep(NA_real_, chunk)
  )
  LDprefix <- paste(chr,"_",POP_NAME,"_",sep="")
  
  # For loop over wpc windows that passed threshold 
  for (i in 1:wpc)
  {
    # 1. SET WINDOW BOUNDS, GET SITES, INCLUDE BOTH BOUNDS
    window_start<-(window[i]-1)*window_width+1
    window_end<-window_start+(window_width-1)
    
    # HOW MANY SITES WITHIN THIS WINDOW
    temp<-(POS>=window_start) & (POS<window_end)
    
    # CHECK FOR >20 SNPS BEFORE FURTHER PROCESSING
    if (!is.na(sum(temp)) && sum(temp) >= sites_per_window){
      # PULL SUBSET WITHIN THE WINDOW
      working_POS<-POS[temp]
      working_genotype<-genotype[, temp]
      
      # 0 AND 1 ARE HOMOZYGOTES
      # 2 IS HET
      # ? IS MISSING
      y<-working_genotype
      y[working_genotype==0]<-0
      y[working_genotype==2]<-1
      y[working_genotype==1]<-2
      y[working_genotype<0]<-'?'
      
      # PASTE STRINGS
      f<-function(x) {paste(x, collapse='')}
      y<-apply(y, 1, f)
      
      # 2. WRITE _loc AND _site FILES
      # WRITE THE site FILE
      site_filename<-paste(c("../data/input/", LDprefix, 'site_', window_start/1000, '_', window_end/1000, '.txt'), collapse='')
      site_first_line<-paste(c(n, length(working_POS), 2), collapse=' ')
      write(site_first_line, site_filename)
      # FOR THE FIRST 50 INDIVIDUALS
      for (j in 1:n)
      {
        write(paste(c('>', sample_id[j]), collapse=''), site_filename, append=T)
        write(y[j], site_filename, append=T)
      }
      # WRITE THE loc FILE
      loc_filename<-paste(c("../data/input/", LDprefix,'loc_', window_start/1000, '_', window_end/1000, '.txt'), collapse='')
      loc_first_line<-paste(length(working_POS), (window_end-window_start+1)/1000, 'L', collapse=' ')
      write(loc_first_line, loc_filename)
      for (j in 1:length(working_POS))
      {
        write((working_POS[j]-window_start)/1000, loc_filename, append=T)
      }
      # OUTPUT FILENAME
      output_prefix<-paste(c("../data/output/", LDprefix, window_start/1000, '_', window_end/1000, 'kb'), collapse='')
      # 4. RUN LD HAT, ONLY IF THERE ARE SOME SNPS
      #print(paste("running LDhat on window ", window_start/1000, "-", window_end/1000, "kb", sep=""))
      command<-paste("./pairwise -seq", site_filename, "-loc", loc_filename, "-prefix", output_prefix,"-lk lk_n100_t0.001")
      print(window_start)
      try(system(command, intern=T, ignore.stdout=T, ignore.stderr=T))
      
      #extract 4Ner and store in list with other params
      outfile_name <- paste(output_prefix, "outfile.txt", sep="")
      
      if(file.exists(outfile_name)){
        #Extract value of rho
        temp<-readLines(outfile_name)[5]
        x1<-gregexpr('= ', temp)[[1]][1]
        x2<-gregexpr(' :', temp)[[1]][1]
        rho <- as.numeric(substr(temp, x1+2, x2))
        
        #add values to their respective columns in the iterator results df
        iterator_results_df$POP_NAME[i] <- POP_NAME
        iterator_results_df$chr[i] <- chr
        iterator_results_df$window_start[i] <- window_start
        iterator_results_df$window_end[i] <- window_end
        iterator_results_df$rho[i] <- rho
      }
    }
  }
  chunk_results <- na.omit(iterator_results_df) #remove empty columns from iterator results
  attr(chunk_results, "na.action") <- NULL #remove attribute left by na.omit
  #attr(chunk_results, "names") <- NULL #remove indexes
  return(chunk_results);
}

#wrapper function for LDhat run all windows
wrapper_LDhat_runall <- function(window_width=window_width, sites_per_window=sites_per_window, chunk=chunk, windows_per_chunk=windows_per_chunk, POP_NAME=POP_NAME, randsamp=randsamp, chr) {
  chr <<- chr
  parameters <- list(window_width=window_width, sites_per_window=sites_per_window, chunk=chunk, windows_per_chunk=windows_per_chunk, n=n, randsamp=randsamp, POP_NAME=POP_NAME)
  out_object_name <- paste0("LDhat_res_", chr, "_", POP_NAME)
  num_window<-ceiling(max(POS)/window_width) #total num of windows (measured from last variable site)
  chunks <- 1:floor(num_window/chunk) #num of function calls - round down for now
  ########################################################### FILTER for SITES_PER_WINDOW FUNCTION ###
  scan_windows <- function(POS, window_width, windows_per_chunk){
    window_lb<-seq(1, max(POS), window_width) #generate all window starts
    window_ub<-window_lb+window_width-1 #generate all window ends
    all_windows<-rep(0, length(window_lb)) #generate empty vector length of total amount of windows
    
    bin<-POS%/%window_width+1 #total number of windows available
    
    num_sites<-by(POS, bin, length) #count number of polymorphic sites in window
    sites_idx<-as.numeric(names(num_sites)) #get index of window
    sites_value<-as.numeric(num_sites) #get amnt of polymorphic sites in window with 
    
    all_windows[sites_idx]<-sites_value
    all_windows<-data.frame(window_ID=seq(1, length(window_lb), 1), window_lb=window_lb, window_ub=window_ub, poly_sites=all_windows)
    above_spw <- subset(all_windows, poly_sites >= 20)
    return(above_spw)
  }
  
  above_spw <- scan_windows(POS, window_width, windows_per_chunk)
  above_spw_IDs <<- above_spw$window_ID
  ########################################################### FILTER for SITES_PER_WINDOW FUNCTION ###
  
  #print parameters
  print(paste("POP_NAME: ", POP_NAME, "; chr: ", chr, sep=""))
  print(paste("window_width: ", window_width, "; sites_per_window: ", sites_per_window, "; chunk: ", chunk, "; numcores: ", numcores, sep=""))
  print(paste("Randomly sampling ", windows_per_chunk, " ", window_width, "bp windows from ", chunk*window_width, "bp chunks (", max(chunks), " total chunks)", sep=""))
  #run LDhat_runall on all chunks, save each out df in a list
  results_list <- mclapply(chunks, safe_LDhat_runall, mc.cores = numcores)
  
  # Check and handle errors
  errors <- sapply(results_list, function(res) is.list(res) && !is.null(res$error) && res$error)
  
  # Print error messages if any
  if (any(errors)) {
    cat("Errors occurred in the following chunks:\n")
    error_messages <- lapply(results_list[errors], function(res) paste("Chunk:", res$iterator, "- Error:", res$message))
    print(do.call(rbind, error_messages))
    combined_df <- NULL
  } else {
    # Combine the results if no errors
    combined_df <- do.call(rbind, results_list[!errors])
  }
  #rbind all dfs together and give sample specific name
  #combined_df <- do.call(rbind, results_list)
  return(list(above_spw = above_spw, combined_df = combined_df))
}

#Error catching wrapper for LDhat_runall
safe_LDhat_runall <- function(iterator) {
  tryCatch({
    LDhat_runall(iterator)
  }, error = function(e) {
    return(list(error = TRUE, message = e$message, iterator = iterator))
  })
}

################# chromosome 3R ###################
# LOAD IN DATA 
load(chr_3R)
res_3R <- wrapper_LDhat_runall(window_width = window_width, sites_per_window = sites_per_window, chunk = chunk, windows_per_chunk = windows_per_chunk, POP_NAME = POP_NAME, chr = "3R", randsamp = randsamp)

# Clean up large objects
rm(POS, genotype, sample_id, chr)
gc()

################# chromosome 3L ###################
# LOAD IN DATA 
load(chr_3L)
res_3L <- wrapper_LDhat_runall(window_width = window_width, sites_per_window = sites_per_window, chunk = chunk, windows_per_chunk = windows_per_chunk, POP_NAME = POP_NAME, chr = "3L", randsamp = randsamp)

# Clean up large objects
rm(POS, genotype, sample_id, chr)
gc()

################# chromosome 2R ###################
# LOAD IN DATA
load(chr_2R)
res_2R <- wrapper_LDhat_runall(window_width = window_width, sites_per_window = sites_per_window, chunk = chunk, windows_per_chunk = windows_per_chunk, POP_NAME = POP_NAME, chr = "2R", randsamp = randsamp)

# Clean up large objects
rm(POS, genotype, sample_id, chr)
gc()

################# chromosome 2L ###################
# LOAD IN DATA
load(chr_2L)
res_2L <- wrapper_LDhat_runall(window_width = window_width, sites_per_window = sites_per_window, chunk = chunk, windows_per_chunk = windows_per_chunk, POP_NAME = POP_NAME, chr = "2L", randsamp = randsamp)

# Clean up large objects
rm(POS, genotype, sample_id, chr)
gc()

################# chromosome X ###################
# LOAD IN DATA
load(chr_X)
res_X <- wrapper_LDhat_runall(window_width = window_width, sites_per_window = sites_per_window, chunk = chunk, windows_per_chunk = windows_per_chunk, POP_NAME = POP_NAME, chr = "X", randsamp = randsamp)

# Clean up large objects
rm(POS, genotype, sample_id, chr)
gc()

allfile_name <- paste("../results/LDhat_out_", POP_NAME, ".RData", sep="")
results_all <- rbind(res_3R=res_3R$combined_df, res_3L=res_3L$combined_df, res_2R=res_2R$combined_df, res_2L=res_2L$combined_df, res_X=res_X$combined_df)
above_spw_all <- list("3R"=res_3R$above_spw, "3L"=res_3L$above_spw, "2R"=res_2R$above_spw, "2L"=res_2L$above_spw, "X"=res_X$above_spw)
save(results_all, parameters, above_spw_all, file = allfile_name) #add parameters object eventually to save the input parameters too.
