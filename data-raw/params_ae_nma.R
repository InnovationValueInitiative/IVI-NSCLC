rm(list = ls())
library("data.table")
library("readxl")
library("iviNSCLC")

# Load data
aes <- c("alt", "ast", "diarrhea", "dry_skin", 
         "eye_problems", "paronychia", "pneumonitis", "pruritus",
         "rash", "stomatitis")
n_aes <- length(aes)

## Check that files exisst
for (i in 1:n_aes){
  file <- paste0("ae_nma/", aes[i], ".xlsx")
  if (!file.exists(file)){
    stop("The file ", paste0("ae_nma/", aes[i], ".xlsx"), " does not exist. ",
         "You must conduct a NMA and produced the required file.")
  }
}

## Create posterior probabilities by treatment
posterior <- vector(mode = "list", length = n_aes)
names(posterior) <- aes
for (i in 1:n_aes){
  file <- paste0("ae_nma/", aes[i], ".xlsx")
  print(paste0("Loading file ", file))
  lookup_i <- data.table(read_excel(file, sheet = "Lookup"))
  posterior_i <- as.matrix(read_excel(file, sheet = "Random Effects Coda"))
  prob_cols <- colnames(posterior_i)[grep("T\\[", colnames(posterior_i))]
  posterior_i <- posterior_i[, prob_cols]
  for (j in 1:nrow(lookup_i)){
    old_var <- paste0("T[", lookup_i$num[j], "]")
    new_var <- paste0("prob_", lookup_i$abb)[j]
    colnames(posterior_i)[colnames(posterior_i) == old_var] <- new_var
  }
  posterior[[i]] <- posterior_i
  attr(posterior[[i]], "lookup") <- lookup_i[, .(name, abb)]
}

## Impute missing results
## Possible 1st and second line treatments
first_line_name <- tx_1L()
first_line_abb <- treatments$tx_abb[match(first_line_name, treatments$tx_name)]
second_line_name <- list()
for (i in 1:length(first_line_name)){
  second_line_name[[i]] <- unlist(tx_2L(first_line_name[i]))
}
second_line_name <- unique(unlist(second_line_name))
second_line_abb <- treatments$tx_abb[match(second_line_name, treatments$tx_name)]

## Imputation
params_ae_nma <- vector(mode = "list", length = n_aes)
names(params_ae_nma) <- aes
params_ae_nma_cols <- paste0("prob_", unique(c(first_line_abb, second_line_abb)))
tki_cols <- paste0("prob_", first_line_abb)
for (i in 1:length(params_ae_nma)){
  params_ae_nma[[i]] <- matrix(NA, nrow = nrow(posterior[[i]]), 
                          ncol = length(params_ae_nma_cols))
  colnames(params_ae_nma[[i]]) <- params_ae_nma_cols
  
  ### Add attributes for imputed 
  attr(params_ae_nma[[i]], "imputed") <- rep(0, length(params_ae_nma_cols))
  
  ### TKI with most adverse events
  tki_cols_subset <- tki_cols[tki_cols %in% colnames(posterior[[i]])]
  tki_means <- apply(posterior[[i]][, tki_cols_subset], 2, mean)
  tki_max_col <- names(which.max(tki_means))
  
  ### Do the imputation
  for (j in 1:ncol(params_ae_nma[[i]])){
    col_name <- colnames(params_ae_nma[[i]])[j]
    if (col_name %in% colnames(posterior[[i]])){
     params_ae_nma[[i]][, col_name] <- posterior[[i]][, col_name] 
    } else{
      attr(params_ae_nma[[i]], "imputed")[j] <- 1
      if (col_name %in% tki_cols){
        params_ae_nma[[i]][, col_name] <- posterior[[i]][, tki_max_col] 
      } else if (col_name %in% c("prob_pbdc", "prob_pbdc_bev")) {
        if ("prob_pbdc_bev" %in% colnames(posterior[[i]])){
          params_ae_nma[[i]][, col_name] <- posterior[[i]][, "prob_pbdc"]
        } else{
          params_ae_nma[[i]][, col_name] <- posterior[[i]][, tki_max_col] 
        }
      }
    }
  } # End loop through columns
} # End loop through adverse events

save(params_ae_nma, file = "../data/params_ae_nma.rda", compress = "bzip2")

