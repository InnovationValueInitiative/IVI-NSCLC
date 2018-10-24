rm(list = ls())
library("data.table")
library("readxl")

# Load data
aes <- c("alt_increase", "ast_increase", "diarrhea", "dry_skin", 
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
params_ae_nma <- vector(mode = "list", length = n_aes)
names(params_ae_nma) <- aes
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
  params_ae_nma[[i]] <- posterior_i
  attr(params_ae_nma[[i]], "lookup") <- lookup_i[, .(name, abb)]
}

save(params_ae_nma, file = "../data/params_ae_nma.rda", compress = "bzip2")

