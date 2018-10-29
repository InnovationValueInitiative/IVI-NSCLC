rm(list = ls())
ae <- rbind(c("Elevated alanine transaminase", "alt"),
            c("Elevated aspartate transaminase", "ast"),
            c("Diarrhea", "diarrhea"),
            c("Dry skin", "dry_skin"),
            c("Eye problems", "eye_problems"),
            c("Paronychia", "paronychia"),
            c("Pneumonitis", "pneumonitis"),
            c("Pruritus", "pruritus"),
            c("Rash", "rash"),
            c("Stomatitis", "stomatitis"))
colnames(ae) <- c("ae_name", "ae_abb")            
adverse_events <- data.table(ae)
save(adverse_events, file = "../data/adverse_events.rda", compress = "bzip2")