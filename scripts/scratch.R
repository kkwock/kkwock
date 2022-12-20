find_p <- function(f, n, power){
  prob <- 0.5
  while ( pbinom(f, n, prob) > 0.8 ){
    prob = prob + 0.0005
  }
  return(prob)
}

power_calc <- function(runs, situation){
  limits <- exactci::binom.exact(runs, runs, p=0.95, "greater")
  probt <- ifelse(situation=="pass", find_p(runs-1, runs, 0.8), 1- find_p(runs-1, runs, 0.8))
  powert <- find_p(runs-1, runs, 0.8)
  return(c("Probability:Acceptance=",limits$conf.int[1],
           "Probability:Rejection=", probt, 
           "Power Level=", pbinom(runs-1, runs, powert)))
}

## Example of original second paragraph, out of 9 runs the fail situation was considered.
power_calc(6, 'pass')

[1] "Probability:Acceptance=" "0.606962231002917"       
"Probability:Rejection="  "0.764999999999971"       
"Power Level="            "0.79956728426678"

## Example of original second paragraph, out of 9 runs the fail situation was considered.
power_calc(9, 'fail')

[1] "Probability:Acceptance=" "0.716871164436886"       
"Probability:Rejection="  "0.163500000000037"       
"Power Level="            "0.799463464521446"








# CCTD Workshop

library(tidyverse) 
library(plotly)

# Find my flowcell

library(RCurl)
library(jsonlite)
myjson <- RCurl::getURL("https://wtf.ghdna.io/flowcell/output/FCID", ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
jsonlite::fromJSON(myjson)

# Read the tsv file
data_gb01 <- read.csv('/Users/kkwock/OneDrive - Guardant Health/Data/Guardbanding/H2YKTBGXK_autoqc_sample_qc.hdr.tsv')


find_input <- function(output_link){
  split <- unlist(strsplit(output_link , '/'))
  id <- split[length(split)-1]
  link <- paste0("bifs.ghdna.io/ghds/apps/rundeck/Reagent_QC/DEV/G360/EIO_iQC/", id, "/input/")
  return(link)
}
# csv file for percent reads
csv <- "/ghds/groups/reagent_qc/g360/eio_qc/output/1872917/HMVN2BGXM_0402392554/HMVN2BGXM_0402392554_g360-eioqc_perc-reads.csv"
perc_heat <- function(csv){
  read.csv(csv)
}

samples_required_for_testing <- function(n_failures_accepted = c(0,1,2),
                                         total_number_of_tests = c(1, 3, 15),
                                         reliabilities = c(0.9,0.95,0.99),
                                         confidence = c(0.90, 0.95)){  
  
  # -------------------------------------------
  # n_failures_accepted: Number of failures to be tolerated to obtain desired confidence/reliability
  # total_number_of_tests: Total number of evaluations to be performed...used for multiplicity correction
  # reliabilities: vector of desired reliability claims
  # confidence: vector of desired confidence claims
  # -------------------------------------------
  
  
  require(tidyverse)
  require(dplyr)
  
  get_number_of_samples_to_test <- function(n_builds_per_lot = 1000,
                                            n_lots = 3,
                                            reliability = 0.99,
                                            prior_shape1 = 199,
                                            prior_shape2 = 1,
                                            confidence = 0.95,
                                            n_failures_tolerated = 0,
                                            bonferroni_correction_factor = 1, 
                                            iter_value = 1E-6, 
                                            standardize_to_2_obs = TRUE){
    
    
    total <- prior_shape1 + prior_shape2
    if(standardize_to_2_obs){
      prior_shape1 <- 2*prior_shape1 / total
      prior_shape2 <- 2*prior_shape2 / total
    }
    alpha <- 1 - confidence
    posterior_density <- function(n){
      return(pbeta(reliability, 
                   shape1 = prior_shape1 + n - n_failures_tolerated, 
                   shape2 = n_failures_tolerated + prior_shape2, 
                   lower.tail = FALSE))
    }
    
    out_n <- seq_along(1:n_builds_per_lot)
    out_post <- sapply(out_n, posterior_density)
    return(out_n[out_post >= (1 - alpha/bonferroni_correction_factor)][1])
  }
  
  grid <- expand.grid(n_failures_accepted, total_number_of_tests, reliabilities, 1, confidence)
  colnames(grid) <- c("n_failures_accepted", "total_number_of_tests", "reliability", "prior_successes", "confidence")
  
  samples_required <- apply(grid, 1, function(i){
    get_number_of_samples_to_test(prior_shape1 = as.numeric(i[4]), 
                                  prior_shape2 = 1,
                                  n_builds_per_lot = 1000, 
                                  reliability = as.numeric(i[3]), 
                                  confidence = i[5], 
                                  n_failures_tolerated = as.numeric(i[1]),
                                  n_lots = 3,
                                  bonferroni_correction_factor = as.numeric(i[2]))
  })
  
  grid <- grid %>% mutate(samples_required = samples_required) %>%
    select(-prior_successes) %>%
    arrange(confidence, reliability, n_failures_accepted)
  return(grid)  
}

# run the function
samples_required_for_testing()