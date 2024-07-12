# A series of functions defined to modularize code and make cleaner

# In a nutshell #####

# Fetch the "in a nutshell" text of an algorithm
fetch_in_a_nutshell = function (algorithm) {
  source("Resources/algorithm_descriptions/in_a_nutshell.R")
  return(desc_list[[algorithm]])
}

# Fetch the citation of an algorithm
fetch_citation = function (algorithm) {
  source("Resources/algorithm_descriptions/citations.R")
  return(citations[[algorithm]])
}

# Generate hyperparameters report text
generate_hyperparams_text = function (algorithm, hyperparameters) {
  source("Resources/algorithm_descriptions/hyperparameter_texts.R")
  return(hyperparameters_text)
}