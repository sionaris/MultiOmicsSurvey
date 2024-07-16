# Automated hyperparameter texts for algorithms

# SNF
if (algorithm == "SNF") {
  # \mu is used in-text instead of \sigma for consistence with description text
  hyperparameter_text = paste0("There are three hyperparameters in SNF: ",
                               "the number of neighbors ($N_i$),", 
                               " the regularization parameter $\sigma$ and ",
                               "the number of iterations.",
                               "Here, we used a range of number of neighbors",
                               " $\[", hyperparameters$num_neighbors_min,
                               ", ", hyperparameters$num_neighbors_max, "\]$,",
                               " $\mu$ values ranging from ",
                               hyperparameters$sigma_min,
                               " to ", hyperparameters$sigma_max, " with a step size of ",
                               hyperparameters$sigma_step, " and a number of ",
                               "iterations $t = `n_iterations`$.",
                               "Results for larger $t$ are usually identical,",
                               " because the algorithm usually converges after roughly 20 iterations.",
                               "\n",
                               conclusion, "\n",
                               "The optimal hyperparameters determined based on ",
                               "`r params$criterion are: $N_i = ",
                               hyperparameters$optimal_N, "$, ",
                               "$\mu = ",
                               hyperparameters$optimal_sigma, "$.")
}
