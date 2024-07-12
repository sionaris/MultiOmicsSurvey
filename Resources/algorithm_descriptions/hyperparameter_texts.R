# Automated hyperparameter texts for algorithms

# SNF
if (algorithm == "SNF") {
  hyperparameter_text = paste0("There are three hyperparameters in SNF: ",
                                          "the number of neighbors ($N_i$),", 
                                          " the regularization parameter $\sigma$ and ",
                                          "the number of iterations.",
                                          "Here, we used a range of number of neighbors",
                                          " $\[", hyperparameters$num_neighbors_min,
                                          ", ", hyperparameters$num_neighbors_max, "\]$,",
                                          " $\sigma$ values ranging from ",
                                          hyperparameters$sigma_min,
                                          " to ", hyperparameters$sigma_max, " with a step size of ",
                                          hyperparameters$sigma_step, " and a number of ",
                                          "iterations $t$ ranging from ",
                                          hyperparameters$iter_min, " to ",
                                          hyperparameters$iter_max, " with a step size of ",
                                          hyperparameters$iter_step, 
                                          ". The optimal hyperparameters determined based on ",
                                          "`r params$criterion are: $N_i = ",
                                          hyperparameters$optimal_N, "$, ",
                                          "$\sigma = ",
                                          hyperparameters$optimal_sigma, "$ and $t = ",
                                          hyperparameters$optimal_iter, "$.")
}
