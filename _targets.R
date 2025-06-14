library(targets)
library(tarchetypes)
library(here)
library(furrr)

#-------------------------------------------------------------------------------
#                                     SETUP
#-------------------------------------------------------------------------------
tar_source(files = "./src") # compile all the script in the src directory

plan(multisession, workers = parallel::detectCores() - 1)

# Package's loading
tar_option_set(packages = c("graphicalExtremes",
                            "fMultivar",
                            "tibble",
                            "tidyr",
                            "dplyr",
                            "ggplot2",
                            "here",
                            "rmdformats",
                            "hrbrthemes",
                            "gridExtra",
                            "igraph",
                            "tidygraph",
                            "ggraph",
                            "evd",
                            "pracma",
                            "furrr",
                            "cowplot"))

#-------------------------------------------------------------------------------
#                                   PIPELINE
#-------------------------------------------------------------------------------
list(
  #------------------------ Simulation study trivariate ------------------------
  # Computation of a the trivariate coefficient in the Husler-Reiss case
  tar_target(
    Gamma_cond_indep_12,
    random_Gamma12(5)          # Random valid HR variogram matrix 3x3
  ),

  tar_target(
    chi_trivariate_example,
    chi_trivariate_HR(Gamma_cond_indep_12)    # Computation of the coefficient
  ),
  # not zero in this case but it was obvious unfortunately...

  # Repetition of the variogram simulation 2000 times (= batches * reps)
  tar_rep(
    Gamma_cond_indep_12_replicate,
    random_Gamma12(25),
    batches = 1,
    reps = 2000,
    iteration = "vector"               # the output of the replicates is a list
  ),

  # Computation of the trivariate coefficient and all chi for all
  # the simulations
  tar_target(
    chi_trivariate_simulation_results,
    data.frame(num_sim = seq_len(nrow(Gamma_cond_indep_12_replicate)),
               Gamma_cond_indep_12_replicate,
               row.names = NULL) |>
      rowwise() |>
      mutate(tri_chi = chi_trivariate_HR(c(X1, X2, X3)),
             chi_12 = 2 - theta(X1),
             chi_13 = 2 - theta(X2),
             chi_23 = 2 - theta(X3),
             prod_chi = (2 - theta(X2)) * (2 - theta(X3)))
  ),

  # Computation of the conditional coefficient for all the simulation
  tar_target(
    cond_trivariate_simulation_results,
    data.frame(num_sim = seq_len(nrow(Gamma_cond_indep_12_replicate)),
               Gamma_cond_indep_12_replicate,
               row.names = NULL) |>
      rowwise() |>
      mutate(tri_chi = cond_trivariate_HR(c(X1, X2, X3), 3))
  ),

  # All the possible coefficient for the matrix)
  tar_target(
    all_matrix_coeff,
    expand.grid(b = seq(0.1, 10, 0.05), c = seq(0.1, 10, 0.05)) |>
      as_tibble() |>
      mutate(a = b + c) |>
      rowwise() |>
      mutate(tri_chi = chi_trivariate_HR(c(a, b, c)),
             chi_12 = 2 - theta(a),
             chi_13 = 2 - theta(b),
             chi_23 = 2 - theta(c),
             prod_chi = (2 - theta(b)) * (2 - theta(c))) |>
      mutate(a10 = b + 10 * c) |>
      rowwise() |>
      mutate(tri_chi10 = chi_trivariate_HR(c(a10, b, c)),
             chi10_12 = 2 - theta(a10))
  ),

  # Resulting plots
  # Plot of conditional values for several simulation
  tar_target(
    cond_trivariate_plot,
    {
      p1 <- cond_trivariate_simulation_results |>
        ggplot() + aes(x = num_sim, y = tri_chi) +
        geom_line(col = "darkorange") + coord_cartesian(ylim = c(0, 1)) +
        labs(title = "Values", x = "Simulation number",
             y = expression(paste(chi[12]^3, " value"))) +
        theme_ipsum(base_family = "serif")

      p2 <- cond_trivariate_simulation_results |>
        ggplot() + aes(x = tri_chi, y = 0) +
        geom_violin(col = "#474F58", fill = "darkorange") +
        coord_cartesian(xlim = c(0, 1)) +
        labs(title = "Violin plot",  y = " ",
             x = expression(paste(chi[12]^3, " value"))) +
        theme_ipsum(base_family = "serif")

      grid.arrange(grobs = list(p1, p2), ncol = 2)
    }
  ),

  # Plot of trivariate values for several simulation
  tar_target(
    chi_trivariate_plot,
    {
      p1 <- chi_trivariate_simulation_results |>
        ggplot() + aes(x = num_sim, y = tri_chi) +
        geom_line(col = "darkorange") + coord_cartesian(ylim = c(0, 1)) +
        labs(title = "Values", x = "Simulation number",
             y = expression(paste(chi[123], " value"))) +
        theme_ipsum(base_family = "serif")

      p2 <- chi_trivariate_simulation_results |>
        ggplot() + aes(x = tri_chi, y = 0) +
        geom_violin(col = "#474F58", fill = "darkorange") +
        coord_cartesian(xlim = c(0, 1)) +
        labs(title = "Violin plot",  y = " ",
             x = expression(paste(chi[123], " value"))) +
        theme_ipsum(base_family = "serif")

      grid.arrange(grobs = list(p1, p2), ncol = 2)
    }
  ),

  # Plot of quotient values for several simulation
  tar_target(
    quotient_chi_plot,
    chi_trivariate_simulation_results |>
      ggplot() + aes(x = num_sim, y = tri_chi / prod_chi) +
      geom_line(col = "darkorange") +
      geom_hline(yintercept = 1, linetype = 2,   # expected value
                 linewidth = 1, colour = "#474F58") +
      labs(x = "Simulation number",
           y = expression(paste(chi[123], " / ", chi[13], chi[23])),
           title = expression(paste("Quotient of ", chi))) +
      theme_ipsum(base_family = "serif")
  ),
  # Quotient plot and relative error
  tar_target(
    error_triprod_chi_plot,
    chi_trivariate_simulation_results |>
      ggplot() + aes(x = num_sim) +
      geom_line(aes(y = abs(tri_chi - prod_chi) / tri_chi),
                col = "darkorange") +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Simulation number",
           y = "Error",
           title = "Relative error between the coefficients") +
      theme_ipsum(base_family = "serif")
  ),

  # Heatmap of the quotient tri_chi/prod_chi over a grid of all possible
  # value for b and c in the variogram matrix when 1 indep 2 | 3
  tar_target(
    plot_all_matrix_coef,
    all_matrix_coeff |>
      ggplot() + aes(x = b, y = c, fill = tri_chi / prod_chi) +
      geom_tile() +
      scale_fill_continuous(name = expression(paste(chi[123],
                                                    "/", chi[13], chi[23])),
                            type = "viridis") +
      labs(x = expression(Gamma[13]), y = expression(Gamma[23]),
           caption = expression(paste("The conditional independance is
             knowing ",
                                      Y[3], "."))) +
      theme_ipsum(base_family = "serif", grid = FALSE, axis_title_size = 12)
  ),

  # Error of the quotient tri_chi/prod_chi over a grid of all possible value for
  # b and c in the variogram matrix when 1 indep 2 | 3,
  # but the expecting results was a mistake
  tar_target(
    plot_all_matrix_error,
    {
      p1 <- all_matrix_coeff |>
        ggplot() +
        aes(x = b, y = c, fill = abs(tri_chi - prod_chi) / tri_chi) +
        geom_tile() +
        scale_fill_continuous(name = expression(paste("|", chi[123], "-",
                                                      chi[13], chi[23], "|",
                                                      "/",
                                                      chi[13], chi[23])),
                              type = "viridis") +
        labs(x = expression(Gamma[13]), y = expression(Gamma[23]),
             title = "Relative error") +
        theme_ipsum(base_family = "serif", grid = FALSE, axis_title_size = 12)

      p2 <- all_matrix_coeff |>
        ggplot() + aes(x = b, y = c, fill = abs(tri_chi - prod_chi)) +
        geom_tile() +
        scale_fill_viridis_c(name = expression(paste(
          "|", chi[123], "-", chi[13], chi[23], "|"
        )
        ),
        option = "inferno") +
        labs(x = expression(Gamma[13]), y = expression(Gamma[23]),
             title = "Absolute error",
             caption = expression(paste("The conditional independance is 
             knowing ",
                                        Y[3], "."))) +
        theme_ipsum(base_family = "serif", grid = FALSE, axis_title_size = 12)

      grid.arrange(grobs = list(p1, p2), nrow = 2)
    }
  ),

  # Trying to see a relationship between bivariate chi and trivariate chi using
  # the SECO representation formula when 1 indep 2 | 3
  tar_target(
    plot_chi_over_trichi,
    all_matrix_coeff |>
      rowwise() |>
      filter(
        checkGamma(to_matrix(c(a, b, c)), returnBoolean =  TRUE, alert = FALSE)
      ) |>
      pivot_longer(cols = starts_with("chi_"),
                   values_to = "chi",
                   names_to = "Coefficient",
                   names_prefix = "chi_") |>
      ggplot() +
      aes(x = chi, y = tri_chi, group = Coefficient) +
      geom_point(col = "grey78") +
      geom_abline(aes(color = "y=x", slope = 1, intercept = 0),
                  linewidth = 1, show.legend = TRUE) +
      facet_grid(. ~ Coefficient) +
      theme_ipsum(base_family = "serif") +
      scale_color_manual(values = "darkorange") +
      labs(x = expression(chi[ij]), y = expression(chi[123]),
           color = "legend",
           title = " ",
           caption = expression(paste("The conditional independance is
            knowing ",
                                      Y[3], ".")))
  ),

  # Verification if the linear form comes from conditional independence or the
  # parametric representation
  tar_target(
    plot_chi10_over_trichi10,
    all_matrix_coeff |>
      rename(chi10_13 = chi_13,
             chi10_23 = chi_23) |>
      filter(!is.na(tri_chi10)) |>
      pivot_longer(cols = starts_with("chi10_"),
                   values_to = "chi",
                   names_to = "Coefficient",
                   names_prefix = "chi10_") |>
      ggplot() +
      aes(x = chi, y = tri_chi10, group = Coefficient) +
      geom_point(col = "grey78") +
      geom_abline(aes(color = "y=x", slope = 1, intercept = 0),
                  linewidth = 1, show.legend = TRUE) +
      facet_grid(. ~ Coefficient) +
      theme_ipsum(base_family = "serif") +
      scale_color_manual(values = "darkorange") +
      labs(x = expression(chi[ij]), y = expression(chi[123]),
           color = "legend",
           title = " ",
           caption = "Here, we have the relation a = b + 10c")
  ),

  #-------------------------- Spectral representation --------------------------
  # Set the graph and a corresponding variogram matrix for a Husler-Reiss
  # graphical model (we use the package graphicalExtremes in order to do that)
  tar_target(
    graphical_model_parameters,
    {
      # Graph creation
      g <- make_graph(edges = c(1, 2,
                                1, 3,
                                2, 3,
                                3, 4),
                      n = 4,
                      directed = FALSE)

      Gamma <- generate_random_graphical_Gamma(g)     # corresponding variogram
      Theta <- Gamma2Theta(Gamma)                     # corresponding theta
      Sigma_Gamma <- Theta2Sigma(Theta)               # general inverse of theta
      list(
        graph = g,
        Gamma = Gamma,
        Theta = Theta,
        Sigma_Gamma = Sigma_Gamma
      )
    }
  ),

  # Seeking of a matrix whose the inverse get an null entry in 12 (and thus
  # the graphical model has no edge between node 1 and node 2)
  tar_target(
    alternative_graph_research,
    {
      x <- seq(-5, 5, by = 0.2)                        # set the mesh
      S <- graphical_model_parameters$Sigma_Gamma
      expand.grid(a = x, b = x, c = x) |>
        tibble() |>
        rowwise() |>
        mutate(
          # search the last coefficient where the precision matrix is null
          # at 1-2 entry
          d_zero = newtonRaphson(
            \(.) (solve(S + ker_gamma(c(a, b, c, .))))[1, 2], -1
          )$root
        )
    }
  ),
  # Notes : we use solve here because everyone is invertible but we should use
  # pinv() from pracma to compute the Moore-Penrose inverse (these two fit when
  # existence).

  # Number of the previous matrices that are semi-definite positive
  tar_target(
    number_alternative_graph,
    {
      S <- graphical_model_parameters$Sigma_Gamma
      alternative_graph_research |>
        rowwise() |>
        # check if the matrix is semi-def pos
        filter(semi_def(S + ker_gamma(c(a, b, c, d_zero)))) |>
        ungroup() |>
        count(name = "Number of semi definite positive")
    }
  ),

  #------------- Conditional independence on linear latent model ---------------
  # Setting the simulation parameters
  tar_target(
    latent_10_sim_parameters,
    {
      d <- 10
      K <- 23
      coef <- 20 * runif(5 * d) - 10.  # all the no zero coefficients

      mat <- coef[1:5]
      for (i in 2:d) {
        mat <- c(mat, rep(0, K - 3), coef[((i - 1) * 5 + 1):(5 * i)])
      }
      M <- matrix(mat, byrow = TRUE, nc = K)
      list(d = d,
           K = K,
           M = M,
           n = 1e5,
           gamma = 0.5)
    }
  ),

  # Simulation of the linear latent model with "path type" graphical model
  # according to the first version of our conditional independence
  tar_target(
    latent_10_simulation,
    {
      n <- latent_10_sim_parameters$n
      d <- latent_10_sim_parameters$d
      z_raw <- rfrechet(n * latent_10_sim_parameters$K,
                        shape = latent_10_sim_parameters$gamma)
      z_sample <- matrix(z_raw, nrow = latent_10_sim_parameters$K)
      epsilon <- matrix(rnorm(d * n), nrow = d)
      t(latent_10_sim_parameters$M %*% z_sample + epsilon)
    }
  ),

  # Estimation of the "classic" extremal independence using MST algorithm
  # to build the tree graphical model as done in (Engelke and Volgushev 2022)
  tar_target(
             latent_10_graph_estimation,
             {
               n <- latent_10_sim_parameters$n
               k <- 300
               res <- (matrixStats::colRanks(latent_10_simulation,
                                             ties.method = "max") / n) > 1 - k / n

               df <- - log((res %*% t(res)) / k)
               completeGraph <- graph_from_adjacency_matrix(df,
                                                            mode = "undirected",
                                                            weighted = TRUE)

               mst(completeGraph, algorithm = "prim")
             }),

  #---------------------- Variable Clustering for HR models --------------------
  ##================ First simulation optimization ================
  tar_target(
    first_sim_param_cluster,
    {
      # Construction of clusters and R matrix
      R <- matrix(c(0.5, -2,
                    -2, 1), nc = 2)
      clusters <- list(1:4, 5:7)

      # Construction of induced theta and corresponding variogram gamma
      Theta <- build_theta(R, clusters)
      Gamma <- Theta2Gamma(Theta)

      list(
        R = R,
        clusters = clusters,
        Theta = Theta,
        Gamma = Gamma,
        chi = 1,    # zeta
        n = 2e3,
        d = 7
      )
    }
  ),

  tar_target(
    first_sim_clustering,
    rmpareto(n = first_sim_param_cluster$n,
             model = "HR",
             par = first_sim_param_cluster$Gamma)
  ),

  tar_target(
    first_sim_optimisation_no_penalty,
    best_clusters(data = first_sim_clustering,
                  chi = first_sim_param_cluster$chi,
                  l_grid = 0)
  ),

  tar_target(
    first_sim_optimisation_results,
    future_map(seq(0, 1e-2, 1e-4), \(.) {
      best_clusters(data = first_sim_clustering,
                    chi = first_sim_param_cluster$chi,
                    l_grid = .,
                    it_max = 2000)
    })
  ),

  tar_target(
    first_sim_results,
    all_info(first_sim_param_cluster$clusters,
             first_sim_optimisation_results,
             lambda = seq(0, 1e-2, 1e-4), one_sim = TRUE)
  ),

  tar_target(
    first_sim_plots,
    plot_info(first_sim_results, one_sim = TRUE)
  ),

  tar_target(
    weight_graph_first_sim,
    {
      W <- compute_W(first_sim_clustering)

      # Graph creation
      g <- graph_from_adjacency_matrix(W, mode = "undirected",
                                       weighted = TRUE, diag = FALSE)
      ggraph(g, layout = "circle") +  # 'dh' = Davidson-Harel layout
        geom_edge_link(aes(edge_width = weight,
                           edge_alpha = weight, color = weight),
                       show.legend = c(edge_width = FALSE,
                                       edge_alpha = FALSE, color = TRUE)) +
        geom_node_point(size = 10, color = "grey70") +
        geom_node_text(aes(label = 1:7)) +
        scale_edge_color_gradient(low = "white", high = "darkorange2",
                                  name = "Weight values") +
        scale_edge_width(range = 0.9) +
        theme_void()
    }
  ),

  # Hierarchical clustering graph for the first simulation

  tar_target(
    hclust_first_sim,
    gg_cluster(first_sim_optimisation_results)
  ),

  # Replication of the first simulation
  tar_rep(
    first_sim_clustering_replicate,
    rmpareto(n = first_sim_param_cluster$n,
             model = "HR",
             par = first_sim_param_cluster$Gamma),
    batches = 1,
    reps = 500,
    iteration = "vector"               # the output of the replicates is a list
  ),

  # Data formatting to separate the replications
  tar_target(
    first_sim_rep_data,
    {
      rep <- first_sim_clustering_replicate

      n <- first_sim_param_cluster$n

      row.names(rep) <- NULL

      rep |>
        tibble() |>
        mutate(sim = (seq_len(nrow(rep)) + (n - 1)) %/% n)
    }
  ),

  # Computation of optimization results for each replication
  tar_target(
    first_sim_rep_pen_results,
    {
      m <- first_sim_rep_data |> summarise(max = max(sim))
      lambda <- seq(0, 3e-2, 1e-4)
      future_map(1:m$max, function(i) {
        data <- first_sim_rep_data |>
          filter(sim == i) |>
          select(-sim) |>
          as.matrix()

        future_map(lambda, \(.) {
          best_clusters(data,
                        chi = first_sim_param_cluster$chi,
                        l_grid = .)
        }
        )
      }
      )
    }
  ),

  tar_target(
    first_sim_rep_nopen_results,
    {
      m <- first_sim_rep_data |> summarise(max = max(sim))
      future_map(1:m$max, function(i) {
        data <- first_sim_rep_data |>
          filter(sim == i) |>
          select(-sim) |>
          as.matrix()


        best_clusters(data,
                      chi = first_sim_param_cluster$chi,
                      l_grid = 0)

      })
    }
  ),

  # Formatting and plot the results
  tar_target(
    first_sim_rep_results,
    {
      lambda <- seq(0, 3e-2, 1e-4)
      all_info(cluster.init = first_sim_param_cluster$clusters,
               list_res = first_sim_rep_pen_results, lambda = lambda)
    }
  ),

  tar_target(
    hclust_first_sim_rep,
    average_hierarchy(first_sim_rep_pen_results)
  ),

  ##======================= Unbalanced class =====================
  tar_target(
    unbal_sim_param_cluster,
    {
      # Construction of clusters and R matrix
      R <- matrix(c(1, -3,
                    -3, 1), nc = 2)
      clusters <- list(1:2, 3:7)

      # Construction of induced theta and corresponding variogram gamma
      Theta <- build_theta(R, clusters)
      Gamma <- Theta2Gamma(Theta)

      list(
        R = R,
        clusters = clusters,
        Theta = Theta,
        Gamma = Gamma,
        chi = 1,
        n = 1e3,
        d = 7
      )
    }
  ),

  # Replication of the unbalanced simulation
  tar_rep(
    unbal_sim_clustering_replicate,
    rmpareto(n = unbal_sim_param_cluster$n,
             model = "HR",
             par = unbal_sim_param_cluster$Gamma),
    batches = 1,
    reps = 200,
    iteration = "vector"               # the output of the replicates is a list
  ),

  # Data formatting to separate the replications
  tar_target(
    unbal_sim_rep_data,
    {
      rep <- unbal_sim_clustering_replicate

      n <- unbal_sim_param_cluster$n

      row.names(rep) <- NULL

      rep |>
        tibble() |>
        mutate(sim = (seq_len(nrow(rep)) + (n - 1)) %/% n)
    }
  ),

  # Computation of optimization results for each replication
  tar_target(
    unbal_sim_rep_pen_results,
    {
      m <- unbal_sim_rep_data |> summarise(max = max(sim))
      lambda <- seq(0, 4e-2, 5e-4)
      future_map(1:m$max, function(i) {
        data <- unbal_sim_rep_data |>
          filter(sim == i) |>
          select(-sim) |>
          as.matrix()

        future_map(lambda, \(.) {
          best_clusters(data,
                        unbal_sim_param_cluster$chi,
                        l_grid = .)
        })
      })
    }
  ),

  tar_target(
    unbal_sim_rep_nopen_results,
    {
      m <- unbal_sim_rep_data |> summarise(max = max(sim))
      future_map(1:m$max, function(i) {
        data <- unbal_sim_rep_data |>
          filter(sim == i) |>
          select(-sim) |>
          as.matrix()

        best_clusters(data, chi = unbal_sim_param_cluster$chi,
                      l_grid = 0)

      })
    }
  ),

  # Formatting and plot the results
  tar_target(
    unbal_sim_rep_results,
    {
      lambda <- seq(0, 4e-2, 5e-4)
      all_info(cluster.init = unbal_sim_param_cluster$clusters,
               list_res = unbal_sim_rep_pen_results, lambda = lambda)
    }
  ),

  tar_target(
    unbal_sim_hierarchy,
    average_hierarchy(unbal_sim_rep_pen_results)
  ),

  ##======================== With three balanced groups ========================
  tar_target(
    gr3_bal_sim_param_cluster,
    {
      # Construction of clusters and R matrix
      R <- matrix(c(1, -3, 0,
                    -3, 2, -2,
                    0, -2, 1), nc = 3)
      clusters <- list(1:5, 6:10, 11:15)

      # Construction of induced theta and corresponding variogram gamma
      Theta <- build_theta(R, clusters)
      Gamma <- Theta2Gamma(Theta)

      list(
        R = R,
        clusters = clusters,
        Theta = Theta,
        Gamma = Gamma,
        chi = 1,
        n = 1e3,
        d = 15
      )
    }
  ),

  # Replication of the 3 groups balanced simulation
  tar_rep(
    gr3_bal_sim_clustering_replicate,
    rmpareto(n = gr3_bal_sim_param_cluster$n,
             model = "HR",
             par = gr3_bal_sim_param_cluster$Gamma),
    batches = 1,
    reps = 200,
    iteration = "vector"               # the output of the replicates is a list
  ),

  # Data formatting to separate the replications
  tar_target(
    gr3_bal_sim_rep_data,
    {
      rep <- gr3_bal_sim_clustering_replicate

      n <- gr3_bal_sim_param_cluster$n

      row.names(rep) <- NULL

      rep |>
        tibble() |>
        mutate(sim = (seq_len(nrow(rep)) + (n - 1)) %/% n)
    }
  ),

  # Computation of optimization results for each replication
  tar_target(
    gr3_bal_sim_rep_pen_results,
    {
      m <- gr3_bal_sim_rep_data |> summarise(max = max(sim))
      lambda <- seq(0, 1, 1e-1)
      future_map(1:m$max, function(i) {
        data <- gr3_bal_sim_rep_data |>
          filter(sim == i) |>
          select(-sim) |>
          as.matrix()

        future_map(lambda, \(.) {
          best_clusters(data,
                        gr3_bal_sim_param_cluster$chi,
                        l_grid = .)
        })
      })
    }
  ),

  tar_target(
    gr3_bal_sim_rep_nopen_results,
    {
      m <- gr3_bal_sim_rep_data |> summarise(max = max(sim))
      future_map(1:m$max, function(i) {
        data <- gr3_bal_sim_rep_data |>
          filter(sim == i) |>
          select(-sim) |>
          as.matrix()

        best_clusters(data, chi = gr3_bal_sim_param_cluster$chi,
                      l_grid = 0)

      })
    }
  ),

  # Formatting and plot the results
  tar_target(
    gr3_bal_sim_rep_results,
    {
      lambda <- seq(0, 1, 1e-1)
      all_info(cluster.init = gr3_bal_sim_param_cluster$clusters,
               list_res = gr3_bal_sim_rep_pen_results, lambda = lambda)
    }
  ),

  tar_target(
    gr3_bal_hierarchy,
    average_hierarchy(gr3_bal_sim_rep_pen_results)
  ),

  ##=================== With three balanced groups with noise ==================
  tar_rep(
    gr3_noise_sim_clustering_replicate,
    {
      # Error matrix creation
      sim <- matrix(rnorm(n = gr3_bal_sim_param_cluster$d**2, 0.1, 1e-2), nc = gr3_bal_sim_param_cluster$d)
      E <- matrix(0, gr3_bal_sim_param_cluster$d, gr3_bal_sim_param_cluster$d)

      # Keep only the upper trinagular matrix produce symetric matrix
      E[upper.tri(E, diag = TRUE)] <- sim[upper.tri(sim, diag = TRUE)]

      # No error in the diagonal to preserve gamma_ii = 0
      diag(E) <- 0

      # Perturbating matrix
      G_noise <- abs(gr3_bal_sim_param_cluster$Gamma + E + t(E))

      # Simulation with the new matrix
      rmpareto(n = gr3_bal_sim_param_cluster$n,
               model = "HR",
               par = G_noise)
    },
    batches = 1,
    reps = 200,
    iteration = "vector"               # the output of the replicates is a list
  ),

  # Data formatting to separate the replications
  tar_target(
    gr3_noise_sim_rep_data,
    {
      rep <- gr3_noise_sim_clustering_replicate

      n <- gr3_bal_sim_param_cluster$n

      row.names(rep) <- NULL

      rep |>
        tibble() |>
        mutate(sim = (seq_len(nrow(rep)) + (n - 1)) %/% n)
    }
  ),

  # Computation of optimization results for each replication
  tar_target(
    gr3_noise_sim_rep_pen_results,
    {
      m <- gr3_noise_sim_rep_data |> summarise(max = max(sim))
      lambda <- seq(0, 2, 25e-2)
      future_map(1:m$max, function(i) {
        data <- gr3_noise_sim_rep_data |>
          filter(sim == i) |>
          select(-sim) |>
          as.matrix()

        future_map(lambda, \(.) {
          best_clusters(data,
                        gr3_bal_sim_param_cluster$chi,
                        l_grid = .)
        })
      })
    }
  ),

  tar_target(
    gr3_noise_sim_rep_nopen_results,
    {
      m <- gr3_noise_sim_rep_data |> summarise(max = max(sim))
      future_map(1:m$max, function(i) {
        data <- gr3_noise_sim_rep_data |>
          filter(sim == i) |>
          select(-sim) |>
          as.matrix()

        best_clusters(data, chi = gr3_bal_sim_param_cluster$chi,
                      l_grid = 0)

      })
    }
  ),

  # Formatting and plot the results
  tar_target(
    gr3_noise_sim_rep_results,
    {
      lambda <- seq(0, 2, 25e-2)
      all_info(cluster.init = gr3_bal_sim_param_cluster$clusters,
               list_res = gr3_noise_sim_rep_pen_results, lambda = lambda)
    }
  ),

  tar_target(
    gr3_noise_sim_numeric_results,
    {
      # Some numerical results in different situations
      ## When 3 clusters is reached
      nb_3clusters <- gr3_noise_sim_rep_results |>
        filter(nb_cluster == 3) |>
        select(simulation) |>
        unique() |>
        count() |>
        as.numeric()

      average_ARI <- gr3_noise_sim_rep_results |>
        filter(nb_cluster == 3) |>
        group_by(simulation) |>
        summarise(minARI = min(ARI)) |>
        summarise(m = mean(minARI)) |>
        as.numeric()

      ## When 3 clusters is not reached before the lambda max
      over_3clusters <- gr3_noise_sim_rep_results |>
        filter(l == 2, nb_cluster > 3) |>
        count() |>
        as.numeric()

      over3_average_ARI <- gr3_noise_sim_rep_results |>
        filter(l == 2, nb_cluster > 3) |>
        summarise(m = mean(ARI)) |>
        as.numeric()

      # When we finally end with 2 clusters or less
      under_3clusters <- gr3_noise_sim_rep_results |>
        filter(l == 2, nb_cluster < 3) |>
        count() |>
        as.numeric()

      under3_average_ARI <- gr3_noise_sim_rep_results |>
        filter(l == 2, nb_cluster < 3) |>
        summarise(m = mean(ARI)) |>
        as.numeric()

      list(
        with_3 = c(nb_3clusters, average_ARI),
        other = data.frame(under3 = c(under_3clusters, under3_average_ARI),
                           over3 = c(over_3clusters, over3_average_ARI))
      )
    }
  ),

  tar_target(
    gr3_noise_hierarchy,
    average_hierarchy(gr3_noise_sim_rep_pen_results)
  ),

  tar_target(
    average_weight_graph_noise_sim,
    {
      data <- future_map(1:200, function(i) {
        gr3_noise_sim_rep_data |>
          filter(sim == i) |>
          select(-sim) |>
          as.matrix()
      }
      )

      W <- (Reduce("+", lapply(data, compute_W))) / 200

      # Graph creation
      g <- graph_from_adjacency_matrix(W, mode = "undirected",
                                       weighted = TRUE, diag = FALSE)
      ggraph(g, layout = "circle") +  # 'dh' = Davidson-Harel layout
        geom_edge_link(aes(edge_width = weight,
                           edge_alpha = weight, color = weight),
                       show.legend = c(edge_width = FALSE,
                                       edge_alpha = FALSE, color = TRUE)) +
        geom_node_point(size = 10, color = "grey70") +
        geom_node_text(aes(label = 1:7)) +
        scale_edge_color_gradient(low = "white", high = "darkorange2",
                                  name = "Weight values") +
        scale_edge_width(range = 0.9) +
        theme_void()
    }
  ),


  ##==================== With chi matrix innitialization =======================
  tar_target(
    chi_init_sim_param_cluster,
    {
      # Construction of clusters and R matrix
      R_chi <- matrix(c(0.2, 0.25, 0.32,
                        0.25, 0.6, 0.45,
                        0.32, 0.45, 0.8), nc = 3)
      clusters <- list(1:5, 6:10, 11:15)

      # Construction of induced theta and corresponding variogram gamma
      Chi <- build_theta(R_chi, clusters)
      diag(Chi) <- 1
      Gamma <- ChiToGamma(Chi)

      list(
        R_chi = R_chi,
        clusters = clusters,
        Chi = Chi,
        Gamma = Gamma,
        zeta = 1,
        n = 1e3,
        d = 15
      )
    }
  ),

  # Replication of the 3 groups balanced simulation
  tar_rep(
    chi_init_sim_clustering_replicate,
    rmpareto(n = chi_init_sim_param_cluster$n,
             model = "HR",
             par = chi_init_sim_param_cluster$Gamma),
    batches = 1,
    reps = 200,
    iteration = "vector"               # the output of the replicates is a list
  ),

  # Data formatting to separate the replications
  tar_target(
    chi_init_sim_rep_data,
    {
      rep <- chi_init_sim_clustering_replicate

      n <- chi_init_sim_param_cluster$n

      row.names(rep) <- NULL

      rep |>
        tibble() |>
        mutate(sim = (seq_len(nrow(rep)) + (n - 1)) %/% n)
    }
  ),

  # Computation of optimization results for each replication
  tar_target(
    chi_init_sim_rep_pen_results,
    {
      m <- chi_init_sim_rep_data |> summarise(max = max(sim))
      lambda <- seq(0, 1, 1e-1)
      n_sim <- 1:m$max
      future_map(n_sim[-c(50, 86, 137)], function(i) {
        data <- chi_init_sim_rep_data |>
          filter(sim == i) |>
          select(-sim) |>
          as.matrix()

        future_map(lambda, \(.) {
          best_clusters(data,
                        chi_init_sim_param_cluster$zeta,
                        l_grid = .)
        })
      })
    }
  ),

  tar_target(
    chi_init_sim_rep_nopen_results,
    {
      m <- chi_init_sim_rep_data |> summarise(max = max(sim))
      future_map(1:m$max, function(i) {
        data <- chi_init_sim_rep_data |>
          filter(sim == i) |>
          select(-sim) |>
          as.matrix()

        best_clusters(data, chi = chi_init_sim_param_cluster$zeta,
                      l_grid = 0)

      })
    }
  ),

  # Formatting and plot the results
  tar_target(
    chi_init_sim_rep_results,
    {
      lambda <- seq(0, 1, 1e-1)
      all_info(cluster.init = chi_init_sim_param_cluster$clusters,
               list_res = chi_init_sim_rep_pen_results, lambda = lambda)
    }
  ),

  tar_target(
    chi_init_hierarchy,
    average_hierarchy(chi_init_sim_rep_pen_results)
  ),

  #----------------------------- Export document -------------------------------
  # Trivariate coefficient documents
  ## HTML
  tar_render(
    trivariate_html,
    path = "./src/trivariate.Rmd",
    output_format = "readthedown",
    output_file = here("public", "trivariate.html")
  ),

  ## PDF
  tar_render(
    trivariate_pdf, error = "continue",
    path = "./src/trivariate.Rmd",
    output_format = "pdf_document",
    output_file = here("report", "trivariate.pdf")
  ),

  # Spectral representation for graphical models documents
  ## HTML
  tar_render(
    spectral_html,
    path = "./src/spectral.Rmd",
    output_format = "readthedown",
    output_file = here("public", "spectral.html")
  ),

  ## PDF
  tar_render(
    spectral_pdf, error = "continue",
    path = "./src/spectral.Rmd",
    output_format = "pdf_document",
    output_file = here("report", "spectral.pdf")
  ),

  # Conditional independance on extremal linear latent model documents
  ## HTML
  tar_render(
    latent_html,
    path = "./src/latent.Rmd",
    output_format = "readthedown",
    output_file = here("public", "latent.html")
  ),

  ## PDF
  tar_render(
    latent_pdf, error = "continue",
    path = "./src/latent.Rmd",
    output_format = "pdf_document",
    output_file = here("report", "latent.pdf")
  ),

  # Variable clustering for Husler-Reiss graphical models
  ## HTML
  tar_render(
    cluster_html,
    path = "./src/cluster.Rmd",
    output_format = "readthedown",
    output_file = here("public", "cluster.html")
  ),

  ## PDF
  tar_render(
    cluster_pdf, error = "continue",
    path = "./src/cluster.Rmd",
    output_format = "pdf_document",
    output_file = here("report", "cluster.pdf")
  )
)

#-------------------------------------------------------------------------------
#                           CORRECTION AND VISUALISATION
#-------------------------------------------------------------------------------
#tar_prune()           clean the obsolete files
#tar_visnetwork()      Viex of the dependancies between datas and functions
