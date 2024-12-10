library(targets)
library(tarchetypes)
library(here)

#-------------------------------------------------------------------------------
#                                     SETUP 
#-------------------------------------------------------------------------------
tar_source(files = "./src") # compile all the script in the src directory

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
                            "gridExtra"))

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
  
  # Computation of the trivariate coefficient for all the simulation
  tar_target(
    chi_trivariate_simulation_results,
    data.frame(num_sim = 1:nrow(Gamma_cond_indep_12_replicate),
                 Gamma_cond_indep_12_replicate, 
               row.names = NULL) |> 
      rowwise() |> 
      mutate(tri_chi = chi_trivariate_HR(c(X1, X2, X3)),
             prod_chi = (2 - theta(X2)) * (2 - theta(X3)))
  ),

  # Computation of the conditional coefficient for all the simulation
  tar_target(
    cond_trivariate_simulation_results,
    data.frame(num_sim = 1:nrow(Gamma_cond_indep_12_replicate),
               Gamma_cond_indep_12_replicate, 
               row.names = NULL) |> 
      rowwise() |> 
      mutate(tri_chi = cond_trivariate_HR(c(X1, X2, X3), 3))
    
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
             y = expression(paste(chi[12]^3 , " value"))) +
        theme_ipsum(base_family = "serif")
      
      p2 <- cond_trivariate_simulation_results |> 
        ggplot() + aes(x = tri_chi, y = 0) + 
        geom_violin(col = "#474F58", fill = "darkorange") + 
        coord_cartesian(xlim = c(0, 1)) +
        labs(title = "Violin plot",  y = " ", 
             x = expression(paste(chi[12]^3 , " value"))) +
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
      geom_line(aes(y = abs(tri_chi - prod_chi) / tri_chi), col = "darkorange") +
      coord_cartesian(ylim = c(0,1)) +
      labs(x = "Simulation number",
           y = "Error",
           title = "Relative error between the coefficients") +
      theme_ipsum(base_family = "serif")
  ),
  
  
  tar_target(
    all_matrix_coeff,
    expand.grid(b = seq(0.1,10, 0.05), c = seq(0.1,10, 0.05)) |>
      as_tibble() |>
      mutate(a = b + c) |> 
      rowwise() |> 
      mutate(tri_chi = chi_trivariate_HR(c(a, b, c)),
             prod_chi = (2 - theta(b)) * (2 - theta(c))) 
  ),
  
  tar_target(
    plot_all_matrix_coef,
    all_matrix_coeff |> 
      ggplot() + aes(x = b, y = c, fill = tri_chi / prod_chi) +
      geom_tile() + 
      scale_fill_continuous(name = expression(paste(chi[123], "/",chi[13],chi[23])),
                            type = "viridis")+
      labs(x = expression(Gamma[13]), y = expression(Gamma[23]),
           caption = expression(paste("The conditional independance is knowing ",
                                      Y[3], ".")))+
      theme_ipsum(base_family = "serif", grid = FALSE, axis_title_size = 12)
  ),
  
  
  tar_target(
    plot_all_matrix_error,
    {
      p1 <- all_matrix_coeff |> 
              ggplot() + aes(x = b, y = c, fill = abs(tri_chi - prod_chi) / tri_chi) +
              geom_tile() + 
              scale_fill_continuous(name = expression(paste("|", chi[123], "-",
                                                            chi[13],chi[23], "|", "/",
                                                            chi[13],chi[23])),
                                  type = "viridis")+
              labs(x = expression(Gamma[13]), y = expression(Gamma[23]),
                   title = "Relative error")+
              theme_ipsum(base_family = "serif", grid = FALSE, axis_title_size = 12)
      
      p2 <- all_matrix_coeff |> 
        ggplot() + aes(x = b, y = c, fill = abs(tri_chi - prod_chi)) +
        geom_tile() + 
        scale_fill_viridis_c(name = expression(paste("|", chi[123], "-",chi[13],chi[23], "|")),
                             option = "inferno") +
        labs(x = expression(Gamma[13]), y = expression(Gamma[23]),
             title = "Absolute error",
             caption = expression(paste("The conditional independance is knowing ",
                                        Y[3], ".")))+
        theme_ipsum(base_family = "serif", grid = FALSE, axis_title_size = 12)
      
      grid.arrange(grobs = list(p1, p2), nrow = 2)
      }
  ),
  
  #----------------------------- Export document -------------------------------
  # Trivariate coefficient documents
  ## HTML
  tar_render(
    trivariate_html, 
    path = "./src/trivariate.Rmd",
    output_format = "readthedown",
    output_file = here("public","trivariate.html")
    ),
  
  ## PDF
  tar_render(
    trivariate_pdf, error = "continue",  
    path = "./src/trivariate.Rmd",
    output_format = "pdf_document",
    output_file = here("report","trivariate.pdf")
    )
)

#-------------------------------------------------------------------------------
#                           CORRECTION AND VISUALISATION
#-------------------------------------------------------------------------------
#tar_prune()          # clean the obsolete files
#tar_visnetwork()     # Viex of the dependancies between datas and functions
