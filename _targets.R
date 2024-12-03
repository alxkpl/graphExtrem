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
                            "tidyr",
                            "dplyr",
                            "ggplot2",
                            "here",
                            "rmdformats"))

#-------------------------------------------------------------------------------
#                                   PIPELINE
#-------------------------------------------------------------------------------
list(
  #------------------------------- Calculation ---------------------------------
  # Computation of a the trivariate coefficient in the Husler-Reiss case
  tar_target(
    Gamma_cond_indep_12,
    random_Gamma12()          # Random valid HR variogram matrix 3x3 
    ),
  
  tar_target(
    chi_trivariate_example,
    chi_trivariate_HR(Gamma_cond_indep_12)    # Computation of the coefficient
    ),
  # not zero in this case but it was obvious unfortunately...
  
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
    trivariate_pdf,
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
