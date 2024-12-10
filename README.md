# Extremal Graphical Models

The main purpose of this git repository is to share my results and work on graphical models applied to Extreme Value Theory.

## Organization

### `_targets.R` file
The code utilizes the `targets` package from `R`, which is useful for building pipelines and ensuring reproducibility of my results.

To compile the `_targets.R` file, you first need to install the `targets` and `tarchetypes` packages in your `R` environment. You can do this by running the following commands in the `R` terminal:

```R
install.packages("targets")
install.packages("tarchetypes")
```
(You may need additional packages to run the pipeline. All dependencies are listed in the requirement.txt file.)

To execute the pipeline, simply run the following command:

```R
targets::tar_make()
```

If you wish to inspect data, results, plots, or other outputs, you can load specific objects using the command:

```R
targets::tar_load(object_name) # Replace 'object_name' with the first argument of each tar_target() in _targets.R
```

The object will then be loaded into your environment as object_name.

### `src` directory
This folder contains all the raw files related to my research, including `.Rmd` files used for reports and `.R` files that store functions utilized in the pipeline.

### `report` directory
In this directory, you will find `.pdf` documents such as reports, summaries, and other materials presenting the results.

### `public` directory
This folder contains `.html` files that correspond to the `.pdf` documents in the report directory. These files can be opened directly in a web browser.

## Works

My initial goal was to identify coefficients that summarize all conditional independence relationships between the components of the MGPD (Multivariate Generalized Pareto Distribution). The results of this work can be found in the `trivariate` output located in the `src` or `public` directories.

## Materials and references

- `targets` [documentation](https://books.ropensci.org/targets/).

- `graphicalExtreme` [git](https://github.com/sebastian-engelke/graphicalExtremes).

- Engelke, Sebastian, and Adrien S. Hitz. 2020. “Graphical Models for Extremes.” Journal of the Royal Statistical Society: Series B (Statistical Methodology) 82 (4): 871–932. https://doi.org/10.1111/rssb.12355.
  
- Hentschel, M., Engelke, S., and Segers, J. (2022). Statistical inference for Hüsler-Reiss graphical models through matrix completions. arXiv. https://doi.org/10.48550/ARXIV.2210.14292

- Ferreira, H. 2011. “Dependence Between Two Multivariate Extremes.” Statistics & Probability Letters 81 (5): 586–91. https://doi.org/10.1016/j.spl.2011.01.014.

- Strokorb, Kirstin. 2020. “Extremal Independence Old and New.” arXiv. https://doi.org/10.48550/arXiv.2002.07808.


