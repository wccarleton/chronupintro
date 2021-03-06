# Chronological Uncertainty Propagation in REC Models
## Overview
This repo contains the data and code used for the study presented in the following pre-print:

[*Improved parameter estimation and uncertainty propagation in Bayesian Radiocarbon-dated Event Count [REC] models*](https://doi.org/10.31219/osf.io/7p9vx)

## Abstract
Data about the past contain chronological uncertainty that needs to be accounted for in statistical models. Recently a method called Radiocarbon-dated Event Count (REC) modelling has been explored as a way to improve the handling of chronological uncertainty in the context of statistical regression. REC modelling has so far employed a Bayesian hierarchical framework for parameter estimation to account for chronological uncertainty in count series of radiocarbon-dates. This approach, however, suffers from a couple of limitations. It is computationally inefficient, which limits the amount of chronological uncertainty that can be accounted for, and the hierarchical framework can produce biased, but highly precise parameter estimates. Here we report the results of an investigation in which we compared hierarchical REC models to an alternative with simulated data and a new R package called "chronup". Our results indicate that the hierarchical framework can produce correct high-precision estimates given enough data, but it is susceptible to sampling bias and has an inflated Type I error rate. In contrast, the alternative better handles small samples and fully propagates uncertainty into parameter estimates. In light of these results, we think the alternative method is more generally suitable for Palaeo Science applications.

## Software
The R scripts contained in this repository are intended for replication efforts and to improve the transparency of research. They are, of course, provided without warranty or technical support. That said, questions about the code can be directed to me, Chris Carleton, at ccarleton@protonmail.com.

### R
This analysis described in the associated manuscript was performed in R. Thus, you may need to download the latest version of [R](https://www.r-project.org/) in order to make use of the scripts described below.

### Nimble
This project made use of a Bayesian Analysis package called [Nimble](https://r-nimble.org/). See the Nimble website for documentation and a tutorial. Then, refer to the R scripts in this repo.

## Contact

[ORCID](https://orcid.org/0000-0001-7463-8638) |
[Google Scholar](https://scholar.google.com/citations?hl=en&user=0ZG-6CsAAAAJ) |
[Website](https://wccarleton.me)

## License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
