
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mteucb：Uniform Confidence Band for the Marginal Treatment Function.

<!-- badges: start -->
<!-- badges: end -->

The **mteucb** package cmputes uniform confidence bands for the MTE
function proposed in Tsuda, Jin and Okui (2024+). A uniform confidence
band covers the true function with a prespecified probability. These
bands are used to make statistical inferences on the shape of the MTE
function. Uniform confidence bands are particularly useful because they
can visualize the statistical uncertainty behind the estimated MTE
function.

## Installation

You can install the development version of mteucb from GitHub:

``` r
# install.packages("devtools") if necessary
devtools::install_github("TXT872/mteucb",build_vignettes = TRUE)
```

## Vignette

For more detalis, see the package vignette with:

``` r
# Getting Started with the mteucb Package Vignette

vignette("p_hut_gen") # Generates Estimated Propensity Scores

vignette("beta_gen") # Estimates a Parametric Part of a MTE Function 

vignette("unif_gen") # Provides the Uniform Confidence Band for a MTE Function 
#Given a Parametric Part of a MTE Function is Estimated

vignette("uniform_con_gen") # Provides the Uniform Confidence Band for a MTE Function from Data.
```

## References

- Tsuda, T, Jin, Y., & Okui, R. (2024+). Uniform Confidence Band for
  Marginal Treatment Effect Function. will be available on arXiv.
