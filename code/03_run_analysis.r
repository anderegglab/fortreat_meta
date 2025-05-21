##---------------------------------------------------------------
## Run complete analysis
##
## author: Jacob Levine; contact: jacob.levine@utah.edu
##---------------------------------------------------------------

## I know this file seems weird, but its just here so that we can quickly regenerate tables and figures
## after making a change to the underlying data.
source("code/01_data_cleaning.r")

## maps
source("code/figures/figure1.r")

## overall treatment effects
source("code/figures/figure2.r") ## survivorship
source("code/figures/figure3.r") ## carbon

## treatment severity
source("code/figures/figure4.r")

## biome
source("code/figures/figure5.r") ## survivorship
source("code/figures/figure6.r") ## carbon

## climate
source("code/figures/figure7.r")

## publication bias
source("code/03_publication_bias.r")
