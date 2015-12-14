# R4Med

A collection of R scripts for analysis of fisheries data in the Mediterranean.
It has been put together to assist the analyses performed during the STECF SGMED meetings.

The repository was started in Brussels, December 2013

## Before each meeting 

To make your life easier, before each meeting you should:

1. Install the latest version of R (32-bit)
2. Install the important dependency packages, using install.packages():

    * triangle
    * copula
    * ggplot2
    * gridExtra 

3. Install the latest FLR packages, using install.packages(repos="http://flr-project.org/R", pkg=PKGNAME):

    * FLCore
    * FLAssess
    * FLash
    * FLXSA
    * FLa4a
    * ggplotFL
    * FLBRP
    * FLEDA

## MEDITS indices preparation

Indices can be prepared using the R/medits_index_preparation/LFD.R script.
This interrogates the TA, TB and TC files (as CSV) and generates CSV files of standardised length-based data.
The length-based data can then be sliced using the slicing_medits.R script.

## Landings data preparation

See the slicing/slicing_dcf.R script.

## Forecasts

The forecasts folder includes scripts for short term forecast, medium term forecast, a simple multifleet forecast and a simple MSE.
