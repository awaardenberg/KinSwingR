#' KinSwingR: A package for predicting kinase activity
#'
#' This package provides functionality for kinase-subtrate prediction, and integration with phosphopeptide
#' fold change and signficance to assess the local connectivity (swing) of kinase-substrate networks. The
#' final output of KinSwingR is a score that is normalised and weighted for prediction of kinase activity.
#'
#' This is a development version of KinSwingR Contact a.waardenberg@gmail.com for questions relating to
#' functionality.
#'
#' @section build.pwm function:
#' Builds PWMs for kinases from a table of kinases and known substrate sequences.
#'
#' @section scores.sequences function:
#' Score kinase PWMs matches against a set of peptide seqeuences.
#'
#' @section swing function:
#' Integrates kinase PWMs matches against peptide seqeuences and directionality as well as significance of
#' peptides for prediction of kinase activity.
#'
#' @section swing.master function:
#' One master function for computing the swing statistic in one command. Run clean.annotation() first, if
#' peptides need to be extracted from the annotations.
#'
#' @section clean.annotation function:
#' Function for extracting peptides from multimapped data
#'
#' @docType package
#' @name KinSwingR
#'
### functions used:
#' @importFrom data.table data.table
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam
#' @importFrom sqldf sqldf
#' @importFrom stats setNames
#' @importFrom stats sd
#' 
NULL
