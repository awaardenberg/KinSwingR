#' KinSwingR: A package for predicting kinase activity
#'
#' This package provides functionality for kinase-subtrate prediction, and 
#' integration with phosphopeptide fold change and signficance to assess the 
#' local connectivity (swing) of kinase-substrate networks. The final output of 
#' KinSwingR is a score that is normalised and weighted for prediction of kinase
#'  activity.
#'
#' Contact a.waardenberg@gmail.com for questions relating to functionality.
#'
#' @section buildPWM function:
#' Builds PWMs for kinases from a table of kinases and known substrate 
#' sequences.
#'
#' @section scoreSequences function:
#' Score kinase PWMs matches against a set of peptide seqeuences.
#'
#' @section swing function:
#' Integrates kinase PWMs matches against peptide seqeuences and directionality 
#' as well as significance of peptides for prediction of kinase activity.
#'
#' @section cleanAnnotation function:
#' Function for extracting peptides from multimapped data
#'
#' @docType package
#' @name KinSwingR
#'
NULL
