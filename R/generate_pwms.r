#' Generate Position Weight Matrices (PWMs)
#'
#' @description Generate Position Weight Matrices (PWMs) for a table containing 
#' centered substrate peptide sequences for a list of kinases. The output of 
#' this function is to be used for scoring PWM matches to peptides via 
#' scoreSequences()
#' @param kinase_table A data.frame of substrate sequences and kinase names. 
#' Format of data must be as follows: column 1 - kinase/kinase family 
#' name/GeneID, column 2 - centered peptide seqeuence.
#' @param wild_card Letter to describe sequences that are outside of the protein
#'  after centering on the phosphosite (e.g ___MERSTRELCLNF). Default: "_".
#' @param substrate_length Full length of substrate sequence (default is 15). 
#' Will be trimmed automatically or report error if sequences in kinase_table 
#' are not long enough.
#' @param substrates_n Number of sequences used to build a PWM model. Low 
#' sequence counts will produce poor representative PWM models. Default: "10"
#' @param pseudo Small number to add to values for PWM log transformation to 
#' prevent log transformation of zero. Default = 0.01
#' @param remove_center Remove all peptide seqeuences with the central amino 
#' acid matching a character (e.g. "y"). Default = FALSE
#' @param verbose Print progress to screen. Default=FALSE
#'
#' @examples
#' ## Build PWM models from phosphositeplus data with default of minimum
#' ## of 10 substrate sequences for building a PWM model.
#'
#' data(phosphositeplus_human)
#'
#' ##randomly sample 1000 substrates for demonstration.
#' set.seed(1)
#' sample_pwm <- phosphositeplus_human[sample(nrow(phosphositeplus_human), 
#' 1000),]
#' pwms <- buildPWM(sample_pwm)
#'
#' ## Data frame of models built and number of sequences used to build each
#' ## PWM model:
#' head(pwms$kinase)
#'
#' @return Output is a list containing two tables, "pwm" and "kinase". To access
#'  PWMs: pwms$pwm and Table of Kinase and sequence counts: 
#'  pwms$kinase
#'
#' @export buildPWM

buildPWM <- function(kinase_table = NULL,
                     wild_card = "_",
                     substrate_length = 15,
                     substrates_n = 10,
                     pseudo = 0.01,
                     remove_center = FALSE,
                     verbose = FALSE) {
  #----------------------------------------------
  #format checks:
  if (is.null(kinase_table))
    stop("kinase_table not provided; you must provide an input table")
  if (!is.matrix(kinase_table))
    stop("kinase_table is not a table; you must provide an
         input table")
  #remove NA's
  kinase_table <- kinase_table[!is.na(kinase_table[, 2]),]
  #remove problematic characters (creates row/column name issues)
  kinase_table[,1] <- gsub("/", "_", kinase_table[,1])
  
  if (substrate_length < 3)
    stop(
      "substrate_length needs to be greater than 2; increase the
      length of substrate window size"
    )
  if ((lapply(substrate_length, "%%", 2) == 0) == TRUE)
    stop("substrate_length must be an odd number! I.e. centered sequence
         of ++X++")
  #----------------------------------------------
  #call trim_seq function and trim sequences:
  kinase_table[, 2] <-
    trimSeqs(kinase_table[, 2], 
             seq_length = substrate_length, 
             verbose = verbose)
  kinase_table[, 2] <- toupper(kinase_table[, 2])
  
  #remove peptides with a certain center character:
  if (remove_center != FALSE) {
    half_window <- (substrate_length + 1) / 2
    center_aa <-
      substr(kinase_table[, 2], half_window, half_window)
    if (verbose) {
      message("You have selected to remove",
              length(which(center_aa == toupper(remove_center))),
              "peptide sequences for building PWMs that contain a centered 
              letter of:",
              toupper(remove_center)
      )
    }
    kinase_table <-
      kinase_table[center_aa != toupper(remove_center),]
  }
  
  kinase_table <- unique(kinase_table)
  #----------------------------------------------
  
  if (verbose) {
    message(nrow(kinase_table), "unique kinase:substrate sequences in 
            table provided."
    )
  }
  
  #initialise table for kinase count data
  kinase_summary <-
    data.frame("kinase" = unique(as.character(kinase_table[, 1])), "n" = NA)
  
  kinase_summary$n <-
    c(sapply(seq_len(nrow(kinase_summary)), function(i)
      length(kinase_table[, 1][kinase_table[, 1] == kinase_summary[i, 1]])))
  
  #filter the list here which is used to build pwms:
  kinase_summary <-
    kinase_summary[kinase_summary$n >= substrates_n, ]
  
  #generate PWM list:
  pwm_list <-
    c(lapply(seq_len(nrow(kinase_summary)), function(i)
      scorePWM(
        as.character(kinase_table[, 2][kinase_table[, 1] == 
                                         kinase_summary[i, 1]]),
        substrate_length = substrate_length,
        wild_card = wild_card,
        pseudo = pseudo
      )))
  
  return(list("pwm" = pwm_list, "kinase" = kinase_summary))
  }

# this helper function performs the calculations for building the PWMs

scorePWM <- function(input_data,
                     substrate_length,
                     wild_card = "_",
                     pseudo = 0.01) {
  #reformat input_data; build a matrix of the sequences (split into individual
  # AA's) = substrate_length
  input_data <-
    data.frame(matrix(unlist(strsplit(input_data , "")) , 
                      ncol = substrate_length , byrow = TRUE))
  uniq_AA <- c(wild_card, "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
               "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  # 1. Build PFM
  pwm <-
    t(cbind(sapply(seq_along(uniq_AA), function(i)
      (
        sapply(seq_len(ncol(input_data)), function(j)
          sum(
            ifelse(as.character(input_data[, j]) == 
                     as.character(uniq_AA[i]), 1, 0)
          ))
      ))))
  rownames(pwm) <- uniq_AA
  colnames(pwm) <- paste("p", seq_len(ncol(pwm)), sep = "")
  # remove wildcard from motif scoring
  pwm <- pwm[!rownames(pwm) %in% c(wild_card),]
  # 2. PPM calculation
  pwm_rownames <- rownames(pwm)
  pwm_colnames <- colnames(pwm)
  # addition of pseduo count to avoid log zero    
  pwm <- pwm + pseudo
  col_counts <- apply(pwm, 2, sum, na.rm = TRUE)
  pwm <- sapply(seq_len(ncol(pwm)), function(i) 
    pwm[,i] / col_counts[i])
  # carry labels forward
  rownames(pwm) <- pwm_rownames
  colnames(pwm) <- pwm_colnames
  # 3. Calculate Position Weight Matrix
  pwm <- log(pwm / (1 / nrow(pwm)))
  return(pwm)
}
