#' Generate Position Weight Matrices (PWMs)
#'
#' @description Generate Position Weight Matrices (PWMs) for a table containing 
#' centered substrate peptide sequences for a list of kinases. The output of 
#' this function is to be used for scoring PWM matches to peptides via 
#' scoreSequences()
#' @param kinase.table A data.frame of substrate sequences and kinase names. 
#' Format of data must be as follows: column 1 - kinase/kinase family 
#' name/GeneID, column 2 - centered peptide seqeuence.
#' @param wild.card Letter to describe sequences that are outside of the protein
#'  after centering on the phosphosite (e.g ___MERSTRELCLNF). Default: "_".
#' @param substrate.length Full length of substrate sequence (default is 15). 
#' Will be trimmed automatically or report error if sequences in kinase.table 
#' are not long enough.
#' @param substrates.n Number of sequences used to build a PWM model. Low 
#' sequence counts will produce poor representative PWM models. Default: "10"
#' @param pseudo Small number to add to values for PWM log transformation to 
#' prevent log transformation of zero. Default = 0.01
#' @param remove.center Remove all peptide seqeuences with the central amino 
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

buildPWM <- function(kinase.table = NULL,
                      wild.card = "_",
                      substrate.length = 15,
                      substrates.n = 10,
                      pseudo = 0.01,
                      remove.center = FALSE,
                      verbose = FALSE) {
  #----------------------------------------------
  #format checks:
  if (is.null(kinase.table))
    stop("kinase.table not provided; you must provide an input table")
  if (!is.matrix(kinase.table))
    stop("kinase.table is not a table; you must provide an
         input table")
  #remove NA's
  kinase.table <- kinase.table[!is.na(kinase.table[, 2]),]
  
  if (substrate.length < 3)
    stop(
      "substrate.length needs to be greater than 2; increase the
      length of substrate window size"
    )
  if ((lapply(substrate.length, "%%", 2) == 0) == TRUE)
    stop("substrate.length must be an odd number! I.e. centered sequence
         of ++X++")
  #----------------------------------------------
  #call trim.seq function and trim sequences:
  kinase.table[, 2] <-
    trimSeqs(kinase.table[, 2], 
              seq.length = substrate.length, 
              verbose = verbose)
  kinase.table[, 2] <- toupper(kinase.table[, 2])
  
  #remove peptides with a certain center character:
  if (remove.center != FALSE) {
    half.window <- (substrate.length + 1) / 2
    center.aa <-
      substr(kinase.table[, 2], half.window, half.window)
    if (verbose) {
      print(
        paste(
          "You have selected to remove",
          length(which(center.aa == toupper(remove.center))),
          "peptide sequences for building PWMs that contain a centered letter 
          of:",
          toupper(remove.center),
          sep = " "
        )
      )
    }
    kinase.table <-
      kinase.table[center.aa != toupper(remove.center),]
  }
  
  kinase.table <- unique(kinase.table)
  #----------------------------------------------
  
  if (verbose) {
    cat(paste(
      nrow(kinase.table),
      "unique kinase:substrate sequences in table provided.",
      sep = " "
    ))
  }
  
  #initialise table for kinase count data
  kinase.summary <-
    data.frame("kinase" = unique(as.character(kinase.table[, 1])), "n" = NA)

  kinase.summary$n <-
    c(sapply(seq_len(nrow(kinase.summary)), function(i)
      length(kinase.table[, 1][kinase.table[, 1] == kinase.summary[i, 1]])))
  
  #filter the list here which is used to build pwms:
  kinase.summary <-
    kinase.summary[kinase.summary$n >= substrates.n, ]
  
  #generate PWM list:
  pwm.list <-
    c(lapply(seq_len(nrow(kinase.summary)), function(i)
      scorePWM(
        as.character(kinase.table[, 2][kinase.table[, 1] == 
                                         kinase.summary[i, 1]]),
        substrate.length = substrate.length,
        wild.card = wild.card,
        pseudo = pseudo
      )))
  
  return(list("pwm" = pwm.list, "kinase" = kinase.summary))
}

# this helper function performs the calculations for building the PWMs

scorePWM <-
  function(input.data,
           substrate.length,
           wild.card = "_",
           pseudo = 0.01) {
    #reformat input.data; build a matrix of the sequences (split into individual
    # AA's) = substrate.length
    input.data <-
      data.frame(matrix(unlist(strsplit(input.data , "")) , 
                        ncol = substrate.length , byrow = TRUE))
    uniq.AA <-
      c(
        wild.card,
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y"
      )
    #1. Generate PFM:
    pwm <-
      t(cbind(sapply(seq_len(length(uniq.AA)), function(i)
        (
          sapply(seq_len(ncol(input.data)), function(j)
            sum(
              ifelse(as.character(input.data[, j]) == 
                       as.character(uniq.AA[i]), 1, 0)
            ))
        ))))
    rownames(pwm) <- uniq.AA
    colnames(pwm) <- paste("p", seq_len(ncol(pwm)), sep = "")
    #remove "_"
    #SET AS AN OPTION:
    pwm <- pwm[!rownames(pwm) %in% c(wild.card),]
    #2. PPM calculation
    pwm <- (pwm + pseudo) / (apply(pwm, 2, sum, na.rm = TRUE))
    #3. Generate Position Weight Matrix
    pwm <- log(pwm / (1 / nrow(pwm)))
    return(pwm)
  }
