#' Score substrate sequences for matches to kinase Position Weight Matrices 
#' (PWMs)
#'
#' @description Scores each input sequence for a match against all PWMs provided
#'  from buildPWM() and generates p-values for scores. The output of this 
#'  function is to be used for building the swing metric, the predicted activity
#'   of kinases.
#' @param input_data A data.frame of phoshopeptide data. Must contain 4 columns
#'  and the following format must be adhered to. Column 1 - Annotation, Column 2
#'   - centered peptide sequence, Column 3 - Fold Change [-ve to +ve], Column 4 
#'   - p-value [0-1]
#' @param pwm_in List of PWMs created using buildPWM()
#' @param background Option to provide a data.frame of peptides to use as 
#' background. If providing a background as a table, this must contain two 
#' columns; Column 1 - Annotation, Column 2 - centered peptide sequence. These 
#' must be centered. OR generate a random background for PWM scoring from the
#' input list - background = random. Default: "random"
#' @param n Number of permutations to perform for generating background. 
#' Default: "1000"
#' @param force_trim This function will detect if a peptide sequence is of 
#' different length to the PWM models generated (provided in pwm_in) and trim 
#' the input sequences to the same length as the PWM models. If a background is 
#' provided, this will also be trimmed to the same width as the PWM models. 
#' Options are: "TRUE, FALSE". Default = FALSE
#' @param verbose Turn verbosity on/off. To turn on, verbose=TRUE. Options are: 
#' "TRUE, FALSE". Default = FALSE
#'
#' @examples
#' ## import data
#' data(example_phosphoproteome)
#' data(phosphositeplus_human)
#'
#' ## clean up the annotations
#' ## sample 100 data points for demonstration
#' sample_data <- head(example_phosphoproteome, 100)
#' annotated_data <- cleanAnnotation(input_data = sample_data)
#'
#' ## build the PWM models:
#' set.seed(1234)
#' sample_pwm <- phosphositeplus_human[sample(nrow(phosphositeplus_human), 
#' 1000),]
#' pwms <- buildPWM(sample_pwm)
#'
#' ## score the PWM - substrate matches
#' ## Using a "random" background, to calculate the p-value of the matches
#' ## Using n=10 for demonstration
#' ## set.seed for reproducibility
#' set.seed(1234)
#' substrate_scores <- scoreSequences(input_data = annotated_data,
#'                                    pwm_in = pwms,
#'                                    background = "random",
#'                                    n = 10)
#'
#' @return A list with 3 elements: 1) PWM-substrate scores: 
#' substrate_scores$peptide_scores, 2) PWM-substrate p-values: 
#' substrate_scores$peptide_p 3) Background used for reproducibility: 
#' substrate_scores$background 4) input_data is returned in the case that it was
#' trimmed.
#'
#' @export scoreSequences
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam


scoreSequences <- function(input_data = NULL,
                            pwm_in = NULL,
                            background = "random",
                            n = 1000,
                            force_trim = FALSE,
                            verbose = FALSE) {
  if (verbose) {
    message("Scoring,", nrow(input_data),
        "sequences provided, against", length(pwm_in[[1]]),
        "PWM models")
    message("Parameters selected:\n", 
            "n =", n, "\n",
            "force_trim =", force_trim, "\n",
            "verbose =", verbose, "\n"
           )
  }
  
  #----------------------------------------------
  #format checks:
  if (is.null(input_data))
    stop("input_data not provided; you must provide an input table")
  if (!is.data.frame(input_data))
    stop("input_data is not a data.frame; you must provide an
         input table")
  if (is.null(pwm_in))
    stop(
      "pwm_in not provided; you must provide an input table containing
      computed position weight matrices using buildPWM()"
    )
  if (!is.list(pwm_in))
    stop(
      "pwm_in is not a list format; something has gone wrong. Make sure
      you compute the position weight matrices using buildPWM()"
    )
  #min peptide sequence length of input data
  min_peptide_seq <- min(nchar(as.character(input_data[, 2])))
  #length of PWM models
  min_pwm_length <- min(sapply(seq_len(length(pwm_in[[1]])), function(X)
    ncol(data.frame(pwm_in[[1]][X]))))
  if ((min_peptide_seq > min_pwm_length) & force_trim == FALSE)
    stop(
      "Centered peptide sequence is greater than length of position weight 
       matrix. Check data and/or consider using force_trim=TRUE."
    )
  if ((min_peptide_seq < min_pwm_length) & force_trim == FALSE)
    stop(
      paste(
        "Centered peptide sequence is less than the length of PWMs.
        Rebuild PWM with length of,",
        min_peptide_seq,
        ",and rerun.",
        sep = ""
      )
    )
  if ((min_peptide_seq > min_pwm_length) & force_trim == TRUE) {
    if (verbose) {
      message("trimming input_data sequences to minimum PWM length, which is,",
              min_pwm_length
              )
    }
    #trim  seqs to the min. sequence length in PWMs
    input_data[, 2] <-
      trimSeqs(
        seqs_to_trim = as.character(input_data[, 2]),
        seq_length = min_pwm_length,
        verbose = verbose
      )
  }
  if (background != "random" && !is.data.frame(background)) {
    stop(
      "Background is not in data.frame() format. Use backgound=random or
      provide a table of peptides as background in data.frame format with
      two columns, Column1 = Annotation, Column2 = Centered peptide"
    )
  }
  if (background != "random" && is.data.frame(background)) {
    #remove NA's
    background <- background[!is.na(background[, 2]),]
    #check peptide sequence length of background data if provided:
    min_peptide_seq_bg <- min(nchar(as.character(background[, 2])))
    if ((min_peptide_seq_bg > min_pwm_length) & force_trim == FALSE)
      stop(
        "Centered peptide sequence of BACKGROUND is greater than width of
        position weight matrix. Set option to force trim to continue."
      )
    if ((min_peptide_seq_bg < min_pwm_length) & force_trim == FALSE)
      stop(
        "Centered peptide sequence if BACKGROUND is less than length of
        position weight matrix. Rebuild PWM with length of X and rerun."
      )
    if (n > nrow(background))
      stop("Background provided is smaller than n selected. Reduce n")
    
    if ((min_peptide_seq_bg > min_pwm_length) & force_trim == FALSE) {
      if (verbose) {
        message("trimming BACKGROUND sequences provided to minimum PWM length,
            which is,",
            min_pwm_length
            )
      }
      #trim seqs to the min. sequence length inPWMs
      background[, 2] <-
        trimSeqs(
          seqs_to_trim = as.character(background[, 2]),
          seq_length = min_pwm_length,
          verbose = verbose
        )
    }
    
  }
  #check "n" not smaller than input dataset
  if (background == "random" && n > nrow(input_data)) {
    stop("n is larger than dataset. Reduce n random sampling")
  }
  
  #----------------------------------------------
  #1. for each peptide, score against all kinase PWMs provided:
  if (verbose) {
    message("[Step1/3] : Scoring peptide sequences against PWMs")
  }
  
  peptide_scores <-
    bplapply(seq_along(pwm_in[[2]]$kinase), function(i)
      sapply(seq_len(nrow(input_data)), function(j)
        seqScore(as.character(input_data[j, 2]), pwm_in[[1]][[i]])))
  
  names(peptide_scores) <- pwm_in[[2]]$kinase
  peptide_scores <-
    cbind("annotation" = input_data[, 1],
          "peptide" = input_data[, 2],
          data.frame(peptide_scores))
  
  #2. extract background set of size "n"
  # if no background is provided:
  if (background == "random" && !is.data.frame(background)) {
    if (verbose) {
      message("[Step2/3] : Random background generated from ", n, " peptides")
    }
    rand_bg <- sample(rownames(peptide_scores), n)
    background_scores <-
      peptide_scores[rownames(peptide_scores) %in% rand_bg, ]
  }
  # or if background is provided
  if (background != "random" && is.data.frame(background)) {
    if (verbose) {
      message("[Step2/3] : Background provided; calculating PWM scores against 
              background generated from ", n, " peptides"
              )
    }
    rand_bg <- sample(rownames(background), n)
    rand_bg <- background[rownames(background) %in% rand_bg, ]
    
    background_scores <-
      bplapply(seq_along(pwm_in[[2]]$kinase), function(i)
        sapply(seq_len(nrow(rand_bg)), function(j)
          seqScore(as.character(rand_bg[j, 2]), pwm_in[[1]][[i]])))
    
    names(background_scores) <- pwm_in[[2]]$kinase
    background_scores <- cbind(
      "annotation" = rand_bg[, 1],
      "peptide" = rand_bg[, 2],
      data.frame(background_scores)
    )
    
  }
  
  #3. compute p-values:
  if (verbose) {
    message("[Step3/3] : Computing p-values for peptide scores:PWM scores")
  }
  
  # multi-core
  p_table <- bplapply(3:ncol(peptide_scores), function(i)
    sapply(seq_len(length(peptide_scores[, i])), function(j)
      (sum(
        ifelse(
          as.numeric(background_scores[, i]) >
            as.numeric(peptide_scores[, i][j]),
          1, 0 )) + 1) / (length(as.numeric(background_scores[, i])) + 1)))
  
  names(p_table) <- colnames(3:ncol(peptide_scores))
  p_table <- data.frame(peptide_scores[1:2], data.frame(p_table))
  colnames(p_table) <- gsub("score_", "p_", colnames(peptide_scores))
  rownames(p_table) <- rownames(peptide_scores)
  
  if (verbose) {
    message("[Finished]")
  }
  return(
    list(
      "peptide_scores" = peptide_scores,
      "peptide_p" = p_table,
      "background" = background_scores,
      "input_data" = input_data
    )
  )
}

# helper function to calculate PWM match scores (vectorised)
seqScore <- function(input_seq, pwm) {
  # vector of letters of sequence to score
  input_seq <- unlist(strsplit(input_seq, ""))
  idx <- input_seq != "_"
  ridx <- match(input_seq[idx], rownames(pwm))
  cidx <- seq_along(input_seq)[idx]
  # return sum of scores for each matched AA
  # create a row/column index (which row matches each AA in string)
  # extracts the corresponding scores and then sums for the total score
  sum(pwm[cbind(ridx, cidx)])
}