#' Score substrate sequences for matches to kinase Position Weight Matrices 
#' (PWMs)
#'
#' @description Scores each input sequence for a match against all PWMs provided
#'  from build.pwm() and generates p-values for scores. The output of this 
#'  function is to be used for building the swing metric, the predicted activity
#'   of kinases.
#' @param input.data A data.frame of phoshopeptide data. Must contain 4 columns
#'  and the following format must be adhered to. Column 1 - Annotation, Column 2
#'   - centered peptide sequence, Column 3 - Fold Change [-ve to +ve], Column 4 
#'   - p-value [0-1]
#' @param pwm.in List of PWMs created using build.pwm()
#' @param background Option to provide a data.frame of peptides to use as 
#' background. If providing a background as a table, this must contain two 
#' columns; Column 1 - Annotation, Column 2 - centered peptide sequence. These 
#' must be centered. OR generate a random background for PWM scoring from the
#' input list - background = random. Default: "random"
#' @param n Number of permutations to perform for generating background. 
#' Default: "1000"
#' @param force.trim This function will detect if a peptide sequence is of 
#' different length to the PWM models generated (provided in pwm.in) and trim 
#' the input sequences to the same length as the PWM models. If a background is 
#' provided, this will also be trimmed to the same width as the PWM models. 
#' Options are: "TRUE, FALSE". Default = FALSE
#' @param seed This is for reproducible results for permutation testing. To not 
#' use, set seed = NULL. Default: "1234"
#' @param verbose Turn verbosity on/off. To turn on, verbose=TRUE. Options are: 
#' "TRUE, FALSE". Default = FALSE
#' @param threads Number of processing cores to use. Default: "1"
#'
#' @examples
#' ## import data
#' data(example_phosphoproteome)
#' data(phosphositeplus_human)
#'
#' ## clean up the annotations
#' ## sample 100 data points for demonstration
#' sample.data <- head(example_phosphoproteome, 100)
#' annotated.data <- clean.annotation(input.data=sample.data)
#'
#' ## build the PWM models:
#' set.seed(1)
#' sample.pwm <- phosphositeplus_human[sample(nrow(phosphositeplus_human), 
#' 1000),]
#' kinase.pwm.models <- build.pwm(sample.pwm)
#'
#' ## score the PWM - substrate matches
#' ## Using a "random" background, to calculate the p-value of the matches
#' ## Using n=10 for demonstration
#'
#' pwm.substrate.scores <- score.sequences(input.data = annotated.data,
#'                                        pwm.in = kinase.pwm.models,
#'                                        background = "random",
#'                                        n = 10,
#'                                        threads = 4)
#'
#' @return A list with 3 elements: 1) PWM-substrate scores: 
#' pwm.substrate.scores$peptide.scores, 2) PWM-substrate p-values: 
#' pwm.substrate.scores$peptide.p and 3) Background used for reproducibility: 
#' pwm.substrate.scores$background
#'
#' @export score.sequences
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam

score.sequences <- function(input.data = NULL,
                            pwm.in = NULL,
                            background = "random",
                            n = 1000,
                            force.trim = FALSE,
                            seed = 1234,
                            verbose = FALSE,
                            threads = 1) {
  if (verbose) {
    cat(
      paste(
        "Scoring,",
        nrow(input.data),
        "sequences provided, against",
        length(pwm.in[[1]]),
        "PWM models\n",
        sep = " "
      )
    )
    cat(
      paste(
        "Parameters selected:\n",
        "n=",
        n,
        "force.trim=",
        force.trim,
        "seed=",
        seed,
        "verbose=",
        verbose,
        "\n",
        sep = ","
      )
    )
  }
  
  set.seed(seed)
  #----------------------------------------------
  #format checks:
  if (is.null(input.data))
    stop("input.data not provided; you must provide an input table")
  if (!is.data.frame(input.data))
    stop("input.data is not a data.frame; you must provide an
         input table")
  if (is.null(pwm.in))
    stop(
      "pwm.in not provided; you must provide an input table containing
      computed position weight matrices using build.pwm()"
    )
  if (!is.list(pwm.in))
    stop(
      "pwm.in is not a list format; something has gone wrong. Make sure
      you compute the position weight matrices using build.pwm()"
    )
  #min peptide sequence length of input data
  min.peptide.seq <- min(nchar(as.character(input.data[, 2])))
  #length of PWM models
  min.pwm.length <- min(sapply(seq_len(length(pwm.in[[1]])), function(X)
    ncol(data.frame(pwm.in[[1]][X]))))
  if ((min.peptide.seq > min.pwm.length) & force.trim == FALSE)
    stop(
      "Centered peptide sequence is greater than length of position weight 
       matrix. Check data and/or consider using force.trim=TRUE."
    )
  if ((min.peptide.seq < min.pwm.length) & force.trim == FALSE)
    stop(
      paste(
        "Centered peptide sequence is less than the length of PWMs.
        Rebuild PWM with length of,",
        min.peptide.seq,
        ",and rerun.",
        sep = ""
      )
    )
  if ((min.peptide.seq > min.pwm.length) & force.trim == TRUE) {
    if (verbose) {
      cat(
        paste(
          "trimming input.data sequences to minimum PWM length, which is,",
          min.pwm.length,
          "\n",
          sep = " "
        )
      )
    }
    #trim  seqs to the min. sequence length in PWMs
    input.data[, 2] <-
      trim.seqs(
        seqs.to.trim = as.character(input.data[, 2]),
        seq.length = min.pwm.length,
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
    min.peptide.seq.bg <- min(nchar(as.character(background[, 2])))
    if ((min.peptide.seq.bg > min.pwm.length) & force.trim == FALSE)
      stop(
        "Centered peptide sequence of BACKGROUND is greater than width of
        position weight matrix. Set option to force trim to continue."
      )
    if ((min.peptide.seq.bg < min.pwm.length) & force.trim == FALSE)
      stop(
        "Centered peptide sequence if BACKGROUND is less than length of
        position weight matrix. Rebuild PWM with length of X and rerun."
      )
    if (n > nrow(background))
      stop("Background provided is smaller than n selected. Reduce n")
    
    if ((min.peptide.seq.bg > min.pwm.length) & force.trim == FALSE) {
      if (verbose) {
        cat(
          paste(
            "trimming BACKGROUND sequences provided to minimum PWM length,
            which is,",
            min.pwm.length,
            "\n",
            sep = " "
          )
        )
      }
      #trim seqs to the min. sequence length inPWMs
      background[, 2] <-
        trim.seqs(
          seqs.to.trim = as.character(background[, 2]),
          seq.length = min.pwm.length,
          verbose = verbose
        )
    }
    
  }
  #check "n" not smaller than input dataset
  if (background == "random" && n > nrow(input.data)) {
    stop("n is larger than dataset. Reduce n random sampling")
  }
  
  #----------------------------------------------
  #1. for each peptide, score against all kinase PWMs provided:
  if (verbose) {
    cat("[Step1/3] : Scoring peptide sequences against PWMs\n")
  }
  peptide.scores <-
    bplapply(seq_len(length(pwm.in[[2]]$kinase)), function(i)
      sapply(seq_len(nrow(input.data)), function(j)
        seq.score(as.character(input.data[j, 2]), pwm.in[[1]][[i]]))
      , BPPARAM = MulticoreParam(workers = threads))
  names(peptide.scores) <- pwm.in[[2]]$kinase
  peptide.scores <-
    cbind("annotation" = input.data[, 1],
          "peptide" = input.data[, 2],
          data.frame(peptide.scores))
  
  #2. extract background set of size "n"
  # if no background is provided:
  if (background == "random" && !is.data.frame(background)) {
    if (verbose) {
      cat(paste(
        "[Step2/3] : Random background generated from ",
        n,
        " peptides\n",
        sep = ""
      ))
    }
    rand.bg <- sample(rownames(peptide.scores), n)
    background.scores <-
      peptide.scores[rownames(peptide.scores) %in% rand.bg, ]
  }
  # or if background is provided
  if (background != "random" && is.data.frame(background)) {
    if (verbose) {
      cat(
        paste(
          "[Step2/3] : Background provided;
          calculating PWM scores against background generated from ",
          n ,
          " peptides",
          "\n",
          sep = ""
        )
      )
    }
    rand.bg <- sample(rownames(background), n)
    rand.bg <- background[rownames(background) %in% rand.bg, ]
    background.scores <-
      bplapply(seq_len(length(pwm.in[[2]]$kinase)), function(i)
        sapply(seq_len(nrow(rand.bg)), function(j)
          seq.score(as.character(rand.bg[j, 2]), pwm.in[[1]][[i]]))
        , BPPARAM = MulticoreParam(workers = threads))
    
    names(background.scores) <- pwm.in[[2]]$kinase
    background.scores <- cbind(
      "annotation" = rand.bg[, 1],
      "peptide" = rand.bg[, 2],
      data.frame(background.scores)
    )
    
  }
  
  #3. compute p-values:
  if (verbose) {
    cat("[Step3/3] : Computing p-values for peptide scores:PWM scores\n")
  }
  
  #multi-core:
  p.table <- bplapply(3:ncol(peptide.scores), function(i)
    sapply(seq_len(length(peptide.scores[, i])), function(j)
      (sum(
        ifelse(
          as.numeric(background.scores[, i]) >
            as.numeric(peptide.scores[, i][j]),
          1,
          0
        )
      ) + 1) /
        (length(
          as.numeric(background.scores[, i])
        ) + 1))
    , BPPARAM = MulticoreParam(workers = threads))
  
  names(p.table) <- colnames(3:ncol(peptide.scores))
  p.table <- data.frame(peptide.scores[1:2], data.frame(p.table))
  colnames(p.table) <- gsub("score.", "p.", colnames(peptide.scores))
  rownames(p.table) <- rownames(peptide.scores)
  
  if (verbose) {
    cat("[Finished]\n")
  }
  return(
    list(
      "peptide.scores" = peptide.scores,
      "peptide.p" = p.table,
      "background" = background.scores
    )
  )
  }

# This helper function calculates PWM matching scores
seq.score <- function(input.seq, pwm) {
  #vector of letters of sequence to score.
  input.seq <- unlist(strsplit(input.seq, ""))
  #function to score ther matching AA to PWM position:
  score.AA <- function(i) {
    if (input.seq[i] == "_") {
      score = 0
    }
    if (input.seq[i] != "_") {
      score <- unlist(pwm[rownames(pwm) == input.seq[i], ][i])
    }
    #return score for each matched AA
    return(score)
  }
  return(sum(sapply(seq_len(length(input.seq)), function(X)
    score.AA(X))))
}
