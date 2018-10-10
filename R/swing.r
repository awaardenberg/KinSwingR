#' Swing statistic
#'
#' @description This function integrates the kinase-substrate predictions,
#' directionality of phosphopeptide fold change and signficance to assess
#' local connectivity (swing) of kinase-substrate networks. The final score
#' is a normalised and weighted score of predicted kinase activity. If 
#' permutations are selected, network node:edges are permutated. P-values will 
#' be calculated for both ends of the distribution of swing scores (positive and
#'  negative swing scores).
#' @param input.data A data.frame of phoshopeptide data. Must contain 4 columns
#'  and the following format must be adhered to. Column 1 - Annotation, Column 2
#'   - centered peptide sequence, Column 3 - Fold Change [-ve to +ve], Column 4 
#'   - p-value [0-1]. This must be the same dataframe used in scoreSequences()
#' @param pwm.in List of PWMs created using buildPWM()
#' @param pwm.scores List of PWM-substrate scores created using
#' score.seqeuences()
#' @param pseudo.count Pseudo-count for avoiding log-zero transformations.
#' Default: "1"
#' @param p.cut.pwm Significance level for determining a significant
#' kinase-substrate enrichment. Default: "0.05"
#' @param p.cut.fc Significance level for determining a significant level of
#' Fold-change in the phosphoproteomics data. Default: "0.05"
#' @param permutations Number of permutations to perform. This will shuffle the 
#' kinase-subtrate edges of the network n times. To not perform permutations and
#'  only generate the scores, set permutations=1 or permutations=FALSE. Default:
#'   "1000"
#' @param verbose Turn verbosity on/off. To turn on, verbose=TRUE. Options are: 
#' "TRUE, FALSE". Default=FALSE
#' @param threads Number of processing cores to use. Default: "1"
#'
#' @examples
#' ## import data
#' data(example_phosphoproteome)
#' data(phosphositeplus_human)
#'
#' ## clean up the annotations
#' ## sample 100 data points for demonstration
#' sample_data <- head(example_phosphoproteome, 100)
#' annotated_data <- cleanAnnotation(input.data = sample_data)
#'
#' ## build the PWM models:
#' set.seed(1234)
#' sample_pwm <- phosphositeplus_human[sample(nrow(phosphositeplus_human), 
#' 1000),]
#' pwms <- buildPWM(sample_pwm)
#'
#' ## score the PWM - substrate matches
#' ## Using a "random" background, to calculate the p-value of the matches
#' ## Using n = 10 for demonstration
#' ## set.seed for reproducibility
#' set.seed(1234)
#' substrate_scores <- scoreSequences(input.data = annotated_data,
#'                                    pwm.in = pwms,
#'                                    background = "random",
#'                                    n = 10,
#'                                    threads = 4)
#'
#' ## Use the pwm.scores and annotated data to predict kinase activity.
#' ## This will permute the network node and edges 10 times for demonstration.
#' ## set.seed for reproducibility
#' set.seed(1234)
#' swing_output <- swing(input.data = annotated_data,
#'                       pwm.in = pwms,
#'                       pwm.scores = substrate_scores,
#'                       permutations = 10,
#'                       threads = 4)
#'
#' @return A data.table of swing scores
#'
#' @export swing
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam
#' @importFrom stats setNames
#' @importFrom stats sd
#' @importFrom data.table melt.data.table

swing <-
  function(input.data = NULL,
           pwm.in = NULL,
           pwm.scores = NULL,
           pseudo.count = 1,
           p.cut.pwm = 0.05,
           p.cut.fc = 0.05,
           permutations = 1000,
           verbose = FALSE,
           threads = 1) {
    #----------------------------------------------
    #format checks:
    if (is.null(input.data))
      stop("input.data not provided; you must provide an input table")
    if (!is.data.frame(input.data))
      stop("input.data is not a data.frame; you must provide an input table")
    if (is.null(pwm.in))
      stop(
        "pwm.in not provided; you must provide an input table containing
        computed position weight matrices using buildPWM()"
      )
    if (!is.list(pwm.in))
      stop(
        "pwm.in is not a list format; something has gone wrong. Make sure you
        compute the position weight matrices using buildPWM()"
      )
    if (is.null(pwm.scores))
      stop(
        "pwm.scores not provided; you must provide an input table containing
        computed scores using scores.sequences()"
      )
    if (!is.list(pwm.scores))
      stop(
        "pwm.scores is not a list format; something has gone wrong. Make sure
        PWM-substrate matches are scored using scoreSequences()"
      )
    if (p.cut.pwm >= 1)
      stop("p.cut.pwm needs to be less than 1; make sure your p-values are not
           log transformed")
    if (p.cut.fc >= 1)
      stop("p.cut.fc needs to be less than 1; make sure your p-values are not
           log transformed")
    if (permutations != FALSE &
        (!is.numeric(permutations) | permutations < 1))
      stop(
        "permutations needs to be a numerical number. Either 1 or FALSE (for no 
        permutation) or greater than 1 to set the number of pemutations"
      )
    
    #----------------------------------------------
    #1. binarise the p-values and the fold change:
    if (verbose) {
      start_time <- Sys.time()
      cat(paste("Start: ", start_time, "\n", sep = ""))
      cat("[Step1/3] : Calculating Swing Scores\n")
    }
    
    pwm.pval <- pwm.scores[[2]]
    #binarise the p-values:
    pwm.pval[, 3:ncol(pwm.pval)] <-
      ifelse(pwm.pval[, 3:ncol(pwm.pval)] > p.cut.pwm, 0, 1)
    #binarise FC:
    input.data[, 3] <- ifelse(as.numeric(input.data[, 3]) > 0, 1, -1)
    #binarise p-value:
    input.data[, 4] <-
      ifelse(as.numeric(input.data[, 4]) > p.cut.fc, 0, 1)
    #compute Sipk statistic (Sipk - see paper):
    input.data$Sipk <-
      as.numeric(input.data[, 3]) * as.numeric(input.data[, 4])
    
    #2.merge the tables together, summarise counts and unique
    input.data$annotation <-
      paste(input.data$annotation, input.data$peptide, sep = "::")
    pwm.pval$annotation <-
      paste(pwm.pval$annotation, pwm.pval$peptide, sep = "::")
    data.merge <-
      unique(merge(input.data, pwm.pval, by = "annotation"))
    
    #3. Swing scores of REAL data.
    swing.out <- swingScore(
      data.merge = data.merge,
      pwm.in = pwm.in,
      permute = FALSE,
      pseudo.count = pseudo.count
    )
    
    #4. permute data and calculate p-values
    if (permutations != FALSE && permutations > 1) {
      if (verbose) {
        cat("[Step2/3] : Permuting Network\n")
      }
      n.permute <- lapply(seq_len(permutations), function(i)
        sample(as.character(colnames(data.merge))[7:ncol(data.merge)],
               length(colnames(data.merge[, 7:ncol(data.merge)])),
               replace = FALSE))
      # not really required
      names(n.permute) <- paste("rand", seq_len(permutations), sep = "")
      # calculate swing scores and return in a single table
      swing.permute <- list(bplapply(seq_len(permutations), function(i)
        swingScore(
          data.merge = setNames(data.frame(data.merge),
                                c(colnames(data.merge)[1:6],
                                  n.permute[[i]])),
          pwm.in = pwm.in,
          permute = TRUE,
          pseudo.count = pseudo.count,
          n = i
        ),
        BPPARAM = MulticoreParam(workers =
                                   threads)))
      
      # uses data.table library to merge:
      swing.permute <-
        data.frame(Reduce(merge, swing.permute[[1]]))
      
      if (verbose) {
        cat("[Step3/3] : Calculating p-values\n")
      }
      #obtain p-values two sided, test independently...
      swing.out$p.greater <-
        unlist(bplapply(seq_len(nrow(swing.out)), function(i)
          (sum(
            ifelse(
              as.numeric(swing.permute[swing.permute$kinase ==
                         swing.out$kinase[i], ][2:ncol(swing.permute)]) >
                         as.numeric(swing.out$swing.raw[i]), 1, 0), 
                         na.rm = TRUE) + 1)
                         / (as.numeric(permutations) + 1),
                         BPPARAM = MulticoreParam(workers =
                         threads)))
      
      swing.out$p.less <-
        unlist(bplapply(seq_len(nrow(swing.out)), function(i)
          (sum(
            ifelse(
              as.numeric(swing.permute[swing.permute$kinase ==
                         swing.out$kinase[i], ][2:ncol(swing.permute)]) <
                         as.numeric(swing.out$swing.raw[i]), 1, 0
                         ), na.rm = TRUE) + 1)
                         / (as.numeric(permutations) + 1),
                         BPPARAM = MulticoreParam(workers =
                         threads)))

    }
    
    if (verbose) {
      cat("[FINISHED]\n")
      end_time <- Sys.time() - start_time
      cat(paste("Finish: ", Sys.time(), "\n", sep = ""))
      end_time#print(paste(end_time, "\n", sep=""))
    }
    
    #5. return interaction table - long
    network <- data.merge
    #convert to long table:
    network <- data.table::melt(network,
                               id.vars = c("annotation"),
                               measure.vars = 
                               colnames(network[7:ncol(network)]))
    # keep "significant" edges
    network <- network[network$value == 1, ]
    network <- data.frame(
      "source" = as.character(network$variable),
      "target" = as.character(network$annotation)
    )
    return(list("scores" = swing.out, "network" = network))
  }

# describeIn swing This helper function performs the swing score calculation
# no verbose in this helper
# permute - for setting right table format (using data.table) for merging.

swingScore <-
  function(data.merge,
           pwm.in,
           pseudo.count = 1,
           permute = FALSE,
           n = 1) {
    #----------------------------------------------
    #count significant positive
    p.k <-
      sapply(7:ncol(data.merge), function (i)
        sum(data.merge$Sipk * data.merge[i] == 1))
    #count significant negative
    n.k <-
      sapply(7:ncol(data.merge), function (i)
        sum(data.merge$Sipk * data.merge[i] == -1))
    #all counts (both positive and negative).
    t.k <- p.k + n.k
    #----------------------------------------------
    #compute swing statistics:
    p.n.all <-
      data.frame(
        "kinase" = colnames(data.merge[7:ncol(data.merge)]),
        "pos" = p.k,
        "neg" = n.k,
        "all" = t.k
      )
    #compute stats:
    #proportion positive/negative:
    p.n.all$pk <- p.n.all$pos / p.n.all$all
    p.n.all$nk <- p.n.all$neg / p.n.all$all
    p.n.all$swing.raw <-
      (p.n.all$pk + pseudo.count) / (p.n.all$nk + pseudo.count)
    #include the count numbers:
    p.n.all <- merge(p.n.all, pwm.in$kinase, by = "kinase")
    #weighted by number of substrates:
    p.n.all$swing.raw <-
      log2(p.n.all$swing.raw) * log2(p.n.all$n + pseudo.count)
    #weighted by size of kinase-substrate network:
    p.n.all$swing.raw <-
      p.n.all$swing.raw * log2(p.n.all$all + pseudo.count)
    #z-score transform
    p.n.all$swing <-
      (p.n.all$swing.raw - mean(p.n.all$swing.raw, na.rm = TRUE)) / 
      sd(p.n.all$swing.raw, na.rm = TRUE)
    #order results and return table::
    p.n.all <- p.n.all[order(p.n.all$swing, decreasing = TRUE),]
    
    # for permuted tables, return the scores in
    # data.table format for merge function!
    if (permute == "TRUE") {
      p.n.all <- data.frame(p.n.all$kinase, p.n.all$swing.raw)
      colnames(p.n.all) <- c("kinase", paste("swing.raw.", n, sep = ""))
      p.n.all = data.table(p.n.all, key = 'kinase')
    }
    
    return(p.n.all)
  }
