#' Swing statistic
#'
#' @description This function integrates the kinase-substrate predictions,
#' directionality of phosphopeptide fold change and signficance to assess
#' local connectivity (swing) of kinase-substrate networks. The final score
#' is a normalised and weighted score of predicted kinase activity. If 
#' permutations are selected, network node:edges are permutated. P-values will 
#' be calculated for both ends of the distribution of swing scores (positive and
#'  negative swing scores).
#' @param input_data A data.frame of phoshopeptide data. Must contain 4 columns
#'  and the following format must be adhered to. Column 1 - Annotation, Column 2
#'   - centered peptide sequence, Column 3 - Fold Change [-ve to +ve], Column 4 
#'   - p-value [0-1]. This must be the same dataframe used in scoreSequences()
#' @param pwm_in List of PWMs created using buildPWM()
#' @param pwm_scores List of PWM-substrate scores created using
#' scoreSequences()
#' @param pseudo_count Pseudo-count acts at two levels. 1) It adds a small
#' number to the counts to avoid zero divisions, which also 2) avoids log-zero 
#' transformations. Note that this means that pos, neg and all values in the
#' output table include the addition of the pseudo-count. Default: "1"
#' @param p_cut_pwm Significance level for determining a significant
#' kinase-substrate enrichment. Default: "0.05"
#' @param p_cut_fc Significance level for determining a significant level of
#' Fold-change in the phosphoproteomics data. Default: "0.05"
#' @param permutations Number of permutations to perform. This will shuffle the 
#' kinase-subtrate edges of the network n times. To not perform permutations and
#'  only generate the scores, set permutations=1 or permutations=FALSE. Default:
#'   "1000"
#' @param return_network Option to return an interaction network for visualising
#' in cystoscape. Default = FALSE
#' @param verbose Turn verbosity on/off. To turn on, verbose=TRUE. Options are: 
#' "TRUE, FALSE". Default=FALSE
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
#' ## Using n = 100 for demonstration
#' ## set.seed for reproducibility
#' set.seed(1234)
#' substrate_scores <- scoreSequences(input_data = annotated_data,
#'                                    pwm_in = pwms,
#'                                    background = "random",
#'                                    n = 100)
#'
#' ## Use substrate_scores and annotated_data data to predict kinase activity.
#' ## This will permute the network node and edges 10 times for demonstration.
#' ## set.seed for reproducibility
#' set.seed(1234)
#' swing_output <- swing(input_data = annotated_data,
#'                       pwm_in = pwms,
#'                       pwm_scores = substrate_scores,
#'                       permutations = 10)
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
  function(input_data = NULL,
           pwm_in = NULL,
           pwm_scores = NULL,
           pseudo_count = 1,
           p_cut_pwm = 0.05,
           p_cut_fc = 0.05,
           permutations = 1000,
           return_network = FALSE,
           verbose = FALSE) {
    #----------------------------------------------
    #format checks:
    if (is.null(input_data))
      stop("input_data not provided; you must provide an input table")
    if (!is.data.frame(input_data))
      stop("input_data is not a data.frame; you must provide an input table")
    if (is.null(pwm_in))
      stop(
        "pwm_in not provided; you must provide an input table containing
        computed position weight matrices using buildPWM()"
      )
    if (!is.list(pwm_in))
      stop(
        "pwm_in is not a list format; something has gone wrong. Make sure you
        compute the position weight matrices using buildPWM()"
      )
    if (is.null(pwm_scores))
      stop(
        "pwm_scores not provided; you must provide an input table containing
        computed scores using scoreSequences()"
      )
    if (!is.list(pwm_scores))
      stop(
        "pwm_scores is not a list format; something has gone wrong. Make sure
        PWM-substrate matches are scored using scoreSequences()"
      )
    if (p_cut_pwm >= 1)
      stop("p_cut_pwm needs to be less than 1; make sure your p-values are not
           log transformed")
    if (p_cut_fc >= 1)
      stop("p_cut_fc needs to be less than 1; make sure your p-values are not
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
      message("Start: ", start_time)
      message("[Step1/3] : Calculating Swing Scores")
    }
    
    pwm_pval <- pwm_scores[[2]]
    #binarise the p-values:
    pwm_pval[, 3:ncol(pwm_pval)] <-
      ifelse(pwm_pval[, 3:ncol(pwm_pval)] > p_cut_pwm, 0, 1)
    #binarise FC:
    input_data[, 3] <- ifelse(as.numeric(input_data[, 3]) > 0, 1, -1)
    #binarise p-value:
    input_data[, 4] <-
      ifelse(as.numeric(input_data[, 4]) > p_cut_fc, 0, 1)
    #compute Sipk statistic (Sipk - see paper):
    input_data$Sipk <-
      as.numeric(input_data[, 3]) * as.numeric(input_data[, 4])
    
    #2.merge the tables together, summarise counts and unique
    input_data$annotation <-
      paste(input_data$annotation, input_data$peptide, sep = "::")
    pwm_pval$annotation <-
      paste(pwm_pval$annotation, pwm_pval$peptide, sep = "::")
    data_merge <-
      unique(merge(input_data, pwm_pval, by = "annotation"))
    
    #3. Swing scores of REAL data.
    swing_out <- swingScore(
      data_merge = data_merge,
      pwm_in = pwm_in,
      permute = FALSE,
      pseudo_count = pseudo_count
    )
    
    #4. permute data and calculate p-values
    if (permutations != FALSE && permutations > 1) {
      if (verbose) {
        message("[Step2/3] : Permuting Network")
      }
      
      # obtain permuted column names
      swing_names <- colnames(data_merge)[7:ncol(data_merge)]
      n_permute <- lapply(seq_len(permutations), function(i)
        sample(as.character(swing_names), length(swing_names),replace = FALSE))
      names(n_permute) <- paste("rand", seq_len(permutations), sep = "")
      
      # calculate swing scores (with permuted names) and return as a vector
      swing_permute <- list(bplapply(seq_len(permutations), function(i)
        swingScore(
          data_merge = setNames(data.frame(data_merge),
                                c(colnames(data_merge)[1:6],
                                  n_permute[[i]])),
          pwm_in = pwm_in,
          permute = TRUE,
          pseudo_count = pseudo_count,
          n = i
        )))
      
      # returns ordered dataframe; order names and merge
      swing_permute_names <- swing_names[order(swing_names)]
      swing_permute <- data.frame("kinase" = swing_permute_names,
                                  swing_permute)
      colnames(swing_permute) <- c("kinase", names(n_permute))
      
      if (verbose) {
        message("[Step3/3] : Calculating p-values")
      }
      #obtain p-values, two sided
      swing_out$p_greater <-
        unlist(bplapply(seq_len(nrow(swing_out)), function(i)
          (sum(
            ifelse(
              as.numeric(swing_permute[swing_permute$kinase ==
                         swing_out$kinase[i], ][2:ncol(swing_permute)]) >
                         as.numeric(swing_out$swing_raw[i]), 1, 0), 
                         na.rm = TRUE) + 1)
                         / (as.numeric(permutations) + 1)))
      
      swing_out$p_less <-
        unlist(bplapply(seq_len(nrow(swing_out)), function(i)
          (sum(
            ifelse(
              as.numeric(swing_permute[swing_permute$kinase ==
                         swing_out$kinase[i], ][2:ncol(swing_permute)]) <
                         as.numeric(swing_out$swing_raw[i]), 1, 0
                         ), na.rm = TRUE) + 1)
                         / (as.numeric(permutations) + 1)))
    }
    # order by p-value
    swing_out <- swing_out[order(swing_out$p_greater),]
    
    if (verbose) {
      message("[FINISHED]\n")
      end_time <- Sys.time() - start_time
      message("Finish: ", Sys.time())
    }
    
    # set option to return a network
    network <- NULL
    if (return_network == "TRUE"){
      network <- data.table::melt(data_merge,
                                 id.vars = c("annotation"),
                                 measure.vars = 
                                 colnames(network[7:ncol(network)]))
      # keep "significant" edges
      network <- network[network$value == 1, ]
      network <- data.frame(
        "source" = as.character(network$variable),
        "target" = as.character(network$annotation)
      )
    }
return(list("scores" = swing_out, "network" = network))
}

# describeIn swing This helper function performs the swing score calculation
# no verbose in this helper
# permute - for setting right table format (using data.table) for merging.

swingScore <-
  function(data_merge,
           pwm_in,
           pseudo_count = 1,
           permute = FALSE,
           n = 1) {
    #----------------------------------------------
    data_merge <- data_merge$Sipk * data_merge[7:ncol(data_merge)]
    # propogation of NaN here due to small number division, add psuedo-count
    # count significant positive
    p_k <- colSums(data_merge == 1, na.rm = TRUE) + pseudo_count
    # count significant negative
    n_k <- colSums(data_merge == -1, na.rm = TRUE) + pseudo_count
    #all counts (both positive and negative).
    t_k <- p_k + n_k
    #----------------------------------------------
    # compute swing statistics
    pk_ratio <- (p_k)/(t_k)
    nk_ratio <- (n_k)/(t_k)
    swing_raw <- (pk_ratio) / (nk_ratio)
    p_n_all <-
      data.frame(
        "kinase" = colnames(data_merge),
        "pos" = p_k,
        "neg" = n_k,
        "all" = t_k,
        "pk" = pk_ratio,
        "nk" = nk_ratio,
        "swing_raw" = swing_raw
      )
    # include the count numbers:
    p_n_all <- merge(p_n_all, pwm_in$kinase, by = "kinase", sort = TRUE)
    # weighted by number of substrates
    # not possible to have zero - must have substrates to build model
    p_n_all$swing_raw <-
      log2(p_n_all$swing_raw) * log2(p_n_all$n)
    # weighted by size of kinase-substrate network:
    p_n_all$swing_raw <-
      p_n_all$swing_raw * log2(p_n_all$all)
    # z-score transform
    p_n_all$swing <-
      (p_n_all$swing_raw - mean(p_n_all$swing_raw, na.rm = TRUE)) / 
      sd(p_n_all$swing_raw, na.rm = TRUE)

    if (permute == "TRUE") {
      p_n_all <- p_n_all$swing_raw
    }
return(p_n_all)
}
