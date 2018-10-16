#' Function for extracting peptide sequences from multimapped or complex 
#' annotated data
#'
#' @description This function extracts unique peptide:annotation combinations 
#' from complex annotated data and formats for further analysis using KinSwingR.
#' For instance, example input annotation may be: 
#' "A0A096MIX2|Ddx17|494|RSRYRTTSSANNPN". This function will extract the peptide
#'  sequence into a second column and associate it all annotations. See vignette
#'   for more details.
#' @param input_data A data.frame of phosphopeptide data. Must contain 4 columns
#'  and the following format must be adhered to. Column 1 - Annotation, Column 2
#'   - centered peptide sequence, Column 3 - Fold Change [-ve to +ve], Column 4 
#'   - p-value [0-1]. This will extract the peptide sequences from Column1 and 
#'   replace all values in Column2 to be used in scoreSequences(). Where 
#'   peptide sequences have not been extracted from the annotation, leave 
#'   Column2 as NA's.
#' @param annotation_delimiter The character used to delimit annotations. 
#' Default="|"
#' @param multi_protein_delimiter The character used to delimit multi-protein 
#' assignments. Default=":". E.g. Ddx17:Ddx2
#' @param multi_site_delimiter The character used to delimit multi-site 
#' assignments. Default=";". E.g. 494;492
#' @param seq_number The annotation frame that contains the sequence after 
#' delimitation. E.g. The sequence "RSRYRTTSSANNPN" is contained in the 4th 
#' annotation frame of the following annotation: 
#' "A0A096MIX2|Ddx17|494|RSRYRTTSSANNPN" and would therefore set seq_number=4.
#' Default=4
#' @param replace Replace a letter that describes sequences outside of the 
#' protein after centering on the phosphosite (e.g X in XXXMERSTRELCLNF). Use in
#'  combination with replace_search and replace_with to replace amino acids. 
#'  Options are "TRUE" or "FALSE". Default="FALSE".
#' @param replace_search Amino Acid to search for when replacing sequences. 
#' Default="X"
#' @param replace_with Amino Acid to replace with when replacing sequences. 
#' Default="_"
#' @param verbose Print progress to screen. Default=FALSE
#'
#' @examples
#' ## Extract peptide sequences from annotation data:
#' data(example_phosphoproteome)
#'
#' ## A0A096MJ61|NA|89|PRRVRNLSAVLAART
#' ## The following will extract all the uniquely annotated peptide
#' ## sequences from the "annotation" column and place these in the
#' ## "peptide" column. Where multi-mapped peptide sequences are input,
#' ## these are placed on a new line.
#' ##
#' ## Here, sequences with a "X" and also replaced with a "_". This is ensure 
#' ## that PWMs are built correctly.
#'
#' ## Sample data for demonstration:
#' sample_data <- head(example_phosphoproteome)
#' annotated_data <- cleanAnnotation(input_data = sample_data,
#'                                    annotation_delimiter = "|",
#'                                    multi_protein_delimiter = ":",
#'                                    multi_site_delimiter = ";",
#'                                    seq_number = 4,
#'                                    replace = TRUE,
#'                                    replace_search = "X",
#'                                    replace_with = "_")
#'
#'## Return the annotated data with extracted peptides:
#' head(annotated_data)
#'
#' @return A data.table with the peptides extracted from the annotation column
#'
#' @export cleanAnnotation
#' @importFrom data.table data.table
#' @importFrom sqldf sqldf

cleanAnnotation <- function(input_data = NULL,
                             annotation_delimiter = "|",
                             multi_protein_delimiter = ":",
                             multi_site_delimiter = ";",
                             seq_number = 4,
                             replace = FALSE,
                             replace_search = "X",
                             replace_with = "_",
                             verbose = FALSE) {
  #----------------------------------------------
  #format checks:
  if (is.null(input_data))
    stop("input_data not provided; you must provide an input table")
  if (!is.data.frame(input_data))
    stop("input_data is not a data.frame; you must provide an
         input table")
  if (replace == TRUE) {
    if (!is.character(replace_search) && !is.character(replace_with)) {
      stop("replace_search AND replace_with MUST both be characters")
    }
  }
  #----------------------------------------------
  colnames(input_data) <- c("annotation", "peptide", "fc", "pval")
  #extract the annotation from the table:
  data_annotation <- as.character(input_data[, 1])
  #1. obtain sequences of interest
  seq_list <-
    strsplit(data_annotation, annotation_delimiter, fixed = TRUE)
  #get the nth element (seq_number) from the list (that contains the sequences)
  seqs <- sapply(seq_list, `[`, seq_number)
  #2. multi-site assignment
  seq_list <-
    strsplit(as.character(seqs), multi_site_delimiter, fixed = TRUE)
  seqs <- unlist(seq_list)
  #3. multi-gene assignment
  seq_list <-
    strsplit(as.character(seqs), multi_protein_delimiter, fixed = TRUE)
  seqs <- unlist(seq_list)
  #4. keep unique sequences for remapping
  seqs <- unique(seqs)
  #5. merge (i.e. remap) with input.data
  dt_in <- data.table(input_data)
  dt_seqs <- data.table(seqs)
  data_out <- data.table(
    sqldf(
      "select *
      from dt_in a inner join dt_seqs b
      on a.annotation like '%' || b.seqs || '%'",
      drv = "SQLite"
    )
  )
  data_out <-
    data.frame(
      "annotation" = data_out$annotation,
      "peptide" = data_out$seqs,
      "fc" = data_out$fc,
      "pval" = data_out$pval
    )
  
  #6. option for reassignment of letters, incase wild-card differs
  if (replace == TRUE) {
    if (verbose) {
      message("Replacing all",
        replace_search,
        "with",
        replace_with
        )
    }
    data_out$peptide <-
      gsub(replace_search, replace_with, data_out$peptide)
  }
  return(unique(data_out))
}

# Helper function for trimming centered sequences
#
# @description Helper function for trimming sequences
# @param seqs_to_trim Vector of sequences to trim
# @param seq_length Full length of substrate sequence (default is 15). Will be 
# trimmed automatically or report error if sequences in kinase_table are not 
# long enough.
# @param verbose Print progress to screen. Default = FALSE
#
# @return Trimmed sequences or HALT

trimSeqs <- function(seqs_to_trim = NULL,
                      seq_length = NULL,
                      verbose = FALSE) {
  #check minimum sequence length is compatible
  min_seq <- min(nchar(seqs_to_trim))
  
  if (min_seq < seq_length) {
    stop("You have selected to trim to", seq_length,
        ", but this is less than minimum sequence length, which is",
        min_seq, ". Trimming will not be performed. Reduce window size to",
        min_seq, "or lower to perform trimming"
        )
  }
  
  #if minimum sequence length is compatible - trim and/or leave same:
  if (min_seq >= seq_length) {
    if (verbose) {
      message("Trimming sequences to", seq_length, "total AA width")
    }
    
    seq_diff <- (min_seq - seq_length) / 2
    trimmed_seqs <-
      substr(seqs_to_trim, 1 + seq_diff, min_seq - seq_diff)
  }
  return(trimmed_seqs)
}
