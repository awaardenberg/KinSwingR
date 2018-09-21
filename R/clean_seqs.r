#' Function for extracting peptide sequences from multimapped or complex annotated
#'  data.
#'
#' @description This function extracts unique peptide:annotation combinations from complex
#' annotated data and formats for further analysis using KinSwingR. For instance, example input
#' annotation may be: "A0A096MIX2|Ddx17|494|RSRYRTTSSANNPN". This function will extract the
#' peptide sequence into a second column and associate it all annotations. See vignette for more
#' details.
#' @param input.data A data.frame of phosphopeptide data. Must contain 4 columns
#'  and the following format must be adhered to. Column 1 - Annotation, Column 2 -
#'  centered peptide sequence, Column 3 - Fold Change [-ve to +ve], Column 4 - p-value [0-1].
#'  This will extract the peptide sequences from Column1 and replace all values in Column2 to
#'   be used in score.sequences(). Where peptide sequences have not been extracted from the
#'   annotation, leave Column2 as NA's.
#' @param annotation.delimiter The character used to delimit annotations. Default="|"
#' @param multi.protein.delimiter The character used to delimit multi-protein assignments.
#' Default=":". E.g. Ddx17:Ddx2
#' @param multi.site.delimiter The character used to delimit multi-site assignments.
#' Default=";". E.g. 494;492
#' @param seq.number The annotation frame that contains the sequence after delimitation. E.g.
#' The sequence "RSRYRTTSSANNPN" is contained in the 4th annotation frame of the following
#' annotation: "A0A096MIX2|Ddx17|494|RSRYRTTSSANNPN" and would therefore set seq.number=4.
#' Default=4
#' @param replace Replace a letter that describes sequences outside of the protein after
#' centering on the phosphosite (e.g X in XXXMERSTRELCLNF). Use in combination with replace.search
#' and replace.with to replace amino acids. Options are "TRUE" or "FALSE". Default="FALSE".
#' @param replace.search Amino Acid to search for when replacing sequences. Default="X"
#' @param replace.with Amino Acid to replace with when replacing sequences. Default="_"
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
#' ## Here, sequences with a "X" and also replaced with a "_". This is ensure that
#' ## PWMs are built correctly.
#' 
#' ## Sample data for demonstration:
#' sample.data <- head(example_phosphoproteome)
#' annotated.data <- clean.annotation(input.data = sample.data, 
#'                                    annotation.delimiter = "|",
#'                                    multi.protein.delimiter = ":", 
#'                                    multi.site.delimiter = ";",
#'                                    seq.number = 4, 
#'                                    replace = TRUE, 
#'                                    replace.search = "X",
#'                                    replace.with = "_")
#'
#'## Return the annotated data with extracted peptides:
#'head(annotated.data) 
#' 
#' @return A data.table with the peptides extracted from the annotation column
#'
#' @export clean.annotation

clean.annotation <- function(input.data = NULL, 
                             annotation.delimiter = "|", 
                             multi.protein.delimiter = ":",
                             multi.site.delimiter = ";", 
                             seq.number = 4, 
                             replace = FALSE,
                             replace.search = "X", 
                             replace.with = "_", 
                             verbose = FALSE
                             ){
  #----------------------------------------------
  #format checks:
  if (is.null(input.data)) stop("input.data not provided; you must provide an input table")
  if (!is.data.frame(input.data)) stop("input.data is not a data.frame; you must provide an
                                       input table")
  if(replace == TRUE){
    if (!is.character(replace.search) && !is.character(replace.with)){
      stop("replace.search AND replace.with MUST both be characters")
    }
  }

  #----------------------------------------------
  colnames(input.data) <- c("annotation", "peptide", "fc", "pval")
  #extract the annotation from the table:
  data.annotation <- as.character(input.data[,1])
  #1. obtain sequences of interest
  seq.list <- strsplit(data.annotation, annotation.delimiter, fixed=T)
  #get the nth element (seq.number) from the list (which contains the sequences):
  seqs <- sapply(seq.list, `[`, seq.number)
  #2. multi-site assignment
  seq.list <- strsplit(as.character(seqs), multi.site.delimiter, fixed=T)
  seqs <- unlist(seq.list)
  #3. multi-gene assignment
  seq.list <- strsplit(as.character(seqs), multi.protein.delimiter, fixed=T)
  seqs <- unlist(seq.list)
  #4. keep unique sequences for remapping
  seqs <- unique(seqs)
  #5. merge (i.e. remap) with input.data
  dt_in <- data.table(input.data)
  dt_seqs <- data.table(seqs)
  data.out <- data.table(sqldf("select *
                            from dt_in a inner join dt_seqs b
                            on a.annotation like '%' || b.seqs || '%'", drv = "SQLite" ))
  data.out <- data.frame("annotation"=data.out$annotation, "peptide"=data.out$seqs,
                         "fc"=data.out$fc, "pval"=data.out$pval)

  #6. option for reassignment of letters, incase wild-card differs
  if(replace==TRUE){
    if(verbose){
      print(paste("Replacing all", replace.search, "with", replace.with, sep=" "))
  }
  data.out$peptide <- gsub(replace.search, replace.with, data.out$peptide)
  }
return(unique(data.out))
}

# Helper function for trimming centered sequences.
#
# @description Helper function for trimming sequences
# @param seqs.to.trim Vector of sequences to trime.
# @param seq.length Full length of substrate sequence (default is 15). Will be trimmed
#  automatically or report error if sequences in kinase.table are not long enough.
# @param verbose Print progress to screen. Default = FALSE
#
# @return Trimmed sequences or HALT

trim.seqs <- function(seqs.to.trim = NULL, 
                      seq.length = NULL, 
                      verbose=FALSE){
  #check minimum sequence length is compatible
  min.seq <- min(nchar(seqs.to.trim))
  
  if (min.seq < seq.length) {
    stop(cat(paste(
      "You have selected to trim to", seq.length, ", but this is less than
      minimum sequence length, which is", min.seq, ". Trimming will not be
      performed. Reduce window size to", min.seq, "or lower to perform
      trimming", sep=" ")))
  }
  
  #if minimum sequence length is compatible - trim and/or leave same:
  if(min.seq >= seq.length) {
    if(verbose){
      cat(paste("Trimming sequences to", seq.length, "total AA width", sep=" "))
    }
    
    seq.diff <- (min.seq - seq.length) / 2
    trimmed.seqs <-
      substr(seqs.to.trim, 1 + seq.diff, min.seq - seq.diff)
  }
  return(trimmed.seqs)
}

