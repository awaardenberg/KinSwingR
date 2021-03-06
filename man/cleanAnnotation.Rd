% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_seqs.r
\name{cleanAnnotation}
\alias{cleanAnnotation}
\title{Function for extracting peptide sequences from multimapped or complex 
annotated data}
\usage{
cleanAnnotation(input_data = NULL, annotation_delimiter = "|",
  multi_protein_delimiter = ":", multi_site_delimiter = ";",
  seq_number = 4, replace = FALSE, replace_search = "X",
  replace_with = "_", verbose = FALSE)
}
\arguments{
\item{input_data}{A data.frame of phosphopeptide data. Must contain 4 columns
and the following format must be adhered to. Column 1 - Annotation, Column 2
 - centered peptide sequence, Column 3 - Fold Change [-ve to +ve], Column 4 
 - p-value [0-1]. This will extract the peptide sequences from Column1 and 
 replace all values in Column2 to be used in scoreSequences(). Where 
 peptide sequences have not been extracted from the annotation, leave 
 Column2 as NA's.}

\item{annotation_delimiter}{The character used to delimit annotations. 
Default="|"}

\item{multi_protein_delimiter}{The character used to delimit multi-protein 
assignments. Default=":". E.g. Ddx17:Ddx2}

\item{multi_site_delimiter}{The character used to delimit multi-site 
assignments. Default=";". E.g. 494;492}

\item{seq_number}{The annotation frame that contains the sequence after 
delimitation. E.g. The sequence "RSRYRTTSSANNPN" is contained in the 4th 
annotation frame of the following annotation: 
"A0A096MIX2|Ddx17|494|RSRYRTTSSANNPN" and would therefore set seq_number=4.
Default=4}

\item{replace}{Replace a letter that describes sequences outside of the 
protein after centering on the phosphosite (e.g X in XXXMERSTRELCLNF). Use in
 combination with replace_search and replace_with to replace amino acids. 
 Options are "TRUE" or "FALSE". Default="FALSE".}

\item{replace_search}{Amino Acid to search for when replacing sequences. 
Default="X"}

\item{replace_with}{Amino Acid to replace with when replacing sequences. 
Default="_"}

\item{verbose}{Print progress to screen. Default=FALSE}
}
\value{
A data.table with the peptides extracted from the annotation column
}
\description{
This function extracts unique peptide:annotation combinations 
from complex annotated data and formats for further analysis using KinSwingR.
For instance, example input annotation may be: 
"A0A096MIX2|Ddx17|494|RSRYRTTSSANNPN". This function will extract the peptide
 sequence into a second column and associate it all annotations. See vignette
  for more details.
}
\examples{
## Extract peptide sequences from annotation data:
data(example_phosphoproteome)

## A0A096MJ61|NA|89|PRRVRNLSAVLAART
## The following will extract all the uniquely annotated peptide
## sequences from the "annotation" column and place these in the
## "peptide" column. Where multi-mapped peptide sequences are input,
## these are placed on a new line.
##
## Here, sequences with a "X" and also replaced with a "_". This is ensure 
## that PWMs are built correctly.

## Sample data for demonstration:
sample_data <- head(example_phosphoproteome)
annotated_data <- cleanAnnotation(input_data = sample_data,
                                   annotation_delimiter = "|",
                                   multi_protein_delimiter = ":",
                                   multi_site_delimiter = ";",
                                   seq_number = 4,
                                   replace = TRUE,
                                   replace_search = "X",
                                   replace_with = "_")

## Return the annotated data with extracted peptides:
head(annotated_data)

}
