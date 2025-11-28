#### Install and load the required packages ####

install.packages("data.table")
library(data.table)
install.packages("plyr")
library(plyr)
install.packages("dplyr")
library(dplyr)
install.packages("stringr")
library(stringr)
install.packages("tidyr")
library(tidyr)
install.packages("writexl")
library(writexl)

# Install Biostrings from Bioconductor if not already available
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Biostrings")
}

library(Biostrings)

### NB. theres tw versions of this script, the first one if for those who want to control how the reconstructions unfold ###
### and the theres a quick version of the script at the end####

################################################################################
################################################################################

################# Version one, the long way around #############################

################################################################################
################################################################################
### Set your working directory (this is where the evidence file should be located) ###

setwd("/path/to/directory/")

# Read the evidence file into a data frame
# Ensure the file contains headers such as "MS.MS.scan.number", "Experiment", "Score", "Sequence"
# These columns are useful for identifying peptide sequences and validating spectra
# Note: fileEncoding = "UTF-16LE" may be required depending on your file format

data <- read.table("evidence_ed.txt", header = TRUE, sep = "\t")

#### Extract peptides for specific genes and save them as text files ####

# Filter the rows that match the gene name "ALB"
alb_data <- data[grep("ALB", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
alb_data1 <- alb_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

alb_data2 <- unite(alb_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/alb_peptides.txt"  
# Replace with your desired file path
write.table(alb_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)

# Filter the rows that match the gene name "AMELX"
amelx_data <- data[grep("AMELX", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
amelx_data1 <- amelx_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

amelx_data2 <- unite(amelx_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/amelx_peptides.txt"  
# Replace with your desired file path
write.table(amelx_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)

# Filter the rows that match the gene name "AMELY"
amely_data <- data[grep("AMELY", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
amely_data1 <- amely_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

amely_data2 <- unite(amely_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/amely_peptides.txt"  
# Replace with your desired file path
write.table(amely_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)


# Filter the rows that match the gene name "AMBN"
ambn_data <- data[grep("AMBN", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
ambn_data1 <- ambn_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

ambn_data2 <- unite(ambn_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/ambn_peptides.txt"  
# Replace with your desired file path
write.table(ambn_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)


# Filter the rows that match the gene name "AHSG"
ahsg_data <- data[grep("AHSG", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
ahsg_data1 <- ahsg_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

ahsg_data2 <- unite(ahsg_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/ahsg_peptides.txt"  
# Replace with your desired file path
write.table(ahsg_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)


# Filter the rows that match the gene name "AMTN"
amtn_data <- data[grep("AMTN", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
amtn_data1 <- amtn_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

amtn_data2 <- unite(amtn_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/amtn_peptides.txt"  
# Replace with your desired file path
write.table(amtn_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)


# Filter the rows that match the gene name "COL17A1"
col17a1_data <- data[grep("COL17A1", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
col17a1_data1 <- col17a1_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

col17a1_data2 <- unite(col17a1_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/col17a1_peptides.txt"  
# Replace with your desired file path
write.table(col17a1_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)

# Filter the rows that match the gene name "COL1A1"
col1a1_data <- data[grep("COL1A1", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
col1a1_data1 <- col1a1_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

col1a1_data2 <- unite(col1a1_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/col1a1_peptides.txt"  
# Replace with your desired file path
write.table(col1a1_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)

# Filter the rows that match the gene name "COL1A2"
col1a2_data <- data[grep("COL1A2", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
col1a2_data1 <- col1a2_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

col1a2_data2 <- unite(col1a2_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/col1a2_peptides.txt"  
# Replace with your desired file path
write.table(col1a2_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)

# Filter the rows that match the gene name "ENAM"
enam_data <- data[grep("ENAM", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
enam_data1 <- enam_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

enam_data2 <- unite(enam_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/enam_peptides.txt"  
# Replace with your desired file path
write.table(enam_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)


# Filter the rows that match the gene name "ODAM"
odam_data <- data[grep("ODAM", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
odam_data1 <- odam_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

odam_data2 <- unite(odam_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/odam_peptides.txt"  
# Replace with your desired file path
write.table(odam_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)


# Filter the rows that match the gene name "MMP20"
mmp20_data <- data[grep("MMP20", data$Proteins), ]

# Select specific columns from filtered_data and create a new data frame
mmp20_data1 <- mmp20_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]

mmp20_data2 <- unite(mmp20_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")

file_path <- "/path/to/directory/1_PSMforFasta/mmp20_peptides.txt"  
# Replace with your desired file path
write.table(mmp20_data2, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)


#### Convert peptide sequences into FASTA format ####

#create a AA stringsetobjcet for MMP20
mmp20sequences <- AAStringSet(mmp20_data2$Sequence)

#create a vector for my sequence names headers

mmp20headers <- mmp20_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(mmp20sequences, mmp20headers, output_file) {
  if (length(mmp20sequences) != length(mmp20headers)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(mmp20sequences))
  for (i in seq_along(mmp20sequences)) {
    fasta_lines[i] <- paste0(">", mmp20headers[i], "\n", as.character(mmp20sequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/mmp20_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(mmp20sequences, mmp20headers, output_fasta_file)

#create a AA stringsetobjcet for ahsg
ahsgsequences <- AAStringSet(ahsg_data2$Sequence)

#create a vector for my sequence names headers

ahsgheaders <- ahsg_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(ahsgsequences, ahsgheaders, output_file) {
  if (length(ahsgsequences) != length(ahsgheaders)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(ahsgsequences))
  for (i in seq_along(ahsgsequences)) {
    fasta_lines[i] <- paste0(">", ahsgheaders[i], "\n", as.character(ahsgsequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/ahsg_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(ahsgsequences, ahsgheaders, output_fasta_file)

#create a AA stringsetobjcet for alb
albsequences <- AAStringSet(alb_data2$Sequence)

#create a vector for my sequence names headers

albheaders <- alb_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(albsequences, albheaders, output_file) {
  if (length(albsequences) != length(albheaders)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(albsequences))
  for (i in seq_along(albsequences)) {
    fasta_lines[i] <- paste0(">", albheaders[i], "\n", as.character(albsequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/alb_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(albsequences, albheaders, output_fasta_file)

#create a AA stringsetobjcet for ambn
ambnsequences <- AAStringSet(ambn_data2$Sequence)

#create a vector for my sequence names headers

ambnheaders <- ambn_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(ambnsequences, ambnheaders, output_file) {
  if (length(ambnsequences) != length(ambnheaders)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(ambnsequences))
  for (i in seq_along(ambnsequences)) {
    fasta_lines[i] <- paste0(">", ambnheaders[i], "\n", as.character(ambnsequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/ambn_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(ambnsequences, ambnheaders, output_fasta_file)

#create a AA stringsetobjcet for amelx
amelxsequences <- AAStringSet(amelx_data2$Sequence)

#create a vector for my sequence names headers

amelxheaders <- amelx_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(amelxsequences, amelxheaders, output_file) {
  if (length(amelxsequences) != length(amelxheaders)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(amelxsequences))
  for (i in seq_along(amelxsequences)) {
    fasta_lines[i] <- paste0(">", amelxheaders[i], "\n", as.character(amelxsequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/amelx_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(amelxsequences, amelxheaders, output_fasta_file)

#create a AA stringsetobjcet for amely
amelysequences <- AAStringSet(amely_data2$Sequence)

#create a vector for my sequence names headers

amelyheaders <- amely_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(amelysequences, amelyheaders, output_file) {
  if (length(amelysequences) != length(amelyheaders)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(amelysequences))
  for (i in seq_along(amelysequences)) {
    fasta_lines[i] <- paste0(">", amelyheaders[i], "\n", as.character(amelysequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/amely_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(amelysequences, amelyheaders, output_fasta_file)

#create a AA stringsetobjcet for amtn
amtnsequences <- AAStringSet(amtn_data2$Sequence)

#create a vector for my sequence names headers

amtnheaders <- amtn_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(amtnsequences, amtnheaders, output_file) {
  if (length(amtnsequences) != length(amtnheaders)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(amtnsequences))
  for (i in seq_along(amtnsequences)) {
    fasta_lines[i] <- paste0(">", amtnheaders[i], "\n", as.character(amtnsequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/amtn_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(amtnsequences, amtnheaders, output_fasta_file)

#create a AA stringsetobjcet for col1a1
col1a1sequences <- AAStringSet(col1a1_data2$Sequence)

#create a vector for my sequence names headers

col1a1headers <- col1a1_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(col1a1sequences, col1a1headers, output_file) {
  if (length(col1a1sequences) != length(col1a1headers)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(col1a1sequences))
  for (i in seq_along(col1a1sequences)) {
    fasta_lines[i] <- paste0(">", col1a1headers[i], "\n", as.character(col1a1sequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/col1a1_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(col1a1sequences, col1a1headers, output_fasta_file)

#create a AA stringsetobjcet for col17a1
col17a1sequences <- AAStringSet(col17a1_data2$Sequence)

#create a vector for my sequence names headers

col17a1headers <- col17a1_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(col17a1sequences, col17a1headers, output_file) {
  if (length(col17a1sequences) != length(col17a1headers)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(col17a1sequences))
  for (i in seq_along(col17a1sequences)) {
    fasta_lines[i] <- paste0(">", col17a1headers[i], "\n", as.character(col17a1sequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/col17a1_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(col17a1sequences, col17a1headers, output_fasta_file)

#create a AA stringsetobjcet for col1a2
col1a2sequences <- AAStringSet(col1a2_data2$Sequence)

#create a vector for my sequence names headers

col1a2headers <- col1a2_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(col1a2sequences, col1a2headers, output_file) {
  if (length(col1a2sequences) != length(col1a2headers)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(col1a2sequences))
  for (i in seq_along(col1a2sequences)) {
    fasta_lines[i] <- paste0(">", col1a2headers[i], "\n", as.character(col1a2sequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/col1a2_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(col1a2sequences, col1a2headers, output_fasta_file)

#create a AA stringsetobjcet for enam
enamsequences <- AAStringSet(enam_data2$Sequence)

#create a vector for my sequence names headers

enamheaders <- enam_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(enamsequences, enamheaders, output_file) {
  if (length(enamsequences) != length(enamheaders)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(enamsequences))
  for (i in seq_along(enamsequences)) {
    fasta_lines[i] <- paste0(">", enamheaders[i], "\n", as.character(enamsequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/enam_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(enamsequences, enamheaders, output_fasta_file)

#create a AA stringsetobjcet for odam
odamsequences <- AAStringSet(odam_data2$Sequence)

#create a vector for my sequence names headers

odamheaders <- odam_data2$id

#Create a function to write the amino acid sequences and their headers to a FASTA file

writeFasta <- function(odamsequences, odamheaders, output_file) {
  if (length(odamsequences) != length(odamheaders)) {
    stop("Number of sequences and headers must be the same.")
  }
  
  fasta_lines <- character(length(odamsequences))
  for (i in seq_along(odamsequences)) {
    fasta_lines[i] <- paste0(">", odamheaders[i], "\n", as.character(odamsequences[i]))
  }
  
  writeLines(fasta_lines, con = output_file)
}

# Specify the output FASTA file
output_fasta_file <- "/path/to/directory/2_PSMtoFasta/odam_peptides.fasta"

# Write the sequences to the FASTA file
writeFasta(odamsequences, odamheaders, output_fasta_file)

#### Align peptides with reference sequences using MAFFT ####

# MAFFT must be installed on your machine: https://mafft.cbrc.jp/alignment/software/
# The following commands align peptide FASTA files with modern reference alignments
#alb
system("einsi --leavegappyregion --addlong /path/to/directory/2_PSMtoFasta/alb_peptides.fasta /path/to/directory/Ref_ALB_aln.fa > /path/to/directory/3_alignedPSMs/alb_peptides_aln.fasta")

#ahsg
system("einsi --leavegappyregion --addlong  /path/to/directory/2_PSMtoFasta/ahsg_peptides.fasta /path/to/directory/Ref_AHSG_aln.fa > /path/to/directory/3_alignedPSMs/ahsg_peptides_aln.fasta")

#ambn
system("einsi --leavegappyregion --addlong /path/to/directory/2_PSMtoFasta/ambn_peptides.fasta /path/to/directory/Ref_AMBN_aln.fa > /path/to/directory/3_alignedPSMs/ambn_peptides_aln.fasta")

#amelx
system("einsi --leavegappyregion --addlong /path/to/directory/2_PSMtoFasta/amelx_peptides.fasta /path/to/directory/Ref_AMELX_aln.fa > /path/to/directory/3_alignedPSMs/amelx_peptides_aln.fasta")

#amely
system("einsi --leavegappyregion --addlong /path/to/directory/2_PSMtoFasta/amely_peptides.fasta /path/to/directory/Ref_AMELY_aln.fa > /path/to/directory/3_alignedPSMs/amely_peptides_aln.fasta")

#amtn
system("einsi --leavegappyregion --addlong /path/to/directory/2_PSMtoFasta/amtn_peptides.fasta /path/to/directory/Ref_AMTN_aln.fa > /path/to/directory/3_alignedPSMs/amtn_peptides_aln.fasta")

#col1a1
system("einsi --leavegappyregion --addlong /path/to/directory/2_PSMtoFasta/col1a1_peptides.fasta /path/to/directory/Ref_COL1A1_aln.fa > /path/to/directory/3_alignedPSMs/col1a1_peptides_aln.fasta")

#col1a2
system("einsi --leavegappyregion --addlong /path/to/directory/2_PSMtoFasta/col1a2_peptides.fasta /path/to/directory/Ref_COL1A2_aln.fa > /path/to/directory/3_alignedPSMs/col1a2_peptides_aln.fasta")

#col17a1
system("einsi --leavegappyregion --addlong /path/to/directory/2_PSMtoFasta/col17a1_peptides.fasta /path/to/directory/Ref_COL17A1_aln.fa > /path/to/directory/3_alignedPSMs/col17a1_peptides_aln.fasta")


#enam
system("einsi --leavegappyregion --addlong /path/to/directory/2_PSMtoFasta/enam_peptides.fasta /path/to/directory/Ref_ENAM_aln.fa > /path/to/directory/3_alignedPSMs/enam_peptides_aln.fasta")

#mmp20
system("einsi --leavegappyregion --addlong /path/to/directory/2_PSMtoFasta/mmp20_peptides.fasta /path/to/directory/Ref_MMP20_aln.fa > /path/to/directory/3_alignedPSMs/mmp20_peptides_aln.fasta")

#col
system("einsi --leavegappyregion --addlong /path/to/directory/2_PSMtoFasta/odam_peptides.fasta /path/to/directory/Ref_ODAM_aln.fa > /path/to/directory/3_alignedPSMs/odam_peptides_aln.fasta")







################################################################################
################################################################################

################# Version two, the short way around#############################

################################################################################
################################################################################



#### Install and load the required packages ####
install.packages("data.table")
library(data.table)

install.packages("plyr")
library(plyr)

install.packages("dplyr")
library(dplyr)

install.packages("stringr")
library(stringr)

install.packages("tidyr")
library(tidyr)

install.packages("writexl")
library(writexl)

# Install Biostrings if not already available
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
library(Biostrings)

### Set working directory ###
setwd("/path/to/directory/")

# Read the evidence file
data <- read.table("evidence_ed.txt", header = TRUE, sep = "\t")

# Define the list of genes you want to process
genes <- c("ALB","AMELX","AMELY","AMBN","AHSG","AMTN",
           "COL17A1","COL1A1","COL1A2","ENAM","ODAM","MMP20")

# Function to filter, export peptides, and create FASTA
process_gene <- function(gene, data) {
  # Filter rows for the gene
  gene_data <- data[grep(gene, data$Proteins), ]
  
  # Select relevant columns
  gene_data1 <- gene_data[, c("MS.MS.scan.number", "Experiment", "Score", "Sequence")]
  
  # Create unique ID
  gene_data2 <- unite(gene_data1, "MS.MS.scan.number", "Experiment", "Score", col = "id", sep = "_")
  
  # Write peptide table
  txt_file <- paste0("/path/to/directory/1_PSMforFasta/", tolower(gene), "_peptides.txt")
  write.table(gene_data2, file = txt_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Create AAStringSet and headers
  sequences <- AAStringSet(gene_data2$Sequence)
  headers <- gene_data2$id
  
  # Write FASTA
  fasta_file <- paste0("/path/to/directory/2_PSMtoFasta/", tolower(gene), "_peptides.fasta")
  writeFasta <- function(sequences, headers, output_file) {
    if (length(sequences) != length(headers)) {
      stop("Number of sequences and headers must be the same.")
    }
    fasta_lines <- character(length(sequences))
    for (i in seq_along(sequences)) {
      fasta_lines[i] <- paste0(">", headers[i], "\n", as.character(sequences[i]))
    }
    writeLines(fasta_lines, con = output_file)
  }
  writeFasta(sequences, headers, fasta_file)
  
  # Run MAFFT alignment (adjust paths to your reference alignments)
  ref_file <- paste0("/path/to/directory/Ref_", gene, "_aln.fa")
  aln_file <- paste0("/path/to/directory/3_alignedPSMs/", tolower(gene), "_peptides_aln.fasta")
  cmd <- paste("einsi --leavegappyregion --addlong", fasta_file, ref_file, ">", aln_file)
  system(cmd)
}

# Iterate through all genes
for (g in genes) {
  process_gene(g, data)
}

