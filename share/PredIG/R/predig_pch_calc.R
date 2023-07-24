### PREDIG 

### GOAL: COMPUTE ALL VARIABLES TO TRAIN PREDIG

### LIBRARIES
library(Peptides, quietly = T)
library(dplyr, quietly = T)
library(stringr, quietly = T)
library(seqinr, quietly = T)
library(argparser, quietly = T)

### ARGUMENTS
# Create a parser
p <- arg_parser("PREDIG PCH FEATURE CALCULATION")

# Add command line arguments
#p <- add_argument(p, "--help", help= 'This script calculates the physicochemical features requires by PredIG. Input requires a CSV with a column named "peptide". Ouput is a CSV at the same path of input ended with "_pch.csv"', type="character")
p <- add_argument(p, "--input", help= "Path to CSV with a required column named 'peptide'", type="character", default = ".")
p <- add_argument(p, "--seed", help= "Seed number", type = "numeric",  default= 123)

# Parse the command line arguments
argv <- parse_args(p)

### UPLOAD RAW DATA
pch_input <- read.csv(argv$input)

### COMPUTE PHYSICOCHEMICAL PROPERTIES

## CONVERT TO VECTOR
pch_input_vector <- as.vector(pch_input$peptide)
length(pch_input_vector)

## EXTRACT TCR CONTACT REGION (RESIDUES 4-W-1)
tcr_contact <- str_sub(pch_input_vector,4,-2)
tcr_contact <- as.vector(tcr_contact)
glimpse(tcr_contact)

## APPEND TCR CONTACT REGION TO DATAFRAME
pch_input <- cbind(pch_input,tcr_contact)

## Bulkiness or Molecular Weight
# Calculate the molecular weight of the Mut peptide # Error due to presence of PTMs and X as aminoacid. Remove both.
mw_peptide <- mw(pch_input_vector)

# Calculate the molecular weight of the TCR contact region (surrogate of bulkiness)
mw_tcr_contact <- mw(tcr_contact)

## HYDROPHOBICITY
# Calculate the hydrophbicity of the Mut peptide
hydroph_peptide <- hydrophobicity(pch_input_vector)

# Calculate the hydrophbicity of the  Mut contact region (surrogate of bulkyness)
hydroph_tcr_contact <- hydrophobicity(tcr_contact)

## NETCHARGE
# Calculate the net charge of the Mut peptide
charge_peptide <- charge(pch_input_vector)

# Calculate the net charge of the  Mut contact region (surrogate of bulkyness)
charge_tcr_contact <- charge(tcr_contact)

# PEPTIDE STABILITY
# Calculate the stability index of the Mut peptide
stab_peptide <- instaIndex(pch_input_vector)

## MAKES NO SENSE TO CALCULATE STABILITY FOR TCR CONTACT REGION SINCE IT IS CONTAINED WITHIN THE PEPTIDE.

# EXPORT DATAFRAME
# Incorporate properties.
pch_input <- cbind(pch_input, mw_peptide, hydroph_peptide, charge_peptide, stab_peptide, mw_tcr_contact, hydroph_tcr_contact, charge_tcr_contact)
write.csv(pch_input, file = paste(gsub(".csv", "", argv$input),"_pch.csv", sep = ""), quote = F, row.names = F)
