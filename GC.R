# Set the working directory (modify as needed)
setwd("C:/Users/mollear/Desktop")

# Load required libraries
library(tidyverse)
library(readr)
library(dplyr)
library(stringr)

# Function to classify molecules based on their names
# Returns a functional group or class for a given molecule name
classify_molecule <- function(molecule_name) {
  # Convert to lowercase and trim whitespace for matching
  name_lower <- tolower(trimws(molecule_name))
  original_name <- trimws(molecule_name)

  # List of amino acids for classification
  amino_acids <- c("alanine", "arginine", "asparagine", "aspartic acid", "aspartate",
                   "cysteine", "glutamic acid", "glutamate", "glutamine", "glycine",
                   "histidine", "isoleucine", "leucine", "lysine", "methionine",
                   "phenylalanine", "proline", "serine", "threonine", "tryptophan",
                   "tyrosine", "valine", "hydroxyproline", "selenocysteine", "pyrrolysine",
                   "ornithine", "citrulline", "homocysteine", "taurine", "beta-alanine")

  # Various chemical pattern lists for classification
  steroid_patterns <- c("sterol", "steroid", "cholesterol", "cortisol", "testosterone", 
                        "estradiol", "progesterone", "androgen", "estrogen", "cortisone",
                        "pregnenolone", "dehydroepiandrosterone", "aldosterone")
  sugar_patterns <- c("glucose", "fructose", "galactose", "mannose", "ribose", "xylose",
                      "arabinose", "sucrose", "lactose", "maltose", "cellobiose", "trehalose",
                      "glycerol", "sorbitol", "mannitol", "xylitol", "erythritol")
  vitamin_patterns <- c("vitamin", "retinol", "thiamine", "riboflavin", "niacin", "pyridoxine",
                        "biotin", "folate", "cobalamin", "ascorbic acid", "tocopherol", "phylloquinone",
                        "cholecalciferol", "ergocalciferol", "pantothenic acid")
  nucleotide_patterns <- c("adenosine", "guanosine", "cytidine", "uridine", "thymidine",
                           "atp", "adp", "amp", "gtp", "gdp", "gmp", "ctp", "cdp", "cmp",
                           "utp", "udp", "ump", "ttp", "tdp", "tmp")
  fatty_acid_patterns <- c("palmitic", "stearic", "oleic", "linoleic", "linolenic", "arachidonic",
                           "myristic", "lauric", "capric", "caprylic", "caproic", "butyric", "acetic")

  # Classification checks: order matters (more specific checks first)
  # Each check returns as soon as a match is found

  # Amino acids
  if (any(sapply(amino_acids, function(aa) str_detect(name_lower, paste0("\\b", aa, "\\b"))))) {
    return("Amino Acid")
  }

  # Heterocyclic compounds (prefixes)
  if (str_detect(original_name, "^Pyrazine,|^Pyridine,|^Pyrazole,|^Furan")) return("Heterocyclic Compound")
  if (str_detect(name_lower, "^pyrazine|^pyridine|^pyrazole|^furan|^thiophene|^imidazole|^oxazole|^thiazole|^indole|^quinoline|^isoquinoline|^purine|^pyrimidine")) return("Heterocyclic Compound")

  # Polyaromatic compounds
  if (str_detect(original_name, "^Naphthalene,|^Anthracene,")) return("Polyaromatic Compound")
  if (str_detect(name_lower, "naphthalene|anthracene|phenanthrene|pyrene|chrysene|benzo\\[a\\]pyrene|fluoranthene|coronene")) return("Polyaromatic Compound")

  # Steroids
  if (any(sapply(steroid_patterns, function(pat) str_detect(name_lower, pat)))) {
    return("Steroid")
  }

  # Sugars and carbohydrates
  if (any(sapply(sugar_patterns, function(pat) str_detect(name_lower, pat)))) {
    return("Carbohydrate/Sugar")
  }
  if (str_detect(name_lower, "saccharide|hexose|pentose|triose|tetrose|heptose")) return("Carbohydrate/Sugar")

  # Vitamins
  if (any(sapply(vitamin_patterns, function(pat) str_detect(name_lower, pat)))) {
    return("Vitamin")
  }

  # Nucleotides and nucleosides
  if (any(sapply(nucleotide_patterns, function(pat) str_detect(name_lower, pat)))) {
    return("Nucleotide/Nucleoside")
  }

  # Fatty acids and lipids
  if (any(sapply(fatty_acid_patterns, function(pat) str_detect(name_lower, paste0(pat, ".*acid"))))) {
    return("Fatty Acid")
  }
  if (str_detect(name_lower, "phospholipid|triglyceride|diglyceride|monoglyceride|sphingolipid|ceramide")) return("Lipid")

  # Alkaloids
  if (str_detect(name_lower, "caffeine|nicotine|morphine|codeine|quinine|strychnine|atropine|cocaine|ephedrine|mescaline")) return("Alkaloid")

  # Terpenes
  if (str_detect(name_lower, "limonene|pinene|camphor|menthol|geraniol|linalool|myrcene|caryophyllene")) return("Terpene")
  if (str_detect(name_lower, "terpene|terpenoid|monoterpene|sesquiterpene|diterpene|triterpene")) return("Terpene")

  # Peptides and proteins
  if (str_detect(name_lower, "peptide|protein|insulin|hemoglobin|albumin|globulin|collagen|elastin")) return("Peptide/Protein")

  # Enzyme names
  if (str_detect(name_lower, "ase$|kinase|phosphatase|dehydrogenase|oxidase|reductase|transferase|hydrolase|lyase|isomerase|ligase")) return("Enzyme")

  # Non-conventional naming (esters)
  if (str_detect(name_lower, "ethyl ester|methyl ester|propyl ester|butyl ester|ester$")) return("Ester")
  if (str_detect(name_lower, "\\w+ ester")) return("Ester")

  # Amides
  if (str_detect(name_lower, "\\w+ amide|amide$")) return("Amide")

  # Carboxylic acids and derivatives
  if (str_detect(name_lower, "oic acid")) return("Organic Acid")
  if (str_detect(name_lower, "ic acid$")) return("Carboxylic Acid")
  if (str_detect(name_lower, "oate$|ate$")) return("Ester")
  if (str_detect(name_lower, "amide$")) return("Amide")
  if (str_detect(name_lower, "nitrile$|cyanide")) return("Nitrile")

  # Aldehydes and ketones
  if (str_detect(name_lower, "aldehyde|al$")) return("Aldehyde")
  if (str_detect(name_lower, "ketone|one$")) return("Ketone")

  # Alcohols and phenols
  if (str_detect(name_lower, "phenol")) return("Phenol")
  if (str_detect(name_lower, "triol$")) return("Triol")
  if (str_detect(name_lower, "diol$|glycol")) return("Diol")
  if (str_detect(name_lower, "ol$")) return("Alcohol")

  # Ethers
  if (str_detect(name_lower, "ether$")) return("Ether")
  if (str_detect(name_lower, "methoxy|ethoxy|propoxy|butoxy")) return("Ether")

  # Amines
  if (str_detect(name_lower, "amine$")) return("Amine")
  if (str_detect(name_lower, "amino")) return("Amine")

  # Alkenes and alkynes
  if (str_detect(name_lower, "triene$")) return("Triene")
  if (str_detect(name_lower, "diene$")) return("Diene")
  if (str_detect(name_lower, "yne$")) return("Alkyne")
  if (str_detect(name_lower, "ene$")) return("Alkene")

  # Alkanes
  if (str_detect(name_lower, "ane$")) return("Alkane")

  # Aromatic compounds
  if (str_detect(name_lower, "benzene|phenyl|toluene|xylene|styrene")) return("Aromatic Hydrocarbon")

  # Halogenated compounds
  if (str_detect(name_lower, "fluoro|chloro|bromo|iodo")) return("Halogenated Compound")
  if (str_detect(name_lower, "fluoride|chloride|bromide|iodide")) return("Halide")

  # Sulfur compounds
  if (str_detect(name_lower, "thiol$|mercaptan")) return("Thiol")
  if (str_detect(name_lower, "sulfide$")) return("Sulfide")
  if (str_detect(name_lower, "sulfate$")) return("Sulfate")
  if (str_detect(name_lower, "sulfonic")) return("Sulfonic Acid")
  if (str_detect(name_lower, "sulfonamide")) return("Sulfonamide")

  # Phosphorus compounds
  if (str_detect(name_lower, "phosphate$")) return("Phosphate")
  if (str_detect(name_lower, "phosphonic")) return("Phosphonic Acid")
  if (str_detect(name_lower, "phosphine")) return("Phosphine")

  # Nitro compounds
  if (str_detect(name_lower, "nitro")) return("Nitro Compound")
  if (str_detect(name_lower, "nitrate")) return("Nitrate")
  if (str_detect(name_lower, "nitrite")) return("Nitrite")

  # Cyclic compounds
  if (str_detect(name_lower, "cyclo")) return("Cyclic Compound")

  # Organometallic compounds
  if (str_detect(name_lower, "lithium|sodium|potassium|magnesium|calcium|zinc|iron|copper|silver|gold|mercury|lead|tin") && 
      !str_detect(name_lower, "chloride|fluoride|bromide|iodide|oxide|sulfate|nitrate|phosphate")) return("Organometallic Compound")

  # If no match is found
  return("Unknown/Other")
}

# This function processes a CSV file of molecule names,
# classifies each molecule, and inserts a new row with classifications as the second row
process_molecule_file_with_insertion <- function(file_path) {
  # Read the original CSV file (no headers)
  original_data <- read_csv(file_path, col_names = FALSE)

  # The first row contains molecule names
  molecule_names <- as.character(original_data[1, ])

  # Classify each molecule, keeping empty cells empty
  classifications <- sapply(molecule_names, function(name) {
    if (is.na(name) || name == "") {
      return("")
    } else {
      return(classify_molecule(name))
    }
  })

  # Prepare the new classification row
  classification_row <- data.frame(t(classifications), stringsAsFactors = FALSE)
  names(classification_row) <- names(original_data)

  # Insert the classification row as the second row
  if (nrow(original_data) >= 2) {
    result <- rbind(
      original_data[1, ],          # First row (molecule names)
      classification_row,          # New second row (classifications)
      original_data[2:nrow(original_data), ]  # Rest of original data
    )
  } else {
    # Only one row in original file
    result <- rbind(
      original_data[1, ],
      classification_row
    )
  }

  return(result)
}

# ----------- Script entry point -----------

# Set your CSV file path here
file_path <- "GC.csv"  # Change to your file path if needed

# Check if the file exists before proceeding
if (file.exists(file_path)) {
  # Process the file and add classifications
  result_with_classifications <- process_molecule_file_with_insertion(file_path)

  # Display results
  cat("File processed successfully!\n")
  cat("Original structure maintained with classifications inserted as second row.\n\n")
  cat("First few rows of the result:\n")
  print(head(result_with_classifications, n = 5))

  # Print a summary of the detected functional groups (excluding empty cells)
  classifications_only <- as.character(result_with_classifications[2, ])
  classifications_only <- classifications_only[classifications_only != "" & !is.na(classifications_only)]
  cat("\nSummary of Functional Groups:\n")
  print(table(classifications_only))

  # Save the resulting data to a new CSV file
  output_file <- "GC_with_classifications.csv"
  write_csv(result_with_classifications, output_file, col_names = FALSE)
  cat("\nResults saved to:", output_file, "\n")
  cat("The file maintains the original structure with chemical classifications added as the second row.\n")

} else {
  cat("File not found. Please check the file path:", file_path, "\n")
}

# Utility function to add classifications to any CSV file
add_classifications_to_csv <- function(input_file, output_file = NULL) {
  if (is.null(output_file)) {
    output_file <- paste0(tools::file_path_sans_ext(input_file), "_with_classifications.csv")
  }

  result <- process_molecule_file_with_insertion(input_file)
  write_csv(result, output_file, col_names = FALSE)

  cat("Classifications added to", input_file, "\n")
  cat("Output saved as:", output_file, "\n")

  return(result)
}

# Quick calculation (for interactive use)
1+1
