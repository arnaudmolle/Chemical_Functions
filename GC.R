setwd("C:/Users/mollear/Desktop")
library(tidyverse)
library(readr)
library(dplyr)
library(stringr)

# Function to classify molecules based on their suffixes/prefixes
classify_molecule <- function(molecule_name) {
  # Convert to lowercase for pattern matching
  name_lower <- tolower(trimws(molecule_name))
  original_name <- trimws(molecule_name)
  
  # Define amino acids list (common 20 + some additional ones)
  amino_acids <- c("alanine", "arginine", "asparagine", "aspartic acid", "aspartate",
                   "cysteine", "glutamic acid", "glutamate", "glutamine", "glycine",
                   "histidine", "isoleucine", "leucine", "lysine", "methionine",
                   "phenylalanine", "proline", "serine", "threonine", "tryptophan",
                   "tyrosine", "valine", "hydroxyproline", "selenocysteine", "pyrrolysine",
                   "ornithine", "citrulline", "homocysteine", "taurine", "beta-alanine")
  
  # Steroid patterns
  steroid_patterns <- c("sterol", "steroid", "cholesterol", "cortisol", "testosterone", 
                        "estradiol", "progesterone", "androgen", "estrogen", "cortisone",
                        "pregnenolone", "dehydroepiandrosterone", "aldosterone")
  
  # Sugar patterns
  sugar_patterns <- c("glucose", "fructose", "galactose", "mannose", "ribose", "xylose",
                      "arabinose", "sucrose", "lactose", "maltose", "cellobiose", "trehalose",
                      "glycerol", "sorbitol", "mannitol", "xylitol", "erythritol")
  
  # Vitamin patterns
  vitamin_patterns <- c("vitamin", "retinol", "thiamine", "riboflavin", "niacin", "pyridoxine",
                        "biotin", "folate", "cobalamin", "ascorbic acid", "tocopherol", "phylloquinone",
                        "cholecalciferol", "ergocalciferol", "pantothenic acid")
  
  # Nucleotide/nucleoside patterns
  nucleotide_patterns <- c("adenosine", "guanosine", "cytidine", "uridine", "thymidine",
                           "atp", "adp", "amp", "gtp", "gdp", "gmp", "ctp", "cdp", "cmp",
                           "utp", "udp", "ump", "ttp", "tdp", "tmp")
  
  # Fatty acid patterns
  fatty_acid_patterns <- c("palmitic", "stearic", "oleic", "linoleic", "linolenic", "arachidonic",
                           "myristic", "lauric", "capric", "caprylic", "caproic", "butyric", "acetic")
  
  # Define functional group patterns (ordered by specificity)
  # Most specific patterns first to avoid misclassification
  
  # Check for specific compound classes first
  
  # Amino acids
  if (any(sapply(amino_acids, function(aa) str_detect(name_lower, paste0("\\b", aa, "\\b"))))) {
    return("Amino Acid")
  }
  
  # Heterocyclic compounds (specific prefixes)
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
  
  # Handle non-conventional naming patterns first (e.g., "2-Butenoic acid, 3-methyl-, ethyl ester")
  # Check for ester patterns anywhere in the name
  if (str_detect(name_lower, "ethyl ester|methyl ester|propyl ester|butyl ester|ester$")) return("Ester")
  if (str_detect(name_lower, "\\w+ ester")) return("Ester")
  
  # Check for amide patterns
  if (str_detect(name_lower, "\\w+ amide|amide$")) return("Amide")
  
  # Carboxylic acids and derivatives - "oic acid" specifically indicates organic acids
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
  
  # Alkenes and alkynes (unsaturated hydrocarbons)
  if (str_detect(name_lower, "triene$")) return("Triene")
  if (str_detect(name_lower, "diene$")) return("Diene")
  if (str_detect(name_lower, "yne$")) return("Alkyne")
  if (str_detect(name_lower, "ene$")) return("Alkene")
  
  # Alkanes (saturated hydrocarbons)
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
  
  # If no pattern matches
  return("Unknown/Other")
}

# Modified function to process file and insert classifications as second row
process_molecule_file_with_insertion <- function(file_path) {
  # Read the entire original file
  original_data <- read_csv(file_path, col_names = FALSE)
  
  # Get the first row (molecule names)
  molecule_names <- as.character(original_data[1, ])
  
  # Remove empty or NA values for classification, but keep original positions
  classifications <- sapply(molecule_names, function(name) {
    if (is.na(name) || name == "") {
      return("")  # Keep empty cells empty
    } else {
      return(classify_molecule(name))
    }
  })
  
  # Create the classification row as a data frame
  classification_row <- data.frame(t(classifications), stringsAsFactors = FALSE)
  names(classification_row) <- names(original_data)
  
  # Insert the classification row as the second row
  if (nrow(original_data) >= 2) {
    # Insert between first row and rest of the data
    result <- rbind(
      original_data[1, ],          # First row (molecule names)
      classification_row,          # New second row (classifications)
      original_data[2:nrow(original_data), ]  # Rest of original data
    )
  } else {
    # Only one row in original file
    result <- rbind(
      original_data[1, ],          # First row (molecule names)
      classification_row           # New second row (classifications)
    )
  }
  
  return(result)
}

# Main execution
file_path <- "GC.csv"  # Change this to your file path

# Check if file exists
if (file.exists(file_path)) {
  # Process the file
  result_with_classifications <- process_molecule_file_with_insertion(file_path)
  
  # Display results
  cat("File processed successfully!\n")
  cat("Original structure maintained with classifications inserted as second row.\n\n")
  
  # Show first few rows
  cat("First few rows of the result:\n")
  print(head(result_with_classifications, n = 5))
  
  # Summary of functional groups (excluding empty cells)
  classifications_only <- as.character(result_with_classifications[2, ])
  classifications_only <- classifications_only[classifications_only != "" & !is.na(classifications_only)]
  
  cat("\nSummary of Functional Groups:\n")
  print(table(classifications_only))
  
  # Save results to a new CSV file
  output_file <- "GC_with_classifications.csv"
  write_csv(result_with_classifications, output_file, col_names = FALSE)
  cat("\nResults saved to:", output_file, "\n")
  cat("The file maintains the original structure with chemical classifications added as the second row.\n")
  
} else {
  cat("File not found. Please check the file path:", file_path, "\n")
}

# Function to add classifications to any CSV file
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
1+1
