# Chemical Functions Classifier

This project provides an R script to identify the main chemical function of a molecule based on its name. It is designed to classify molecule names from a CSV file and insert the corresponding functional group as a new row in the data.

## Features

- Detects various chemical classes (amino acids, steroids, sugars, vitamins, nucleotides, fatty acids, etc.)
- Handles both conventional and non-conventional molecule naming
- Processes CSV files and inserts functional group classifications as the second row
- Provides summary statistics for the detected groups

## Getting Started

### Prerequisites

- R (tested on recent versions)
- R packages: `tidyverse`, `readr`, `dplyr`, `stringr`
  
You can install the required packages with:

```r
install.packages(c("tidyverse", "readr", "dplyr", "stringr"))
```

### Usage

1. Place your CSV file (e.g., `GC.csv`) in the working directory. The first row should contain molecule names.
2. Edit `GC.R` to set the correct file path if needed.
3. Run the script in R.

#### Example

```
Rscript GC.R
```

The script will:

- Add a new row after the first row in your CSV with the predicted functional group for each molecule.
- Save the new file as `GC_with_classifications.csv`.
- Print a summary of functional groups detected.

### Functions

- `classify_molecule(molecule_name)`: Classifies a molecule based on its name.
- `process_molecule_file_with_insertion(file_path)`: Processes a CSV, adds classifications as the second row.
- `add_classifications_to_csv(input_file, output_file = NULL)`: Adds classifications to any CSV.

### Example Output

The output CSV will look like:

| Molecule1 | Molecule2 | ... |
|-----------|-----------|-----|
| Alcohol   | Alkane    | ... |

## Contributing

Contributions, suggestions, and bug reports are welcome! Please open an issue or submit a pull request.

## License

_No license specified yet. Consider adding one (e.g., MIT, GPL) to let others know how they can use your code._

## Author

[arnaudmolle](https://github.com/arnaudmolle)
