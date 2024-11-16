# GEO-Data-Processor

This Python script automates the process of downloading, parsing, and processing gene expression data from GEO datasets. It supports merging information from GPL and GSM tables, mapping gene symbols, and generating a comprehensive CSV file for downstream analysis.

---

## Features

- Downloads GEO datasets using **GEOparse**.
- Processes GPL tables to retain specific columns and map gene symbols using **mygene.info**.
- Processes GSM tables to include sample metadata and expression values.
- Merges all tables into a single CSV file, ready for analysis.

---

## Requirements

Install the required Python libraries before running the script:

```bash
pip install GEOparse pandas requests functools
```

---

## Usage

1. Run the python script

2. Enter the GEO Series Number (e.g., `GSE12345`) when prompted.

3. Follow the prompts to select relevant columns and enter metadata information.

4. The output CSV file will be saved in the current directory with the name `<GSE_No>.csv`.

---

## Functions

### 1. **get_file(GSE_No)**
Downloads the GEO dataset. If the automated download fails, the script provides instructions for manual download.

### 2. **trim_unigene(string)**
Trims unnecessary text from unigene annotations.

### 3. **get_gene_symbol(input_list)**
Maps gene symbols using the **mygene.info** API.

### 4. **process_gpl(gse)**
Processes GPL tables to retain relevant columns and map gene symbols.

### 5. **process_gsm(gse)**
Processes GSM tables to include sample metadata and expression data.

### 6. **merge_dataframes(dfs)**
Merges processed GPL and GSM tables into a single CSV file.

---

## Example Output

The output is a merged CSV file containing:

- `ID`: Unique probe or gene identifiers.
- Gene expression data for all samples in the dataset.
- Additional metadata columns based on user input.

---

## Notes

- Ensure a stable internet connection for downloading GEO datasets and mapping gene symbols.
- If the script fails to fetch data automatically, place the required `.soft.gz` file in the `./Data` directory as instructed.
- The script assumes human species for mapping gene symbols. Update the parameters in the `get_gene_symbol` function if analyzing data from other species.

---


## Contributing

Feel free to submit issues or contribute to this repository by creating pull requests.

---

## Acknowledgments

- [GEOparse](https://github.com/guma44/GEOparse) for handling GEO data.
- [mygene.info](http://mygene.info/) for gene symbol mapping.
