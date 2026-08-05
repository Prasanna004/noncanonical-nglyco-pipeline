# NonCanonical-NGlyco Pipeline

## Expanding the N-glycoproteome: Identification of novel non-canonical N-glycosylation sites in existing datasets using a bioinformatics approach

This repository contains the Python scripts, reference files, and example data used in our study for the identification of canonical and non-canonical N-glycosylation sites using **pGlyco3.1**. The workflow includes generation of modified protein databases, entrapment database construction, and post-processing of pGlyco3.1 search results to identify and validate canonical (**N-X-S/T**) and non-canonical (**N-X-C** and **N-X-V**) N-glycosylation sites.

# Repository Contents

```
NonCanonical-NGlyco-Pipeline
│
├── create_JXV_database.py
├── create_entrapment_database.py
├── filter_pglyco_results.py
│
├── human_data_uniprot.xlsx
├── unmodified_proteome.fasta
│
├── Example.csv
├── Example.xlsx
│
├── README.md
└── requirements.txt
```

---

# Workflow

The complete workflow implemented in this repository is illustrated below.

```
Human UniProt Reference Proteome
               │
               ▼
      Create J-X-V Database
(create_JXV_database.py)
               │
               ▼
      Modified FASTA Database
               │
               ▼
      pGlyco3.1 Database Search
               │
               ▼
     Post-processing of Results
 (filter_pglyco_results.py)
               │
               ▼
Canonical (N-X-S/T)
Non-canonical (N-X-C)
Non-canonical (N-X-V)
               │
               ▼
UniProt Annotation
Manual Validation
```

### Entrapment Analysis

```
Human UniProt Reference Proteome
               │
               ▼
 Create Entrapment Database
(create_entrapment_database.py)
               │
               ▼
 Modified Entrapment FASTA
               │
               ▼
      pGlyco3.1 Database Search
               │
               ▼
      Entrapment Analysis
```

---

# Repository Files

## 1. create_JXV_database.py

This script generates a modified UniProt human proteome by replacing every **N-X-V** motif with **J-X-V**.

### Purpose

The modified database enables pGlyco3.1 to search for potential **N-X-V glycosylation sites**.

### Input

- Human UniProt proteome FASTA

### Output

- `modified_proteome_JXV.fasta`

---

## 2. create_entrapment_database.py

Constructs an entrapment database by replacing biologically implausible motifs with placeholder residues.

### Purpose

The entrapment database is used to estimate the false-positive identification rate of non-canonical glycopeptide assignments.

### Input

- Human UniProt proteome FASTA

### Outputs

- `modified_proteome_entrapment.fasta`
- `changes_log_entrapment.txt`

---

## 3. filter_pglyco_results.py

Processes pGlyco3.1 search results.

The script performs:

- Extraction of UniProt accessions
- Identification of canonical glycopeptides (N-X-S/T)
- Identification of non-canonical N-X-C glycopeptides
- Identification of non-canonical N-X-V glycopeptides
- Duplicate removal
- UniProt glycosylation site annotation
- Export of processed Excel reports

---

# Input Files

| File | Description |
|------|-------------|
| unmodified_proteome.fasta | Human UniProt reference proteome |
| human_data_uniprot.xlsx | UniProt glycosylation annotation reference |
| Example.csv | Example pGlyco3.1 output |

---

# Output Files

The scripts generate:

- Modified J-X-V FASTA database
- Modified entrapment FASTA database
- Entrapment modification log
- Processed Excel reports
- Canonical glycopeptide identifications
- Non-canonical N-X-C glycopeptide identifications
- Non-canonical N-X-V glycopeptide identifications

---

# Software Requirements

- Python 3.10 or later

Required Python packages

```
pandas
openpyxl
```

Install using

```bash
pip install pandas openpyxl
```

---

# Example Usage

### Generate the modified J-X-V database

```bash
python create_JXV_database.py
```

### Generate the entrapment database

```bash
python create_entrapment_database.py
```

### Process pGlyco3.1 results

```bash
python filter_pglyco_results.py
```

---

# Example Data

This repository includes example input and output files to demonstrate the workflow.

| File | Description |
|------|-------------|
| Example.csv | Example pGlyco3.1 search result |
| Example.xlsx | Example processed output |


---

## Associated Manuscript

**Expanding the N-glycoproteome: Identification of novel non-canonical N-glycosylation sites in existing datasets using a bioinformatics approach**

**Journal:** *Clinical Proteomics*

**Status:** Under Review

---

# Reproducibility

The scripts in this repository reproduce the computational workflow described in the associated manuscript, including:

- Modified proteome generation
- Entrapment database construction
- pGlyco3.1 result processing
- Identification of canonical and non-canonical N-glycopeptides
- UniProt-based annotation

The complete pGlyco3.1 search outputs, processed datasets, and supplementary tables are provided with the manuscript.

---


# Citation

If you use this repository in your research, please cite the associated manuscript:

> Garapati K, Selvam PK, Ghose V, *et al.* **Expanding the N-glycoproteome: Identification of novel non-canonical N-glycosylation sites in existing datasets using a bioinformatics approach.** *Clinical Proteomics*. Under review.

---

# License

This repository is released under the **MIT License**.

---

# Contact

For questions, suggestions, or bug reports, please open a GitHub Issue or contact the corresponding authors through the details provided in the associated manuscript.
