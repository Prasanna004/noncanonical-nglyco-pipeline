import glob
import os
import re
import sys
import pandas as pd


def extract_accession(protein_str):
    """Extract UniProt accession numbers from FASTA header formatted string."""
    if pd.isna(protein_str):
        return ""
    entries = str(protein_str).split(";")
    result = []
    for entry in entries:
        parts = entry.split("|")
        if len(parts) == 3:
            result.append(parts[1])
        else:
            result.append(entry.strip())
    return " ".join(result)


def check_prosite_reported(prosite, glycosylation):
    """Check if the specific ProSite number exists in the UniProt glycosylation list."""
    if pd.notna(prosite) and pd.notna(glycosylation):
        prosite_str = str(prosite).strip()
        glyco_numbers = re.findall(r"\d+", str(glycosylation))
        if prosite_str in glyco_numbers:
            return "reported"
    return ""


def split_prosites(df):
    """Expand rows with multiple ProSites (separated by ';') into distinct rows."""
    if "ProSites" not in df.columns:
        return df

    rows = []
    for _, row in df.iterrows():
        prosite_values = str(row["ProSites"]).split(";")
        for prosite in prosite_values:
            new_row = row.copy()
            new_row["ProSites"] = prosite.strip()
            rows.append(new_row)
    return pd.DataFrame(rows)


def load_uniprot_reference(uniprot_file="human_data_uniprot.xlsx"):
    """Load UniProt reference data into a dictionary once."""
    if os.path.exists(uniprot_file):
        print(f"Loading UniProt reference data from '{uniprot_file}'...")
        human_data = pd.read_excel(uniprot_file)
        human_data = human_data.rename(columns={"Entry": "Proteins"})
        return human_data.set_index("Proteins").to_dict(orient="index")
    else:
        print(
            f"Warning: '{uniprot_file}' not found. Glycosylation mapping will be skipped."
        )
        return {}


def process_glyco_file(csv_file, human_dict, output_file=None):
    """Process a single CSV file and save to Excel."""
    if output_file is None:
        base_name = os.path.splitext(os.path.basename(csv_file))[0]
        output_file = f"{base_name}.xlsx"

    print(f"\n--- Processing '{csv_file}' ---> '{output_file}' ---")

    # Step 1: Load Initial CSV
    data = pd.read_csv(csv_file)

    if "Proteins" in data.columns:
        data["Proteins"] = data["Proteins"].apply(extract_accession)

    def match_multiple_proteins(protein_str):
        proteins = str(protein_str).split()
        glyco_list = []
        for protein in proteins:
            if protein in human_dict:
                glyco = human_dict[protein].get("Glycosylation", "")
                if glyco and isinstance(glyco, str):
                    glyco_list.append(glyco)
        return "; ".join(glyco_list) if glyco_list else ""

    # Step 2: Define Subsets based on Motif Patterns
    sheets = {"Main": data.copy()}

    if "Peptide" in data.columns:
        amino_acids = "ACDEFGHIKLMNQRSTVWY"
        cond_st = data["Peptide"].str.contains(
            rf"J[{amino_acids}][ST]", case=False, regex=True, na=False
        )
        cond_c = data["Peptide"].str.contains(
            rf"J[{amino_acids}]C", case=False, regex=True, na=False
        )
        cond_v = data["Peptide"].str.contains(
            rf"J[{amino_acids}]V", case=False, regex=True, na=False
        )

        sheets["Peptide_J_any_ST"] = data[cond_st].copy()
        sheets["Peptide_J_any_C"] = data[cond_c].copy()
        sheets["Peptide_J_any_V"] = data[cond_v].copy()

    # Step 3: Process, Annotate, and Deduplicate
    columns_to_dedup = ["Peptide", "Glycan(H,N,A,F)", "Proteins", "ProSites"]

    with pd.ExcelWriter(output_file) as writer:
        for sheet_name, df in sheets.items():
            if (
                sheet_name != "Main"
                and "ProSites" in df.columns
                and "Proteins" in df.columns
            ):
                df = split_prosites(df)

                if human_dict:
                    df["Glycosylation"] = df["Proteins"].apply(
                        match_multiple_proteins
                    )
                    df["ProSites Reported"] = df.apply(
                        lambda row: check_prosite_reported(
                            row["ProSites"], row["Glycosylation"]
                        ),
                        axis=1,
                    )

                    cols = df.columns.tolist()
                    if "Glycosylation" in cols and "ProSites Reported" in cols:
                        glyco_idx = cols.index("Glycosylation")
                        cols.insert(
                            glyco_idx + 1,
                            cols.pop(cols.index("ProSites Reported")),
                        )
                        df = df[cols]

            existing_cols = [c for c in columns_to_dedup if c in df.columns]
            if existing_cols:
                df = df.drop_duplicates(subset=existing_cols)

            df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"Done: '{output_file}' generated.")


def process_all_csvs(directory="."):
    """Find and process all CSV files in the target directory."""
    csv_files = glob.glob(os.path.join(directory, "*.csv"))

    if not csv_files:
        print(f"No CSV files found in '{directory}'.")
        return

    print(f"Found {len(csv_files)} CSV file(s) to process.")

    # Load UniProt lookup once to avoid re-reading reference data for every CSV file
    human_dict = load_uniprot_reference()

    for csv_file in csv_files:
        process_glyco_file(csv_file, human_dict)

    print("\nAll files processed successfully!")


if __name__ == "__main__":
    # Optional: Pass target directory via command-line (e.g., `python script.py /path/to/folder`)
    # Defaults to current directory if not provided
    target_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    process_all_csvs(target_dir)