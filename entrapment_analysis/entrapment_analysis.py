import pandas as pd
import re
import csv
from collections import Counter
from pathlib import Path

# Define directories
INPUT_DIR = Path(".")
OUTPUT_DIR = Path("Results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Deduplicate ONLY on these columns
REQ_COLS = [
    "Peptide",
    "Glycan(H,N,A,F)",
    "Proteins",
    "ProSites"
]

# All standard amino acids
ALL_AA = list("ACDEFGHIKLMNPQRSTVWY")
CANONICAL_AA = ["S", "T", "C"]
ENTRAPMENT_AA = [aa for aa in ALL_AA if aa not in CANONICAL_AA]


# ============================================================
# HELPER FUNCTIONS
# ============================================================
def read_input(file_path):
    """
    Read csv/tsv with an explicit delimiter based on extension.
    """
    file_path = Path(file_path)
    if file_path.suffix.lower() == '.csv':
        delim = ','
    else:
        delim = '\t'

    print("=" * 70)
    print(f"Reading file: {file_path.name} | Delimiter: {repr(delim)}")
    print("=" * 70)

    return pd.read_csv(
        file_path,
        sep=delim,
        dtype=str,
        low_memory=False
    )
    
def clean_accessions(protein_string):
    """
    sp|P02763|A1AG1_HUMAN becomes P02763
    Multiple proteins remain separated by ;
    """
    if pd.isna(protein_string):
        return protein_string

    proteins = []
    for item in str(protein_string).split(";"):
        item = item.strip()
        fields = item.split("|")
        if len(fields) >= 2:
            proteins.append(fields[1])
        else:
            proteins.append(item)

    return ";".join(proteins)


def expand_prosites(df):
    """
    Expand 52;56;118 into separate rows
    """
    expanded = []
    for _, row in df.iterrows():
        value = row.get("ProSites")
        if pd.isna(value):
            expanded.append(row.copy())
            continue

        values = [x.strip() for x in str(value).split(";")]
        for site in values:
            newrow = row.copy()
            newrow["ProSites"] = site
            expanded.append(newrow)

    return pd.DataFrame(expanded)


def extract_j_motifs(peptide):
    if pd.isna(peptide):
        return []
    return J_PATTERN.findall(str(peptide).upper())

def motif_last_aa(motif):
    return motif[2]

def classify_motif_list(motifs):
    canonical = []
    entrapment = []
    for motif in motifs:
        aa = motif_last_aa(motif)
        if aa in CANONICAL_AA:
            canonical.append(motif)
        elif aa in ENTRAPMENT_AA:
            entrapment.append(motif)
    return canonical, entrapment

def unique_terminal_aa(motif_string):
    if motif_string == "":
        return ""
    motifs = motif_string.split("; ")
    aa = sorted(set(m[2] for m in motifs))
    return "; ".join(aa)

def count_unique_terminal_aa(x):
    if x == "":
        return 0
    return len(set(x.split("; ")))

def classify_peptide(row):
    total = row["Total_J_Motifs"]
    canon = row["Canonical_Count"]
    entrap = row["Entrapment_Count"]
    j_count = row["J_Count"]

    if total == 0:
        return "No_Motif"

    if total == 1:
        if canon == 1:
            return "Single_Canonical"
        if entrap == 1:
            if j_count > 1:
                return "Multiple_Motifs"
            return "Single_Entrapment"

    if canon > 0 and entrap > 0:
        return "Mixed_Multiple"
    if canon == total:
        return "Multiple_Canonical"
    if entrap == total:
        return "Multiple_Entrapment"

    return "Mixed"


input_files = sorted(INPUT_DIR.glob("*"))
input_files = [f for f in input_files if f.suffix.lower() in ['.csv', '.txt']]

print(f"Found {len(input_files)} input file(s) to process.")

J_PATTERN = re.compile(r"J[A-Z][A-Z]", re.IGNORECASE)


for input_file in input_files:
    print(f"\nProcessing: {input_file}")
    
    project_id = input_file.name.split("_")[0]
    output_filename = OUTPUT_DIR / f"{project_id}_Entrapment_Analysis.xlsx"

    # --- STEP 1: Read Dataset ---
    df = read_input(input_file)
    print(f"Rows loaded : {len(df):,}")

    # --- STEP 2: Clean Protein Accessions ---
    df["Proteins"] = df["Proteins"].apply(clean_accessions)

    # --- STEP 3: Expand Multi-ProSites ---
    df = expand_prosites(df)
    print(f"Rows after ProSite expansion : {len(df):,}")

    # --- STEP 4: Deduplicate ---
    dedup_cols = [c for c in REQ_COLS if c in df.columns]
    df = df.drop_duplicates(subset=dedup_cols).copy()
    print(f"Rows after deduplication : {len(df):,}")

    # --- PART 2: Detect ALL J Motifs ---
    df["Motif_List"] = df["Peptide"].apply(extract_j_motifs)
    df["Total_J_Motifs"] = df["Motif_List"].apply(len)

    df["Canonical_List"] = df["Motif_List"].apply(lambda x: classify_motif_list(x)[0])
    df["Entrapment_List"] = df["Motif_List"].apply(lambda x: classify_motif_list(x)[1])

    df["Canonical_Count"] = df["Canonical_List"].apply(len)
    df["Entrapment_Count"] = df["Entrapment_List"].apply(len)

    df["J_Count"] = df["Peptide"].fillna("").str.upper().str.count("J")

    for column in ["Motif_List", "Canonical_List", "Entrapment_List"]:
        df[column] = df[column].apply(lambda x: "; ".join(x) if len(x) else "")

    df["Unique_Motif_AA"] = df["Motif_List"].apply(unique_terminal_aa)

    # --- PART 3: Peptide Classification ---
    df["Classification"] = df.apply(classify_peptide, axis=1)
    df["Has_Canonical"] = df["Canonical_Count"] > 0
    df["Has_Entrapment"] = df["Entrapment_Count"] > 0
    df["Multiple_J"] = df["J_Count"] > 1
    df["Multiple_Motifs"] = df["Total_J_Motifs"] > 1
    df["Unique_Motif_AA_Count"] = df["Unique_Motif_AA"].apply(count_unique_terminal_aa)

    # --- PART 4: Generate Motif-Specific DataFrames & Summaries ---
    motif_dataframes = {}
    summary_rows = []

    for aa in ALL_AA:
        motif = f"J[A-Z]{aa}"
        df_motif = df[
            df["Motif_List"].str.contains(motif, regex=True, case=False, na=False)
        ].copy()

        motif_dataframes[f"Motif_{aa}"] = df_motif
        canonical_motif_count = (df_motif["Canonical_Count"] > 0).sum()

        summary_rows.append({
            "AA": aa,
            "Motif": f"J-X-{aa}",
            "Type": "Canonical" if aa in CANONICAL_AA else "Entrapment",
            "Presence of Canonical motif": int(canonical_motif_count),
            "Total_Entries": len(df_motif),
            "Unique_Peptides": df_motif["Peptide"].nunique(),
            "Unique_Proteins": df_motif["Proteins"].nunique(),
            "Unique_ProSites": df_motif["ProSites"].nunique(),
            "Single_Motif": (df_motif["Total_J_Motifs"] == 1).sum(),
            "Multiple_Motifs": (df_motif["Total_J_Motifs"] > 1).sum(),
            "Single_Canonical": (df_motif["Classification"] == "Single_Canonical").sum(),
            "Single_Entrapment": (df_motif["Classification"] == "Single_Entrapment").sum(),
            "Multiple_Canonical": (df_motif["Classification"] == "Multiple_Canonical").sum(),
            "Multiple_Entrapment": (df_motif["Classification"] == "Multiple_Entrapment").sum(),
            "Mixed_Multiple": (df_motif["Classification"] == "Mixed_Multiple").sum()
        })

    AA_summary = pd.DataFrame(summary_rows).sort_values(by=["Type", "AA"]).reset_index(drop=True)

    overall_summary = pd.DataFrame({
        "Metric": [
            "Rows after deduplication", "Unique peptides", "Unique proteins",
            "Unique ProSites", "Peptides with ≥1 motif", "Canonical motif entries",
            "Entrapment motif entries", "Total canonical motifs", "Total entrapment motifs"
        ],
        "Value": [
            len(df), df["Peptide"].nunique(), df["Proteins"].nunique(),
            df["ProSites"].nunique(), (df["Total_J_Motifs"] > 0).sum(),
            (df["Canonical_Count"] > 0).sum(), (df["Entrapment_Count"] > 0).sum(),
            df["Canonical_Count"].sum(), df["Entrapment_Count"].sum()
        ]
    })

    # --- PART 5: Save Results to Excel ---
    cols_to_drop = [
        "Total_J_Motifs", "Canonical_Count", "Entrapment_Count", "J_Count", 
        "Unique_Motif_AA", "Classification", "Has_Canonical", "Has_Entrapment", 
        "Multiple_J", "Multiple_Motifs", "Unique_Motif_AA_Count"
    ]

    df_cleaned = df.drop(columns=[c for c in cols_to_drop if c in df.columns])

    df_single_entrapment_all = df[df["Classification"] == "Single_Entrapment"].copy()
    df_single_entrapment_cleaned = df_single_entrapment_all.drop(columns=[c for c in cols_to_drop if c in df_single_entrapment_all.columns])

    df_multiple_motifs_all = df[df["Classification"] == "Multiple_Motifs"].copy()
    df_multiple_motifs_cleaned = df_multiple_motifs_all.drop(columns=[c for c in cols_to_drop if c in df_multiple_motifs_all.columns])

    with pd.ExcelWriter(output_filename, engine="openpyxl") as writer:
        df_cleaned.to_excel(writer, sheet_name="All_Processed_Data", index=False)
        AA_summary.to_excel(writer, sheet_name="AA_Summary", index=False)
        overall_summary.to_excel(writer, sheet_name="Overall_Summary", index=False)
        
        df_single_entrapment_cleaned.to_excel(writer, sheet_name="Single_Entrapment_All", index=False)
        df_multiple_motifs_cleaned.to_excel(writer, sheet_name="Multiple_Motifs", index=False)
        
        for sheet_name, motif_df in motif_dataframes.items():
            safe_sheet_name = sheet_name[:31]
            motif_cleaned = motif_df.drop(columns=[c for c in cols_to_drop if c in motif_df.columns])
            motif_cleaned.to_excel(writer, sheet_name=safe_sheet_name, index=False)

    print(f"Successfully saved: {output_filename}")

print("\n=" * 70)
print("BATCH PROCESSING COMPLETED FOR ALL FILES.")
print("=" * 70)
