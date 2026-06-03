#!/usr/bin/env python

import pandas as pd
import argparse
import subprocess
import re

def get_species_keywords(species): #species is "genus species"
    """For each specie match important tissues/developemental stages."""
    #matches dictionary, important values extracted form table (description column)
    species_keyword={ #genus: keywords
    "Plasmodium": ["sporozoite", "merozoite", "ring", "trophozoite", "schizont", "gametocyte", "ookinete"],
    "Eimeria": ["oocyst", "sporozoite", "merozoite", "microgamete", "macrogamete"],
    "Toxoplasma": ["tachyzoite", "bradyzoite", "sporozoite", "oocyst"],
    "Sarcocystis": ["sarcocyst", "bradyzoite", "merozoite", "sporocyst"],
    "Cryptosporidium": ["oocyst", "sporozoite", "merozoite", "gamont"],
    "Paramecium": ["vegetative", "conjugation", "autogamy"],
    "Tetrahymena": ["vegetative", "starvation", "conjugation"]}

    for g in species_keyword: #g is "Plasmodium"
        if g.lower() in str(species).lower(): #species is "Plasmodium vivax"
            return species_keyword[g] #return important tissues/dev stage for that specie
    return [] #Specie is not in dictionary

def keyword_extraction(description, keyword_list):
    """Parses the description of a SRA and returns found keywords"""
    #search for keywords in description, store them in a list
    found_words=[]
    for word in keyword_list:
        if word in description.lower():
            found_words.append(word.capitalize())
    
    #if keywords are found
    if len(found_words)>0:
        return found_words
    
    #no keywords are found
    #special cases
    if "hpi" in description or "hrs post infection" in description:
        return ["Timecourse_HPI"]
    if "mixed" in description:
        return ["Mixed stages"]
    #none match
    return ["Unspecified"]

def get_size_gb(srr_id):
    """Returns the SRR file size in GB using vdb-dump, or None if it can't be determined."""
    try:
        result = subprocess.run(['vdb-dump', srr_id, '--info'], capture_output=True, text=True, timeout=15)
        for line in result.stdout.split('\n'):
            if line.strip().startswith('size'):
                # Extract everything after the colon
                size_str = line.split(':')[-1].strip()
                # Remove any commas or spaces just in case
                clean_size = re.sub(r'[, ]', '', size_str)

                if clean_size.isdigit():
                    size_bytes = int(clean_size)
                    size_gb = size_bytes / (1024**3)
                    print(f"Size for {srr_id}: {size_gb} GB")
                    return size_gb
    except Exception as e:
        print(f"Warning: Could not determine size for {srr_id} ({e}).")
        return None
    return None



def log_no_sra(species_label, candidates, error_file):
    """Append [Species, SRA, size] rows to the error file when no SRA survived the filters."""
    with open(error_file, "a") as fh:
        if candidates.empty:
            fh.write(f"{species_label}\tNA\tNA\n")
            return
        for _, row in candidates.iterrows():
            size_str=f"{row['Size_GB']:.2f}" if pd.notna(row["Size_GB"]) else "NA"
            fh.write(f"{row['Species']}\t{row['SRA_id']}\t{size_str}\n")

def sra_for_annotation(df, topReads=1, max_gb=9999999999999999.0, error_file=None):
    """Takes SRA with n top most reads per each Platform, tissue and developement stage of a particular species"""
    out_cols=["SRA_id", "Description", "Tissue_stage", "TaxonID", "Lineage", "Species", "Source", "Strategy", "Platform", "Read_count", "Date"]

    #species label for error logging, captured before any filtering can empty the frame
    species_label=str(df["Species"].iloc[0]) if not df.empty else "Unknown"

    #Source selection, single cell tends to be worse
    df=df[df["Source"].str.contains("TRANSCRIPTOMIC", case=False, na=False)].reset_index(drop=True)

    #failsafe: no TRANSCRIPTOMIC SRA, log the species (no candidates to size) and return empty
    if df.empty:
        print("No TRANSCRIPTOMIC SRA found. Returning empty selection.")
        if error_file:
            log_no_sra(species_label, pd.DataFrame(columns=["Species", "SRA_id", "Size_GB"]), error_file)
        return pd.DataFrame(columns=out_cols)

    #clean SRA id(keep only ID)
    df[["Exp_id", "SRA_id"]]=df["SRA_id"].str.split(":", n=1, expand=True)

    #remove srrs with huge filesize
    print(f"Evaluating {len(df)} SRA ids against the {max_gb} GB size limit.")
    df["Size_GB"]=df["SRA_id"].apply(get_size_gb) #evaluate srrs
    candidates=df[["Species", "SRA_id", "Size_GB"]].copy() #keep for error logging
    valid_size_mask=df["Size_GB"].apply(lambda s: s is not None and s <= max_gb)
    df=df[valid_size_mask].reset_index(drop=True) #select
    print(f"Candidates remaining after size filtering: {len(df)}")

    #failsafe: no SRA survived the filters, log candidates and return empty so empty files are written
    if df.empty:
        print("No SRA passed the filters. Returning empty selection.")
        if error_file:
            log_no_sra(species_label, candidates, error_file)
        return pd.DataFrame(columns=out_cols)

    #get keywords
    current_specie=df["Species"].iloc[0]
    keywords=get_species_keywords(current_specie)

    #create column for keywords per SRA(summarized description)
    df["Tissue_stage"]=df["Description"].apply(lambda description: keyword_extraction(description, keywords))
    #join list elements(remove list bracketing)
    df["Tissue_stage"]=df["Tissue_stage"].apply(lambda x: ", ".join(map(str, x)))
    
    #sort by platform and tissue/stage, keeping SRA with most reads on top
    sorted_df=df.sort_values(by=["Platform", "Tissue_stage", "Read_count"], ascending=[True, True, False])
    #select 1rst SRA(most reads) for platform and each tissue/stage
    best=sorted_df.groupby(["Platform", "Tissue_stage"]).head(topReads).reset_index(drop=True) #top reads are adjustable
    
    best=best[out_cols]
    
    return best

def main():
    #argsparse
    parser=argparse.ArgumentParser(description="Select the best SRA experiments for LyRic annotation.")

    parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input TSV file containing species specific SRA data.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to save the output TSV file.")
    parser.add_argument("-s", "--srr_id", type=str, help="Path to save the SRR indexes to a TSV file.")
    parser.add_argument("-t", "--topReads", type=int, default=2, help="Optional: Number of top SRA runs to select per group (default is 2).")
    parser.add_argument("-m", "--max_size", type=float, default=9999999999999999.0, help="Optional: Maximum file size in GB (if none is provided all are accepted).")
    parser.add_argument("-e", "--error_file", type=str, help="Optional: Path to append [Species, SRA, size] rows when no SRA pass the filters.")

    args=parser.parse_args()

    #load pandas data
    colnames=["SRA_id", "Description", "TaxonID", "Lineage", "Species", "Source", "Strategy", "Platform", "Read_count", "Date"]
    data=pd.read_csv(args.input, sep="\t", header=None, names=colnames)

    #get best sra
    best_sra=sra_for_annotation(data, topReads=args.topReads, max_gb=args.max_size, error_file=args.error_file)
    best_sra.to_csv(args.output, sep="\t", index=False, header=False)
    print(f"Saved SRA experiments({best_sra.shape[0]}) at {args.output}.")

    if args.srr_id:
        best_sra["SRA_id"].to_csv(args.srr_id, sep="\t", index=False, header=False)
        print(f"Saved SRR IDs at {args.srr_id}")

if __name__ == "__main__":
    main()
