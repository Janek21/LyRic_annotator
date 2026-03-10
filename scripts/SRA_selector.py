import pandas as pd
import argparse

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
    



def sra_for_annotation(df, topReads=1):
    """Takes SRA with n top most reads per each Platform, tissue and developement stage of a particular species"""
    #Source selection, single cell tends to be worse
    df=df[df["Source"]=="TRANSCRIPTOMIC"].reset_index(drop=True)

    #clean SRA id(keep only ID)
    df[["Exp_id", "SRA_id"]]=df["SRA_id"].str.split(":", n=1, expand=True)

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
    
    best=best[["SRA_id", "Description", "Tissue_stage", "TaxonID", "Lineage", "Species", "Source", "Strategy", "Platform", "Read_count", "Date"]]
    
    return best

def main():
    #argsparse
    parser=argparse.ArgumentParser(description="Select the best SRA experiments for LyRic annotation.")

    parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input TSV file containing species specific SRA data.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to save the output TSV file.")
    parser.add_argument("-s", "--srr_id", type=str, help="Path to save the SRR indexes to a TSV file.")
    parser.add_argument("-t", "--topReads", type=int, default=2, help="Optional: Number of top SRA runs to select per group (default is 1).")

    args=parser.parse_args()

    #load pandas data
    colnames=["SRA_id", "Description", "TaxonID", "Lineage", "Species", "Source", "Strategy", "Platform", "Read_count", "Date"]
    data=pd.read_csv(args.input, sep="\t", header=None, names=colnames)

    #get best sra
    best_sra=sra_for_annotation(data, topReads=args.topReads)
    best_sra.to_csv(args.output, sep="\t", index=False, header=False)
    print(f"Saved SRA experiments({best_sra.shape[0]}) at {args.output}.")

    if args.srr_id:
        best_sra["SRA_id"].to_csv(args.srr_id, sep="\t", index=False, header=False)
        print(f"Saved SRR IDs at {args.srr_id}")

if __name__ == "__main__":
    main()