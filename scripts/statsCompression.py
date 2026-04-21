import os, sys, argparse
import pandas as pd

def gffstats_2data(filepath):
    """Parses a gffcompare .stats file using basic string splits."""
    # Extract species name from filename
    filename=os.path.basename(filepath)
    species=filename.replace(".stats", "").replace(".cmp", "").replace(".txt", "")
    
    data={"Species": species}
    
    with open(filepath, "r") as f:
        for line in f:
            line=line.strip()
            
            # Skip empty lines or header comments
            if not line or line.startswith("#"):
                continue
            
            # 1. Parse Sensitivity and Precision (these lines contain "|")
            # Example: "Base level:    88.9     |    86.7    |"
            if "|" in line and "level:" in line:
                parts=line.split("|")
                left_side=parts[0].split(":") # Splits "Base level" and "88.9"
                
                level_name=left_side[0].strip()
                sensitivity=left_side[1].strip()
                precision=parts[1].strip()
                
                data[f"{level_name} Sensitivity"]=sensitivity
                data[f"{level_name} Precision"]=precision
            
            # 2. Parse Matching counts
            # Example: "Matching intron chains:    4921"
            elif line.startswith("Matching"):
                parts=line.split(":")
                metric=parts[0].strip()
                value=parts[1].strip()
                data[metric]=value
                
            # 3. Parse Missed and Novel features
            # Example: "Missed exons:   11535/62990    ( 18.3%)"
            elif line.startswith("Missed") or line.startswith("Novel"):
                parts=line.split(":")
                metric=parts[0].strip()
                values=parts[1].strip().split() # Splits by whitespace
                print(values)
                if len(values) >= 3:
                    fraction=values[0]
                    # Clean up the percentage string (remove parenthesis and "%")
                    percentage=values[2].replace("(", "").replace(")", "").replace("%", "")
                    
                    data[f"{metric} fraction"]=fraction
                    data[f"{metric} (%)"]=percentage

    return data

def gffcmp_compress(infiles, outfile):

    allData=[]
    #append to list data from each specie stats
    for filepath in infiles:
        print(f"Parsing {filepath}...")
        parsed_dict=gffstats_2data(filepath)
        allData.append(parsed_dict)
        #transform data to dataframe
        df=pd.DataFrame(allData)
        df=df.fillna("NA")
        #write dataframe
        df.to_csv(outfile, index=False, sep="\t")


def main():
    #argument parser
    parser=argparse.ArgumentParser(description="Compile multiple gffcompare .stats files into a table.")

    #get input (-i)
    parser.add_argument("-i", "--input", nargs="+", required=True, 
                        help="Input gffcompare .stats files (e.g., -i file1.stats file2.stats)")
    
    #get output (-o)
    parser.add_argument("-o", "--output", required=True, 
                        help="Output file name (e.g., -o summary_table.tsv)")
    
    #parse
    args=parser.parse_args()

    all_data=[]
    # Loop through the files provided to the -i argument
    for filepath in args.input:
        print(f"Parsing {filepath}...")
        parsed_dict=gffstats_2data(filepath)
        all_data.append(parsed_dict)
        
    if not all_data:
        print("No valid data parsed.")
        sys.exit(1)
        
    # Create a Pandas DataFrame from the list of dictionaries
    df=pd.DataFrame(all_data)
    
    # Optional: Fill missing values with None or "NA" if some files are missing metrics
    df=df.fillna("NA")
    
    # Save the DataFrame using the filename provided to the -o argument
    output_file=args.output
    df.to_csv(output_file, index=False, sep="\t")
    
    print(f"\nSuccess! Compiled metrics of {len(args.input)} stat files into {output_file} using Pandas.")

if __name__ == "__main__":
    main()
