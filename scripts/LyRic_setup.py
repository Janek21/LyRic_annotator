#!/usr/bin/env python

import os
import argparse
import sys
import shutil
import gzip

def get_shortname(species_name):
    """Helper function to generate the shortname."""
    parts=species_name.split("_")
    if len(parts) < 2:
        raise ValueError("Species name must be in 'Genus_species' format (e.g., Plasmodium_vivax).")
    return parts[0][0] + parts[1]


### Step 1: Prepare the Configuration File
def config_prep(species, config_path):
    """Replaces 'Pvivax' with the target shortname in the yaml file."""
    #default out for config_file
    if config_path==None: 
        config_path=f"{species}/config/default.yaml"

    #define shortname
    shortname=get_shortname(species)
    if not os.path.exists(config_path):
        print(f"Warning: {config_path} not found.")
        return

    with open(config_path, "r") as f:
        content=f.read()

    new_content=content.replace("shortname", shortname)

    with open(config_path, "w") as f:
        f.write(new_content)
    print(f"[Step 1] Updated {config_path} with shortname: '{shortname}'.")
    
### Step 2: On bash, speceis/genome >data/fasta/shortname.fa.gz
def filetransfer(species, file_origin, file_out): #for annotation and fasta
    """Copies a file (can also compresses it) from the remote folder to a local one"""

    #define shortname
    shortname=get_shortname(species)

    #create folder (up to last "/" as onwards its file)
    res_folder=file_out[:file_out.rfind("/")]
    os.makedirs(res_folder, exist_ok=True)

    #input path exists
    if not os.path.exists(file_origin):
        print(f"Warning: {file_origin} not found.")
        return

    #compress file
    if file_out[-3:]==".gz":

        gz_out=f"{file_out}.alt" #.alt for identifying the commpressed
        with open(file_origin, "rb") as f_in: #read straight from the origin, never write plain bytes to file_out
            with gzip.open(gz_out, "wb") as f_comp:
                shutil.copyfileobj(f_in, f_comp)
        os.rename(gz_out, file_out) #only the fully compressed file ever lands on the .gz name
    elif file_origin[-3:]==".gz":
        #decompress a gzipped source onto the plain output name
        plain_out=f"{file_out}.alt" #.alt while decompressing, never expose a partial plain file
        with gzip.open(file_origin, "rb") as f_in: #read the gz source, write plain bytes
            with open(plain_out, "wb") as f_plain:
                shutil.copyfileobj(f_in, f_plain)
        os.rename(plain_out, file_out) #only the fully decompressed file lands on the final name
    else:
        #copy genome fasta from remot location to local
        shutil.copy(file_origin, file_out)

    print(f"[Step 2] Created {file_origin} at: {file_out}")

### Step 3: Fill Sample Annotations
def sampleAnn_filler(species, srr_list_path, output_file):
    """Reads SRR IDs from TSV and generates data/sample_annotations.tsv."""

    #set up default output_file
    if output_file==None: 
        output_file=f"{species}/data/sample_annotations.tsv"
    
    #define shortname
    shortname=get_shortname(species)

    header="sample_name\tuse_matched_HiSeq\tfilter_SJ_Qscore\tseqPlatform\tcaptureDesign\tcellLine\tsubProject\tuse_dhs_peaks\tuse_cage_peaks\tuse_repeats\tlibraryPrep\n"
    
    #create out folder if it does not exist
    output_folder=output_file[:output_file.rfind("/")]#out folder is out_path until last /
    os.makedirs(output_folder, exist_ok=True)

    if not os.path.exists(srr_list_path):
        print(f"Error: SRR list file {srr_list_path} not found. Skipping annotations.")
        return
        
    count=0
    with open(srr_list_path, "r") as infile, open(output_file, "w") as outfile:
        outfile.write(header)
        
        for line in infile:
            srr_id=line.strip()
            #ensure we only grab valid SRA accessions
            if srr_id.startswith(("SRR", "ERR", "DRR")):
                row=f"ont_HpreCap_0+_{srr_id}\tFALSE\t10\tONT\tHpreCap\tSB210\t{shortname}\tFALSE\tFALSE\tTRUE\tTotalRNA\n"
                outfile.write(row)
                count+=1

    print(f"[Step 3] Generated {output_file} with {count} sample entries.")


def main():
    parser=argparse.ArgumentParser(description="Toolkit for setting up genomic project configurations and data.")
    subparsers=parser.add_subparsers(dest="command", help="Available commands to run separately")

    #Auxiliary/test shortname
    parser_sn=subparsers.add_parser("shortname", help="Just generate and print the shortname")
    parser_sn.add_argument("-s", "--species", required=True, help="Species (e.g., Plasmodium_vivax)")

    # Command 1: config file
    parser_config=subparsers.add_parser("config", help="Update the default.yaml config file")
    parser_config.add_argument("-s", "--species", required=True, help="Species name (e.g., Plasmodium_vivax)")
    parser_config.add_argument("-o", "--out_path", required=False, help="Path to output the config file")

    # Command 2: fasta file compression
    parser_transfer=subparsers.add_parser("file_transfer", help="Copy files, can also compress them")
    parser_transfer.add_argument("-s", "--species", required=True, help="Species name (e.g., Plasmodium_vivax)")
    parser_transfer.add_argument("-i", "--input", required=True, help="Path to input file")
    parser_transfer.add_argument("-o", "--out_path", required=False, help="Path to output file (.gz for compression)")

    # Command 3: sample_annotation config file
    parser_ann=subparsers.add_parser("annotate_config", help="Create config file sample_annotations.tsv")
    parser_ann.add_argument("-s", "--species", required=True, help="Species name (e.g., Plasmodium_vivax)")
    parser_ann.add_argument("-i", "--srr_list", required=True, help="Path to input SRR list")
    parser_ann.add_argument("-o", "--out_path", required=False, help="Path to output sample configs")

    # Command 4: all (convenience method)
    parser_all=subparsers.add_parser("all", help="Run all three steps sequentially")
    parser_all.add_argument("-s", "--species", required=True, help="Species name (e.g., Plasmodium_vivax)")
    parser_all.add_argument("-g", "--genome", required=True, help="Path to remote genome fasta")
    parser_all.add_argument("-a", "--annotation", required=True, help="Path to remote annotation fasta")
    parser_all.add_argument("-i", "--srr_list", required=True, help="Path to input SRR list")
    parser_all.add_argument("-oc", "--out_path_cf", required=False, help="Path to output the config file")
    parser_all.add_argument("-os", "--out_path_sa", required=False, help="Path to output sample configs")
    parser_all.add_argument("-of", "--out_path_fasta", required=False, help="Path to output the fasta file")
    parser_all.add_argument("-oa", "--out_path_annotation", required=False, help="Path to output the annotation file")

    #Parse the bash arguments
    args=parser.parse_args()

    # If no command is provided, show the help menu
    if args.command is None:
        parser.print_help()
        sys.exit(1)

    try:
        # Generate shortname centrally to ensure consistency
        shortname=get_shortname(args.species) #true input is species name, it gets shortened here and goes to functions

        if args.command=="shortname": #debugging function
            print(shortname)

        # Route to the requested function
        elif args.command=="config":
            config_prep(args.species, args.out_path)
            
        elif args.command=="file_transfer":
            filetransfer(args.species, args.input, args.out_path)
            
        elif args.command=="annotate_config":
            sampleAnn_filler(args.species, args.srr_list, args.out_path)
            
        elif args.command=="all":
            print(f"--- Starting Full Setup for {args.species} (Shortname: {shortname}) ---")
            config_prep(args.species, args.out_path_cf)
            filetransfer(args.species, args.genome, args.out_path_fasta)
            filetransfer(args.species, args.annotation, args.out_path_annotation)
            sampleAnn_filler(args.species, args.srr_list, args.out_path_sa)
            print("--- Setup Complete! ---")

    except ValueError as e:
        print(f"Execution Error: {e}")

if __name__=="__main__":
    main()
