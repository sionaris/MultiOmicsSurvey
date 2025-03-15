#!/usr/bin/env python3
import json
from pathlib import Path
import argparse

def _get_pathway_child(pathway_data, key):
    """
    Given the pathway JSON data (a list of dictionaries), returns the value 
    for the first dictionary that contains the specified key.
    """
    for d in pathway_data:
        if key in d:
            return d[key]
    return None

def extract_uniprot_ids_from_file(cx_file):
    """
    Loads a CX file (JSON format), extracts the nodes, and returns a set of UniProt IDs.
    
    For each node, it checks for the 'alias' field and/or 'r' field.
    It splits each alias by ":" and, if present, takes the second part as the UniProt ID.
    """
    uniprot_ids = set()
    try:
        with open(cx_file, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {cx_file}: {e}")
        return uniprot_ids

    nodes = _get_pathway_child(data, 'nodes')
    if not nodes:
        return uniprot_ids

    for node in nodes:
        # First try the 'alias' field
        if 'alias' in node:
            aliases = node['alias']
            if not isinstance(aliases, list):
                aliases = [aliases]
            for a in aliases:
                parts = a.split(":")
                if len(parts) > 1 and parts[1].strip():
                    uniprot_ids.add(parts[1].strip())
        # Also check if there is an 'r' field (which may be a string or a list)
        if 'r' in node:
            r_val = node['r']
            if isinstance(r_val, list):
                for a in r_val:
                    parts = a.split(":")
                    if len(parts) > 1 and parts[1].strip():
                        uniprot_ids.add(parts[1].strip())
            elif isinstance(r_val, str):
                parts = r_val.split(":")
                if len(parts) > 1 and parts[1].strip():
                    uniprot_ids.add(parts[1].strip())
    return uniprot_ids

def main():
    parser = argparse.ArgumentParser(
        description="Extract unique UniProt IDs from CX pathway files."
    )
    parser.add_argument(
        "--cx-dir",
        type=str,
        default="data/cx",
        help="Directory containing CX files (default: data/cx)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="unique_uniprot_ids.txt",
        help="Output text file (default: unique_uniprot_ids.txt)"
    )
    args = parser.parse_args()

    cx_dir = Path(args.cx_dir)
    if not cx_dir.exists() or not cx_dir.is_dir():
        print(f"Error: {cx_dir} is not a valid directory.")
        return

    all_ids = set()
    # Iterate through all .cx files in the directory
    for file in cx_dir.glob("*.cx"):
        print(f"Processing file: {file}")
        ids = extract_uniprot_ids_from_file(file)
        print(f"  Found {len(ids)} unique UniProt IDs in this file.")
        all_ids.update(ids)
        
    print(f"\nTotal unique UniProt IDs found across all files: {len(all_ids)}")

    # Write the unique IDs to the output file, one per line.
    output_path = Path(args.output)
    with open(output_path, 'w') as out_f:
        for uid in sorted(all_ids):
            out_f.write(uid + "\n")
    print(f"Unique UniProt IDs written to {output_path}")

if __name__ == "__main__":
    main()
