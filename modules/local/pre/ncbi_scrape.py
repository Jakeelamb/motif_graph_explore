import pandas as pd
import numpy as np
import requests
import json
import time
import random
import os
import xml.etree.ElementTree as ET # Import for XML parsing
import argparse # Import argparse
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Base URL for NCBI E-utils
EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# Function to get assembly summaries in batches
def get_assembly_summaries(assembly_ids, api_key, batch_size=500):
    all_summaries = []
    for i in range(0, len(assembly_ids), batch_size):
        batch_ids = assembly_ids[i:i+batch_size]
        ids_str = ",".join(batch_ids)
        summary_url = f"{EUTILS_BASE_URL}esummary.fcgi?db=assembly&id={ids_str}&retmode=json&api_key={api_key}"
        print(f"Fetching summaries for batch {i//batch_size + 1}...")
        try:
            response = requests.get(summary_url)
            response.raise_for_status() # Raise HTTPError for bad responses (4xx or 5xx)
            data = response.json() # Use response.json() directly
            result_data = data.get('result', {})
            uids_in_batch = result_data.get('uids', [])
            batch_summaries = [result_data[uid] for uid in uids_in_batch if uid in result_data]
            all_summaries.extend(batch_summaries)
        except requests.exceptions.RequestException as e:
            print(f"Error fetching summary batch: {e}")
        except json.JSONDecodeError:
            print(f"Error decoding JSON response for summary batch: {response.text}")
        # Respect NCBI rate limits
        time.sleep(0.2) # Adjust sleep time as needed (0.1 for API key, 0.34 without)

    return all_summaries

# Function to get taxonomy lineage using efetch and parse XML
def get_taxonomy_lineages(taxids, api_key, batch_size=200):
    print(f"\nFetching taxonomy lineage for {len(taxids)} unique taxids...")
    lineage_map = {}
    ranks_to_extract = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus'] # Target ranks

    for i in range(0, len(taxids), batch_size):
        batch_taxids = taxids[i:i+batch_size]
        ids_str = ",".join(map(str, batch_taxids)) # Ensure taxids are strings
        efetch_url = f"{EUTILS_BASE_URL}efetch.fcgi?db=taxonomy&id={ids_str}&retmode=xml&api_key={api_key}"
        print(f"Fetching lineage batch {i//batch_size + 1}...")
        try:
            response = requests.get(efetch_url)
            response.raise_for_status()
            xml_root = ET.fromstring(response.content)

            for taxon in xml_root.findall('./Taxon'):
                taxid = taxon.find('TaxId').text
                lineage = {rank: None for rank in ranks_to_extract} # Initialize with None
                lineage_ex = taxon.find('./LineageEx')
                if lineage_ex is not None:
                    for lineage_taxon in lineage_ex.findall('./Taxon'):
                        rank = lineage_taxon.find('Rank').text
                        if rank in ranks_to_extract:
                            lineage[rank] = lineage_taxon.find('ScientificName').text
                # Add the genus from the main taxon if not found in lineage (sometimes happens)
                if lineage.get('genus') is None and taxon.find('Rank').text == 'species':
                     genus_name = taxon.find('ScientificName').text.split()[0]
                     lineage['genus'] = genus_name

                lineage_map[taxid] = lineage

        except requests.exceptions.RequestException as e:
            print(f"Error fetching taxonomy batch: {e}")
        except ET.ParseError as e:
            print(f"Error parsing XML response for taxonomy batch: {e}")
            print(f"Response text: {response.text[:500]}...") # Print partial response on error
        except Exception as e:
            print(f"An unexpected error occurred during taxonomy fetching/parsing: {e}")

        time.sleep(0.2) # Rate limiting

    print("Finished fetching taxonomy lineages.")
    return lineage_map

# Updated function to get all metazoan chromosome-level assemblies
def get_metazoan_chromosome_assemblies(api_key, max_ids=10000):
    print("Searching for Metazoan chromosome-level assemblies...")
    search_term = '("Metazoa"[Organism]) AND "Chromosome"[Assembly Level]'
    search_url = f"{EUTILS_BASE_URL}esearch.fcgi?db=assembly&term={search_term}&retmax={max_ids}&retmode=json&api_key={api_key}"

    assembly_ids = []
    try:
        response = requests.get(search_url)
        response.raise_for_status()
        search_data = response.json()
        assembly_ids = search_data.get('esearchresult', {}).get('idlist', [])
        count = int(search_data.get('esearchresult', {}).get('count', 0))
        print(f"Found {len(assembly_ids)} assembly IDs (out of {count} total).")
        if len(assembly_ids) >= max_ids:
            print(f"Warning: Number of results might be truncated at {max_ids}. Consider increasing max_ids or implementing pagination.")

    except requests.exceptions.RequestException as e:
        print(f"Error during ESearch: {e}")
        return pd.DataFrame() # Return empty DataFrame on search error
    except json.JSONDecodeError:
        print(f"Error decoding JSON response for ESearch: {response.text}")
        return pd.DataFrame()

    if not assembly_ids:
        print("No assembly IDs found.")
        return pd.DataFrame()

    # Get summaries for the found IDs
    print("Fetching assembly summaries...")
    assembly_summaries = get_assembly_summaries(assembly_ids, api_key)

    # Create the DataFrame
    if assembly_summaries:
        df = pd.DataFrame(assembly_summaries)
        return df
    else:
        print("Failed to retrieve any assembly summaries.")
        return pd.DataFrame()

# --- Main execution part ---

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Fetch NCBI Metazoan chromosome-level assembly data or load from CSV.")
    parser.add_argument("--generate", action="store_true",
                        help="Force fetch data from NCBI and generate/overwrite the CSV file.")
    parser.add_argument("--csv_file", default="metazoan_chromosome_assemblies_with_lineage.csv",
                        help="Path to the CSV file for loading/saving data.")
    args = parser.parse_args()

    df = pd.DataFrame() # Initialize DataFrame
    csv_filename = args.csv_file

    if args.generate:
        print("--- Generate Mode: Fetching data from NCBI ---")
        # get the api key from the .env file
        api_key = os.getenv("NCBI_API_KEY")
        if not api_key:
            print("Error: NCBI_API_KEY not found in environment variables. Cannot generate data.")
            return # Exit if no API key for generation

        # use the function to get the metazoan chromosome assemblies
        df = get_metazoan_chromosome_assemblies(api_key)

        if not df.empty and 'taxid' in df.columns:
            # Convert taxid column to string for consistency before fetching/merging
            # Handle potential NaN values before converting to string
            df['taxid'] = df['taxid'].fillna(-1).astype(int).astype(str) # Or handle as needed
            unique_taxids = df[df['taxid'] != '-1']['taxid'].unique().tolist()

            if unique_taxids:
                # Get lineage information
                lineage_info = get_taxonomy_lineages(unique_taxids, api_key)

                # Convert lineage map to DataFrame
                lineage_df = pd.DataFrame.from_dict(lineage_info, orient='index')
                lineage_df = lineage_df.reset_index().rename(columns={'index': 'taxid'})
                # Ensure taxid in lineage_df is also string for merging
                lineage_df['taxid'] = lineage_df['taxid'].astype(str)

                # Merge lineage data into the main DataFrame
                print("\nMerging taxonomy lineage information...")
                df = pd.merge(df, lineage_df, on='taxid', how='left')
            else:
                print("\nNo valid unique taxids found to fetch lineage information.")

        # Save the generated DataFrame to CSV (only in generate mode)
        if not df.empty:
            try:
                df.to_csv(csv_filename, index=False)
                print(f"\nDataFrame successfully saved to {csv_filename}")
            except Exception as e:
                print(f"\nError saving DataFrame to CSV: {e}")
        else:
             print("\nGenerated DataFrame is empty. CSV file not saved.")

    else:
        print(f"--- Load Mode: Attempting to load data from {csv_filename} ---")
        if os.path.exists(csv_filename):
            try:
                df = pd.read_csv(csv_filename)
                # Ensure 'taxid' is string after loading, handling potential NaNs read as float
                if 'taxid' in df.columns:
                    df['taxid'] = df['taxid'].apply(lambda x: str(int(x)) if pd.notna(x) else None)
                print("DataFrame successfully loaded from CSV.")
            except FileNotFoundError:
                 print(f"Error: CSV file not found at {csv_filename}.")
                 print("Use the --generate flag to create it.")
                 df = pd.DataFrame() # Ensure df is empty
            except Exception as e:
                print(f"\nError loading DataFrame from CSV: {e}")
                df = pd.DataFrame() # Ensure df is empty
        else:
            print(f"CSV file not found: {csv_filename}")
            print("Please run the script with the --generate flag first.")
            df = pd.DataFrame() # Ensure df is empty

    # --- Taxonomic Diversity Analysis ---
    print("\n--- Taxonomic Diversity Analysis ---")
    if not df.empty:
        diversity_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'organism']
        print(f"Total assemblies: {len(df)}")
        for rank in diversity_ranks:
            if rank in df.columns:
                unique_count = df[rank].nunique()
                print(f"Unique {rank.capitalize()}s: {unique_count}")
                # Optional: Print top 5 most common entries for each rank
                # print(f"  Top 5 {rank.capitalize()}s:")
                # print(df[rank].value_counts().head().to_string())
            else:
                print(f"Column '{rank}' not found for diversity analysis.")

        print("\nGoal: Sample for maximal taxonomic coverage with minimal species.")
        print("Next steps could involve implementing a sampling strategy, such as:")
        print("- A greedy approach: iteratively select species that add the most new higher taxa.")
        print("- Using phylogenetic diversity metrics if a phylogeny is available/constructible.")

    else:
        print("DataFrame is empty, skipping diversity analysis.")

    # Display information from the loaded/generated DataFrame
    print("\n--- Displaying Sample Assembly and Taxonomy Information ---") # Modified title
    if not df.empty:
        # Expanded columns list
        columns_to_print = [
            'uid', 'taxid', 'organism', 'assemblyname', 'assemblystatus',
            'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus'
        ]
        # Check if columns exist before printing
        existing_columns = [col for col in columns_to_print if col in df.columns]
        if existing_columns:
            print(f"Displaying {len(df)} entries:")
            # Optionally limit printing if DataFrame is very large
            print(df[existing_columns].head()) # Print head to avoid huge output
        else:
            print(f"Requested columns ({columns_to_print}) not found.")
            print("Available DataFrame columns:", df.columns.tolist())
    else:
        print("Final DataFrame is empty.")

# --- Script entry point ---
if __name__ == "__main__":
    main()