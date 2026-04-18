#!/usr/bin/env python3
"""
Script to retrieve assembly information from NCBI for a given organism or search term.
This script uses BioPython to query NCBI's Assembly database and fetch metadata,
then processes taxonomic information and outputs to Excel files.

Usage: python get_assembly_info_20210819.py --input "Organism Name"
"""

from Bio import Entrez
from datetime import datetime
from bs4 import BeautifulSoup
from pandas import DataFrame, concat
from math import ceil
import pandas as pd
import argparse
from tqdm import tqdm

# Set your email for NCBI Entrez (required by NCBI policy)
Entrez.email = 'A.N.Other@example.com'

# Commented out custom progress bar function (using tqdm instead)
# def ProcBar(percent, StartStr='', EndStr='', TotalLength=50):
#     bar = ''.join(["\033[0;42m%s\033[0m"%' '] * int(percent * TotalLength)) + ''
#     bar = '\r' + StartStr + bar.ljust(TotalLength) + ' {:0>4.1f}%'.format(percent*100) + EndStr
#     print(bar, end='', flush=True)

def get_assembly_info(id):
    """
    Retrieve detailed information for a specific assembly ID from NCBI.

    Args:
        id (str): NCBI Assembly ID

    Returns:
        list: Assembly metadata including accession, name, species, etc.
    """
    # Fetch raw XML data from NCBI Assembly database
    raw_info_id = Entrez.esummary(db="assembly", id=id)
    record_id = raw_info_id.read()
    id_info = BeautifulSoup(record_id, "lxml")

    # Check if assembly status is suitable for download
    AssemblyStatus = id_info.find_all("assemblystatus")[0].string
    if str(AssemblyStatus) in ['Contig', 'Scaffold', 'Chromosome', 'Complete Genome']:
        # Extract key metadata fields
        AssemblyAccession = id_info.find_all("assemblyaccession")[0].string
        AssemblyName = id_info.find_all("assemblyname")[0].string
        SpeciesName = id_info.find_all("speciesname")[0].string
        SpeciesTaxid = int(id_info.find_all("speciestaxid")[0].string)
        Coverage = id_info.find_all("coverage")[0].string
        ScaffoldN50 = id_info.find_all("scaffoldn50")[0].string

        # Find total genome length from statistics
        for stat in id_info.find_all("stat"):
            if str(stat['category']) == 'total_length':
                total_length = stat.string
                break

        # Construct FTP download URL for genomic FASTA
        FtpPath = id_info.find_all("ftppath_genbank")[0].string
        Download_Ftp = FtpPath + '/' + str(FtpPath.rsplit('/', 1)[1]) + '_genomic.fna.gz'

        # Return assembled metadata as list
        return [AssemblyAccession, AssemblyName, SpeciesName, SpeciesTaxid, AssemblyStatus,
                total_length, Coverage, ScaffoldN50, Download_Ftp]
    else:
        # Skip assemblies that are not suitable
        return None

def get_taxinfo_df(taxid_list):
    """
    Fetch taxonomic lineage information for a list of TaxIDs from NCBI Taxonomy.

    Args:
        taxid_list (list): List of NCBI Taxonomy IDs

    Returns:
        DataFrame: Taxonomic information for each species
    """
    # Initialize DataFrame for taxonomic data
    taxinfo_tab = DataFrame(columns=['SpeciesTaxid', 'Phylum', 'SubPhylum', 'Class', 'Order',
                                     'SubOrder', 'InfraOrder', 'SuperFamily', 'Family', 'SubFamily', 'ScientificName'])

    total_length = len(taxid_list)
    # Process in batches of 10,000 to comply with NCBI limits
    stp = [x * 10000 for x in range(ceil(total_length / 10000))]

    for i in stp:
        print(f' Fetching taxonomic data from NCBI: {i} - {i + 10000}')
        # Fetch taxonomy records in batch
        TaxInfo = Entrez.read(Entrez.efetch(db="taxonomy", retstart=i, id=taxid_list))

        print(' Processing taxonomic information...')
        for num_id in tqdm(range(len(TaxInfo)), colour='blue'):
            taxinfo_dict = {}
            tax_id = int(TaxInfo[num_id]['TaxId'])
            taxinfo_dict['SpeciesTaxid'] = tax_id
            taxinfo_dict['ScientificName'] = TaxInfo[num_id]['ScientificName']

            # Extract specific taxonomic ranks
            for rank in ['Phylum', 'SubPhylum', 'Class', 'Order', 'SubOrder', 'InfraOrder',
                        'SuperFamily', 'Family', 'SubFamily']:
                try:
                    info_rank = [item['ScientificName'] for item in TaxInfo[num_id]['LineageEx']
                                if item['Rank'].lower() == rank.lower()][0]
                except IndexError:
                    info_rank = None
                taxinfo_dict[rank] = info_rank

            taxinfo_tab.loc[tax_id] = taxinfo_dict
        print(' ')

    print(f'Final count: {len(taxinfo_tab)} species processed.')
    return taxinfo_tab

def esearch_workfolw(Find_What):
    """
    Main workflow: Search NCBI Assembly database, fetch assembly info, get taxonomy, and output results.

    Args:
        Find_What (str): Search term for organism (e.g., "Coleoptera")
    """
    print(f"0-Starting ESearch for: {Find_What}")

    # Perform search in NCBI Assembly database
    FullRecords_id = Entrez.read(Entrez.esearch(db="assembly", RetMax=10000000, term=Find_What + "[Organism]"))
    idlist = FullRecords_id['IdList']
    print(f"1-Total assemblies found: {len(idlist)}")

    print("2-Retrieving assembly information:")
    print(' Fetching assembly details...')

    # Initialize DataFrame for assembly data
    df_assembly_info = DataFrame(columns=['AssemblyAccession', 'AssemblyName', 'SpeciesName', 'SpeciesTaxid',
                                          'AssemblyStatus', 'total_length', 'Coverage', 'ScaffoldN50', 'Download_Ftp'])

    # Fetch info for each assembly ID
    for num_id in tqdm(range(len(idlist)), colour='green'):
        try:
            id_sum = get_assembly_info(idlist[num_id])
            if id_sum:
                df_assembly_info.loc[idlist[num_id]] = id_sum
        except Exception as e:
            print(f"Error processing ID {idlist[num_id]}: {e}")
            continue

    print("")

    print("3-Retrieving taxonomic information:")
    # Get unique TaxIDs and fetch taxonomy
    taxid_list = list(set(df_assembly_info['SpeciesTaxid'].dropna()))
    print(f"Unique TaxIDs: {taxid_list}")
    df_taxinfo = get_taxinfo_df(taxid_list)

    print("4-Merging taxonomic and assembly data")
    # Merge DataFrames on TaxID
    Summarydf = pd.merge(df_taxinfo, df_assembly_info, how='right', right_on="SpeciesTaxid", left_index=True)

    print("5-Outputting summary table")
    # Sort and filter results
    Assembly_Sort_table = Summarydf.sort_values(by=['AssemblyStatus', 'total_length'], ascending=[True, False])
    raw_table_rows = len(Assembly_Sort_table.index)
    print(f" Original data count: {raw_table_rows}")

    # Remove duplicates, keeping the longest assembly per species
    Assembly_final_table = Assembly_Sort_table.drop_duplicates(subset=['SpeciesTaxid'], keep='first')
    filtered_table_rows = len(Assembly_final_table.index)
    print("     >Removed duplicate species, keeping longest assembly.")
    print(f" Filtered data count: {filtered_table_rows}")

    # Generate timestamp for output files
    data = datetime.now().strftime('%Y%m%d')
    prefix = Find_What

    # Output raw and filtered data to Excel
    Assembly_Sort_table.to_excel(f'0-rawdata_{prefix}_{data}.xlsx')
    Assembly_final_table.to_excel(f"1-{prefix}_Assembly_{data}.xlsx")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch NCBI assembly information for a given organism.")
    parser.add_argument('--input', '-i',
                        type=str,
                        required=True,
                        help='Search term or file containing search terms line by line, formatted for NCBI.')
    args = parser.parse_args()
    esearch_workfolw(args.input)
    # Example usage: esearch_workfolw('Coleoptera')
    # Another example: Staphylinidae

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i',
                        type=str,
                        help='Search content or can be a file containing search content line by line.'
                             'Format requirements according to NCBI')
    args = parser.parse_args()
    esearch_workfolw(args.input)
    #esearch_workfolw('Coleoptera')
    #Staphylinidae

