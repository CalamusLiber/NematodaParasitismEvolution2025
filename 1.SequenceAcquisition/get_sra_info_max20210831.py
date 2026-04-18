#!/usr/bin/env python3
"""
Script to retrieve SRA (Sequence Read Archive) information from NCBI for a given organism or search term.
This script fetches metadata for sequencing runs, filters for paired-end data, and outputs to Excel files
with taxonomic information for RNA-seq and WGS data.

Usage: python get_sra_info_max20210831.py --input "Organism Name"
"""

from Bio import Entrez
from datetime import datetime
from bs4 import BeautifulSoup
from pandas import DataFrame
from math import ceil
import pandas as pd
import argparse

# Set your email for NCBI Entrez (required by NCBI policy)
Entrez.email = 'A.N.Other@example.com'

def ProcBar(percent, StartStr='', EndStr='', TotalLength=50):
    """
    Custom progress bar function for displaying download/analysis progress.

    Args:
        percent (float): Completion percentage (0-1)
        StartStr (str): String to prepend
        EndStr (str): String to append
        TotalLength (int): Total length of progress bar
    """
    bar = ''.join(["\033[0;42m%s\033[0m" % ' '] * int(percent * TotalLength)) + ''
    bar = '\r' + StartStr + bar.ljust(TotalLength) + ' {:0>4.1f}%'.format(percent * 100) + EndStr
    print(bar, end='', flush=True)

def get_SRAinfo_df(id_list):
    """
    Fetch SRA metadata for a list of IDs from NCBI SRA database.

    Args:
        id_list (list): List of SRA IDs

    Returns:
        DataFrame: SRA metadata for paired-end sequencing runs
    """
    total_length = len(id_list)
    # Process in batches of 10,000 to comply with NCBI limits
    stp = [x * 10000 for x in range(ceil(total_length / 10000))]
    df_SRAinfo = DataFrame(columns=['Accession', 'SpeciesTaxid', 'Organism', 'Run_num', 'Library_strategy',
                                    'Library_source', 'Total_bases', 'Download_address', 'Org'])

    for i in stp:
        print(f' Fetching SRA data from NCBI: {i} - {i + 10000}')
        # Fetch XML data from SRA database
        raw_info_id = Entrez.efetch(db="sra", retstart=i, id=id_list, retmode="xml")
        print(' Loading data...')
        record_id = raw_info_id.read()
        list_info = BeautifulSoup(record_id, "lxml")
        print(' Analyzing data...')

        num = 1
        for id_info in list_info.find_all("experiment_package"):
            length = len(list_info.find_all("experiment_package"))
            ProcBar(num / length, StartStr='>>>', EndStr='|' + str(length) + ' records')
            num += 1

            # Check for paired-end layout
            library_layout = str(id_info.find_all("library_layout")[0].contents[0])
            if library_layout.endswith('</paired>'):
                # Extract metadata for paired-end runs
                taxid = int(id_info.find_all('taxon_id')[0].string)
                for srainfo in id_info.find_all("alternatives"):
                    org = srainfo['org']
                    download_address = srainfo['url']
                    accession = str(id_info.find_all("member")[0]['accession'])
                    organism = id_info.find_all("member")[0]['organism']
                    library_strategy = id_info.find_all("library_strategy")[0].string
                    library_source = id_info.find_all("library_source")[0].string
                    run_num = id_info.find_all("run")[0]['accession']
                    total_bases = int(id_info.find_all("run")[0]['total_bases'])

                # Store metadata in DataFrame
                id_sum = [accession, taxid, organism, run_num, library_strategy, library_source,
                         total_bases, download_address, org]
                df_SRAinfo.loc[accession] = id_sum

        print(' ')
    print(f'Final count: {len(df_SRAinfo)} paired-end SRA records from NCBI.')
    return df_SRAinfo

def get_taxinfo_df(taxid_list):
    """
    Fetch taxonomic lineage information for a list of TaxIDs from NCBI Taxonomy.

    Args:
        taxid_list (list): List of NCBI Taxonomy IDs

    Returns:
        DataFrame: Taxonomic information for each species
    """
    taxinfo_tab = DataFrame(columns=['tax_id', 'Phylum', 'SubPhylum', 'Class', 'Order', 'SubOrder',
                                     'InfraOrder', 'SuperFamily', 'Family', 'SubFamily', 'ScientificName'])
    total_length = len(taxid_list)
    stp = [x * 10000 for x in range(ceil(total_length / 10000))]

    for i in stp:
        print(f' Fetching taxonomic data from NCBI: {i} - {i + 10000}')
        TaxInfo = Entrez.read(Entrez.efetch(db="taxonomy", retstart=i, id=taxid_list))
        print(' Analyzing taxonomic information...')

        for num_id in range(len(TaxInfo)):
            ProcBar(num_id / len(TaxInfo), StartStr='>>>', EndStr='|' + str(len(TaxInfo)) + ' records')
            taxinfo_dict = {}
            tax_id = int(TaxInfo[num_id]['TaxId'])
            taxinfo_dict['tax_id'] = tax_id
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
    print(f'Final count: {len(taxinfo_tab)} taxonomic records.')
    return taxinfo_tab

def esearch_workfolw(Find_What):
    """
    Main workflow: Search NCBI SRA, fetch run info, get taxonomy, filter and output results.

    Args:
        Find_What (str): Search term for organism (e.g., "Coleoptera")
    """
    print(f"0-Starting ESearch for: {Find_What}")

    # Perform search in NCBI SRA database
    FullRecords_id = Entrez.read(Entrez.esearch(db="sra", RetMax=10000000, term=Find_What + "[Organism]"))
    idlist = FullRecords_id['IdList']
    print(f"1-Total SRA records found: {len(idlist)}")

    print("2-Retrieving SRA information:")
    df_SRAinfo = get_SRAinfo_df(idlist)

    print("3-Retrieving taxonomic information:")
    taxid_list = list(set(df_SRAinfo["SpeciesTaxid"].dropna()))
    df_taxinfo = get_taxinfo_df(taxid_list)

    print("4-Merging taxonomic and SRA data")
    # Merge DataFrames on TaxID
    Summarydf = pd.merge(df_taxinfo, df_SRAinfo, how='right', right_on="SpeciesTaxid", left_index=True)

    print("5-Outputting summary tables")

    # Process RNA-seq data
    print("     >Filtering for RNA-seq data and removing duplicates by organism.")
    RNA_Sort_table = Summarydf.loc[Summarydf["Library_source"] == "TRANSCRIPTOMIC"].sort_values(by=['Total_bases'], ascending=False)
    raw_rna_table_rows = len(RNA_Sort_table.index)
    print(f" Original RNA-seq data count: {raw_rna_table_rows}")

    RNA_final_table = RNA_Sort_table.drop_duplicates(subset=['Organism'], keep='first')
    filtered_table_rows = len(RNA_final_table.index)
    print(f" Filtered RNA-seq data count: {filtered_table_rows}")

    # Process WGS data
    WGS_Sort_table = Summarydf.loc[Summarydf["Library_source"] == 'GENOMIC'].sort_values(by=['Total_bases'], ascending=False)
    raw_wgs_table_rows = len(WGS_Sort_table.index)
    print(f" Original WGS data count: {raw_wgs_table_rows}")

    # Filter WGS data by minimum total bases (>5Gb)
    WGS_filter_table = WGS_Sort_table[WGS_Sort_table['Total_bases'] > 5000000000]
    WGS_final_table = WGS_filter_table.drop_duplicates(subset=['Organism'], keep='first')
    print("     >Filtered out runs with Total_bases < 5,000,000,000 and removed duplicates by organism.")
    filtered_table_rows = len(WGS_final_table.index)
    print(f" Filtered WGS data count: {filtered_table_rows}")

    # Generate timestamp for output files
    data = datetime.now().strftime('%Y%m%d')
    prefix = Find_What

    # Output data to Excel files
    Summarydf.to_excel(f'0-rawdata_{prefix}_{data}.xlsx')
    WGS_Sort_table.to_excel(f'1-WGS_all_{prefix}_{data}.xlsx')
    RNA_Sort_table.to_excel(f'1-RNA_all_{prefix}_{data}.xlsx')
    RNA_final_table.to_excel(f"2-{prefix}_RNA-Seq_{data}.xlsx")
    WGS_final_table.to_excel(f"2-{prefix}_WGS_{data}.xlsx")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch NCBI SRA information for a given organism.")
    parser.add_argument('--input', '-i',
                        type=str,
                        required=True,
                        help='Search term or file containing search terms line by line, formatted for NCBI.')
    args = parser.parse_args()
    esearch_workfolw(args.input)
    # Example usage: esearch_workfolw('Mecoptera')
