#-----------------------------------------------------------------------------------------------------
# Cleaning
#-----------------------------------------------------------------------------------------------------

# Import libraries
from pathlib import Path
import pandas as pd
import numpy as np

def cleaning_dataset():
    """
    Cleans raw hERG IC50 data from ChEMBL.
    Returns a deduplicated DataFrame with pIC50 values.
    """
    
    # Import dataset from csv file
    input_path = Path(__file__).parent.parent / 'data' / 'raw' / 'herg_raw.csv'
    df = pd.read_csv(input_path)
    
    # Dropping columns
    columns_to_drop = ['action_type', 'assay_variant_accession', 'assay_variant_mutation', 'data_validity_comment', 'data_validity_description', 
                   'modality', 'molecule_pref_name', 'standard_text_value', 'standard_upper_value', 'text_value', 'toid', 'upper_value',
                   'ligand_efficiency', 'activity_comment', 'qudt_units', 'uo_units', 'relation', 'units', 'value', 'document_journal', 
                   'document_year', 'pchembl_value']
    
    df = df.drop(columns = columns_to_drop)
    
    # Filtering on the units and relation against standard of the tested molecules
    df_filtered = df.query("standard_units == 'nM' and standard_relation == '='")
    df_filtered_2 = df_filtered.query("canonical_smiles.notna()").copy()
    
    # Computation of the pIC50
    df_filtered_2['pIC50'] = -np.log10(df_filtered_2['standard_value'] * 1e-9)
    
    # Investigation of the presence of duplicated and filter/aggregation on them
    duplicate = df_filtered_2.groupby(['molecule_chembl_id']).agg(
        count = ('pIC50', 'count'),
        mean = ('pIC50', 'mean'),
        std = ('pIC50', 'std')
        )
    duplicate_filter = duplicate.query("std >= 1")
    ids_duplicated = duplicate_filter.index
    df_filtered_2 = df_filtered_2.query("molecule_chembl_id not in @ids_duplicated").copy()
    
    # Aggregation of duplicate and computation to obtain the mean pIC50 on them
    columns_list = ['activity_id', 'activity_properties', 'assay_chembl_id', 'assay_description', 'assay_type', 'bao_endpoint', 'bao_format', 'bao_label', 'canonical_smiles', 'document_chembl_id',
                'parent_molecule_chembl_id', 'potential_duplicate', 'record_id', 'src_id', 'standard_flag', 'standard_type', 'target_chembl_id', 'target_organism', 'target_pref_name', 'target_tax_id', 'type']
    df_filtered_3 = df_filtered_2.groupby('molecule_chembl_id').agg(pIC50=('pIC50', 'mean'),
                                                                    **{column: (column,'first') for column in columns_list})

    # Dropping new useless columns
    columns_to_drop_2 = ['target_organism', 'standard_flag', 'type', 'standard_type']
    df_filtered_3 = df_filtered_3.drop(columns=columns_to_drop_2).copy()

    # Loading the raw data
    output_path = Path(__file__).parent.parent / 'data' / 'raw' / 'herg_cleaned.csv'
    df_filtered_3.to_csv(output_path)
    
    return df_filtered_3

if __name__ == "__main__":
    cleaning_dataset()