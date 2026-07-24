#-----------------------------------------------------------------------------------------------------
# Feature engineering
#-----------------------------------------------------------------------------------------------------

# Import libraries
from pathlib import Path
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import MACCSkeys
from rdkit.Chem import DataStructs

def features_engineering():
    """
    Computes molecular features from cleaned hERG dataset.
    Returns a DataFrame with Lipinski descriptors (MW, LogP, HBD, HBA, TPSA),
    additional molecular descriptors (aromaticity, flexibility, charge),
    Morgan fingerprints (ECFP4, 2048-bit) and MACCS keys (167-bit) as additional columns.
    """
    
    # Import dataset from csv file
    input_path = Path(__file__).parent.parent / 'data' / 'raw' / 'herg_cleaned.csv'
    df = pd.read_csv(input_path, index_col='molecule_chembl_id')
    
    # Parsing SMILES into RDKit Mol objects and computing Lipinski descriptors
    df['smiles']    = df['canonical_smiles'].apply(lambda x: Chem.MolFromSmiles(x))
    df['MW']        = df['smiles'].apply(lambda x: Descriptors.MolWt(x))
    df['LogP']      = df['smiles'].apply(lambda x: Descriptors.MolLogP(x))
    df['HBD']       = df['smiles'].apply(lambda x: Descriptors.NumHDonors(x))
    df['HBA']       = df['smiles'].apply(lambda x: Descriptors.NumHAcceptors(x))
    df['TPSA']      = df['smiles'].apply(lambda x: Descriptors.TPSA(x))
    
    # Computing additional molecular descriptors relevant for hERG binding (aromaticity, flexibility, charge)
    df['RotBonds']      = df['smiles'].apply(lambda x: Descriptors.NumRotatableBonds(x))
    df['AromaticRings'] = df['smiles'].apply(lambda x: Descriptors.NumAromaticRings(x))
    df['RingCount']     = df['smiles'].apply(lambda x: Descriptors.RingCount(x))
    df['FractionCSP3']  = df['smiles'].apply(lambda x: Descriptors.FractionCSP3(x))
    df['HetAtoms']      = df['smiles'].apply(lambda x: Descriptors.NumHeteroatoms(x))
    df['MaxCharge']     = df['smiles'].apply(lambda x: Descriptors.MaxPartialCharge(x))
    df['MinCharge']     = df['smiles'].apply(lambda x: Descriptors.MinPartialCharge(x))
    df['MolRefractivity']  = df['smiles'].apply(lambda x: Descriptors.MolMR(x))
    df['Complexity']  = df['smiles'].apply(lambda x: Descriptors.BertzCT(x))
    df['KierAlpha']  = df['smiles'].apply(lambda x: Descriptors.HallKierAlpha(x))
    df['RadicalElectrons']  = df['smiles'].apply(lambda x: Descriptors.NumRadicalElectrons(x))
    
    # Computing Morgan fingerprints (ECFP4, 2048-bit) using MorganGenerator
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    df['Morgan_FP'] = df['smiles'].apply(lambda x: generator.GetFingerprint(x))
    
    # MACCS keys
    df['MACCSkeys'] = df['smiles'].apply(lambda x: MACCSkeys.GenMACCSKeys(x))
    
    # Converting RDKit fingerprint objects to numpy arrays
    def fp_to_array(fp, size): 
        arr = np.zeros(size, dtype=np.int8) 
        DataStructs.ConvertToNumpyArray(fp, arr) 
        return arr    

    # Expanding Morgan fingerprints into 2048 binary columns
    Morgan_FP_matrice = pd.DataFrame(df['Morgan_FP'].apply(lambda x: fp_to_array(x, 2048)), columns = [f'Morgan_{n}' for n in range(2048)])
    df = pd.concat([df, Morgan_FP_matrice], axis = 1)
    
    # Expanding MACCS keys into 167 binary columns
    MACCS_keys_matrice = pd.DataFrame(df['MACCSkeys'].apply(lambda x: fp_to_array(x, 167)), columns = [f'MACCSkey_{n}' for n in range(167)])
    df = pd.concat([df, MACCS_keys_matrice], axis = 1)
    
    # Exporting featured dataset
    output_path = Path(__file__).parent.parent / 'data' / 'raw' / 'herg_featured.csv'
    df.to_csv(output_path)
    print("export done!")
    return df

if __name__ == "__main__":
    features_engineering()