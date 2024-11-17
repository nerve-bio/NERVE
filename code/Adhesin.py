#!/usr/local/bin/python
"""Run Adhesin protein preprocessing and Adhesin prediction"""
import os
import numpy as np
from tensorflow import keras
from code.Utils import bashCmdMethod

def extract_features(list_of_proteins, NERVE_dir, iFeature_dir, working_dir, proteome1) -> list:
    """Extract features from protein sequences with iFeature and store them as numpy array in Protein.model_raw_data for subsequent preprocessing and model prediction. These features will be used by adhesin and virulent factor predictos.
    param: working_dir: str
    param: iFeauture_dir: str
    param: proteome1: str
    output: list_of_proteins: pupulates model_raw_data in every Protein element and returns updated list_of_proteins
    """
    #
    features = ["AAC", "DPC", "CTDC", "CTDT", "CTDD"]
    extension = ".out"
    
    # run iFeauture
    for feature in features:
        bashCmdMethod(f"python3 {os.path.join(iFeature_dir, 'iFeature.py')} --file {proteome1} --type {feature}\
        --out {os.path.join(working_dir, feature+'.out')}")
    # parse files and update Protein entires
    for i in range(len(features)):
        with open(os.path.join(working_dir, features[i]+extension)) as f:
            lines = f.readlines()[1:]
            for line in lines:
                information = line.split('\t')
                # Append to the correct protein
                for protein in list_of_proteins:
                    if information[0] in protein.id:
                        protein.model_raw_data.append(np.array([float(el) for el in information[1:]]))
                        break

    # delete files after computation
    for file in features:
        os.remove(os.path.join(working_dir, file+extension))
    # standardization
    for protein in list_of_proteins:
        # concatenate
        protein.model_raw_data = np.concatenate(protein.model_raw_data)
    
    return list_of_proteins
    
def adhesin_predict(list_of_proteins, NERVE_dir)->list:
    """Run Adhesin prediction: str
    param: list_of_proteins: str
    param: NERVE_dir: str
    output: list_of_proteins: pupulates p_ad in every Protein element and returns updated list_of_proteins
    """
    
    model_dir = os.path.join(NERVE_dir, 'models_data/')
    std_devs = np.load(os.path.join(model_dir, 'std_devs_adhesin.npy'))
    means = np.load(os.path.join(model_dir, 'means_adhesin.npy'))
    projection_matrix = np.load(os.path.join(model_dir, 'projection_matrix_adhesin.npy'))
    model = keras.models.load_model(os.path.join(model_dir, 'adhesin_classifier.h5'))
    for protein in list_of_proteins:
        data = protein.standardize(means, std_devs, projection_matrix)
        data = np.array(data)[None, ...]
        protein.p_ad = float(model.predict(data))
        
    return list_of_proteins
        
