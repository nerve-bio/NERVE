#!/usr/local/bin/python
"""Run Adhesin protein preprocessing and Adhesin prediction"""
import os
import numpy as np
from tensorflow import keras
from code.Utils import Progress, run_command, build_id_index, match_protein

_ADHESIN_CHUNK = 2048


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
        output_path = os.path.join(working_dir, feature + extension)
        if os.path.exists(output_path):
            os.remove(output_path)
        run_command(['python3', os.path.join(iFeature_dir, 'iFeature.py'),
                     '--file', proteome1, '--type', feature, '--out', output_path])
        if not os.path.isfile(output_path) or os.path.getsize(output_path) == 0:
            raise RuntimeError(f'iFeature {feature} did not produce {output_path}')
    # parse files and update Protein entires
    index = build_id_index(list_of_proteins)
    for i in range(len(features)):
        with open(os.path.join(working_dir, features[i]+extension)) as f:
            next(f, None)
            matched = set()
            for line in f:
                information = line.split('\t')
                # Append to the correct protein.
                protein = match_protein(information[0], index, list_of_proteins)
                if protein is not None:
                    if protein.accession in matched:
                        raise ValueError(f'iFeature {features[i]} duplicated {protein.accession}')
                    matched.add(protein.accession)
                    protein.model_raw_data.append(np.array([float(el) for el in information[1:]]))
            if len(matched) != len(list_of_proteins):
                raise ValueError(f'iFeature {features[i]} returned {len(matched)}/{len(list_of_proteins)} proteins')

    # delete files after computation
    for file in features:
        os.remove(os.path.join(working_dir, file+extension))
    # standardization
    for protein in list_of_proteins:
        # concatenate
        protein.model_raw_data = np.concatenate(protein.model_raw_data)
    
    return list_of_proteins
    
def adhesin_predict(list_of_proteins, NERVE_dir, inference_batch_size=1)->list:
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
    total = len(list_of_proteins)
    progress = Progress('adhesin', total)
    for start in range(0, total, _ADHESIN_CHUNK):
        proteins = list_of_proteins[start:start + _ADHESIN_CHUNK]
        data = np.array([protein.standardize(means, std_devs, projection_matrix)
                         for protein in proteins])
        predictions = model.predict(data, batch_size=inference_batch_size, verbose=0)
        for protein, prediction in zip(proteins, predictions):
            protein.p_ad = float(prediction)
        progress.update(start + len(proteins))
        
    return list_of_proteins
        
