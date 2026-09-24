#!/usr/local/bin/python
"""Runs virulent factor prediction"""
import os
import numpy as np
from tensorflow import keras

_VIRULENCE_CHUNK = 2048


def virulent_factor_predict(list_of_proteins, NERVE_dir, inference_batch_size=1) -> list:
    """Run Adhesin prediction
    param: list_of_proteins: list
    param: NERVE_dir: str
    output: list_of_proteins: pupulates p_vir in every Protein element and returns updated list_of_proteins
    """

    model_dir = os.path.join(NERVE_dir, 'models_data/')
    std_devs = np.load(os.path.join(model_dir, 'std_devs_virulent.npy'))
    means = np.load(os.path.join(model_dir, 'means_virulent.npy'))
    projection_matrix = np.load(os.path.join(model_dir, 'projection_matrix_virulent.npy'))
    model = keras.models.load_model(os.path.join(model_dir, 'virulent_classifier.h5'))
    
    for start in range(0, len(list_of_proteins), _VIRULENCE_CHUNK):
        proteins = list_of_proteins[start:start + _VIRULENCE_CHUNK]
        data = np.array([protein.standardize(means, std_devs, projection_matrix)
                         for protein in proteins])
        predictions = model.predict(data, batch_size=inference_batch_size, verbose=0)
        for protein, prediction in zip(proteins, predictions):
            protein.p_vir = float(prediction)
        
    return list_of_proteins
