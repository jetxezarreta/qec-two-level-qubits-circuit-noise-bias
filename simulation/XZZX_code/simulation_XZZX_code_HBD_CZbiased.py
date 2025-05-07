import stim
import sinter
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../qec-two-level-qubits-circuit-noise-bias'))
# Here the circuit file at the circuits folder must be selected and its associated functions to generate the circuits of interest
# 1- XZZX_surface_code_HybridBiasCLN --> HBD with no residual bias, CNOTs add depolarizing noise
# 2- XZZX_surface_code_HybridBiasCLN_ResidualCNOT --> HBD with residual bias in CNOTs
# 3- CZcompilation_XZZX_surface_code_HybridBiasCLN --> XZZX with CZ preserving bias and compilation with only CZ gates
from circuits.XZZX_surface_code_HybridBiasCLN_ResidualCNOT import create_rotated_XZZX_surface_code_circuit_biased_noise, CircuitGenParametersXZZXBiased
 
"""
Rotated XZZX surface code simulation with HBD and HBD with remainder bias noise (CNOTS) models:
Resets and measurements are followed by flips. Similarly each data qubit is followed by a biased channel before each round. Hadamards are followed by depolarizing errors.
CZ gates are followed by a biased channel (assumed to be bias preserving). CNOTS are followed either by depolarizing (HBD) or remainder biased channel.
"""

# Generates surface code circuit tasks using Stim's circuit generation.
def generate_example_tasks(is_memory_H=False):
    # Select biases to simulate
    etas = [0.5]
    # Select remainder bias for CNOTs if such circuit is simulated
    etasCNOT = [0.5]
    # Select error rates to be simulated
    probabilities = [0.001]
    # Select distances to be simulated
    distances = [5]
    for idx, eta in enumerate(etas):
        for p in probabilities:
            for d in distances:            
                rounds = 3 * d   # number of rounds for simulation
                # This setup is for the HBD noise model with remainder errors in the CNOTs
                params = CircuitGenParametersXZZXBiased(
                                    rounds=rounds,
                                    distance=d,
                                    after_clifford_depolarization = p, # this refers to the depolarizing errors of the Hadamard and CNOT
                                    before_round_data_bias_probability = (p, eta), # biased error for idling noise
                                    before_measure_flip_probability = p, # probability of flipped measurement
                                    after_reset_flip_probability = p, # probability of preparing wrong state
                                    afterCZ_bias_probability = (p,eta), # biased error after CZ gate
                                    afterCNOT_bias_probability= (p,etasCNOT[idx]) # biased error after CNOT gates (remainder  bias)
                                )
                # This setup is fo the generic HBD model, i.e. CNOTs add depolarizing noise (they do not have any residual bias)
                # This setup is also for the compilation only using CZ gates
                # params = CircuitGenParametersXZZXBiased(
                #                                     rounds=rounds,
                #                                     distance=d,
                #                                     after_clifford_depolarization = p,
                #                                     before_round_data_bias_probability=(p, eta),
                #                                     before_measure_flip_probability = p,
                #                                     after_reset_flip_probability = p,
                #                                     afterCZ_bias_probability= (p, eta), # biased error after CZ gate
                #                                 )
                circuit = create_rotated_XZZX_surface_code_circuit_biased_noise(params, is_memory_H=is_memory_H)
                
                yield sinter.Task(
                    circuit=circuit,
                    decoder=None,
                    # detector_error_model=decoder_dem,
                    json_metadata={
                        'p': p,
                        'd': d,
                        "eta": eta,
                        "params": params.__dict__,
                        "memory": ["V", "H"][is_memory_H]}       
                        )               

def main():
    # SELECT THE NAME OF THE FOLDER FOR THE RESULTS
    filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test")
    # Collect the samples (takes a few minutes).
    samples = sinter.collect(
        num_workers=multiprocessing.cpu_count()-1,
        max_shots=20_000_000,
        max_errors=200000, # now I changed this to test!
        tasks=[task for task in generate_example_tasks(is_memory_H=True)] + [task for task in generate_example_tasks(is_memory_H=True)],
        decoders=["pymatching"],
        #count_detection_events=True,
        save_resume_filepath= os.path.join(filepath, "results.csv")
    )
    # Print samples as CSV data.
    print(sinter.CSV_HEADER)
    for sample in samples:
        print(sample.to_csv_line())


# NOTE: This is actually necessary! If the code inside 'main()' was at the
# module level, the multiprocessing children spawned by sinter.collect would
# also attempt to run that code.
if __name__ == '__main__':
    main()