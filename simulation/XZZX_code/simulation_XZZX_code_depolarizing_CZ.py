import stim
import sinter
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../qec-two-level-qubits-circuit-noise-bias'))
# Here the circuit file at the circuits folder must be selected and its associated functions to generate the circuits of interest
from circuits.XZZX_surface_code_depolarizing_CZ import create_rotated_XZZX_surface_code_circuit_depolarizing_CZ, CircuitGenParametersXZZXDepolCZ
 
"""
Rotated XZZX surface code simulation with biased data qubits and depolarizing CZ and CNOT gates:
Resets and measurements are followed by flips. Similarly each data qubit is followed by a biased channel before each round. Hadamards are followed by depolarizing errors.
CNOT and CZ gates are followed by a depolarizing channel.
"""

# Generates surface code circuit tasks using Stim's circuit generation.
def generate_example_tasks(is_memory_H=False):
    # Select the bias values to simulated (here only for idling)
    etas = [0.5, 1, 10, 100, 1000, 10000]
    # Select the error rates to be simulated
    probabilities = list(np.round(np.arange(0.006, 0.0095, 0.0001), 6))
    # Select the code distances to be simulated
    distances = [5, 9, 13, 17]
    for idx,eta in enumerate(etas):
        for p in probabilities:
            for d in distances:            
                rounds = 3 * d   # number of rounds for simulation
                params = CircuitGenParametersXZZXDepolCZ(
                                                    rounds=rounds,
                                                    distance=d,
                                                    after_clifford_depolarization = p,
                                                    before_round_data_bias_probability=(p, eta),
                                                    before_measure_flip_probability = p,
                                                    after_reset_flip_probability = p,
                                                    after_CZ_depolarization=p,
                                                )
                circuit = create_rotated_XZZX_surface_code_circuit_depolarizing_CZ(params, is_memory_H=is_memory_H)
                
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
        max_shots=200_000_000,
        max_errors=200000,
        tasks=[task for task in generate_example_tasks(is_memory_H=False)] + [task for task in generate_example_tasks(is_memory_H=True)],
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