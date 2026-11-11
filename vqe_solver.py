import numpy as np

from qiskit.circuit.library import real_amplitudes
from qiskit.primitives import StatevectorEstimator
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA
from qiskit.quantum_info import SparsePauliOp

def run_vqe(hamiltonian, reps=2):

    ansatz = real_amplitudes(
        num_qubits=hamiltonian.num_qubits,
        reps=reps
    )

    #record the energy during optimisation
    energy_history = []

    def store_intermediate(eval_count, parameters, mean, metadata):
        energy_history.append(float(np.real(mean)))

    estimator = StatevectorEstimator()
    optimizer = COBYLA(maxiter=500)

    #VQE configuration
    vqe = VQE(
        estimator=estimator,
        ansatz=ansatz,
        optimizer=optimizer,
        initial_point=np.full(ansatz.num_parameters, 0.1),
        callback=store_intermediate
    )


    #run the calculation
    result = vqe.compute_minimum_eigenvalue(hamiltonian)

    return result.eigenvalue.real, energy_history



if __name__ == "__main__":

    H = SparsePauliOp.from_list([
        ("ZI", -1.0),
        ("IZ", -1.0),
        ("XX", 0.5)
    ])

    energy, history = run_vqe(H, reps=2)

    print("VQE energy:", energy)
    print("Energy evaluations:", len(history))