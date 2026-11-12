
import numpy as np

from qiskit.circuit.library import real_amplitudes
from qiskit.primitives import StatevectorEstimator
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA
from qiskit.quantum_info import SparsePauliOp


def run_vqe(hamiltonian, reps=2, ansatz=None, initial_point=None, maxiter=500, use_aer=False):

    # Create a default ansatz if none is provided
    if ansatz is None:
        ansatz = real_amplitudes(
            num_qubits=hamiltonian.num_qubits,
            reps=reps
        )

    # Check that the ansatz matches the Hamiltonian
    if ansatz.num_qubits != hamiltonian.num_qubits:
        raise ValueError(
            "Ansatz and Hamiltonian must have the same number of qubits"
        )

    # Set initial parameters for the optimiser
    if initial_point is None:
        initial_point = np.full(ansatz.num_parameters, 0.1)

    # Record the energy during optimisation
    energy_history = []

    def store_intermediate(eval_count, parameters, mean, metadata):
        energy = float(np.real(mean))
        energy_history.append(energy)

        print(
            f"Evaluation {eval_count}: Energy = {energy:.8f}",
            flush=True
        )

    # Choose between the original estimator and Qiskit Aer
    if use_aer:
        from qiskit import transpile
        from qiskit_aer import AerSimulator
        from qiskit_aer.primitives import EstimatorV2

        simulator = AerSimulator(method="statevector")

        # Convert the circuit into gates the simulator supports
        ansatz = transpile(
            ansatz,
            simulator,
            optimization_level=1
        )

        estimator = EstimatorV2(
            options={
                "backend_options": {"method": "statevector"}
            }
        )

    else:
        estimator = StatevectorEstimator()

    # Configure the classical optimiser
    optimizer = COBYLA(maxiter=maxiter)

    # Construct the VQE
    vqe = VQE(
        estimator=estimator,
        ansatz=ansatz,
        optimizer=optimizer,
        initial_point=initial_point,
        callback=store_intermediate
    )

    # Run the calculation
    result = vqe.compute_minimum_eigenvalue(hamiltonian)

    print("\nOptimal parameters:", result.optimal_point, flush=True)

    return float(np.real(result.eigenvalue)), energy_history


# Test the solver independently
if __name__ == "__main__":

    H = SparsePauliOp.from_list([
        ("ZI", -1.0),
        ("IZ", -1.0),
        ("XX", 0.5)
    ])

    energy, history = run_vqe(H, reps=2)

    print("VQE energy:", energy)
    print("Energy evaluations:", len(history))
