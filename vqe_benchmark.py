import numpy as np
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import real_amplitudes
from qiskit.primitives import StatevectorEstimator
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA
import matplotlib.pyplot as plt

#define a two qunit hamiltonian
H = SparsePauliOp.from_list([
    ("ZI", -1.0),
    ("IZ", -1.0),
    ("XX", 0.5)
])

#comver the Hamiltonian into a matrix
H_matrix = H.to_matrix()

#calculate its eigenvalues
eigenvalues = np.linalg.eigvalsh(H_matrix)

#lowest eigenvalue is its ground state energy
exact_energy = eigenvalues[0]

print("Hamiltonian:")
print(H_matrix)

print("\nEigenvalues:")
print(eigenvalues)

print("\nExact ground-state energy:")
print(exact_energy)


#create a parameterised quantum circuit
ansatz = real_amplitudes(
    num_qubits=2,
    reps=2
)

print("\nVQE Ansatz:")
print(ansatz.draw())

print("\nNUmber of parameters:")
print(ansatz.num_parameters)

energy_history = []

def store_intermediate_result(eval_count, parameters, mean, metadata):
    energy_history.append(float(np.real(mean)))

#configure the quatum energy estiator
estimator = StatevectorEstimator()

#configure the classical estimator
optimizer = COBYLA(maxiter=500)

#construct the VQE
vqe = VQE(
    estimator=estimator,
    ansatz=ansatz,
    optimizer=optimizer,
    initial_point=np.full(ansatz.num_parameters, 0.1),
    callback=store_intermediate_result
)

#run VQE on the hamiltonian
result = vqe.compute_minimum_eigenvalue(H)
print("\nNumber of energy evaluations")
print(len(energy_history))

print("\nVQE ground state energy:")
print(result.eigenvalue.real)

vqe_energy = result.eigenvalue.real

error = abs(vqe_energy - exact_energy)

print("\n--- Benchmark Results ---")
print(f"Exact energy: {exact_energy:.8f}")
print(f"VQE energy:   {vqe_energy:.8f}")
print(f"Absolute error: {error:.8f}")




plt.figure(figsize=(9, 5))

plt.plot(
    range(1, len(energy_history) + 1),
    energy_history,
    label="VQE energy"
)

plt.axhline(
    y=exact_energy,
    color="red",
    linestyle="--",
    label="Exact Ground-State Energy"
)

plt.xlabel("Function Evalution")
plt.ylabel("Energy")
plt.title("VQE Convergence")
plt.legend()
plt.grid(True)
plt.show()
