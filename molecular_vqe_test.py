
import numpy as np

from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD

from vqe_solver import run_vqe

# Generate a genuine molecular Hamiltonian for H2
problem = PySCFDriver(
    atom="H 0 0 0; H 0 0 0.735",
    basis="sto3g"
).run()

# Map the electronic Hamiltonian to qubits
mapper = ParityMapper(num_particles=problem.num_particles)

hamiltonian = mapper.map(
    problem.hamiltonian.second_q_op()
)

# Prepare the chemistry-specific quantum circuit
hf_state = HartreeFock(
    problem.num_spatial_orbitals,
    problem.num_particles,
    mapper
)

ansatz = UCCSD(
    problem.num_spatial_orbitals,
    problem.num_particles,
    mapper,
    initial_state=hf_state
)

# Run your existing VQE solver
electronic_energy, history = run_vqe(
    hamiltonian,
    ansatz=ansatz,
    initial_point=np.zeros(ansatz.num_parameters)
)

# Add nuclear repulsion to obtain total molecular energy
total_energy = electronic_energy + problem.nuclear_repulsion_energy

print(f"Qubits: {hamiltonian.num_qubits}")
print(f"VQE electronic energy: {electronic_energy:.8f} Ha")
print(f"VQE total energy: {total_energy:.8f} Ha")
print(f"Energy evaluations: {len(history)}")
print(f"Difference from approximate H2 reference: {abs(total_energy + 1.137306):.2e} Ha")
