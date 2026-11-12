
import numpy as np
import time
import sys
import csv
from pathlib import Path

from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCC

from molecular_connection import load_molecular_models
from vqe_solver import run_vqe


# Start timing
start_time = time.perf_counter()

# Load the group's prepared molecular models
print("Loading molecular models...", flush=True)

models = load_molecular_models()

# Select which molecular structure to run
structure = sys.argv[1] if len(sys.argv) > 1 else "1JFF"

model = models[structure]

hamiltonian = model["hamiltonian"]
particles = model["num_particles"]
orbitals = model["num_spatial_orbitals"]

print(f"Selected structure: {structure}")
print("Qubits:", hamiltonian.num_qubits, flush=True)
print("Particles:", particles, flush=True)
print("Spatial orbitals:", orbitals, flush=True)


# Set up the parity mapper
mapper = ParityMapper(num_particles=particles)


# Prepare the Hartree-Fock initial state
print("\nBuilding Hartree-Fock state...", flush=True)

hf_state = HartreeFock(
    orbitals,
    particles,
    mapper
)

print("Hartree-Fock state built successfully.", flush=True)




def selected_excitation(num_spatial_orbitals, num_particles):

    n_alpha, n_beta = num_particles
    n = num_spatial_orbitals

    if not (1 < n_alpha < n and 1 < n_beta < n):
        raise ValueError("Not enough occupied or virtual orbitals")

    return [
        # First excitation: HOMO to LUMO
        (
            (n_alpha - 1, n + n_beta - 1),
            (n_alpha, n + n_beta)
        ),

        # Second excitation: HOMO-1 to LUMO
        (
            (n_alpha - 2, n + n_beta - 2),
            (n_alpha, n + n_beta)
        )
    ]


# Construct a small molecular ansatz
print("\nBuilding UCC ansatz...", flush=True)

ansatz = UCC(
    num_spatial_orbitals=orbitals,
    num_particles=particles,
    excitations=selected_excitation,
    qubit_mapper=mapper,
    initial_state=hf_state
)

print("UCC ansatz built successfully.", flush=True)
print("Selected excitations:", ansatz.excitation_list, flush=True)
print("Number of parameters:", ansatz.num_parameters, flush=True)


# Run VQE
print("\nStarting molecular VQE optimisation...", flush=True)

vqe_start = time.perf_counter()

electronic_energy, history = run_vqe(
    hamiltonian,
    ansatz=ansatz,
    initial_point=np.zeros(ansatz.num_parameters),
    maxiter=40,
    use_aer=True
)

vqe_time = time.perf_counter() - vqe_start

print("\nVQE optimisation completed!", flush=True)


# Add the molecular energy correction
total_energy = electronic_energy + model["energy_offset"]

# Get the classical reference energy
reference = model["reference_total_energy"]

# Calculate the difference
error = total_energy - reference


# Display results
print("\n--- VQE Results ---")

print(f"Selected structure: {structure}", flush=True)
print(f"VQE electronic energy: {electronic_energy:.8f} Ha")
print(f"VQE total energy: {total_energy:.8f} Ha")
print(f"Classical reference: {reference:.8f} Ha")
print(f"Energy difference: {error:+.8f} Ha")
print(f"Energy evaluations: {len(history)}")

print(f"\nVQE optimisation time: {vqe_time:.2f} seconds")
print(f"Total execution time: {time.perf_counter() - start_time:.2f} seconds")



# Save the VQE results to CSV
csv_path = Path("tubulin_vqe_results.csv")

fieldnames = [
    "structure",
    "qubits",
    "pauli_terms",
    "ansatz",
    "maxiter",
    "evaluations",
    "vqe_total_energy_Ha",
    "classical_reference_energy_Ha",
    "optimization_time_s",
    "total_time_s",
    "absolute_error_mHa"
]

row = {
    "structure": structure,
    "qubits": hamiltonian.num_qubits,
    "pauli_terms": len(hamiltonian),
    "ansatz": "UCC (2 selected double excitations)",
    "maxiter": 40,
    "evaluations": len(history),
    "vqe_total_energy_Ha": total_energy,
    "classical_reference_energy_Ha": reference,
    "optimization_time_s": round(vqe_time, 2),
    "total_time_s": round(time.perf_counter() - start_time, 2),
    "absolute_error_mHa": round(abs(total_energy - reference) * 1000, 5)
}

# Append to existing results
with csv_path.open("a", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=fieldnames)

    if file.tell() == 0:
        writer.writeheader()

    writer.writerow(row)

print(f"\nResults saved to {csv_path}")


# Save the convergence history for this structure
history_path = Path(f"{structure}_convergence.csv")

with history_path.open("w", newline="") as file:
    writer = csv.writer(file)

    writer.writerow([
        "evaluation",
        "electronic_energy_Ha",
        "total_energy_Ha"
    ])

    for evaluation, energy in enumerate(history, start=1):
        writer.writerow([
            evaluation,
            energy,
            energy + model["energy_offset"]
        ])

print(f"Convergence history saved to {history_path}")

