
import numpy as np
from qiskit.quantum_info import SparsePauliOp
from vqe_solver import run_vqe
import csv

# Artificial Hamiltonians for testing the screening loop
candidates = {
    "Test candidate A": SparsePauliOp.from_list([
        ("ZI", -1.0),
        ("IZ", -1.0),
        ("XX", 0.5)
    ]),

    "Test candidate B": SparsePauliOp.from_list([
        ("ZI", -1.0),
        ("IZ", -1.0),
        ("XX", 0.7)
    ]),

    "Test candidate C": SparsePauliOp.from_list([
        ("ZI", -1.0),
        ("IZ", -1.0),
        ("XX", 0.3)
    ])
}

# Store the results of each calculation
results = []

# Run VQE for each candidate
for name, hamiltonian in candidates.items():

    print(f"\nRunning VQE for {name}...")

    vqe_energy, history = run_vqe(hamiltonian, reps=2)

    # Exact result for benchmarking these small test systems
    exact_energy = np.linalg.eigvalsh(hamiltonian.to_matrix())[0]

    results.append({
        "candidate": name,
        "vqe_energy": vqe_energy,
        "exact_energy": exact_energy,
        "error": abs(vqe_energy - exact_energy),
        "evaluations": len(history)
    })

# Display the results
print("\n--- Screening Test Results ---")

for result in results:
    print(
        f"{result['candidate']}: "
        f"VQE = {result['vqe_energy']:.6f}, "
        f"Error = {result['error']:.2e}, "
        f"Evaluations = {result['evaluations']}"
    )

# Save screening results to CSV
with open("screening_results.csv", "w", newline="") as file:

    writer = csv.DictWriter(
        file,
        fieldnames=[
            "candidate",
            "vqe_energy",
            "exact_energy",
            "error",
            "evaluations"
        ]
    )

    writer.writeheader()
    writer.writerows(results)

print("\nResults saved to screening_results.csv")