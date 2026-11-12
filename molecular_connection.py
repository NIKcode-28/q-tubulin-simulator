
from pathlib import Path

from qiskit_nature.second_q.mappers import ParityMapper
from new_task0_h2_lih import load_taxol_active_space

DATA_FOLDER = Path("Taxol_Active_Space_Complete")
STRUCTURES = {"9WDA", "9WDB", "1JFF"}


def load_molecular_models():
    """Load molecular Hamiltonians prepared by the group."""

    files = sorted(DATA_FOLDER.rglob("base.FCIDUMP"))

    if not files:
        raise FileNotFoundError(
            f"No base.FCIDUMP files found in {DATA_FOLDER}"
        )

    models = {}

    for file in files:
        structure = file.parent.name

        if structure.upper() not in STRUCTURES:
            continue

        # Load the prepared molecular problem
        problem, report, path, offset = load_taxol_active_space(
            root=file.parent.parent,
            structure=structure,
            space="base"
        )

        # Convert the molecular Hamiltonian to a qubit Hamiltonian
        mapper = ParityMapper(num_particles=problem.num_particles)

        hamiltonian = mapper.map(
            problem.hamiltonian.second_q_op()
        ).simplify(atol=1e-12)

        # Store everything required by our screening pipeline
        models[structure] = {
            "hamiltonian": hamiltonian,
            "energy_offset": offset,
            "reference_total_energy": float(report["CASCI_Eh"]),
            "num_particles": problem.num_particles,
            "num_spatial_orbitals": problem.num_spatial_orbitals
        }

        print(f"\n{structure}")
        print(f"Qubits: {hamiltonian.num_qubits}")
        print(f"Pauli terms: {len(hamiltonian)}")
        print(f"Energy offset: {offset:.6f} Ha")
        print(f"Classical reference: {report['CASCI_Eh']:.6f} Ha")

    if not models:
        raise FileNotFoundError(
            "No matching 9WDA, 9WDB or 1JFF models found."
        )

    return models


if __name__ == "__main__":
    models = load_molecular_models()
    print(f"\nSuccessfully loaded {len(models)} molecular models.")
