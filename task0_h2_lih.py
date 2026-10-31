"""Person 2, Task 0: validate the chemistry -> qubit Hamiltonian chain on H2 and LiH.

Stages (each printed so you can see where you are):
  1. PySCF driver      -> ElectronicStructureProblem (integrals + HF)
  2. Second-quantised  -> fermionic operator (a-dagger / a terms)
  3. Qubit mapping     -> Pauli-string operator (Jordan-Wigner, Parity, Bravyi-Kitaev)
  4. Exact solve       -> lowest eigenvalue (the "right answer" VQE must match)
"""
from qiskit_algorithms import NumPyMinimumEigensolver
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import (
    BravyiKitaevMapper,
    JordanWignerMapper,
    ParityMapper,
)
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer

H2_REFERENCE = -1.137  # hartree, STO-3G at 0.735 angstrom (approximate literature value)


def build_problem(atom: str, basis: str = "sto3g", charge: int = 0, spin: int = 0):
    """Stage 1: run PySCF and wrap the result as a Qiskit Nature problem."""
    driver = PySCFDriver(atom=atom, basis=basis, charge=charge, spin=spin)
    return driver.run()


def describe_problem(name: str, problem) -> None:
    """Stage 2: show what the second-quantised Hamiltonian looks like."""
    fermionic_op = problem.hamiltonian.second_q_op()
    print(f"\n=== {name}: second-quantised Hamiltonian ===")
    print(f"spatial orbitals : {problem.num_spatial_orbitals}")
    print(f"particles (a, b) : {problem.num_particles}")
    print(f"nuclear repulsion: {problem.nuclear_repulsion_energy:.6f} Ha")
    print(f"fermionic terms  : {len(fermionic_op)}")
    # Show the first few terms: '+_0 -_1' means a-dagger on spin-orbital 0, a on 1.
    for label, coeff in list(fermionic_op.items())[:4]:
        print(f"   {label:<18} {coeff:.6f}")


def mappers_for(problem):
    return {
        "JordanWigner": JordanWignerMapper(),
        "Parity (2-qubit reduction)": ParityMapper(num_particles=problem.num_particles),
        "BravyiKitaev": BravyiKitaevMapper(),
    }


def map_and_solve(name: str, problem) -> dict:
    """Stages 3 and 4: map to qubits, then solve exactly with every mapper."""
    print(f"\n=== {name}: qubit mapping and exact solution ===")
    fermionic_op = problem.hamiltonian.second_q_op()
    energies = {}
    for mapper_name, mapper in mappers_for(problem).items():
        qubit_op = mapper.map(fermionic_op)
        solver = GroundStateEigensolver(mapper, NumPyMinimumEigensolver())
        result = solver.solve(problem)
        total = float(result.total_energies[0])
        energies[mapper_name] = total
        print(
            f"{mapper_name:<28} qubits={qubit_op.num_qubits:<3} "
            f"pauli_terms={len(qubit_op):<5} E_total={total:.8f} Ha"
        )
    return energies


def main() -> None:
    # ---- H2 -------------------------------------------------------------
    h2 = build_problem("H 0 0 0; H 0 0 0.735")
    describe_problem("H2 / STO-3G", h2)
    h2_energies = map_and_solve("H2 / STO-3G", h2)

    spread = max(h2_energies.values()) - min(h2_energies.values())
    best = next(iter(h2_energies.values()))
    print(f"\nH2 mapper spread : {spread:.2e} Ha (should be ~1e-8 or smaller)")
    print(f"H2 vs reference  : {best - H2_REFERENCE:+.5f} Ha (should be well under 1e-3)")
    assert spread < 1e-6, "Mappers disagree: something is wrong with the setup."
    assert abs(best - H2_REFERENCE) < 1e-3, "H2 energy is off the known value."

    # ---- LiH with an active space ---------------------------------------
    # Full LiH in STO-3G: 6 spatial orbitals, 4 electrons -> 12 qubits (JW).
    # Active space (2e, 3o): freeze the Li 1s core and drop the highest orbitals.
    lih = build_problem("Li 0 0 0; H 0 0 1.6")
    describe_problem("LiH / STO-3G (full)", lih)
    full_energies = map_and_solve("LiH / STO-3G (full)", lih)

    transformer = ActiveSpaceTransformer(num_electrons=2, num_spatial_orbitals=3)
    lih_active = transformer.transform(lih)
    describe_problem("LiH / STO-3G (2e, 3o active space)", lih_active)
    active_energies = map_and_solve("LiH / STO-3G (2e, 3o active space)", lih_active)

    full = full_energies["JordanWigner"]
    active = active_energies["JordanWigner"]
    print(f"\nLiH full space    : {full:.6f} Ha")
    print(f"LiH active space  : {active:.6f} Ha")
    print(
        f"Active-space error: {(active - full) * 1000:+.3f} mHa "
        "(positive is expected: the active space leaves out correlation; "
        "this is why the orbital choice in Task 2 matters)"
    )
    print("\nTask 0 chain validated.")


if __name__ == "__main__":
    main()
