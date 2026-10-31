"""Person 2, Task 0: validate the chemistry -> qubit Hamiltonian chain on H2 and LiH.

Stages (each printed so you can see where you are):
  1. PySCF driver      -> ElectronicStructureProblem (integrals + HF)
  2. Second-quantised  -> fermionic operator (a-dagger / a terms)
  3. Qubit mapping     -> Pauli-string operator (Jordan-Wigner, Parity, Bravyi-Kitaev)
  4. Exact solve       -> lowest eigenvalue (the "right answer" VQE must match)
"""
import argparse
import json
import math
from pathlib import Path


def load_quantum_dependencies():
    """Keep coordinate inspection usable without installing quantum packages."""
    global NumPyMinimumEigensolver, GroundStateEigensolver, PySCFDriver
    global BravyiKitaevMapper, JordanWignerMapper, ParityMapper, ActiveSpaceTransformer
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


def map_and_solve(name: str, problem, *, filter_criterion=None) -> dict:
    """Stages 3 and 4: map to qubits, then solve exactly with every mapper."""
    print(f"\n=== {name}: qubit mapping and exact solution ===")
    fermionic_op = problem.hamiltonian.second_q_op()
    energies = {}
    for mapper_name, mapper in mappers_for(problem).items():
        qubit_op = mapper.map(fermionic_op)
        solver = GroundStateEigensolver(mapper, NumPyMinimumEigensolver(filter_criterion=filter_criterion))
        result = solver.solve(problem)
        total = float(result.total_energies[0])
        energies[mapper_name] = total
        print(
            f"{mapper_name:<28} qubits={qubit_op.num_qubits:<3} "
            f"pauli_terms={len(qubit_op):<5} E_total={total:.8f} Ha"
        )
    return energies


def main() -> None:
    load_quantum_dependencies()
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


def load_pocket_input(path, region_number=1):
    """Return a checked PySCFDriver atom string, charge, spin and source record.

    Coordinates are in Angstrom. The electron count is for the entire prepared
    molecular region, NOT the downstream electronic active space.
    """
    pocket = json.loads(Path(path).read_text())
    regions = pocket["regions"]
    if not 1 <= region_number <= len(regions):
        raise ValueError(f"Choose a region from 1 to {len(regions)}")
    region = regions[region_number - 1]
    if not region.get("chemistry_checks_passed", False):
        raise ValueError("The selected region failed preparation chemistry checks")
    if region.get("environment") != "isolated_cluster":
        raise ValueError("Only isolated-cluster inputs are supported")
    elements = region["elements"]
    coords = region["coordinates_angstrom"]
    if not elements or len(elements) != len(coords):
        raise ValueError("Element and coordinate counts differ or are empty")
    # Elements present in the supplied taxane models; fail on unfamiliar inputs.
    atomic_numbers = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16}
    if any(e not in atomic_numbers for e in elements):
        raise ValueError("Unsupported element: extend the validated element table")
    for xyz in coords:
        if len(xyz) != 3 or not all(math.isfinite(float(v)) for v in xyz):
            raise ValueError("Invalid coordinates")
    charge, spin = region["net_charge"], region["spin_2s"]
    if type(charge) is not int or type(spin) is not int:
        raise ValueError("Charge and spin must be integers")
    electrons = sum(atomic_numbers[e] for e in elements) - charge
    if electrons != region["electron_count"]:
        raise ValueError("Stored electron count does not match elements and charge")
    if electrons < abs(spin) or (electrons - spin) % 2:
        raise ValueError("Electron count and spin are inconsistent")
    atom = "; ".join(
        f"{e} {float(x):.12f} {float(y):.12f} {float(z):.12f}"
        for e, (x, y, z) in zip(elements, coords)
    )
    return {"atom": atom, "charge": charge, "spin": spin}, region




def pocket_main(args):
    """Use prepared coordinates in the original driver/transform/map/solve chain.

    This does not discover chemically relevant orbitals. The caller supplies
    zero-based spatial MO indices from an independently reviewed selection.
    Full-cluster SCF and integral generation still precede active-space reduction.
    """
    kwargs, region = load_pocket_input(args.pocket_json, args.region)
    print("Prepared region:", region["region_id"])
    print("Atoms:", len(region["elements"]), "total electrons:", region["electron_count"])
    print("Charge:", kwargs["charge"], "spin:", kwargs["spin"], "units: Angstrom")
    print("Protonation:", region.get("protonation", {}))
    print("Basis:", args.basis)
    for warning in region.get("warnings", []):
        print("Preparation warning:", warning)
    if not args.run_pocket:
        print("Input checked; no calculations run. Electronic active space is not selected by this loader.")
        return
    if kwargs["spin"] != 0:
        raise ValueError("This extension supports the supplied closed-shell singlet models only")
    orbitals = args.active_orbitals
    ne = args.active_electrons
    if orbitals is None or ne is None:
        raise ValueError("Running requires --active-electrons and --active-orbitals (zero-based spatial MO indices)")
    if len(orbitals) != len(set(orbitals)) or not orbitals or min(orbitals) < 0:
        raise ValueError("Active orbital indices must be nonnegative and unique")
    if ne <= 0 or ne % 2 or ne > 2 * len(orbitals) or ne > region["electron_count"]:
        raise ValueError("Active electrons must be positive, even and fit the selected orbitals and molecule")
    if 2 * len(orbitals) > args.max_exact_qubits:
        raise ValueError("Active space exceeds --max-exact-qubits; exact diagonalization scales exponentially")
    print("Running full-region SCF/integrals; a small active space does not reduce this initial cost.")
    load_quantum_dependencies()
    from qiskit_nature.units import DistanceUnit
    problem = PySCFDriver(**kwargs, unit=DistanceUnit.ANGSTROM, basis=args.basis).run()
    if sum(problem.num_particles) != region["electron_count"]:
        raise ValueError("Driver electron count differs from the prepared model")
    if max(orbitals) >= problem.num_spatial_orbitals:
        raise ValueError("An active orbital index exceeds the computed orbital count")
    # The original LiH active-space stage, with explicit protein orbital choices.
    transformer = ActiveSpaceTransformer(
        num_electrons=ne, num_spatial_orbitals=len(orbitals), active_orbitals=orbitals)
    active = transformer.transform(problem)
    name = region["region_id"] + " / active space"
    describe_problem(name, active)
    # Restrict protein exact solutions to the intended particle/spin sector.
    energies = map_and_solve(name, active, filter_criterion=active.get_default_filter_criterion())
    if max(energies.values()) - min(energies.values()) >= 1e-6:
        raise ValueError("Protein active-space mapper energies disagree")
    print("Active-space mapper agreement checked; this is not a binding-affinity result.")
    return energies


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pocket-json", type=Path,
                        help="A preparation pocket.json (1JFF, 9WDA or 9WDB); omit for original H2/LiH run")
    parser.add_argument("--region", type=int, default=1,
                        help="One-based region number; alternatives are protonation hypotheses, not rankings")
    parser.add_argument("--basis", default="sto3g", help="Pocket basis; original default STO-3G is a test basis")
    parser.add_argument("--run-pocket", action="store_true", help="Run pocket calculations; otherwise inspect only")
    parser.add_argument("--active-electrons", type=int)
    parser.add_argument("--active-orbitals", type=int, nargs="+", help="Explicit zero-based spatial MO indices")
    parser.add_argument("--max-exact-qubits", type=int, default=12,
                        help="Upper bound for untapered pocket exact solve (default 12)")
    args = parser.parse_args()
    if args.pocket_json is None:
        if args.run_pocket or args.active_electrons is not None or args.active_orbitals is not None:
            parser.error("Pocket calculation options require --pocket-json")
        main()
    else:
        try:
            pocket_main(args)
        except (ValueError, KeyError, OSError) as error:
            parser.error(str(error))


if __name__ == "__main__":
    cli()
