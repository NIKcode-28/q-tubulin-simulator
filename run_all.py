"""Person 2 pipeline: config.yaml -> HF -> (active space) -> Hamiltonian -> qubits -> exact energies.

Usage:  python run_all.py [config.yaml]

For each system in the config it writes, into <output_dir>/<system>/:
  qubit_op.json         Pauli strings + coefficients (real and imaginary parts)
  hamiltonian_meta.json mapper, qubits, Pauli terms, energy shift, particle numbers
  reference_energies.json  HF, CASCI (PySCF) and qubit-operator exact energies
"""
import json
import sys
from pathlib import Path

import numpy as np
import yaml
from pyscf import gto, mcscf, scf
from qiskit_algorithms import NumPyMinimumEigensolver
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import (
    BravyiKitaevMapper,
    JordanWignerMapper,
    ParityMapper,
)
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer


def load_config(path="config.yaml"):
    with open(path) as f:
        return yaml.safe_load(f)


def get_mapper(name, problem):
    if name == "jordan_wigner":
        return JordanWignerMapper()
    if name == "parity":
        return ParityMapper(num_particles=problem.num_particles)
    if name == "bravyi_kitaev":
        return BravyiKitaevMapper()
    raise ValueError(f"Unknown mapper '{name}'")


def atom_spec(system):
    """Atom string from `atom:` or from an `xyz:` file path."""
    if system.get("xyz"):
        lines = Path(system["xyz"]).read_text().splitlines()[2:]
        return "; ".join(l.strip() for l in lines if l.strip())
    return system["atom"]


def pyscf_references(atom, scf_cfg, active):
    """HF and CASCI energies straight from PySCF (independent of Qiskit Nature)."""
    mol = gto.M(atom=atom, basis=scf_cfg["basis"], charge=scf_cfg["charge"],
                spin=scf_cfg["spin"], verbose=0)
    mf = scf.RHF(mol).run()
    out = {"hf": float(mf.e_tot)}
    if active:
        cas = mcscf.CASCI(mf, active["num_spatial_orbitals"], active["num_electrons"]).run()
    else:
        cas = mcscf.CASCI(mf, mol.nao, mol.nelectron).run()
    out["casci"] = float(cas.e_tot)
    out["correlation"] = out["casci"] - out["hf"]
    return out


def run_system(name, system, cfg):
    scf_cfg = cfg["scf"]
    atom = atom_spec(system)
    active = system.get("active_space")

    # 1-2. PySCF: AO integrals, RHF, AO->MO transform, wrapped as a problem
    problem = PySCFDriver(atom=atom, basis=scf_cfg["basis"],
                          charge=scf_cfg["charge"], spin=scf_cfg["spin"]).run()

    # 3. Optional active space (folds frozen orbitals into a constant + effective 1e term)
    if active:
        problem = ActiveSpaceTransformer(
            num_electrons=active["num_electrons"],
            num_spatial_orbitals=active["num_spatial_orbitals"],
        ).transform(problem)

    # 4. Second quantisation -> qubit operator
    mapper = get_mapper(cfg["mapper"], problem)
    qubit_op = mapper.map(problem.hamiltonian.second_q_op())

    # 5. Exact energy of the qubit operator. Constants: nuclear repulsion + frozen-core shift.
    shift = float(sum(problem.hamiltonian.constants.values()))
    e_qubit = float(NumPyMinimumEigensolver().compute_minimum_eigenvalue(qubit_op).eigenvalue.real) + shift

    refs = pyscf_references(atom, scf_cfg, active)
    refs["qubit_exact"] = e_qubit
    refs["qubit_minus_casci"] = e_qubit - refs["casci"]

    # Checks
    tol = cfg["tolerances"]
    ok = abs(refs["qubit_minus_casci"]) < tol["casci_vs_qubit_hartree"]
    ref_val = system.get("reference_energy")
    if ref_val is not None:
        ok = ok and abs(e_qubit - ref_val) < tol["reference_hartree"]

    # Save
    outdir = Path(cfg["output_dir"]) / name
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "qubit_op.json").write_text(json.dumps(
        [{"pauli": p.to_label(), "re": float(c.real), "im": float(c.imag)}
         for p, c in zip(qubit_op.paulis, qubit_op.coeffs)], indent=1))
    meta = {"system": name, "atom": atom, "basis": scf_cfg["basis"], "mapper": cfg["mapper"],
            "num_qubits": qubit_op.num_qubits, "pauli_terms": len(qubit_op),
            "energy_shift_hartree": shift, "num_particles": list(problem.num_particles),
            "num_spatial_orbitals": problem.num_spatial_orbitals, "active_space": active}
    (outdir / "hamiltonian_meta.json").write_text(json.dumps(meta, indent=1))
    (outdir / "reference_energies.json").write_text(json.dumps(refs, indent=1))

    print(f"{name:<11} qubits={meta['num_qubits']:<3} terms={meta['pauli_terms']:<4} "
          f"HF={refs['hf']:.6f}  CASCI={refs['casci']:.6f}  qubit={e_qubit:.6f}  "
          f"diff={refs['qubit_minus_casci']:+.1e}  {'OK' if ok else 'CHECK FAILED'}")
    return ok


def main():
    cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else "config.yaml")
    print(f"mapper={cfg['mapper']}  basis={cfg['scf']['basis']}\n")
    results = [run_system(n, s, cfg) for n, s in cfg["systems"].items()]
    print("\nAll checks passed." if all(results) else "\nSome checks FAILED.")
    sys.exit(0 if all(results) else 1)


if __name__ == "__main__":
    main()
