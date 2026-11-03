# Taxol active spaces: the finished proof of concept

Updated 2 November 2026. Start with **Taxol_Active_Space.ipynb**. This guide explains what the notebook does, the choices it makes, and how to use its results. The matching README.pdf contains the same information.

## 1. What is now complete?

The notebook now builds and calculates real Taxol-protein contact models for **9WDA, 9WDB and 1JFF**. It no longer waits for blank chemistry forms or externally supplied force-field templates for this defined recipe.

Implemented and executed: assembly inspection, contact extraction, sequence correspondence, chemical bond construction, protonation, peptide and quantum-region caps, hydrogen relaxation, environmental partial charges, embedded electronic calculations, automatic orbital checks, a larger active-space calculation, FCIDUMP export, and qubit-mapping verification.

The current results are in **active_space_outputs/completed/**. Older folders from the manual-input version are not used. Only three main documents are needed: the notebook and the two versions of this README. Scientific inputs and results are generated data, not extra setup documents.

## 2. Open and run

The project's existing environment has been prepared. From this folder:

```sh
source .venv-taxol-notebook/bin/activate
python -m jupyterlab Taxol_Active_Space.ipynb
```

Run the notebook from top to bottom. The delivered results are cached. The notebook rebuilds the small chemical models and reuses electronic results only when the exact inputs, numerical versions and output hashes match. It also rereads the portable file's electron counts, orbital counts and scalar constant. Set **force_recalculate=True** to repeat the expensive calculations. Restart the kernel after changing settings or code.

On another machine, make a Python environment and install:

```sh
python -m pip install "numpy>=2" scipy biopython gemmi rdkit
python -m pip install pyscf qiskit qiskit-nature nbformat
python -m pip install jupyterlab ipykernel
```

Use that environment as the notebook kernel. Internet is needed only for installation and any uncached structure or chemical-definition files. Long function cells start collapsed to keep the notebook readable. ChimeraX is optional for viewing the exported CIF/SDF structures.

Full recomputation is more demanding than reading cached results. The calculations are sequential, use a small quantum region and include explicit size limits. PySCF's memory setting is a working-memory target, not a strict cap on total process memory. Keep other heavy applications closed on a small laptop.

## 3. What molecular model did we choose?

A **pocket** is a geometric region. A **quantum region** is a chemically complete selection of atoms. An **active space** is a selection of electrons and orbitals inside the quantum region. These are different things.

The quantum region is **31 atoms**: a Taxol benzamide fragment and a methylimidazole fragment representing the contacting histidine. Both rings and the amide remain intact. Three peripheral carbon-carbon cuts are completed with link hydrogens. The whole quantum region is neutral and has 116 electrons; most remain inactive during the final active-space calculation.

The classical environment contains the remaining Taxol and capped peptide segments at 29 corresponding beta-tubulin positions. We take the union of residues found within 5 angstroms in the three structures, fill one-residue gaps to avoid overlapping caps, and explicitly add the mutation site. There are about 720 environmental charge sites per model.

**Residue numbering matters.** Human His227 is author-numbered His229 in 1JFF. Human Ser239/Cys239 corresponds to author Cys241 in 1JFF. The notebook checks sequence correspondence and saves both numbering systems. 1JFF is a bovine reference, not a matched human mutation control.

Residue 239 is about 10 angstroms from Taxol in the human structures, so a 5 angstrom extraction alone would omit it. The model includes it in the classical environment. This captures a defined local electrostatic environment, not the complete long-range allosteric mechanism.

## 4. Calculated active spaces

| Structure | Base active space | Larger active space | Qubits: base / larger |
|---|---|---|---|
| 9WDA | 12 electrons / 10 orbitals | 12 electrons / 11 orbitals | 20 / 22 |
| 9WDB | 12 electrons / 10 orbitals | 12 electrons / 11 orbitals | 20 / 22 |
| 1JFF | 12 electrons / 9 orbitals | 12 electrons / 11 orbitals | 18 / 22 |

The counts are calculation results, not targets forced to fit a device. A spatial orbital maps to two spin orbitals and therefore two Jordan-Wigner qubits before reductions. The full protein is never sent to the quantum calculation.

The base AVAS threshold is 0.10; the larger-space threshold is 0.08. Both use the same O14 and histidine NE2 valence-p targets. The larger-space calculation measures sensitivity to orbital selection. A change in energy is not a binding-energy change, and this one comparison does not prove convergence.

## 5. Assumptions and why we use them

- **Static finite environment:** keeps the calculation feasible. It omits bulk protein outside the selected fragments, solvent, waters and molecular dynamics. Nearby excluded nonprotein groups are recorded.
- **MMFF94 parameters and charges:** RDKit supplies parameters for the complete prepared organic fragments and verifies coverage. This replaces the unfinished custom OpenMM-template step. It is a practical finite-fragment choice, not a validated full-protein force field or a claim of superior accuracy.
- **Explicit protonation:** Asp/Glu are negative, Lys/Arg positive, histidine is HIE with NE2-H, and cysteine is neutral. These are conventional modelling choices, not a pKa calculation. In particular, HIE is a hypothesis supporting the selected donor geometry.
- **Chemical definitions rather than guessed bonds:** the PDB Chemical Component Dictionary provides connectivity and bond orders. Taxol stereochemistry is checked. The selected real components have their required heavy atoms; no heavy-atom repair was needed. Missing required atoms in a different input cause a clear error.
- **Fixed heavy atoms:** MMFF94 relaxes hydrogens while experimental heavy atoms and heavy cap anchors stay fixed. This retains each deposited pose. A coherent alternate conformation is selected by occupancy.
- **ACE/NME peptide caps and link-H quantum caps:** these complete the chosen fragments without cutting through the quantum model's rings or amide bond. Cap locations and atom mappings are saved.
- **Boundary charge treatment:** remove QM atoms from the charge environment, zero the outside cut endpoints and their hydrogens, then redistribute a correction within the remaining component. Up to six eligible MM heavy atoms receive the correction. The same geometric rule is used across structures; recipients can differ and are recorded. Combined formal charge is conserved. This boundary approximation needs sensitivity testing before quantitative claims.
- **Fixed electrostatic embedding:** includes the field of the assigned charges, with no mutual polarisation. A 1.2 angstrom minimum distance rejects extremely close point charges; it does not establish that boundary artefacts are absent.
- **def2-SVP, closed-shell Hartree-Fock and AVAS:** manageable starting methods for this proof of concept. Hartree-Fock uses density fitting; the active Hamiltonian uses exact integrals in those orbitals. CASCI is exact only within that specified active-space model.

The 1JFF donor-H-acceptor angle is unfavourable and is retained as a warning. Its model must not be described as a validated hydrogen bond. The human models pass that simple geometric screen, which still does not prove binding strength.

## 6. What to send for qubit mapping

Start with **active_space_outputs/completed/9WDA/base.FCIDUMP** and its accompanying **base_validation.json**. Use the equivalent file for 9WDB to work with the mutant structure. The larger-space files use the prefix **expanded**.

Each structure also has:

- **qm.sdf / qm.xyz:** prepared quantum geometry; **environment.sdf:** the full capped parent cluster used to assign charges.
- **model.json:** atom provenance, cuts, protonation, charges and boundary adjustments.
- **base_orbitals.molden / base_orbital_character.csv:** orbitals and automatic character checks. These remain available for visual inspection.
- **base_hamiltonian.npz:** numerical integrals for other software frameworks.
- **base_paulis.csv:** Jordan-Wigner Pauli labels and coefficients in hartree.
- **base_validation.json:** electron counts, qubit count, reference energy, scalar offset, mapping checks and file hash.
- **sensitivity.json:** the change between the base and larger active spaces.

The notebook's final code cell reads the real 9WDA FCIDUMP and constructs the Qiskit qubit Hamiltonian. No new unembedded molecular driver is needed.

**Add the scalar offset exactly once.** FCIDUMP stores a constant containing inactive-electron energy, quantum nuclear repulsion and nuclear/environment electrostatics. Qiskit Nature keeps this separately from the mapped electronic operator. The Pauli CSV excludes it. The exported electronic Hamiltonian does not include MM-MM, mechanical, restraint or van der Waals energies.

## 7. What was checked?

The real models passed chemical validation, MMFF94 parameter coverage, hydrogen relaxation, charge conservation, source-coordinate checks, and quantum-boundary checks. SCF converged and passed internal stability. Active orbitals passed orthonormality, electron-count and cap-weight checks. Orbital files are available for human inspection; the automatic checks are not a claim that an expert visually certified every orbital.

Both active spaces for all three structures were calculated. The exported one- and two-electron integrals and constant were checked. Reloaded FCIDUMP energies matched CASCI within 1e-7 hartree. The mapped ground-state expectation matched the classical reference within the same tolerance; the mapped eigenstate residual and particle-number leakage were also checked.

Mapping was verified by applying Pauli terms in the correct electron-number sector, not by allocating a huge full qubit matrix. A separate embedded LiH test checks frozen-electron offsets. The final notebook was run in order using the verified electronic cache, and all source and output hashes were checked.

These are software and model-consistency checks. They do not validate binding affinity, clinical efficacy, or the adequacy of the chosen fragment, protonation state, basis and embedding.

## 8. What this result can and cannot show

This is a completed **structure-to-Hamiltonian proof of concept**. It demonstrates how specified Taxol-protein contact models become active-space Hamiltonians and qubit operators, including a representation of their local environment.

Do not rank drugs or mutation effects by these total energies. Binding affinity requires appropriate bound/unbound comparisons, solvation and sampling. The experimental structures also differ in geometry. S239C is paclitaxel-sensitising and allosteric; these calculations do not establish its full mechanism.

No VQE optimisation, physical quantum-hardware execution, binding-free-energy calculation or quantum advantage is claimed. These are outside this notebook's active-space extraction scope.

## Sources

- Structures: https://www.rcsb.org/structure/9WDA ; https://www.rcsb.org/structure/9WDB ; https://www.rcsb.org/structure/1JFF
- RDKit force fields: https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html
- PySCF embedding: https://pyscf.org/user/qmmm.html
- PySCF active spaces: https://pyscf.org/user/mcscf.html
- Qiskit FCIDUMP import: https://qiskit-community.github.io/qiskit-nature/stubs/qiskit_nature.second_q.formats.fcidump_to_problem.html
