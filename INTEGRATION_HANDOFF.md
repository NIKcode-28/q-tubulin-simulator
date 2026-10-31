# Integration handoff: Taxol–tubulin molecular-region preparation

Updated 31 October 2026. Intended audience: the next developer/researcher implementing electronic-structure calculations, electronic active-space selection, fermionic Hamiltonians, qubit mapping and subsequent solvers.

## 1. Current state and agreed scope

The working implementation is `human_taxane_pocket_interface.ipynb`. It downloads and verifies selected structures, handles deposited finite assemblies, extracts bound Taxol pockets, prepares capped isolated molecular regions, and returns machine-readable inputs.

The default run processes only:

| ID | Role | Resolution |
|---|---|---:|
| 1JFF | Approved bovine reference; the only nonhuman exception | 3.50 Å |
| 9WDA | Human alpha1A/beta3 tubulin with paclitaxel | 2.35 Å |
| 9WDB | Human alpha1A/beta3 S239C comparator with paclitaxel | 2.39 Å |

Both human structures represent recombinant human tubulin; they are not patient-tumor structures. S239C is a sensitizing experimental comparator in the associated study, not a proven naturally occurring cancer-resistance mutation. Other taxanes, docking, unoccupied-pocket prediction and arbitrary modified residues are not implemented.

**User decision to preserve:** include chemically prepared human regions in the downstream inputs even when structural-density flags are present. Keep those flags; do not silently mark them resolved. Chemistry failures remain excluded.

## 2. What is complete versus what remains

| Stage | Status |
|---|---|
| Three-entry retrieval, caching, metadata/species/TA1 filters | Implemented and exercised on actual files |
| Optional bounded RCSB discovery | Implemented; pagination tested with mocked transport; not an exhaustive live search |
| Assembly operators, bound-ligand enumeration, distance-based residue selection | Implemented for selected deposited assemblies; exercised on the three default entries |
| CCD bond assignment, peptide capping, hydrogen preparation, protonation scenarios | Implemented for supported standard chemistry |
| Geometry/charge/electron checks and structural flags | Implemented as screening checks, not proof of chemical accuracy |
| Molecular output and PySCF-input dictionary adapter | Implemented |
| Classical orbitals and molecular integrals | Not implemented |
| Actual electronic active-space selection | Not implemented; only a budget/policy is returned |
| Active-space integral transformation and fermionic Hamiltonian | Not implemented |
| Qubit mapping, classical eigensolver, VQE/hardware | Not implemented |
| Quantitative binding affinity, resistance prediction, drug optimization | Not established by these outputs |

A molecular region is a set of atoms. An electronic active space is a set of electrons and orbitals. Do not rename the coordinate selection as an already selected CAS.

## 3. Files and startup

The project folder contains the notebook, `requirements.txt`, `README.md`, `validation_summary.json`, and the earlier scaffold backup. Use the current notebook, not the scaffold backup. Some earlier explanatory material may describe the former strict handoff; the current cell 21 and this document describe the updated policy.

The companion `CODE_WALKTHROUGH.md` explains every function, its rationale and limitations.

Install the notebook's dependencies into the Python environment used by Jupyter:

```bash
python -m pip install -r requirements.txt
```

Core scientific dependencies: NumPy, Gemmi and RDKit. JupyterLab/ipykernel support notebook execution. Registering a specially named kernel is optional if the correct environment is already selected. PySCF and Qiskit are not needed to run this preparation notebook; add and compatibility-test the next-stage packages in your own environment before freezing their versions.

The recorded preparation test environment was Python 3.12, NumPy 2.5.4, Gemmi 0.7.5 and RDKit 2026.09.1. The requirements file specifies version ranges, not a complete reproducibility lock.

Run the notebook in order. Cell 17 performs retrieval/preparation. It sets `PROJECT_DIR = Path.cwd()`; verify the printed directory. Cell 19 prints diagnostics and cell 21 builds `next_stage_inputs`. If `result` is already loaded and only the adapter changed, rerun cells 15 and 21.

First execution requires internet access to RCSB. Later executions reuse cached files but repeat extraction and preparation. Cache files are not automatically checked for newer revisions.

## 4. Observed baseline results

These are recorded test results, not guaranteed output for changed files/settings/software.

| Entry | Selected residues including bridges | Histidine scenarios | Prepared atoms | Net charges | Structural review |
|---|---:|---:|---|---|---|
| 1JFF | 23 | HID, HIE, HIP | 586, 586, 587 | −1, −1, 0 | No implemented flags |
| 9WDA | 27 | HID, HIE, HIP | 718, 718, 719 | +2, +2, +3 | ARG 276 and ARG 359 |
| 9WDB | 26 | HID, HIE, HIP | 682, 682, 683 | +1, +1, +2 | ARG 276 |

All nine passed the implemented chemical checks. The manifest retains three `ready_for_integrals=True` regions and two `needs_review` pockets. The updated handoff includes nine inputs, six with `requires_structural_review=True`.

### Why the human cases are flagged

The recorded validation reports contain:

| Entry/residue | Residue inclusion | Q-score |
|---|---:|---:|
| 9WDA b:276 ARG | 0.5455 | 0.609 |
| 9WDA b:359 ARG | 0.4545 | 0.607 |
| 9WDB b:276 ARG | 0.6364 | 0.612 |

The notebook's residue-inclusion cutoff is 0.7; its Q-score cutoff is 0.4. Inclusion triggers these particular warnings. Inclusion is a contour-dependent measure of atoms inside the experimental map, not a probability that the residue is correct. These are project screening thresholds, not universal rejection criteria.

Use the validation report/maps to inspect which atoms are affected and whether they form the contact of interest. The authorized handoff permits computation with these flagged geometries; it does not establish that uncertainty has disappeared. Absence of a flag for 1JFF does not prove it is the superior structure because available validation evidence differs.

## 5. Data contract

### Full result: `result` or `results/manifest.json`

Top-level fields:

- `schema_version`: currently `"1.0"`.
- `configuration`: the stored selection settings; not a complete environment/policy lock.
- `discovery`: optional search response; `None` for the normal approved-ID run.
- `accepted_entries`: metadata, coordinate hashes and validation-report locations/warnings.
- `excluded_entries`: entries failing scientific selection.
- `retrieval_errors`: download/parsing failures.
- `pockets`: successful or failed site/assembly records.
- `coverage`: counts of accepted/excluded entries, errors, pocket records, prepared scenarios, review pockets and ready regions.

Each successful pocket carries IDs, selected/bridging residue keys, contacts, beta instance, local validation metrics, excluded nearby components, warnings, files and `regions`. Assembly-failure records can omit `pocket_id` or `regions`; use `.get("regions", [])` when iterating.

Each prepared region contains:

| Field | Meaning / convention |
|---|---|
| `elements` | Element symbols in the molecular atom order |
| `coordinates_angstrom` | N×3 coordinates, **Å** |
| `bonds` | `(atom_i, atom_j, bond_order)` with **zero-based** indices |
| `formal_charges`, `net_charge` | Per-atom and total formal charge |
| `spin_2s` | Nα − Nβ; currently 0 for supported singlets |
| `electron_count` | Sum of nuclear atomic numbers minus net charge; no ECP used in preparation |
| `source_atom_mapping` | Experimental-atom identity, cap anchors, or parent index for added H |
| `caps`, `protonation`, `protonation_scenario` | Chemical preparation decisions |
| `environment` | Currently `isolated_cluster` |
| `chemistry_checks_passed` | Whether implemented chemical/geometry checks passed |
| `ready_for_integrals` | Whether chemical checks passed and no structural warnings remain |
| `warnings` | Unresolved chemical or structural messages |
| `readiness_scope` | Limits of the readiness claim |
| `orbital_selection_policy` | A downstream budget and analysis-method label |
| `orbital_indices`, `active_electrons` | **None**, because electronic selection has not been done |
| `files` | Prepared SDF and XYZ locations |

The source mapping is needed if you later choose atom-localized orbitals or revise the QM region. It is not duplicated in the compact handoff; retrieve it by region ID from the full manifest.

### Compact integration list: `next_stage_inputs`

Each item includes:

```python
{
    "pdb_id": "9WDA",
    "pocket_id": "...",
    "region_id": "...",
    "molecule_kwargs": {
        "atom": [("C", [x, y, z]), ...],
        "unit": "Angstrom",
        "charge": 2,
        "spin": 0,
        "basis": "def2-svp",
    },
    "orbital_selection_policy": {
        "max_spatial_orbitals": 4,
        "method": "downstream_contact_orbital_analysis",
    },
    "chemistry_checks_passed": True,
    "ready_for_integrals": False,
    "requires_structural_review": True,
    "warnings": ["..."],
    "protonation_scenario": {"...": "HID"},
    "readiness_scope": "...",
}
```

This is a schematic shape, not a complete executable molecule. Use actual coordinates from the generated data.

**Do not re-filter this list by `ready_for_integrals` if you intend to preserve the user's requested human inputs.** Retain review status in every subsequent output. Exclude chemistry failures, which the current notebook already does.

## 6. Loading and transferring the handoff

Cell 21 constructs the list in memory. It is not automatically written as a standalone JSON file. To export it after running the notebook:

```python
from pathlib import Path
import json

handoff_path = Path("results") / "next_stage_inputs.json"
handoff_path.write_text(json.dumps(next_stage_inputs, indent=2))
```

Share the manifest and coordinate files as well; they contain richer provenance. JSON converts tuples to lists, so consumers should accept both as element/coordinate pairs. File paths in manifests may refer to the original computer: copy the actual files and resolve paths relative to the transferred project.

To reconstruct the compact inputs from a saved manifest without rerunning preparation, load the JSON, then use the notebook's `as_pyscf_kwargs` definition:

```python
import json
from pathlib import Path

result = json.loads(Path("results/manifest.json").read_text())
inputs = []
for pocket in result["pockets"]:
    for region in pocket.get("regions", []):
        if not region.get("chemistry_checks_passed", False):
            continue
        inputs.append({
            "pdb_id": pocket["pdb_id"],
            "pocket_id": pocket["pocket_id"],
            "region_id": region["region_id"],
            "molecule_kwargs": as_pyscf_kwargs(
                region, basis="def2-svp", allow_structural_warnings=True),
            "orbital_selection_policy": region["orbital_selection_policy"],
            "chemistry_checks_passed": True,
            "ready_for_integrals": region["ready_for_integrals"],
            "requires_structural_review": not region["ready_for_integrals"],
            "warnings": list(region.get("warnings", [])),
            "protonation_scenario": region.get("protonation_scenario", {}),
            "readiness_scope": region.get("readiness_scope"),
        })
```

This code makes no scientific calculations and does not remove warnings. The adapter's default remains strict; allowing structural warnings is an explicit choice at this boundary.

## 7. Implementing the next stages

### A. Resource assessment and classical molecule

The clusters contain approximately 586–719 atoms. Do not launch nine full dense integral calculations merely because there are only nine inputs. Estimate atomic-orbital count, memory, integral-storage requirements and SCF cost for one model first. A four-spatial-orbital target does not reduce the initial full-cluster integral cost automatically.

If a smaller QM region or embedding is needed, make that a documented model revision with a new identifier and retained source mapping. Do not silently drop atoms to fit a quantum budget. The current adapter supports isolated clusters only.

Build the classical molecule with the supplied unit, charge and spin. Confirm element/coordinate counts and electron populations. Select a basis according to the scientific/resource tradeoff; `def2-svp` in the handoff is a starting specification, not an accuracy guarantee. Record method, basis, convergence, nuclear repulsion and any environment assumptions. A converged SCF solution still needs suitable stability/spin checks.

### B. Actual electronic active-space selection

Compute initial orbitals and inspect their character. Connect orbital selection to the contact atoms using the source mapping. Select relevant occupied and virtual orbitals on both sides of the interaction. Do not automatically use the frontier orbitals if they lie on caps or unrelated groups.

Save exact orbital coefficients/order, selected indices, active electron populations, inactive orbitals and selection rationale. Four spatial orbitals imply eight spin-orbitals/qubits under untapered Jordan–Wigner encoding only if that selection is actually used. No active electron count is currently prescribed.

### C. Integral transformation and Hamiltonian construction

Transform the same one-/two-electron integrals into the selected orbital basis. Keep frozen/inactive-space contributions and constant energies explicitly. Record tensor conventions, units, orbital ordering and the chosen spin sector. Preserve nuclear repulsion and inactive offsets without omitting or counting them twice.

The appropriate comparison is between implementations of the **same Hamiltonian**, not energies obtained from separately optimized or differently ordered orbital spaces.

### D. Classical validation, mapping and solvers

Before VQE, validate a small selected problem against a classical reference such as CASCI/exact diagonalization in the intended particle sector. Check Hermiticity, particle number, energy offsets and mapping size. Only then attach a variational circuit or hardware backend. Propagate structural-review metadata into solver results so those results cannot be mistaken for structurally unqualified predictions.

### E. Scientific interpretation

A complex's total electronic energy alone is not binding affinity. Different atom counts, net charges or protonation scenarios cannot be ranked by raw total energy. Binding calculations require balanced references and appropriate treatment of geometry, environment and thermodynamics. The current workflow does not establish quantum advantage, cancer resistance or an optimized drug.

## 8. Acceptance checks for an integration

- Original preparation results are retained unchanged; downstream results have their own location and identifiers.
- The tested handoff has nine scenarios, of which six carry structural-review flags.
- Every input has matching element/coordinate counts, explicit Å units and preserved charge/spin.
- Non-ready human inputs remain marked for review; chemistry failures cannot enter through the override.
- Protonation alternatives are distinct; the three alternatives are not treated as independent experimental pockets.
- SCF/orbital failures produce explicit results/errors rather than fabricated energies.
- Orbital indices and active electron counts are recorded only after actual selection.
- Hamiltonian constants and mapping conventions are reproducible and tested against a classical reference.
- New output records retain PDB/pocket/region IDs and warnings.
- Resource estimates precede large-scale runs.

## 9. Known limitations and engineering follow-up

**Important scientific limits:** isolated finite assemblies; no reconstructed complete lattice; no solvent/embedding; heavy atoms fixed; pH-7.4 template assumptions rather than pKa prediction; no repair of missing heavy atoms; coarse bond/clash screening; closed-shell spin assumption; incomplete local validation can remain unknown without blocking readiness.

**Identity/assembly limits:** verification uses deposited descriptions and species, not independent sequence alignment. The output stores raw mutation/sequence metadata. Pocket identity uses TA1 proximity to beta-tubulin rather than an independently aligned taxane-site motif. Broader entries with complex alternate conformers, insertions or modified residues need additional tests.

**Software limits:** the notebook is self-contained, not a separately importable library. A useful first refactor is to move definition cells into a module while preserving behavior and keeping the notebook as the demonstration. Add fixture-based tests for broader assembly operators, unusual chemistry and API schema changes. Expand configuration/version logging and structured warning types before a schema revision.

**Filesystem limits:** relative output paths depend on the kernel working directory. Reruns overwrite matching outputs, but stale files from older settings may remain. The current manifest, not a directory glob, is the authoritative index. Cache refresh/revision checking is not automatic. Some export failures can stop execution rather than returning a per-site status.

**Documentation precedence:** the current code and these two documents supersede the original bovine-first blueprint and any text saying human cases are automatically excluded from the handoff. They are included with warnings by explicit user request.

## 10. Reference links

- [1JFF structure record](https://www.rcsb.org/structure/1JFF)
- [9WDA structure record](https://www.rcsb.org/structure/9WDA)
- [9WDB structure record](https://www.rcsb.org/structure/9WDB)
- [Human tubulin/paclitaxel study](https://www.nature.com/articles/s41589-026-02204-2)
- [wwPDB EM validation guide](https://www.wwpdb.org/validation/2017/EMValidationReportHelp)
- [RCSB Search API](https://search.rcsb.org/)
- [PySCF CASCI/CASSCF documentation](https://pyscf.org/user/mcscf.html)
- [Qiskit Nature problem transformers](https://qiskit-community.github.io/qiskit-nature/tutorials/05_problem_transformers.html)

These links identify the structural and software context; consult the documentation corresponding to the versions used when implementing the next module.
