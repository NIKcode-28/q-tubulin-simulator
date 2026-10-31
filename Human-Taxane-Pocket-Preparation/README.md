# Human Taxane Pocket Preparation

Open `human_taxane_pocket_interface.ipynb` and follow its installation section.
The notebook is self-contained; no separate Python module is required.

Install with `python -m pip install -r requirements.txt`, select that environment
as your notebook kernel, and run cells in order. The run cell downloads only
1JFF, 9WDA and 9WDB by default. An optional bounded database-search helper is
provided but is not automatically invoked.

Outputs are saved relative to the kernel's current directory in `data/structures/`
and `results/`. The notebook prints this directory. A run updates the manifest;
use that manifest as the authoritative index rather than leftover files from
older configurations.

## Verified behavior

Tested with real PDB inputs: three accepted entries, three bound pockets, and
nine prepared histidine scenarios. All nine passed implemented chemistry checks.
Human scenarios were flagged by local density screening; only the three 1JFF
scenarios passed all current handoff checks. See `validation_summary.json` for
exact results. Do not clear warnings simply to make a model pass.

The prepared regions are isolated clusters containing complete Taxol and capped
local peptide segments, not final electronic active spaces. They contain hundreds
of atoms. An orbital budget does not remove the cost of their classical integral
calculation or guarantee that the pocket fits accurately in a small active space.

The notebook exports coordinates, bonds, charge, spin, source mappings and
selection policies. It does not run PySCF or Qiskit. For new chemistry, unsupported
protonation or missing atoms it returns review reasons instead of fabricating
ready inputs. Extended microtubule lattice reconstruction is outside scope.
