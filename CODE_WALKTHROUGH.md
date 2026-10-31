# Taxol–tubulin notebook: function-by-function code walkthrough

Updated 31 October 2026. This document describes the current `human_taxane_pocket_interface.ipynb`, including the change that includes chemically prepared human regions in the handoff while preserving structural warnings. It documents implemented behavior, not the earlier proposed design.

## 1. What the code produces

The notebook converts approved experimental structures into isolated molecular clusters containing complete Taxol and nearby capped peptide segments. It exports atom coordinates, bonds, charges, spin, preparation provenance and quality flags. These are inputs for later electronic-structure work.

A **structural pocket** is a collection of residues around a ligand. A **prepared molecular region** is a chemically defined molecule/cluster constructed from those coordinates. An **electronic active space** is a selection of molecular orbitals and electrons. The notebook builds the first two; it only records a policy for the third.

Default entries are 1JFF (the sole bovine exception), 9WDA (human reference) and 9WDB (human S239C comparator). The normal run does not call the database-search function.

## 2. Execution and dependencies

Numbering counts every notebook cell, including Markdown. Cells 5, 7, 9, 11, 13 and 15 load dependencies or define functions. Cell 17 calls the pipeline; cell 19 prints diagnostics; cell 21 constructs `next_stage_inputs`. Defining a function does not execute its body.

| Dependency | Reason it is used |
|---|---|
| `numpy` | Arrays, coordinate transformations, distances and numerical checks |
| `gemmi` | Reading mmCIF categories and writing pocket coordinate files |
| `rdkit.Chem` | Molecular graphs, bond orders, formal charges, stereochemistry and molecular files |
| `rdkit.Chem.AllChem` | UFF parameter checks and hydrogen-position minimization |
| `urllib` | HTTP downloads and the optional RCSB search request |
| `json`, `csv`, `gzip`, `xml.etree.ElementTree` | Output serialization and compressed validation-report parsing |
| `pathlib`, `hashlib` | File locations and source-coordinate hashes |
| `defaultdict`, `product`, `re` | Grouping atoms, enumerating combinations and parsing identifiers/operators |
| `time` | Brief retry delays after temporary download failures |

`APPROVED_IDS` defines the default three entries. `STANDARD_AA` is a set of the 20 standard amino-acid component names. Set membership is used to separate ordinary protein residues from ligands, waters, ions and unsupported modifications.

## 3. Download and entry verification — cell 7

### `download(url, destination, refresh=False)`

**What:** obtains a remote file and returns its local `Path`.

**How:** if the destination already exists and refresh is false, it returns immediately. Otherwise it creates the parent directory, requests the URL with a 60-second timeout, checks that the response is nonempty, writes a temporary file, and replaces the destination only after writing completes. It makes up to three attempts for network/time-out failures and selected HTTP failures. Most HTTP errors below 500 are immediately raised; HTTP 429 can be retried.

**Why:** structure and chemical-definition files are reused across scenarios. Caching reduces requests and allows later runs to reuse downloaded evidence. The temporary-file replacement reduces the chance of using a partially written download.

**Side effects:** network access and cache-file creation/replacement. The caller handles most failures.

**Limits:** an existing file is trusted without checking its revision or checksum. A successful nonempty HTTP response is not itself proof that the contents are a valid CIF; parsing occurs later. Refresh is available on this helper, not as a public pipeline parameter. Simultaneous runs sharing the same cache are not coordinated.

### `search_human_taxane_ids(max_resolution=3.0, page_size=100, max_results=300)`

**What:** optionally discovers candidate entry IDs through the RCSB Search API.

**How:** creates a query requiring human source taxonomy 9606, Taxol component TA1, an eligible resolution and a tubulin description. Requests experimental entries in pages until all reported matches are fetched or the configured cap is reached. Returns `ids`, `total_matches` and `truncated`. An HTTP 204 is treated as an empty result.

**Why:** provides an extension beyond the approved list without making every default run a discovery run.

**Important:** entry-level search conditions may match different entities within the same entry. Results therefore still pass `inspect_entry()`. The search is not proof that the Taxol is bound to the relevant human tubulin chain.

**Side effects/limits:** makes POST requests; errors propagate. The default notebook never calls it. Its pagination/result logic was tested with mocked transport; the full discovery path has not been established by an exhaustive live database run.

#### Nested helper `term(attribute, operator, value)`

Returns the RCSB JSON representation of one search condition. This avoids repeating the same nested dictionary structure for every filter. It makes no request itself.

### `rows(block, category)`

**What:** converts a Gemmi CIF category into a list of row dictionaries.

**How:** Gemmi supplies a dictionary whose values are column lists. The function transposes the lists and pairs values with their column names. An absent category becomes an empty list.

**Why:** the rest of the code can read `row["entity_id"]` instead of manually managing parallel columns. This is an in-memory conversion with no file changes.

### `number(value, default=None)`

**What:** safely converts a value to `float`.

**How:** returns the default for `None`, booleans, or values that cannot be converted.

**Why:** CIF missing-value representations must not become plausible measurements. In particular, `float(False)` would be zero in Python, which could incorrectly look like an excellent resolution. This helper avoids that conversion. It does not, on its own, reject every nonfinite float; callers provide additional checks where implemented.

### `inspect_entry(block, pdb_id, max_resolution)`

**What:** accepts or rejects the entry based on metadata and returns metadata for accepted entries.

**How:** reads entity descriptions and source-organism tables, identifies descriptions containing both “beta” and “tubulin,” and checks the species of every entity whose description contains “tubulin.” Human taxonomy 9606 is required, except that 1JFF requires bovine taxonomy 9913. It reads experimental method and positive resolution values, applies the threshold, and requires TA1 in the nonpolymer table. 1JFF receives an explicit 3.5 Å limit; other entries use the configured human cutoff.

**Returns:** beta entity IDs, entity descriptions, species, resolution, method, exception status, raw entity metadata under `mutations`, and sequence-reference metadata.

**Why:** prevents a human expression partner or an unrelated human protein elsewhere in a complex from satisfying the tubulin-species requirement.

**Limits:** protein identity is checked from deposited descriptions/source metadata, not an independent sequence alignment. Sequence references are retained but not used to establish an isotype or mutation by alignment. The `mutations` field contains raw entity rows, not a normalized mutation list. The function uses the minimum positive resolution across supported categories; mixed-method entries beyond the tested set need scrutiny. Presence of TA1 is an entry-level test; proximity to beta-tubulin is checked later.

## 4. Assembly handling and contacts — cell 9

### `operation_sequences(expression)`

**What:** expands a deposited assembly-operation expression into operator-ID tuples.

**Example:** `(1-2)(3,4)` becomes `[("1", "3"), ("1", "4"), ("2", "3"), ("2", "4")]`.

**How:** extracts groups, expands numerical ranges and uses `itertools.product` to enumerate combinations.

**Why:** some biological assemblies require combinations of rotations/translations rather than one coordinate copy. The returned tuples tell `assembly_atoms()` which matrices to compose. This is a small parser for the supported deposited syntax, not a general mathematical-expression evaluator.

### `assembly_atoms(block, assembly_id="1", model_id="1", max_atoms=300000)`

**What:** returns atom dictionaries for one model and one deposited finite assembly.

**How:** reads `_atom_site`, filters the model, omits H/D atoms, and records element, coordinates, names, source IDs, chain/residue identifiers, occupancy and B-factor. It groups atoms by residue and chooses one alternate-location label using summed occupancy, retaining shared blank-location atoms. It builds homogeneous 4×4 transformation matrices from deposited assembly operators, applies the requested operation sequences to the relevant chains, and removes duplicate transformed source atoms. A maximum-atom check prevents accidental huge expansion.

**Why:** pocket contacts may involve neighboring chains. Performing selection before applying the intended assembly could omit those neighbors or confuse repeated copies.

**Identifier example:** `B:1:227:None:HIS` encodes label chain B, operation 1, sequence position 227, no insertion code, and histidine. The literal `None` here is a missing insertion-code value in the generated key; it is not a missing residue.

**Limits:** the default is assembly 1, not every possible assembly. The code does not reconstruct an infinite helical lattice or predict missing neighboring subunits. Residue-level alternate selection does not guarantee a globally coherent multi-residue conformer. Source/model identifiers are available in metadata; simplified exported pocket coordinates use their own presentation identifiers.

### `extract_pockets(atoms, beta_entities, radius)`

**What:** creates a pocket record for every bound TA1 instance in the supplied assembly.

**How:** groups atoms into residues, identifies TA1 groups, and computes a residue-versus-ligand heavy-atom distance matrix. It records the minimum distance and nearest atom pair for nearby residues/components. Complete standard residues within the cutoff are retained. A TA1 instance with no sufficiently close verified beta-tubulin entity is rejected. A single unselected residue between two selected sequential residues is added as a bridge, preventing duplicate/overlapping cap anchors.

**Returns:** selected residue keys, bridging residues, contact rows, nearest beta instance and temporary ligand/residue lookups used for preparation.

**Why:** this produces an auditable geometric selection while preserving full residue chemistry. The bridge rule favors defensible peptide boundaries over strict adherence to a sphere.

**Limits:** one nearest atom pair is stored per residue, not all atom contacts. The flags `within_4A`/`within_5A` describe contacts actually collected at the requested radius. Distances do not establish hydrogen bonds. Beta-tubulin proximity is not an independently aligned taxane-site motif check. Unoccupied pockets, docking and non-TA1 taxanes are outside scope. The pairwise distance approach is suitable for these bounded inputs, but a spatial index would improve larger searches.

## 5. Chemical definitions and protonation — cell 11

### `component_definition(component, cache_dir)`

**What:** downloads/caches the component CIF and returns its CCD atom and bond tables.

**Why:** coordinate files do not necessarily provide enough information to recover chemical bond orders, formal charges and stereochemistry. The chemical dictionary supplies the expected topology without replacing the experimental bound coordinates.

### `add_component(builder, coordinates, component, cache_dir, mapping, positions)`

**What:** inserts a complete heavy-atom component into an RDKit editable molecule and returns a name-to-index mapping.

**How:** compares observed heavy-atom names with CCD expectations, excluding OXT for standard amino acids. Duplicate, missing or unexpected names cause a failure. It adds atoms with formal charges, makes backbone nitrogen neutral for the subsequent peptide assembly, and adds single/double/triple/aromatic bonds from CCD definitions. It appends experimental coordinates and source-atom provenance. For TA1 it also records expected R/S labels for a later check.

**Why:** explicit topology avoids guessing chemical bonds from atom proximity. OXT handling reflects the intended capped peptide model rather than a free amino-acid terminus.

**Side effects:** mutates `builder`, `mapping` and `positions`. It does not add hydrogens or caps itself.

**Limits:** missing heavy atoms are rejected, not rebuilt. Terminal residues, modified amino acids and unfamiliar chemistry may require different templates. This is intentionally narrower than a general protein-repair tool.

### `apply_protonation(builder, indices, residue_name, choice=None)`

**What:** modifies charges and selected bond/hydrogen representations for an explicit microstate, returning its name.

**Implemented choices:** ASP/GLU deprotonated or neutral; LYS protonated or neutral; ARG protonated only; HIS HID/HIE/HIP; CYS neutral or deprotonated. Other supported standard residues use the neutral choice. Acid carbonyl/carboxylate bonds are made explicit, and arginine's guanidinium charge is assigned consistently with its double bond.

**Why:** electron count and electrostatics depend on protonation. Histidine can differ in both proton location and charge; treating all scenarios as one molecule would hide that uncertainty.

**Limits:** this is a template policy, not a pKa prediction. The helper requires explicit choices for HIS/CYS; the public pipeline provides histidine scenarios and a neutral CYS default. Other pH values are rejected by the public function. It does not detect or implement arbitrary post-translational modifications or general disulfide chemistry.

## 6. Prepared region construction — cell 13

### `prepare_region(pocket, cache_dir, protonation_overrides, max_atoms=1800)`

**What:** creates one chemically prepared region for one protonation scenario. Returns `(rdkit_molecule, region_dictionary)`.

**How, in execution order:**

1. Add full TA1 and all selected standard residues from CCD definitions.
2. Apply protonation choices.
3. Find adjacent residues using chain, operation and sequential label numbering.
4. Verify peptide-boundary C–N distances (between 1.1 and 1.7 Å).
5. Connect retained consecutive residues; cap unconnected boundaries with ACE or NME.
6. Sanitize the molecular graph and attach the experimental heavy-atom coordinates.
7. Assign stereochemistry from 3D coordinates and compare Taxol R/S centers against CCD expectations.
8. Add hydrogens with coordinates; enforce the prepared-atom budget.
9. Verify UFF parameter availability. Optimize hydrogens for up to 2,000 minimization iterations while keeping all heavy atoms fixed; interactions between disconnected components are included.
10. Record provenance for added hydrogens, calculate charge and electron number, and reject radicals or odd electron counts in this closed-shell workflow.
11. Check finite coordinates, severe nonbonded clashes below 1 Å, and bond distances outside 0.65–1.35 times summed covalent radii.
12. Return geometry, topology, charges, spin 0, source mappings, cap descriptions, microstates, warnings and readiness fields.

**Why caps are needed:** removing the protein beyond a selected peptide segment leaves cut covalent bonds. ACE supplies an acetyl-like N-terminal cap; NME supplies an N-methylamide-like C-terminal cap. Their heavy-atom positions are anchored to the excluded neighboring backbone.

**Why only hydrogens move:** unconstrained isolated-cluster optimization could move the model away from the experimental bound arrangement. Fixing heavy atoms preserves that arrangement but also preserves any underlying experimental heavy-atom errors.

**Limits:** closed-shell singlet treatment is assumed, not energetically established against other spin states. UFF convergence is a preparation check, not evidence of quantum accuracy. The bond/clash rules are coarse screens. Entire residues and Taxol remain, so these regions have hundreds of atoms. The returned `orbital_indices` and `active_electrons` are deliberately `None`.

#### Nested `neighbor(key, offset)`

Looks up the preceding/following standard residue using sequence number ±1 within the same transformed chain. This avoids renumbering the extracted fragments. It assumes meaningful integer label sequence IDs.

#### Nested `atom_record(key, name)`

Retrieves exactly one named source atom. Missing neighbors or missing/duplicate cap anchors raise an error rather than generating an arbitrary geometry.

#### Nested `link(left_key, right_key)`

Checks the source C–N distance at a peptide boundary. It validates geometry; the parent function separately adds the corresponding bond.

#### Nested `cap_atom(key, name, element, cap_kind)`

Adds a cap atom at an excluded neighbor's backbone position, recording the source anchor. The atom belongs to a new cap model, so its provenance is distinguished from an unchanged experimental atom.

### `save_pocket_files(atoms, contacts, folder)`

**What:** writes `pocket.cif` and `contacts.csv`.

**Why:** these are human-inspectable structural selections even if chemical preparation later fails. The CIF contains selected heavy atoms before hydrogen/cap preparation; the SDF/XYZ files created later contain prepared regions.

**Limits:** the export uses generated chain identifiers and model 1 for the standalone selection. Detailed original identity should be recovered from JSON provenance and the original mmCIF, not inferred solely from this presentation file. The function does not create the folder; its caller does.

## 7. Pipeline orchestration and adapter — cell 15

### `json_safe(value)`

Converts NumPy arrays to lists, NumPy scalars to Python scalars, and Paths to strings for JSON serialization. Unknown types raise `TypeError`. This preserves numerical data without converting the entire object indiscriminately to text.

### `prepare_human_taxane_regions(...)`

**What:** public entry point; returns the complete result dictionary and writes a manifest.

| Parameter | Default | Meaning |
|---|---|---|
| `pdb_ids` | `APPROVED_IDS` | Explicit list; `None` enables bounded discovery |
| `max_resolution_angstrom` | 3.0 | Human cutoff; 1JFF has its separate exception |
| `contact_radius_angstrom` | 5.0 | Local residue-selection radius; accepted range (0, 8] |
| `assembly_ids` | `("1",)` | Deposited assemblies to process |
| `protonation_overrides` | `None` | Shared residue choices or per-entry choices |
| `ph` | 7.4 | Declared template assumption; other values rejected |
| `max_spatial_orbitals` | 4 | Stored downstream budget, not executed orbital selection |
| `cache_dir` | `data/structures` | Downloads and CCD definitions |
| `output_dir` | `results` | Structural/preparation outputs |
| `max_region_atoms` | 1800 | Cap on the prepared molecule's atom count |

**Execution:** validate settings → obtain IDs → download and verify each entry → read validation XML if available → process all deposited models and requested assemblies → enumerate pockets → save selections → evaluate structural screening flags → enumerate protonation scenarios → prepare each scenario → save SDF/XYZ and JSON → print progress → write `manifest.json` and return the result.

**Structural screening:** occupancy below 0.9, alternate conformations and certain unmodeled nearby components are flagged. Reported local Q-score below 0.4, RSCC below 0.7 or residue inclusion below 0.7 also triggers a flag. These thresholds are project heuristics. Missing validation reports/metrics are recorded but do not automatically block readiness.

**Scenario enumeration:** each unresolved HIS contributes HID/HIE/HIP alternatives. At most two unresolved histidines are supported automatically, giving at most nine scenarios. CYS defaults to neutral unless overridden. Example per-entry override: `{"9WDA": {"B:1:227:None:HIS": "HIE"}}`.

**Status distinction:** preparation sets `chemistry_checks_passed`. Structural warnings are appended afterward, and `ready_for_integrals` becomes false if any warning remains. A pocket is `ready` if any region is ready, otherwise `needs_review`. Failed regions may have reasons but no region object; some assembly-failure records have no pocket ID or regions list.

**Errors:** retrieval problems, rejected metadata, assembly failures and scenario-preparation failures are generally retained separately so other entries/sites can proceed. Not every operation is enclosed in a recovery block: invalid arguments, discovery errors or certain filesystem/export failures can still stop the run.

**Limits:** finite isolated assemblies only; waters/ions excluded from the prepared cluster; no missing-heavy-atom reconstruction. The manifest stores useful configuration and the coordinate hash, but not a complete environment lock, all policy settings, or every cache revision. Broader reproducibility logging is an extension task.

### `as_pyscf_kwargs(region, basis="def2-svp", *, allow_structural_warnings=False)`

**What:** returns a PySCF-compatible molecule dictionary containing `atom`, `unit`, `charge`, `spin` and `basis`.

**Checks:** chemistry must have passed; non-ready regions require explicit `allow_structural_warnings=True`; the environment must be `isolated_cluster`.

**Why:** the adapter is the integration boundary. It retains a strict default while allowing the user's requested inclusion of structurally flagged human models. It never bypasses a chemistry failure.

**Limits:** it trusts the preparation fields; it is not a second independent chemical validator. It does not import PySCF, compute integrals, select orbitals or calculate an energy. Specifying `def2-svp` only names a basis for the next module.

## 8. Code outside functions

### Cell 17: run configuration

`PROJECT_DIR = Path.cwd()` chooses the kernel's current directory, which is not guaranteed to be the notebook's directory. The printed path is the location used for downloads/results. `PROTONATION_OVERRIDES = {}` requests automatic histidine scenarios. Calling the public function triggers the actual workflow and assigns its return to `result`.

### Cell 19: diagnostics

Loops through `result`, printing pocket statuses, scenario atom counts/charges, warnings, retrieval errors and scientific exclusions. This is live output. Markdown descriptions of a previous test run are static text.

### Cell 21: current handoff behavior

Builds `next_stage_inputs` from every region with `chemistry_checks_passed=True`, calling the adapter with `allow_structural_warnings=True`. Each item retains `pdb_id`, `pocket_id`, `region_id`, `molecule_kwargs`, `orbital_selection_policy`, both readiness fields, `requires_structural_review`, warnings, protonation scenario and readiness scope.

The user's requested inclusion does **not** change `ready_for_integrals` in the original records. Expected tested counts remain **three ready regions in the manifest**, but **nine included handoff inputs**, six carrying review flags. `next_stage_inputs` exists in notebook memory; the notebook does not automatically export it as a separate JSON file.

## 9. Practical interpretation

The tested configurations produced one bound pocket per structure: 23 retained residues for 1JFF, 27 for 9WDA and 26 for 9WDB, including bridges. Three histidine scenarios per pocket yielded nine molecular models. These are not nine independent binding pockets.

A passed preparation check does not establish biological accuracy. In particular, the fixed heavy-atom approximation, isolated environment, large model size, unresolved local-density uncertainty and unselected electronic orbitals remain relevant to the next user's calculations. The accompanying `INTEGRATION_HANDOFF.md` describes how to carry those assumptions forward.
