# PyRosetta Notebook Experiments

Exploratory notebooks for teaching and prototyping PyRosetta workflows: membrane
protein modeling, enzyme design, symmetry setup, active-site analysis, structural
alignment, and early ideas for a full online book chapter.

This repository is best read as a curated lab notebook. Some examples are polished
chapter prototypes; others preserve useful experiments from older PyRosetta
bindings so the reasoning and implementation details remain available.

## Featured Work

### Membrane Proteins: Moving Metals

![Membrane protein metal-sampling movie](resources/images/mp.gif)

The membrane-protein notebook demonstrates a visual workflow for opening a
membrane protein in PyRosetta, introducing a metal, and sampling nearby residues
with repacking over a sequence of steps.

- Notebook: [Membrane 1 - Moving Metals.ipynb](membrane_proteins_moving_metals/Membrane%201%20-%20Moving%20Metals.ipynb)
- Supporting files: [1bl8.span](membrane_proteins_moving_metals/1bl8.span),
  [1bl8.nometal.pdb](membrane_proteins_moving_metals/1bl8.nometal.pdb)

### Prototype Book Chapter: Writing Protein Design Algorithms

| Histidine tautomer example | Alternate tautomer example |
| --- | --- |
| ![Histidine structure](prototype/his.png) | ![Alternate histidine structure](prototype/his_d.png) |

The prototype chapter explores how to break down Rosetta packing behavior into
readable data tables, then uses that visibility to build and debug a simple
protein-design algorithm. It is the clearest example of the larger book format
this repository was originally exploring.

- Notebook: [Prototype of PyRosetta Jupyter Notebook.ipynb](prototype/Prototype%20of%20PyRosetta%20Jupyter%20Notebook.ipynb)
- Section notes: [prototype/README.md](prototype/README.md)

## Notebook Catalog

| Section | Status | What it covers |
| --- | --- | --- |
| [membrane_proteins_moving_metals](membrane_proteins_moving_metals/) | Featured, PyRosetta 4 | Membrane setup, metal placement, residue repacking, and animated sampling output. |
| [prototype](prototype/) | Featured prototype | A book-chapter style walkthrough for extracting energies into dataframes and writing a simple design algorithm. |
| [enzymedesign_in_pyrosetta](enzymedesign_in_pyrosetta/) | Working example | Enzyme-design setup in PyRosetta using Bagel/Foldit supporting files. |
| [activesiteenergycalc](activesiteenergycalc/) | Archival example | Active-site RMSD and energy calculations against a reference structure. |
| [poses_and_dataframes](poses_and_dataframes/) | Archival example, older bindings | Accessing pose energies with pandas and comparing Rosetta/PyRosetta scoring. |
| [symmetry](symmetry/) | Archival example, older bindings | Building symmetry setup logic directly in PyRosetta. |
| [getting_tmalign_to_work_w_ligands](getting_tmalign_to_work_w_ligands/) | Archival example, older bindings | Using TMalign-style superposition for proteins with different sequences. |

## Repository Layout

```text
.
|-- README.md
|-- resources/
|   |-- images/                 # Shared README and notebook visuals
|   `-- protein_structures/      # Shared structure inputs
|-- membrane_proteins_moving_metals/
|-- prototype/
|-- enzymedesign_in_pyrosetta/
|-- activesiteenergycalc/
|-- poses_and_dataframes/
|-- symmetry/
|-- getting_tmalign_to_work_w_ligands/
`-- sjb_util.py                  # Helper utilities used by some legacy demos
```

## Running The Notebooks

These notebooks were written across multiple PyRosetta eras. The featured
membrane and prototype notebooks expect PyRosetta 4-era APIs, while several
archival examples use older bindings and may need small updates before running
on a modern installation.

General setup:

1. Install PyRosetta using the license and platform-specific instructions from
   RosettaCommons.
2. Create a Python environment with Jupyter, pandas, matplotlib, and seaborn.
3. Launch Jupyter from the repository root so relative paths to PDB, params,
   span, image, and XML files resolve correctly.

## Relationship To The Official PyRosetta Notebooks

A later RosettaCommons effort produced the broad PyRosetta notebook collection
this project originally hoped would exist. This repository now serves as a
smaller companion archive of experiments, chapter prototypes, and teaching ideas.

- Official notebooks: <https://github.com/RosettaCommons/PyRosetta.notebooks>
- Associated preprint: <https://www.preprints.org/manuscript/202002.0097/v1>

## Notes

- PyRosetta and Rosetta are distributed under RosettaCommons licensing terms;
  this repository does not include PyRosetta itself.
- `sjb_util.py` is retained for examples that depend on it, but it has not been
  cleaned up as a public API.
- Data files are included only where needed to make individual notebooks easier
  to inspect and reproduce.
