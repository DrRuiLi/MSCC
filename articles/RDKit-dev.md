# RDKit-dev

``` r

library(MSCC)
library(reticulate)

# Default conda env for RDKit helpers (before any other reticulate Python init):
ensure_RDKit_python("env_for_r")
# or: options(MSCC.rdkit.condaenv = "env_for_r")
```

## 1. Introduction

MSCC calls RDKit through **reticulate**. This vignette pairs the native
Python API with the MSCC R wrappers in `dev_rdkit.R`.

Requires a Python environment with `rdkit` (and `numpy` / Cairo-enabled
Draw for plotting). Default helpers target the conda env `env_for_r`. On
Windows,
[`ensure_RDKit_python()`](https://drruili.github.io/MSCC/reference/ensure_RDKit_python.md)
adds `Library/bin` to the DLL search path so `rdMolDraw2D` can load;
prefer conda-forge builds of `rdkit` / `cairo` / `harfbuzz`.

Default demo molecule is glycine: `NCC(O)=O`.

## 2. Main functions

### 2.1 Create a molecule

Python:

``` python
from rdkit import Chem

mol = Chem.MolFromSmiles("NCC(O)=O")
mol
```

R (MSCC wrapper):

``` r

Chem <- get_RDKit_Chem()
mol <- rdkit_mol_from_smiles("NCC(O)=O")
mol
```

### 2.2 Draw a molecule

Python:

``` python
from rdkit import Chem
from rdkit.Chem import Draw

mol = Chem.MolFromSmiles("NCC(O)=O")
img = Draw.MolToImage(mol, size=(300, 300))
img  # PIL Image
```

R (MSCC wrapper — in-memory `Draw.MolToImage` → ggplot):

``` r

mol <- rdkit_mol_from_smiles("NCC(O)=O")
rdkit_plot_mol(mol)
```

### 2.3 Molecule properties (formula and exact mass)

Python:

``` python
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

mol = Chem.MolFromSmiles("NCC(O)=O")
rdMolDescriptors.CalcMolFormula(mol)
Descriptors.ExactMolWt(mol)
```

R (MSCC wrappers):

``` r

mol <- rdkit_mol_from_smiles("NCC(O)=O")
rdkit_mol_formula(mol)
rdkit_mol_exact_mass(mol)
```

### 2.4 Substructure match (SMARTS)

Python:

``` python
from rdkit import Chem

mol = Chem.MolFromSmiles("NCC(O)=O")
pat = Chem.MolFromSmarts("[NH2]")
mol.HasSubstructMatch(pat)
```

R (MSCC wrapper):

``` r

mol <- rdkit_mol_from_smiles("NCC(O)=O")
rdkit_has_substruct(mol, "[NH2]")
rdkit_has_substruct(rdkit_mol_from_smiles("CCO"), "[OH]")
```
