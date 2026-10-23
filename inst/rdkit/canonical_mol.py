"""Canonical molecules for MSCC, using the same steps as ICEBERG canonical_root.

This file is executed by the ICEBERG RDKit interpreter. It does not import
ms-pred or torch. Atom indices are the parse order of MolFromSmiles on the
canonical SMILES. Hydrogens are not written into the graph.
"""

import json
import sys
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdmolops import RemoveStereochemistry

RDLogger.DisableLog("rdApp.*")


def _tautomer_enumerator():
    try:
        from rdkit.Chem.MolStandardize.tautomer import (
            TautomerCanonicalizer,
            TautomerTransform,
        )

        transforms = (
            TautomerTransform(
                "1,3 heteroatom H shift",
                "[#7,S,O,Se,Te;!H0]-[#7X2,#6,#15]=[#7,#16,#8,Se,Te]",
            ),
            TautomerTransform(
                "1,3 (thio)keto/enol r",
                "[O,S,Se,Te;X2!H0]-[C]=[C]",
            ),
        )
        return TautomerCanonicalizer(transforms=transforms).canonicalize
    except Exception:
        from rdkit.Chem.MolStandardize.rdMolStandardize import (
            CleanupParameters,
            TautomerEnumerator,
        )

        params = CleanupParameters()
        params.tautomerTransformsFile = str(
            Path(__file__).with_name("tautomer_transforms.txt")
        )
        return TautomerEnumerator(params).Canonicalize


_CANONICALIZE_TAUTOMER = _tautomer_enumerator()


def _clear_maps(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)


def rm_stereo(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    _clear_maps(mol)
    RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)


def smiles_from_inchi(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    inchi = Chem.MolToInchi(mol)
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        return ""
    return Chem.MolToSmiles(_CANONICALIZE_TAUTOMER(mol))


def canonical_root(smiles):
    return smiles_from_inchi(rm_stereo(smiles))


def _mol(smiles):
    canon = canonical_root(smiles)
    if not canon:
        raise ValueError("RDKit could not canonicalize this SMILES.")
    mol = Chem.MolFromSmiles(canon)
    if mol is None:
        raise ValueError("RDKit could not parse the canonical SMILES.")
    return canon, mol


def graph(smiles):
    canon, mol = _mol(smiles)
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    atoms = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        point = conf.GetAtomPosition(idx)
        atoms.append(
            {
                "id": "%s_%d" % (symbol, idx + 1),
                "element": symbol,
                "index": idx,
                "x": float(point.x),
                "y": float(point.y),
                "isotope": int(atom.GetIsotope()),
                "charge": int(atom.GetFormalCharge()),
            }
        )
    bonds = []
    for bond in mol.GetBonds():
        bonds.append(
            {
                "from": atoms[bond.GetBeginAtomIdx()]["id"],
                "to": atoms[bond.GetEndAtomIdx()]["id"],
                "bond_type": int(bond.GetBondTypeAsDouble()),
            }
        )
    return {
        "smiles": canon,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "molblock": Chem.MolToMolBlock(mol),
        "atoms": atoms,
        "bonds": bonds,
    }


def mcs(left, right):
    left_canon, left_mol = _mol(left)
    right_canon, right_mol = _mol(right)
    result = rdFMCS.FindMCS([left_mol, right_mol])
    pairs = []
    if result.numAtoms:
        pattern = Chem.MolFromSmarts(result.smartsString)
        left_match = left_mol.GetSubstructMatch(pattern)
        right_match = right_mol.GetSubstructMatch(pattern)
        for src, dst in zip(left_match, right_match):
            pairs.append(
                {
                    "from": "%s_%d" % (left_mol.GetAtomWithIdx(src).GetSymbol(), src + 1),
                    "to": "%s_%d" % (right_mol.GetAtomWithIdx(dst).GetSymbol(), dst + 1),
                }
            )
    return {"from_smiles": left_canon, "to_smiles": right_canon, "pairs": pairs}


def _mapped_mol(mapped_smiles):
    mol = Chem.MolFromSmiles(mapped_smiles)
    if mol is None:
        raise ValueError("RDKit could not parse a mapped SMILES.")
    numbers = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol, numbers


def rxn_ids(mapped_smiles):
    """Align RXNMapper atom-map numbers onto canonical atom ids."""
    mapped, numbers = _mapped_mol(mapped_smiles)
    plain = Chem.MolToSmiles(mapped)
    canon, canon_mol = _mol(plain)
    match = canon_mol.GetSubstructMatch(mapped)
    pairs = []
    if len(match) == mapped.GetNumAtoms():
        pairs = list(enumerate(match))
    else:
        result = rdFMCS.FindMCS([mapped, canon_mol])
        if result.numAtoms:
            pattern = Chem.MolFromSmarts(result.smartsString)
            left = mapped.GetSubstructMatch(pattern)
            right = canon_mol.GetSubstructMatch(pattern)
            pairs = list(zip(left, right))
    out = {}
    for src, dst in pairs:
        atom = canon_mol.GetAtomWithIdx(dst)
        out["%s_%d" % (atom.GetSymbol(), dst + 1)] = int(numbers[src])
    for atom in canon_mol.GetAtoms():
        key = "%s_%d" % (atom.GetSymbol(), atom.GetIdx() + 1)
        out.setdefault(key, None)
    return {"smiles": canon, "rxn_id": out}


def main():
    request = json.load(sys.stdin)
    task = request["task"]
    if task == "graph":
        payload = graph(request["smiles"])
    elif task == "canonical":
        payload = {"smiles": canonical_root(request["smiles"])}
    elif task == "mcs":
        payload = mcs(request["from"], request["to"])
    elif task == "rxn":
        payload = rxn_ids(request["smiles"])
    else:
        raise ValueError("Unknown task.")
    json.dump(payload, sys.stdout)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.stderr.write(str(exc) + "\n")
        sys.exit(1)
