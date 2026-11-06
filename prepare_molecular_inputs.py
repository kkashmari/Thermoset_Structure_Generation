#!/usr/bin/env python3
"""Resolve MOL2 paths and convert SMILES inputs into explicit-H 3D MOL2 files."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem


def existing_mol2(value: str):
    try:
        path = Path(value).expanduser().resolve()
        return path if path.is_file() and path.suffix.lower() == ".mol2" else None
    except OSError:
        return None


def smiles_to_mol2(smiles: str, destination: Path, name: str):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES input: {smiles}")
    molecule = Chem.AddHs(molecule)
    molecule.SetProp("_Name", name)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 0xF00D
    parameters.useSmallRingTorsions = True
    embedded = AllChem.EmbedMolecule(molecule, parameters)
    if embedded != 0:
        parameters.useRandomCoords = True
        embedded = AllChem.EmbedMolecule(molecule, parameters)
    if embedded != 0:
        raise ValueError(f"RDKit could not generate a 3D conformer for SMILES: {smiles}")
    if AllChem.MMFFHasAllMoleculeParams(molecule):
        AllChem.MMFFOptimizeMolecule(molecule, maxIters=1000)
        optimization = "MMFF94"
    elif AllChem.UFFHasAllMoleculeParams(molecule):
        AllChem.UFFOptimizeMolecule(molecule, maxIters=1000)
        optimization = "UFF"
    else:
        optimization = "RDKit ETKDG coordinates; no complete MMFF/UFF parameters"

    executable = shutil.which("obabel")
    if not executable:
        raise FileNotFoundError(
            "Open Babel executable 'obabel' is required to convert generated 3D "
            "structures to MOL2.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="smiles_3d_", dir=destination.parent) as temp:
        sdf_path = Path(temp) / f"{name}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(molecule)
        writer.close()
        result = subprocess.run(
            [executable, "-isdf", str(sdf_path), "-omol2", "-O", str(destination)],
            capture_output=True, text=True, timeout=180, check=False)
    if result.returncode != 0 or not destination.is_file() or destination.stat().st_size == 0:
        detail = (result.stderr or result.stdout or "no diagnostic output").strip()
        raise ValueError(f"Open Babel SMILES-to-MOL2 conversion failed: {detail}")
    check = Chem.MolFromMol2File(str(destination), sanitize=False, removeHs=False)
    if check is None:
        raise ValueError(f"Generated MOL2 cannot be read by RDKit: {destination}")
    failure = Chem.SanitizeMol(check, catchErrors=True)
    if failure != Chem.SanitizeFlags.SANITIZE_NONE:
        raise ValueError(
            f"Generated MOL2 failed RDKit sanitization at operation {failure}: {destination}")
    return {
        "atom_count": check.GetNumAtoms(),
        "explicit_hydrogen_count": sum(atom.GetAtomicNum() == 1 for atom in check.GetAtoms()),
        "coordinate_method": optimization,
    }


def prepare(inputs, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for index, value in enumerate(inputs, 1):
        path = existing_mol2(value)
        if path:
            records.append({
                "index": index, "kind": "mol2", "original_input": value,
                "resolved_mol2": str(path),
            })
            continue
        if value.strip().lower().endswith(".mol2"):
            raise FileNotFoundError(f"MOL2 input does not exist: {value}")
        name = f"monomer_{index:02d}_from_smiles"
        destination = output_dir / f"{name}.mol2"
        if destination.exists():
            raise FileExistsError(f"Refusing to overwrite generated input: {destination}")
        details = smiles_to_mol2(value.strip(), destination, name)
        records.append({
            "index": index, "kind": "smiles", "original_input": value,
            "canonical_smiles": Chem.MolToSmiles(Chem.MolFromSmiles(value.strip())),
            "resolved_mol2": str(destination.resolve()), **details,
        })
    return {
        "schema_version": 1,
        "rdkit_version": rdBase.rdkitVersion,
        "inputs": records,
        "resolved_mol2": [record["resolved_mol2"] for record in records],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args(argv)
    result = prepare(args.input, args.output_dir.expanduser().resolve())
    manifest = args.manifest.expanduser().resolve()
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for record in result["inputs"]:
        print(
            f"Input {record['index']}: {record['kind']} -> {record['resolved_mol2']}")
    print(f"Wrote {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
