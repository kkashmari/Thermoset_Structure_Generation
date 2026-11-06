#!/usr/bin/env python3
"""Interactive multi-monomer topology editor for MOL2 and LAMMPS data files.

This is intentionally a topology editor, not a force-field parameterizer.
It merges two or more complete monomers, applies user-requested bond/type changes, and
writes the complete product. For LAMMPS data, existing higher-order topology
is retained but new angles/dihedrals/impropers are not inferred.
"""

from __future__ import annotations

import argparse
import copy
import json
import re
import shutil
import shlex
import subprocess
import sys
import webbrowser
from dataclasses import dataclass, field
from pathlib import Path


# ---------------------------------------------------------------------------
# DEFAULT INPUT/OUTPUT PATHS
# Edit these values if you prefer to run the script without command-line paths.
# Any paths supplied on the command line override these defaults.
# All input files in one run must be the same format: all MOL2 or all LAMMPS.
# ---------------------------------------------------------------------------
INPUT_PATHS = []
OUTPUT_PATH = ""

# Interactive graph settings.
GENERATE_GRAPH = True
OPEN_GRAPH_IN_BROWSER = True
# Leave blank to create reaction_graph.html beside the product.
GRAPH_OUTPUT_PATH = ""
# Leave blank to create final_product_graph.html beside the product.
FINAL_GRAPH_OUTPUT_PATH = ""

# Automatically minimize MOL2 products with the installed Open Babel program.
MINIMIZE_WITH_OPENBABEL = True
OPENBABEL_FORCEFIELD = "UFF"
OPENBABEL_STEPS = 1500
# Leave blank to create <product_name>_minimized.mol2 beside the raw product.
MINIMIZED_OUTPUT_PATH = ""


SECTION_NAMES = [
    "Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
    "Dihedral Coeffs", "Improper Coeffs", "BondBond Coeffs",
    "BondAngle Coeffs", "AngleAngle Coeffs", "AngleAngleTorsion Coeffs",
    "EndBondTorsion Coeffs", "MiddleBondTorsion Coeffs",
    "BondBond13 Coeffs", "AngleTorsion Coeffs", "Atoms", "Velocities",
    "Bonds", "Angles", "Dihedrals", "Impropers",
]

TYPE_OFFSET_FOR_SECTION = {
    "Masses": "atom", "Pair Coeffs": "atom", "Bond Coeffs": "bond",
    "Angle Coeffs": "angle", "Dihedral Coeffs": "dihedral",
    "Improper Coeffs": "improper", "BondBond Coeffs": "angle",
    "BondAngle Coeffs": "angle", "AngleAngle Coeffs": "improper",
    "AngleAngleTorsion Coeffs": "dihedral",
    "EndBondTorsion Coeffs": "dihedral",
    "MiddleBondTorsion Coeffs": "dihedral",
    "BondBond13 Coeffs": "dihedral", "AngleTorsion Coeffs": "dihedral",
}


def split_comment(line: str) -> tuple[str, str]:
    if "#" in line:
        body, comment = line.split("#", 1)
        return body.rstrip(), " #" + comment.rstrip()
    return line.rstrip(), ""


def numbered_lines(lines):
    return [x for x in lines if x.strip() and not x.lstrip().startswith("#")]


@dataclass
class Mol2:
    path: Path
    name: str = "product"
    atoms: list[list[str]] = field(default_factory=list)
    bonds: list[list[str]] = field(default_factory=list)

    @classmethod
    def read(cls, path: Path) -> "Mol2":
        text = path.read_text(errors="replace").splitlines()
        blocks: dict[str, list[str]] = {}
        current = None
        for line in text:
            if line.startswith("@<TRIPOS>"):
                current = line[9:]
                blocks.setdefault(current, [])
            elif current is not None:
                blocks[current].append(line)
        if "ATOM" not in blocks or "BOND" not in blocks:
            raise ValueError(f"{path}: MOL2 must contain ATOM and BOND blocks")
        mol = cls(path=path)
        molecule = blocks.get("MOLECULE", [])
        if molecule and molecule[0].strip():
            mol.name = molecule[0].strip()
        for line in numbered_lines(blocks["ATOM"]):
            p = line.split()
            if len(p) < 6:
                raise ValueError(f"{path}: malformed MOL2 atom line: {line}")
            mol.atoms.append(p)
        for line in numbered_lines(blocks["BOND"]):
            p = line.split()
            if len(p) < 4:
                raise ValueError(f"{path}: malformed MOL2 bond line: {line}")
            mol.bonds.append(p)
        return mol

    def atom_rows(self):
        for a in self.atoms:
            charge = a[8] if len(a) > 8 else ""
            yield int(a[0]), a[1], a[5], charge, tuple(map(float, a[2:5]))


@dataclass
class LammpsData:
    path: Path
    sections: dict[str, list[str]] = field(default_factory=dict)
    section_headers: dict[str, str] = field(default_factory=dict)
    counts: dict[str, int] = field(default_factory=dict)
    type_counts: dict[str, int] = field(default_factory=dict)

    @classmethod
    def read(cls, path: Path) -> "LammpsData":
        obj = cls(path=path)
        lines = path.read_text(errors="replace").splitlines()
        current = None
        for line in lines:
            stripped = line.strip()
            matched = next((s for s in SECTION_NAMES
                            if stripped == s or stripped.startswith(s + " #")), None)
            if matched:
                current = matched
                obj.sections.setdefault(current, [])
                obj.section_headers[current] = stripped
                continue
            if current is not None:
                obj.sections[current].append(line)
            else:
                m = re.match(r"\s*(\d+)\s+(atoms|bonds|angles|dihedrals|impropers)\s*$", line)
                if m:
                    obj.counts[m.group(2)] = int(m.group(1))
                m = re.match(r"\s*(\d+)\s+(atom|bond|angle|dihedral|improper) types\s*$", line)
                if m:
                    obj.type_counts[m.group(2)] = int(m.group(1))
        if "Atoms" not in obj.sections or "Bonds" not in obj.sections:
            raise ValueError(f"{path}: expected LAMMPS Atoms and Bonds sections")
        if "full" not in obj.section_headers.get("Atoms", "").lower():
            raise ValueError(f"{path}: this prototype supports 'Atoms # full' only")
        return obj

    def atom_rows(self):
        labels = {}
        for line in numbered_lines(self.sections.get("Masses", [])):
            body, comment = split_comment(line)
            p = body.split()
            labels[int(p[0])] = comment.replace("#", "").strip() or "?"
        for line in numbered_lines(self.sections["Atoms"]):
            body, comment = split_comment(line)
            p = body.split()
            if len(p) < 7:
                raise ValueError(f"{self.path}: malformed full atom line: {line}")
            aid, atype = int(p[0]), int(p[2])
            note = comment.replace("#", "").strip()
            yield aid, labels.get(atype, "?"), str(atype), p[3], tuple(map(float, p[4:7])), note


def detect(path: Path):
    suffix = path.suffix.lower()
    if suffix == ".mol2":
        return Mol2.read(path)
    if suffix in {".data", ".lmp"}:
        return LammpsData.read(path)
    if suffix == ".mol":
        raise ValueError("MDL .mol is not MOL2. Convert it to .mol2, or extend the parser for MDL MOL.")
    raise ValueError(f"Unsupported extension {suffix!r}; use .mol2, .data, or .lmp")


def show_atoms(mol, number: int):
    print(f"\nMonomer {number}: {mol.path}")
    if isinstance(mol, Mol2):
        print(" atom   name   MOL2 type        charge          x          y          z")
        for aid, name, typ, q, xyz in mol.atom_rows():
            print(f" {aid:4d}   {name:<5}  {typ:<14} {q:>9}  {xyz[0]:9.4f} {xyz[1]:9.4f} {xyz[2]:9.4f}")
    else:
        print(" atom   element/label   numeric type     charge          x          y          z   annotation")
        for aid, label, typ, q, xyz, note in mol.atom_rows():
            print(f" {aid:4d}   {label:<14} {typ:>8} {q:>11}  {xyz[0]:9.4f} {xyz[1]:9.4f} {xyz[2]:9.4f}   {note}")


def show_lammps_type_map(molecules):
    print("\nCombined LAMMPS type map used in the product:")
    for title, section, kind in (("atom", "Masses", "atom"), ("bond", "Bond Coeffs", "bond")):
        print(f"  {title} types:")
        offset = 0
        for mi, mol in enumerate(molecules, 1):
            for line in numbered_lines(mol.sections.get(section, [])):
                body, comment = split_comment(line); p = body.split()
                if p and p[0].isdigit():
                    combined = int(p[0]) + offset
                    label = comment.replace("#", "").strip() or "(unlabeled)"
                    print(f"    {combined:4d}  monomer {mi} original {p[0]:>3}  {label}")
            offset += mol.type_counts.get(kind, 0)


def bond_rows(mol):
    """Yield local bond ID, local atom IDs, and local bond type/order."""
    if isinstance(mol, Mol2):
        for bond in mol.bonds:
            yield int(bond[0]), int(bond[1]), int(bond[2]), bond[3]
    else:
        for line in numbered_lines(mol.sections["Bonds"]):
            p = split_comment(line)[0].split()
            if len(p) < 4:
                raise ValueError(f"{mol.path}: malformed bond line: {line}")
            yield int(p[0]), int(p[2]), int(p[3]), p[1]


def bond_catalog(molecules):
    """Return the deterministic pre-reaction combined bond numbering."""
    rows = []
    product_number = 1
    type_offset = 0
    for mi, mol in enumerate(molecules, 1):
        for local_id, atom1, atom2, local_type in bond_rows(mol):
            product_type = local_type
            if isinstance(mol, LammpsData):
                product_type = str(int(local_type) + type_offset)
            rows.append((product_number, mi, local_id, atom1, atom2,
                         local_type, product_type))
            product_number += 1
        if isinstance(mol, LammpsData):
            type_offset += mol.type_counts.get("bond", 0)
    return rows


def show_bonds(molecules):
    print("\nExisting bonds (before reaction):")
    print(" product#   local bond   atom pair       local type   product type")
    for product_no, mi, local_id, atom1, atom2, local_type, product_type in bond_catalog(molecules):
        print(f" {product_no:8d}   {mi}:{local_id:<9d} {mi}:{atom1:<5d} {mi}:{atom2:<7d} "
              f"{str(local_type):>10}   {str(product_type):>12}")


def _element_from_label(label):
    text = str(label).strip().lower()
    if text.startswith("cl"): return "Cl"
    if text.startswith("br"): return "Br"
    if text.startswith("n"): return "N"
    if text.startswith("o"): return "O"
    if text.startswith("h"): return "H"
    if text.startswith("s"): return "S"
    if text.startswith("p"): return "P"
    return "C"


def _mol2_type_from_lammps_label(label):
    """Map the common PCFF/IFF labels in the source files to Tripos MOL2 types."""
    text = str(label).strip().lower()
    if text == "cp": return "C.ar"
    if text in {"ct", "c1"}: return "C.1"
    if text.startswith("c"): return "C.3"
    if text in {"nt", "n1"}: return "N.1"
    if text.startswith("n"): return "N.3"
    if text.startswith("o"): return "O.3"
    if text.startswith("h"): return "H"
    if text.startswith("s"): return "S.3"
    if text.startswith("p"): return "P.3"
    return _element_from_label(label)


def _mol2_order_from_lammps_type(mol, local_type):
    """Infer only unambiguous aromatic bonds; otherwise use a single bond."""
    for line in numbered_lines(mol.sections.get("Bond Coeffs", [])):
        body, comment = split_comment(line)
        fields = body.split()
        if fields and fields[0].isdigit() and int(fields[0]) == int(local_type):
            label = comment.replace("#", "").strip().lower()
            return "ar" if label == "cp-cp" else "1"
    return "1"


def graph_data(molecules):
    """Build browser-ready graph data without reaction-specific assumptions."""
    catalog = bond_catalog(molecules)
    product_by_local = {(mi, local_id): product_no
                        for product_no, mi, local_id, *_ in catalog}
    result = []
    for mi, mol in enumerate(molecules, 1):
        atoms = []
        if isinstance(mol, Mol2):
            rows = [(aid, name, typ, q, xyz) for aid, name, typ, q, xyz in mol.atom_rows()]
        else:
            rows = [(aid, label, typ, q, xyz) for aid, label, typ, q, xyz, _ in mol.atom_rows()]
        ranges = []
        for axis in range(3):
            values = [row[4][axis] for row in rows]
            ranges.append(max(values) - min(values))
        axes = sorted(range(3), key=lambda axis: ranges[axis], reverse=True)[:2]
        axis_names = "xyz"
        for aid, label, typ, charge, xyz in rows:
            element = _element_from_label(label)
            atoms.append({"id": aid, "label": str(label), "type": str(typ),
                          "charge": str(charge), "element": element,
                          "x": round(xyz[axes[0]], 6), "y": round(xyz[axes[1]], 6)})
        bonds = []
        adjacency = {a["id"]: set() for a in atoms}
        for local_id, atom1, atom2, local_type in bond_rows(mol):
            bonds.append({"id": local_id, "a": atom1, "b": atom2,
                          "type": str(local_type),
                          "product": product_by_local[(mi, local_id)]})
            adjacency[atom1].add(atom2); adjacency[atom2].add(atom1)
        result.append({"number": mi, "name": mol.path.name, "atoms": atoms,
                       "bonds": bonds,
                       "projection": axis_names[axes[0]] + axis_names[axes[1]]})
    return result


def write_graph_html(molecules, path: Path):
    """Write a dependency-free interactive molecular graph viewer."""
    data_json = json.dumps(graph_data(molecules), separators=(",", ":")).replace("</", "<\\/")
    template = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Molecular topology editor</title>
<style>
:root{color-scheme:light dark;--bg:#f6f7fb;--panel:#fff;--text:#172033;--muted:#64748b;--line:#94a3b8;--accent:#2563eb}
@media(prefers-color-scheme:dark){:root{--bg:#10141f;--panel:#181e2b;--text:#eef2ff;--muted:#a6b0c3;--line:#64748b;--accent:#60a5fa}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--text);font:14px system-ui,-apple-system,sans-serif}.app{max-width:1280px;margin:auto;padding:18px}.toolbar,.selection{display:flex;gap:10px;align-items:center;flex-wrap:wrap;background:var(--panel);padding:12px 14px;border-radius:10px;margin-bottom:12px}.toolbar h1{font-size:18px;margin:0 auto 0 0}.toolbar select,.toolbar button,.selection button,.selection input{font:inherit;padding:7px 10px}.viewer{background:var(--panel);border-radius:10px;overflow:hidden;position:relative}.viewer svg{display:block;width:100%;height:min(72vh,720px);min-height:480px;touch-action:none}.bond{stroke:var(--line);stroke-width:2}.bond-hit{stroke:transparent;stroke-width:16}.bond-label{fill:var(--muted);font-size:11px;paint-order:stroke;stroke:var(--panel);stroke-width:4px;stroke-linejoin:round;pointer-events:none}.node circle{stroke:var(--panel);stroke-width:3}.node text{fill:var(--text);font-size:12px;font-weight:600;paint-order:stroke;stroke:var(--panel);stroke-width:4px}.node.selected circle{stroke:var(--accent);stroke-width:6}.C{fill:#64748b}.H{fill:#d7dce5}.O{fill:#ef4444}.N{fill:#3b82f6}.S{fill:#eab308}.P{fill:#f59e0b}.Cl,.Br{fill:#22c55e}.selection label{color:var(--muted)}#command{flex:1;min-width:280px}.hint{color:var(--muted);margin:8px 2px}.status{min-height:20px;color:var(--muted)}button{border:1px solid color-mix(in srgb,var(--text) 20%,transparent);border-radius:7px;background:var(--panel);color:var(--text)}button:hover{border-color:var(--accent)}
</style></head><body><main class="app">
<div class="toolbar"><h1>Reaction atom and bond selector</h1><label>Monomer <select id="monomer"></select></label><label><input id="hydrogens" type="checkbox" checked> Show H</label><button id="minus">−</button><button id="fit">Fit</button><button id="plus">+</button><button id="clear">Clear selection</button></div>
<div class="viewer"><svg id="graph" viewBox="0 0 1000 650" aria-label="Interactive molecular connectivity graph"><g id="scene"></g></svg></div>
<p class="hint">Click atoms to build an add/remove command. Click a bond to select it. Drag empty space to pan; use the mouse wheel to zoom. Node text is atom ID · atom type. Edge text is bond ID · bond type.</p>
<div class="selection"><label for="command">Suggested command</label><input id="command" readonly value="Select two atoms or one bond"><button id="copy">Copy</button><span id="status" class="status" aria-live="polite"></span></div>
</main><script>
const molecules=__DATA__;
const svg=document.getElementById('graph'),scene=document.getElementById('scene'),picker=document.getElementById('monomer'),hydrogens=document.getElementById('hydrogens'),command=document.getElementById('command'),status=document.getElementById('status');
let active=0,selected=[],scale=1,panX=0,panY=0,drag=null;
const NS='http://www.w3.org/2000/svg';
function el(name,attrs={},text=''){const n=document.createElementNS(NS,name);for(const[k,v]of Object.entries(attrs))n.setAttribute(k,v);if(text)n.textContent=text;return n}
function fitPoints(atoms){const xs=atoms.map(a=>a.x),ys=atoms.map(a=>a.y),minX=Math.min(...xs),maxX=Math.max(...xs),minY=Math.min(...ys),maxY=Math.max(...ys),w=Math.max(maxX-minX,1),h=Math.max(maxY-minY,1),s=Math.min(850/w,530/h);return atoms.map(a=>({...a,px:75+(a.x-minX)*s,py:590-(a.y-minY)*s}))}
function updateTransform(){scene.setAttribute('transform',`translate(${panX} ${panY}) scale(${scale})`)}
function render(){scene.replaceChildren();const mol=molecules[active],visible=hydrogens.checked?mol.atoms:mol.atoms.filter(a=>a.element!=='H'),atoms=fitPoints(visible),byId=new Map(atoms.map(a=>[a.id,a]));
 for(const b of mol.bonds){const a=byId.get(b.a),c=byId.get(b.b);if(!a||!c)continue;const g=el('g',{'data-bond':b.id});g.append(el('line',{x1:a.px,y1:a.py,x2:c.px,y2:c.py,class:'bond'}));const hit=el('line',{x1:a.px,y1:a.py,x2:c.px,y2:c.py,class:'bond-hit'});hit.addEventListener('click',e=>{e.stopPropagation();selected=[];command.value=`remove-bond ${mol.number}:${b.id}`;status.textContent=`Product bond ${b.product}: atoms ${mol.number}:${b.a}–${mol.number}:${b.b}, type ${b.type}`;render()});g.append(hit);g.append(el('text',{x:(a.px+c.px)/2,y:(a.py+c.py)/2-5,'text-anchor':'middle',class:'bond-label'},`${b.id} · ${b.type}`));scene.append(g)}
 for(const a of atoms){const g=el('g',{class:`node ${a.element}${selected.some(s=>s.m===mol.number&&s.a===a.id)?' selected':''}`,'data-atom':a.id});g.append(el('circle',{cx:a.px,cy:a.py,r:a.element==='H'?12:18}));g.append(el('text',{x:a.px,y:a.py-22,'text-anchor':'middle'},`${a.id} · ${a.type}`));g.addEventListener('click',e=>{e.stopPropagation();chooseAtom(mol,a)});scene.append(g)}updateTransform();status.textContent=status.textContent||`${mol.name}: ${mol.atoms.length} atoms, ${mol.bonds.length} bonds; projection ${mol.projection}`}
function chooseAtom(mol,a){const ref={m:mol.number,a:a.id};const at=selected.findIndex(s=>s.m===ref.m&&s.a===ref.a);if(at>=0)selected.splice(at,1);else{if(selected.length===2)selected=[];selected.push(ref)}if(selected.length===1){command.value=`Selected ${ref.m}:${ref.a}; select one more atom`;status.textContent=`Atom ${ref.m}:${ref.a}, ${a.label}, type ${a.type}, charge ${a.charge}`}else if(selected.length===2){const[x,y]=selected;const xm=molecules[x.m-1];const bond=x.m===y.m?xm.bonds.find(b=>(b.a===x.a&&b.b===y.a)||(b.a===y.a&&b.b===x.a)):null;command.value=bond?`remove ${x.m}:${x.a} ${y.m}:${y.a}`:`add ${x.m}:${x.a} ${y.m}:${y.a} <bond-type>`;status.textContent=bond?`Existing bond ${x.m}:${bond.id}, type ${bond.type}`:'No existing bond between the selected atoms'}render()}
molecules.forEach((m,i)=>{const o=document.createElement('option');o.value=i;o.textContent=`${m.number}: ${m.name}`;picker.append(o)});picker.addEventListener('change',()=>{active=+picker.value;status.textContent='';render()});hydrogens.addEventListener('change',render);
document.getElementById('clear').onclick=()=>{selected=[];command.value='Select two atoms or one bond';status.textContent='';render()};document.getElementById('plus').onclick=()=>{scale*=1.2;updateTransform()};document.getElementById('minus').onclick=()=>{scale/=1.2;updateTransform()};document.getElementById('fit').onclick=()=>{scale=1;panX=panY=0;updateTransform()};document.getElementById('copy').onclick=async()=>{try{await navigator.clipboard.writeText(command.value);status.textContent='Command copied'}catch{command.select();status.textContent='Press Ctrl/Cmd+C to copy'}};
svg.addEventListener('wheel',e=>{e.preventDefault();scale*=e.deltaY<0?1.1:.9;updateTransform()},{passive:false});svg.addEventListener('pointerdown',e=>{if(e.target===svg||e.target===scene){drag={x:e.clientX-panX,y:e.clientY-panY};svg.setPointerCapture(e.pointerId)}});svg.addEventListener('pointermove',e=>{if(drag){panX=e.clientX-drag.x;panY=e.clientY-drag.y;updateTransform()}});svg.addEventListener('pointerup',()=>drag=null);render();
</script></body></html>'''
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(template.replace("__DATA__", data_json), encoding="utf-8")


def pair_from_local_bond(token, molecules):
    m = re.fullmatch(r"(\d+):(\d+)", token)
    if not m:
        raise ValueError(f"Bad bond reference {token!r}; use MONOMER_NUMBER:BOND_ID")
    mi, bond_id = int(m.group(1)), int(m.group(2))
    if not 1 <= mi <= len(molecules):
        raise ValueError(f"Monomer {mi} is unavailable; choose 1..{len(molecules)}")
    matches = [(a1, a2) for bid, a1, a2, _ in bond_rows(molecules[mi - 1]) if bid == bond_id]
    if not matches:
        raise ValueError(f"Monomer {mi} has no bond {bond_id}")
    a1, a2 = matches[0]
    return (mi, a1), (mi, a2)


def pair_from_product_bond(token, molecules):
    try:
        number = int(token)
    except ValueError as exc:
        raise ValueError("Combined bond number must be an integer") from exc
    matches = [(mi, a1, a2) for product_no, mi, _, a1, a2, _, _
               in bond_catalog(molecules) if product_no == number]
    if not matches:
        raise ValueError(f"There is no pre-reaction combined bond number {number}")
    mi, a1, a2 = matches[0]
    return (mi, a1), (mi, a2)


def parse_ref(token: str, molecules) -> tuple[int, int]:
    m = re.fullmatch(r"(\d+):(\d+)", token)
    if not m:
        raise ValueError(f"Bad atom reference {token!r}; use MONOMER_NUMBER:ATOM_ID")
    mi, aid = int(m.group(1)), int(m.group(2))
    if not 1 <= mi <= len(molecules):
        raise ValueError(f"Monomer {mi} is unavailable; choose 1..{len(molecules)}")
    ids = {row[0] for row in molecules[mi - 1].atom_rows()}
    if aid not in ids:
        raise ValueError(f"Monomer {mi} has no atom {aid}")
    return mi, aid


def read_commands(molecules, is_mol2):
    removes, adds, retypes, reorders = [], [], [], []
    type_word = "MOL2 bond order (1, 2, 3, ar, am)" if is_mol2 else "combined numeric bond type"
    new_type_word = "NEW_ORDER (1, 2, 3, ar, am)" if is_mol2 else "NEW_COMBINED_BOND_TYPE"
    print("\nDescribe the reaction one command per line:")
    print("  remove        MONOMER:ATOM MONOMER:ATOM")
    print("  remove-bond   MONOMER:BOND_ID")
    print("  remove-number COMBINED_BOND_NUMBER")
    print(f"  add           MONOMER:ATOM MONOMER:ATOM {type_word}")
    print(f"  add-number    PRODUCT_POSITION MONOMER:ATOM MONOMER:ATOM {type_word}")
    print(f"  order         MONOMER:ATOM MONOMER:ATOM {new_type_word}")
    print(f"  order-bond    MONOMER:BOND_ID {new_type_word}")
    print(f"  order-number  COMBINED_BOND_NUMBER {new_type_word}")
    print("  type   MONOMER:ATOM NEW_ATOM_TYPE")
    print("  show atoms | show bonds")
    print("  done   (write the product)")
    while True:
        try:
            raw = input("reaction> ").strip()
        except EOFError:
            raw = "done"
        if not raw or raw.lower() == "done":
            break
        try:
            p = shlex.split(raw)
            cmd = p[0].lower()
            if cmd == "show" and (len(p) == 1 or (len(p) == 2 and p[1].lower() == "atoms")):
                for number, molecule in enumerate(molecules, 1):
                    show_atoms(molecule, number)
            elif cmd == "show" and len(p) == 2 and p[1].lower() == "bonds":
                show_bonds(molecules)
            elif cmd == "remove" and len(p) == 3:
                removes.append((parse_ref(p[1], molecules), parse_ref(p[2], molecules)))
            elif cmd == "remove-bond" and len(p) == 2:
                removes.append(pair_from_local_bond(p[1], molecules))
            elif cmd == "remove-number" and len(p) == 2:
                removes.append(pair_from_product_bond(p[1], molecules))
            elif cmd == "add" and len(p) == 4:
                if not is_mol2 and (not p[3].isdigit() or int(p[3]) < 1):
                    raise ValueError("LAMMPS bond type must be a positive integer")
                adds.append((parse_ref(p[1], molecules), parse_ref(p[2], molecules), p[3], None))
            elif cmd == "add-number" and len(p) == 5:
                if not p[1].isdigit() or int(p[1]) < 1:
                    raise ValueError("Product bond position must be a positive integer")
                if not is_mol2 and (not p[4].isdigit() or int(p[4]) < 1):
                    raise ValueError("LAMMPS bond type must be a positive integer")
                adds.append((parse_ref(p[2], molecules), parse_ref(p[3], molecules),
                             p[4], int(p[1])))
            elif cmd == "order" and len(p) == 4:
                if not is_mol2 and (not p[3].isdigit() or int(p[3]) < 1):
                    raise ValueError("LAMMPS bond type must be a positive integer")
                reorders.append((parse_ref(p[1], molecules), parse_ref(p[2], molecules), p[3]))
            elif cmd == "order-bond" and len(p) == 3:
                if not is_mol2 and (not p[2].isdigit() or int(p[2]) < 1):
                    raise ValueError("LAMMPS bond type must be a positive integer")
                atom1, atom2 = pair_from_local_bond(p[1], molecules)
                reorders.append((atom1, atom2, p[2]))
            elif cmd == "order-number" and len(p) == 3:
                if not is_mol2 and (not p[2].isdigit() or int(p[2]) < 1):
                    raise ValueError("LAMMPS bond type must be a positive integer")
                atom1, atom2 = pair_from_product_bond(p[1], molecules)
                reorders.append((atom1, atom2, p[2]))
            elif cmd == "type" and len(p) == 3:
                if not is_mol2 and (not p[2].isdigit() or int(p[2]) < 1):
                    raise ValueError("LAMMPS atom type must be a positive integer")
                retypes.append((parse_ref(p[1], molecules), p[2]))
            else:
                raise ValueError("Command not understood")
        except ValueError as exc:
            print(f"  Error: {exc}")
    return removes, adds, retypes, reorders


def x_translations(molecules, gap=5.0):
    """Lay any number of intact monomers side by side along x."""
    shifts = []
    right_edge = None
    for mol in molecules:
        xs = [r[-1][0] if isinstance(mol, Mol2) else r[4][0] for r in mol.atom_rows()]
        if not xs:
            raise ValueError(f"{mol.path}: no atoms")
        shift = 0.0 if right_edge is None else right_edge + gap - min(xs)
        shifts.append(shift)
        right_edge = max(xs) + shift
    return shifts


def merge_mol2(molecules, removes, adds, retypes, reorders, output, shifts):
    atoms, bonds, idmap = [], [], {}
    for mi, mol in enumerate(molecules, 1):
        if isinstance(mol, Mol2):
            source_atoms = [a.copy() for a in mol.atoms]
        else:
            source_atoms = []
            for aid, label, _numeric_type, charge, xyz, _note in mol.atom_rows():
                element = _element_from_label(label)
                source_atoms.append([str(aid), f"{element}{aid}", str(xyz[0]),
                                     str(xyz[1]), str(xyz[2]),
                                     _mol2_type_from_lammps_label(label),
                                     str(mi), f"MON{mi}", charge])
        for atom in source_atoms:
            new = atom.copy(); old = int(atom[0]); new_id = len(atoms) + 1
            idmap[(mi, old)] = new_id; new[0] = str(new_id)
            new[2] = f"{float(new[2]) + shifts[mi - 1]:.6f}"
            if len(new) > 6: new[6] = str(mi)
            if len(new) > 7: new[7] = f"MON{mi}"
            atoms.append(new)
        for local_id, atom1, atom2, local_type in bond_rows(mol):
            # LAMMPS bond type IDs are force-field types, not MOL2 bond orders.
            # Use a single bond on conversion; user-added MOL2 bonds keep the
            # explicitly supplied order.
            order = (local_type if isinstance(mol, Mol2)
                     else _mol2_order_from_lammps_type(mol, local_type))
            bonds.append(["0", str(idmap[(mi, atom1)]), str(idmap[(mi, atom2)]), order])
    remove_keys = {frozenset((idmap[a], idmap[b])) for a, b in removes}
    before = len(bonds)
    bonds = [b for b in bonds if frozenset((int(b[1]), int(b[2]))) not in remove_keys]
    if before - len(bonds) != len(remove_keys):
        raise ValueError("At least one requested bond removal did not match an existing bond")
    for a, b, order in reorders:
        key = frozenset((idmap[a], idmap[b]))
        matches = [bond for bond in bonds
                   if frozenset((int(bond[1]), int(bond[2]))) == key]
        if len(matches) != 1:
            raise ValueError(f"Cannot change order of {a}-{b}: existing bond not found")
        matches[0][3] = order
    for a, b, order, requested_number in adds:
        key = frozenset((idmap[a], idmap[b]))
        if any(frozenset((int(x[1]), int(x[2]))) == key for x in bonds):
            raise ValueError(f"Cannot add {a}-{b}: bond already exists")
        new_bond = ["0", str(idmap[a]), str(idmap[b]), order]
        if requested_number is None:
            bonds.append(new_bond)
        else:
            if requested_number > len(bonds) + 1:
                raise ValueError(
                    f"Requested product bond number {requested_number} is out of range; "
                    f"choose 1..{len(bonds) + 1} after removals")
            bonds.insert(requested_number - 1, new_bond)
    for ref, typ in retypes:
        atoms[idmap[ref] - 1][5] = typ
    for i, b in enumerate(bonds, 1): b[0] = str(i)
    with output.open("w") as f:
        f.write("@<TRIPOS>MOLECULE\nreacted_product\n")
        f.write(f" {len(atoms)} {len(bonds)} 0 0 0\nSMALL\nUSER_CHARGES\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for a in atoms:
            q = a[8] if len(a) > 8 else "0.0"
            f.write(f"{int(a[0]):7d} {a[1]:<8} {float(a[2]):11.6f} {float(a[3]):11.6f} {float(a[4]):11.6f} {a[5]:<8} {a[6] if len(a)>6 else '1':>3} {a[7] if len(a)>7 else 'MOL':<8} {float(q):11.6f}\n")
        f.write("@<TRIPOS>BOND\n")
        for b in bonds:
            f.write(f"{int(b[0]):6d} {int(b[1]):6d} {int(b[2]):6d} {b[3]:>4}\n")


def transform_coeff_line(line, offset):
    body, comment = split_comment(line)
    p = body.split()
    if not p or not p[0].isdigit(): return line
    p[0] = str(int(p[0]) + offset)
    return " ".join(p) + comment


def merge_lammps(molecules, removes, adds, retypes, reorders, output, shifts):
    kinds = ("atom", "bond", "angle", "dihedral", "improper")
    running_types = {k: 0 for k in kinds}
    type_offsets = []
    for mol in molecules:
        type_offsets.append(running_types.copy())
        for kind in kinds:
            running_types[kind] += mol.type_counts.get(kind, 0)
    combined_counts = running_types

    idmap = {}
    next_atom_id = 1
    molecule_offsets = []
    molecule_base = 0
    for mi, mol in enumerate(molecules, 1):
        atom_lines = numbered_lines(mol.sections["Atoms"])
        old_atom_ids = [int(split_comment(line)[0].split()[0]) for line in atom_lines]
        for old_id in old_atom_ids:
            idmap[(mi, old_id)] = next_atom_id
            next_atom_id += 1
        old_molecule_ids = [int(split_comment(line)[0].split()[1]) for line in atom_lines]
        molecule_offsets.append(molecule_base)
        molecule_base += max(old_molecule_ids, default=0)

    merged: dict[str, list[str]] = {}
    all_sections = list(dict.fromkeys(sec for mol in molecules for sec in mol.sections))
    topology = {"Bonds": ("bond", 2), "Angles": ("angle", 3),
                "Dihedrals": ("dihedral", 4), "Impropers": ("improper", 4)}
    for sec in all_sections:
        out = []
        for mi, src in enumerate(molecules, 1):
            lines = numbered_lines(src.sections.get(sec, []))
            for line in lines:
                if sec in TYPE_OFFSET_FOR_SECTION:
                    off = type_offsets[mi - 1][TYPE_OFFSET_FOR_SECTION[sec]]
                    out.append(transform_coeff_line(line, off))
                elif sec == "Atoms":
                    body, comment = split_comment(line); p = body.split()
                    old = int(p[0]); p[0] = str(idmap[(mi, old)])
                    p[1] = str(int(p[1]) + molecule_offsets[mi - 1])
                    p[2] = str(int(p[2]) + type_offsets[mi - 1]["atom"])
                    p[4] = f"{float(p[4]) + shifts[mi - 1]:.10f}"
                    out.append(" ".join(p) + comment)
                elif sec == "Velocities":
                    body, comment = split_comment(line); p = body.split()
                    p[0] = str(idmap[(mi, int(p[0]))]); out.append(" ".join(p) + comment)
                elif sec in topology:
                    kind, nrefs = topology[sec]; body, comment = split_comment(line); p = body.split()
                    p[0] = "0"; p[1] = str(int(p[1]) + type_offsets[mi - 1][kind])
                    for j in range(2, 2 + nrefs): p[j] = str(idmap[(mi, int(p[j]))])
                    out.append(" ".join(p) + comment)
                else:
                    out.append(line)
        merged[sec] = out
    # Apply atom retyping after type namespaces have been merged.
    for ref, typ in retypes:
        if int(typ) > combined_counts["atom"]:
            raise ValueError(f"Atom type {typ} has no mass/coefficient entry; available combined types are 1..{combined_counts['atom']}")
    type_changes = {idmap[ref]: typ for ref, typ in retypes}
    changed_atoms = []
    for line in merged["Atoms"]:
        body, comment = split_comment(line); p = body.split()
        if int(p[0]) in type_changes: p[2] = type_changes[int(p[0])]
        changed_atoms.append(" ".join(p) + comment)
    merged["Atoms"] = changed_atoms
    # Apply bond removals/additions.
    remove_keys = {frozenset((idmap[a], idmap[b])) for a, b in removes}
    old_bonds = merged.get("Bonds", [])
    kept = []
    for line in old_bonds:
        body, comment = split_comment(line); p = body.split()
        if frozenset((int(p[2]), int(p[3]))) not in remove_keys: kept.append(line)
    if len(old_bonds) - len(kept) != len(remove_keys):
        raise ValueError("At least one requested bond removal did not match an existing bond")
    for a, b, typ in reorders:
        if not typ.isdigit() or int(typ) > combined_counts["bond"]:
            raise ValueError(
                f"Bond type {typ} has no coefficient; available combined types are "
                f"1..{combined_counts['bond']}")
        key = frozenset((idmap[a], idmap[b]))
        match_indexes = []
        for index, line in enumerate(kept):
            fields = split_comment(line)[0].split()
            if frozenset((int(fields[2]), int(fields[3]))) == key:
                match_indexes.append(index)
        if len(match_indexes) != 1:
            raise ValueError(f"Cannot change type of {a}-{b}: existing bond not found")
        index = match_indexes[0]
        body, comment = split_comment(kept[index]); fields = body.split()
        fields[1] = typ
        kept[index] = " ".join(fields) + comment
    existing = {frozenset((int(split_comment(x)[0].split()[2]), int(split_comment(x)[0].split()[3]))) for x in kept}
    for a, b, typ, requested_number in adds:
        key = frozenset((idmap[a], idmap[b]))
        if key in existing: raise ValueError(f"Cannot add {a}-{b}: bond already exists")
        if int(typ) > combined_counts["bond"]:
            raise ValueError(f"Bond type {typ} has no coefficient; available combined types are 1..{combined_counts['bond']}")
        new_bond = f"0 {typ} {idmap[a]} {idmap[b]} # added by react_monomers"
        if requested_number is None:
            kept.append(new_bond)
        else:
            if requested_number > len(kept) + 1:
                raise ValueError(
                    f"Requested product bond number {requested_number} is out of range; "
                    f"choose 1..{len(kept) + 1} after removals")
            kept.insert(requested_number - 1, new_bond)
        existing.add(key)
    for i, line in enumerate(kept, 1):
        body, comment = split_comment(line); p = body.split(); p[0] = str(i); kept[i-1] = " ".join(p) + comment
    merged["Bonds"] = kept
    # Renumber all other topology records.
    for sec in ("Angles", "Dihedrals", "Impropers"):
        for i, line in enumerate(merged.get(sec, []), 1):
            body, comment = split_comment(line); p = body.split(); p[0] = str(i); merged[sec][i-1] = " ".join(p) + comment
    counts = {"atoms": len(merged.get("Atoms", [])), "bonds": len(merged.get("Bonds", [])),
              "angles": len(merged.get("Angles", [])), "dihedrals": len(merged.get("Dihedrals", [])),
              "impropers": len(merged.get("Impropers", []))}
    xyz = []
    for line in merged["Atoms"]:
        p = split_comment(line)[0].split(); xyz.append(tuple(map(float, p[4:7])))
    with output.open("w") as f:
        f.write(f"LAMMPS data file: {len(molecules)} complete monomers merged by react_monomers.py\n\n")
        for k in ("atoms", "bonds", "angles", "dihedrals", "impropers"): f.write(f"{counts[k]:8d} {k}\n")
        f.write("\n")
        for k in ("atom", "bond", "angle", "dihedral", "improper"): f.write(f"{combined_counts[k]:8d} {k} types\n")
        f.write("\n")
        for axis, idx in (("x",0),("y",1),("z",2)):
            vals = [v[idx] for v in xyz]; f.write(f"{min(vals)-5:16.8f} {max(vals)+5:16.8f} {axis}lo {axis}hi\n")
        for sec in all_sections:
            if not merged.get(sec): continue
            header = next((mol.section_headers[sec] for mol in molecules if sec in mol.section_headers), sec)
            f.write(f"\n{header}\n\n")
            f.write("\n".join(merged[sec]) + "\n")


def minimize_mol2_openbabel(source: Path, destination: Path, forcefield: str, steps: int):
    """Minimize a MOL2 product with Open Babel and validate its output."""
    executable = shutil.which("obabel")
    if not executable:
        raise ValueError(
            "Open Babel executable 'obabel' was not found. Install Open Babel "
            "or run with --no-minimize.")
    if steps < 1:
        raise ValueError("Open Babel minimization steps must be positive")
    if destination.exists():
        raise ValueError(f"Refusing to overwrite existing minimized file: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    command = [executable, str(source), "-O", str(destination), "--minimize",
               "--ff", forcefield, "--steps", str(steps)]
    try:
        result = subprocess.run(command, capture_output=True, text=True,
                                timeout=300, check=False)
    except subprocess.TimeoutExpired as exc:
        raise ValueError("Open Babel minimization exceeded 300 seconds") from exc
    if result.returncode != 0 or not destination.exists() or destination.stat().st_size == 0:
        details = (result.stderr or result.stdout or "no diagnostic output").strip()
        raise ValueError(f"Open Babel minimization failed: {details}")
    # Ensure the minimized file is readable before using it for visualization.
    detect(destination)
    return result


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("monomers", nargs="*", type=Path,
                    help="one or more molecular files; if omitted, use INPUT_PATHS in the script")
    ap.add_argument("-o", "--output", type=Path,
                    help="product file; overrides OUTPUT_PATH in the script")
    ap.add_argument("--gap", type=float, default=5.0, help="x-axis separation between monomers (angstrom; default 5)")
    ap.add_argument("--graph-output", type=Path,
                    help="interactive graph HTML path; overrides GRAPH_OUTPUT_PATH")
    ap.add_argument("--final-graph-output", type=Path,
                    help="final reacted-product graph path; overrides FINAL_GRAPH_OUTPUT_PATH")
    ap.add_argument("--no-graph", action="store_true", help="do not generate the interactive graph")
    ap.add_argument("--no-open-graph", action="store_true", help="generate the graph without opening a browser")
    ap.add_argument("--no-minimize", action="store_true", help="skip automatic Open Babel minimization")
    ap.add_argument("--minimized-output", type=Path,
                    help="minimized MOL2 path; overrides MINIMIZED_OUTPUT_PATH")
    ap.add_argument("--minimize-steps", type=int, default=OPENBABEL_STEPS,
                    help=f"Open Babel minimization steps (default {OPENBABEL_STEPS})")
    args = ap.parse_args(argv)
    try:
        input_paths = args.monomers or [Path(p) for p in INPUT_PATHS]
        if not input_paths:
            raise ValueError("Provide at least one molecular input path")
        molecules = [detect(path) for path in input_paths]
        if len({type(mol) for mol in molecules}) != 1:
            raise ValueError("All inputs in one run must have the same format")
        for number, molecule in enumerate(molecules, 1):
            show_atoms(molecule, number)
        input_is_mol2 = isinstance(molecules[0], Mol2)
        # MOL2 is the default product even when the source monomers are LAMMPS data.
        ext = ".mol2"
        configured_output = args.output or (Path(OUTPUT_PATH) if OUTPUT_PATH else Path(f"reacted_product{ext}"))
        if configured_output.is_dir() or configured_output.suffix.lower() not in {".mol2", ".data", ".lmp"}:
            output = configured_output / f"reacted_product{ext}"
        else:
            output = configured_output
        output_is_mol2 = output.suffix.lower() == ".mol2"
        if not output_is_mol2 and input_is_mol2:
            raise ValueError("LAMMPS output requires LAMMPS data inputs; use a .mol2 output for MOL2 inputs")
        if not input_is_mol2 and not output_is_mol2:
            show_lammps_type_map(molecules)
        show_bonds(molecules)
        if GENERATE_GRAPH and not args.no_graph:
            graph_path = args.graph_output or (Path(GRAPH_OUTPUT_PATH) if GRAPH_OUTPUT_PATH
                                                else output.parent / "reaction_graph.html")
            write_graph_html(molecules, graph_path)
            print(f"\nInteractive graph: {graph_path.resolve()}")
            if OPEN_GRAPH_IN_BROWSER and not args.no_open_graph:
                if not webbrowser.open(graph_path.resolve().as_uri()):
                    print("The browser did not open automatically; open the graph path manually.")
        removes, adds, retypes, reorders = read_commands(molecules, output_is_mol2)
        if output.exists(): raise ValueError(f"Refusing to overwrite existing file: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        shifts = x_translations(molecules, args.gap)
        if output_is_mol2:
            merge_mol2(molecules, removes, adds, retypes, reorders, output, shifts)
        else:
            merge_lammps(molecules, removes, adds, retypes, reorders, output, shifts)
        visualized_product = output
        minimized_output = None
        if output_is_mol2 and MINIMIZE_WITH_OPENBABEL and not args.no_minimize:
            default_minimized = output.with_name(f"{output.stem}_minimized.mol2")
            configured_minimized = args.minimized_output or (
                Path(MINIMIZED_OUTPUT_PATH) if MINIMIZED_OUTPUT_PATH else default_minimized)
            if configured_minimized.is_dir() or configured_minimized.suffix.lower() != ".mol2":
                minimized_output = configured_minimized / default_minimized.name
            else:
                minimized_output = configured_minimized
            if minimized_output.resolve() == output.resolve():
                raise ValueError("The minimized output path must differ from the raw product path")
            print(f"\nMinimizing product with Open Babel {OPENBABEL_FORCEFIELD} "
                  f"({args.minimize_steps} steps)...")
            minimize_mol2_openbabel(output, minimized_output,
                                    OPENBABEL_FORCEFIELD, args.minimize_steps)
            visualized_product = minimized_output
        if GENERATE_GRAPH and not args.no_graph:
            final_graph_path = args.final_graph_output or (
                Path(FINAL_GRAPH_OUTPUT_PATH) if FINAL_GRAPH_OUTPUT_PATH
                else output.parent / "final_product_graph.html")
            final_product = detect(visualized_product)
            write_graph_html([final_product], final_graph_path)
            print(f"\nFinal product graph: {final_graph_path.resolve()}")
            if OPEN_GRAPH_IN_BROWSER and not args.no_open_graph:
                if not webbrowser.open(final_graph_path.resolve().as_uri()):
                    print("The final-product graph did not open automatically; open the path manually.")
        print(f"\nWrote raw reacted structure: {output.resolve()}")
        if minimized_output is not None:
            print(f"Wrote minimized final structure: {minimized_output.resolve()}")
        for number, shift in enumerate(shifts[1:], 2):
            print(f"Monomer {number} was translated by {shift:+.3f} angstrom along x to avoid overlap.")
        if not output_is_mol2:
            print("WARNING: Existing angles/dihedrals/impropers were preserved, but terms spanning the new bond were not generated.")
            print("Run the product through your force-field typing/topology tool before simulation.")
        elif not input_is_mol2:
            print("NOTE: LAMMPS bonds were converted to MOL2 orders (cp-cp aromatic; otherwise single); no coefficients were written.")
        return 0
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr); return 2


if __name__ == "__main__":
    raise SystemExit(main())
