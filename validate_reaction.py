#!/usr/bin/env python3
"""Validate an edited multi-monomer reaction and write an HTML review report."""

from __future__ import annotations

import argparse
import html
import json
import sys
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import react_monomers as editor


def molecule_atoms(molecule):
    if isinstance(molecule, editor.Mol2):
        return [(aid, editor._element_from_label(atom_type), float(charge or 0.0))
                for aid, _name, atom_type, charge, _xyz in molecule.atom_rows()]
    return [(aid, editor._element_from_label(label), float(charge or 0.0))
            for aid, label, _atom_type, charge, _xyz, _note in molecule.atom_rows()]


def pair_text(pair, references):
    atom1, atom2 = sorted(pair)
    ref1, ref2 = references[atom1], references[atom2]
    return f"{ref1} — {ref2} (product atoms {atom1}–{atom2})"


def iframe_url(report_path, target):
    try:
        return Path(target).resolve().relative_to(report_path.parent.resolve()).as_posix()
    except ValueError:
        return Path(target).resolve().as_uri()


def validate(monomer_paths, product_path):
    molecules = [editor.detect(path) for path in monomer_paths]
    product = editor.Mol2.read(product_path)
    errors, warnings = [], []
    references, elements, pre_charges = {}, {}, {}
    pre_bonds = set()
    pre_bond_orders = {}
    offset = 0
    for monomer_number, molecule in enumerate(molecules, 1):
        atoms = molecule_atoms(molecule)
        local_ids = {atom_id for atom_id, _element, _charge in atoms}
        for atom_id, element, charge in atoms:
            global_id = offset + atom_id
            references[global_id] = f"{monomer_number}:{atom_id}"
            elements[global_id] = element
            pre_charges[global_id] = charge
        for _bond_id, atom1, atom2, bond_type in editor.bond_rows(molecule):
            if atom1 not in local_ids or atom2 not in local_ids:
                errors.append(f"Monomer {monomer_number} has a bond to an unknown atom")
            pair = frozenset((offset + atom1, offset + atom2))
            pre_bonds.add(pair)
            pre_bond_orders[pair] = str(bond_type)
        offset += len(atoms)

    product_ids = [int(atom[0]) for atom in product.atoms]
    if len(product_ids) != offset:
        errors.append(f"Atom count changed: monomers contain {offset}, product contains {len(product_ids)}")
    if len(product_ids) != len(set(product_ids)):
        errors.append("Product contains duplicate atom IDs")
    if set(product_ids) != set(references):
        errors.append(
            "Product atom IDs do not preserve the concatenated monomer atom mapping: "
            f"missing={sorted(set(references) - set(product_ids))}, "
            f"extra={sorted(set(product_ids) - set(references))}"
        )

    post_pairs = []
    post_bond_orders = {}
    for bond in product.bonds:
        atom1, atom2 = int(bond[1]), int(bond[2])
        if atom1 not in references or atom2 not in references:
            errors.append(f"Product bond {bond[0]} references unknown atom {atom1} or {atom2}")
            continue
        if atom1 == atom2:
            errors.append(f"Product bond {bond[0]} is a self-bond on atom {atom1}")
        pair = frozenset((atom1, atom2))
        post_pairs.append(pair)
        post_bond_orders[pair] = str(bond[3])
    duplicate_pairs = [tuple(sorted(pair)) for pair, count in Counter(post_pairs).items() if count > 1]
    if duplicate_pairs:
        errors.append(f"Product contains duplicate bonds: {duplicate_pairs}")
    post_bonds = set(post_pairs)

    removed = sorted(pre_bonds - post_bonds, key=lambda pair: tuple(sorted(pair)))
    added = sorted(post_bonds - pre_bonds, key=lambda pair: tuple(sorted(pair)))
    order_changes = sorted(
        ((pair, pre_bond_orders[pair], post_bond_orders[pair])
         for pair in pre_bonds & post_bonds
         if pre_bond_orders[pair] != post_bond_orders[pair]),
        key=lambda item: tuple(sorted(item[0])),
    )
    if not removed and not added and not order_changes:
        warnings.append("No connectivity or bond-order change was detected")

    degree = Counter(atom for pair in post_bonds for atom in pair)
    pre_degree = Counter(atom for pair in pre_bonds for atom in pair)

    product_charge = sum(float(atom[8]) if len(atom) > 8 else 0.0 for atom in product.atoms)
    pre_charge = sum(pre_charges.values())
    charge_delta = product_charge - pre_charge
    if abs(charge_delta) > 1.0e-4:
        warnings.append(f"Raw-product charge changed by {charge_delta:+.6f} e")

    changed_pairs = list(removed) + list(added) + [item[0] for item in order_changes]
    changed_atoms = sorted(set().union(*changed_pairs)) if changed_pairs else []
    degree_changes = [
        f"{references[atom_id]} ({elements[atom_id]}): "
        f"degree {pre_degree[atom_id]} → {degree[atom_id]}"
        for atom_id in changed_atoms
    ]
    return {
        "status": "invalid" if errors else ("review" if warnings else "valid"),
        "monomers": [str(path) for path in monomer_paths],
        "product": str(product_path),
        "atom_count": len(product_ids),
        "pre_bond_count": len(pre_bonds),
        "post_bond_count": len(post_bonds),
        "removed_bonds": [pair_text(pair, references) for pair in removed],
        "added_bonds": [pair_text(pair, references) for pair in added],
        "bond_order_changes": [
            f"{pair_text(pair, references)}: {before} → {after}"
            for pair, before, after in order_changes
        ],
        "changed_atoms": [references[atom_id] for atom_id in changed_atoms],
        "degree_changes": degree_changes,
        "pre_charge": pre_charge,
        "product_charge": product_charge,
        "charge_delta": charge_delta,
        "errors": errors,
        "warnings": warnings,
    }


def render_report(result, output, before_graph, after_graph):
    def list_items(values, empty="None"):
        return "<p>" + html.escape(empty) + "</p>" if not values else (
            "<ul>" + "".join(f"<li>{html.escape(str(value))}</li>" for value in values) + "</ul>")

    status = result["status"]
    color = {"valid": "#16803c", "review": "#b36b00", "invalid": "#b42318"}[status]
    document = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Reaction validation report</title><style>
body{{font:15px system-ui;background:#eef3f8;color:#17324d;margin:0}}main{{max-width:1400px;margin:auto;padding:22px}}
.head,.card{{background:white;border-radius:12px;padding:18px;margin-bottom:14px;box-shadow:0 1px 4px #0002}}
.status{{display:inline-block;background:{color};color:white;border-radius:999px;padding:6px 12px;font-weight:700}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:14px}}iframe{{width:100%;height:650px;border:1px solid #ccd6e0;border-radius:9px}}
table{{border-collapse:collapse;width:100%}}td,th{{padding:8px;border-bottom:1px solid #dde5ec;text-align:left}}
@media(max-width:900px){{.grid{{grid-template-columns:1fr}}}}</style></head><body><main>
<section class="head"><span class="status">{status.upper()}</span><h1>Reaction validation report</h1>
<p>This report verifies topology and structural consistency. Chemical intent still requires your approval.</p></section>
<section class="card"><h2>Checks</h2><table>
<tr><th>Atoms</th><td>{result['atom_count']}</td></tr>
<tr><th>Bonds before / after</th><td>{result['pre_bond_count']} / {result['post_bond_count']}</td></tr>
<tr><th>Charge before / after</th><td>{result['pre_charge']:+.6f} / {result['product_charge']:+.6f} e</td></tr>
<tr><th>Changed atoms</th><td>{html.escape(', '.join(result['changed_atoms']) or 'None')}</td></tr></table>
<h3>Coordination changes</h3>{list_items(result['degree_changes'])}</section>
<div class="grid"><section class="card"><h2>Removed bonds</h2>{list_items(result['removed_bonds'])}</section>
<section class="card"><h2>Added bonds</h2>{list_items(result['added_bonds'])}</section></div>
<section class="card"><h2>Bond-order changes</h2>{list_items(result['bond_order_changes'])}</section>
<div class="grid"><section class="card"><h2>Errors</h2>{list_items(result['errors'])}</section>
<section class="card"><h2>Warnings</h2>{list_items(result['warnings'])}</section></div>
<div class="grid"><section class="card"><h2>Before reaction</h2><iframe src="{html.escape(iframe_url(output, before_graph))}"></iframe></section>
<section class="card"><h2>Final minimized product</h2><iframe src="{html.escape(iframe_url(output, after_graph))}"></iframe></section></div>
</main></body></html>"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--monomer", action="append", required=True, type=Path)
    parser.add_argument("--product", required=True, type=Path)
    parser.add_argument("--before-graph", required=True, type=Path)
    parser.add_argument("--after-graph", required=True, type=Path)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    args = parser.parse_args(argv)
    result = validate(
        [path.expanduser().resolve() for path in args.monomer],
        args.product.expanduser().resolve(),
    )
    html_output = args.html.expanduser().resolve()
    json_output = args.json.expanduser().resolve()
    render_report(result, html_output, args.before_graph, args.after_graph)
    json_output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Validation status: {result['status']}")
    print(f"Removed bonds: {result['removed_bonds']}")
    print(f"Added bonds: {result['added_bonds']}")
    print(f"Bond-order changes: {result['bond_order_changes']}")
    print(f"Errors: {result['errors']}")
    print(f"Warnings: {result['warnings']}")
    print(f"Wrote validation report: {html_output}")
    print(f"Wrote validation data: {json_output}")
    return 2 if result["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
