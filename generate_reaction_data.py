#!/usr/bin/env python3
"""Extract reduced pre/post LAMMPS data templates from generic graph differences."""

from __future__ import annotations

import argparse
import json
from collections import deque
from pathlib import Path

import lammps_data_io as base


def bfs(adjacency, centers, radius):
    selected = set()
    for center in centers:
        distances = {center: 0}
        queue = deque([center])
        while queue:
            atom = queue.popleft()
            selected.add(atom)
            if distances[atom] >= radius:
                continue
            for neighbor in adjacency.get(atom, ()):
                if neighbor not in distances:
                    distances[neighbor] = distances[atom] + 1
                    queue.append(neighbor)
    return selected


def parse_source_ref(value, source_to_global):
    try:
        monomer, atom = map(int, value.split(":"))
    except ValueError as exc:
        raise ValueError(f"Bad atom reference {value!r}; use MONOMER:ATOM") from exc
    key = (monomer, atom)
    if key not in source_to_global:
        raise ValueError(f"Unknown source atom reference {value}")
    return source_to_global[key]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--monomer", action="append", required=True, type=Path)
    parser.add_argument("--product", required=True, type=Path)
    parser.add_argument("--pre-output", required=True, type=Path)
    parser.add_argument("--post-output", required=True, type=Path)
    parser.add_argument("--atom-map", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--metadata", type=Path,
                        help="optional machine-readable reaction-candidate JSON")
    parser.add_argument("--radius", type=int, default=3)
    parser.add_argument(
        "--initiator", nargs=2, metavar=("MONOMER:ATOM", "MONOMER:ATOM"),
        help="initiator source atoms; needed if more than one new bond is possible",
    )
    args = parser.parse_args(argv)
    if args.radius < 0:
        raise ValueError("Template radius must be non-negative")

    monomers = [base.DataFile(path.expanduser().resolve()) for path in args.monomer]
    product = base.DataFile(args.product.expanduser().resolve())

    source_to_global = {}
    global_to_source = {}
    next_global = 1
    for monomer_number, data in enumerate(monomers, 1):
        for source_atom in data.atom_records():
            source_to_global[(monomer_number, source_atom)] = next_global
            global_to_source[next_global] = (monomer_number, source_atom)
            next_global += 1
    expected_ids = set(range(1, next_global))
    product_ids = set(product.atom_records())
    if product_ids != expected_ids:
        raise ValueError(
            "Parameterized product did not preserve concatenated atom IDs: "
            f"missing={sorted(expected_ids - product_ids)}, extra={sorted(product_ids - expected_ids)}"
        )

    pre_bonds = set()
    pre_adjacency = {atom: set() for atom in expected_ids}
    for monomer_number, data in enumerate(monomers, 1):
        for _record, _bond_type, refs in data.topology("Bonds"):
            atom1 = source_to_global[(monomer_number, refs[0])]
            atom2 = source_to_global[(monomer_number, refs[1])]
            pair = tuple(sorted((atom1, atom2)))
            pre_bonds.add(pair)
            pre_adjacency[atom1].add(atom2)
            pre_adjacency[atom2].add(atom1)

    post_bonds = set()
    post_adjacency = {atom: set() for atom in expected_ids}
    for _record, _bond_type, refs in product.topology("Bonds"):
        atom1, atom2 = refs
        pair = tuple(sorted((atom1, atom2)))
        post_bonds.add(pair)
        post_adjacency[atom1].add(atom2)
        post_adjacency[atom2].add(atom1)
    removed_global = sorted(pre_bonds - post_bonds)
    added_global = sorted(post_bonds - pre_bonds)
    if not removed_global and not added_global:
        raise ValueError("No bond connectivity changes were found between monomers and product")

    changed = set(atom for pair in removed_global + added_global for atom in pair)
    selected = bfs(pre_adjacency, changed, args.radius) | bfs(post_adjacency, changed, args.radius)
    ordered_global = sorted(selected)
    global_to_local = {global_id: index for index, global_id in enumerate(ordered_global, 1)}
    edges = {
        atom for atom in selected
        if any(neighbor not in selected for neighbor in pre_adjacency[atom] | post_adjacency[atom])
    }

    initiator_global = None
    if args.initiator:
        initiator_global = tuple(parse_source_ref(value, source_to_global) for value in args.initiator)
        if initiator_global[0] == initiator_global[1]:
            raise ValueError("Initiator atoms must be different")
        if not all(atom in selected for atom in initiator_global):
            raise ValueError(
                f"Selected initiators {args.initiator} are outside the reduced reaction template")

    categories = ("atom", "bond", "angle", "dihedral", "improper")
    offsets = []
    totals = {kind: 0 for kind in categories}
    pre_coeffs = {section: [] for section in base.COEFF_SECTIONS}
    for data in monomers:
        offsets.append(dict(totals))
        for section in base.COEFF_SECTIONS:
            category = base.SECTION_KIND[section]
            for fields, comment in data.records(section):
                updated = fields.copy()
                updated[0] = str(int(updated[0]) + totals[category])
                pre_coeffs[section].append(base.format_record(updated, comment))
        for kind in totals:
            totals[kind] += data.type_counts.get(kind, 0)

    pre_atom_lines = []
    for global_id in ordered_global:
        monomer_number, source_atom = global_to_source[global_id]
        fields, comment = monomers[monomer_number - 1].atom_records()[source_atom]
        product_fields, _product_comment = product.atom_records()[global_id]
        updated = fields.copy()
        updated[0] = str(global_to_local[global_id])
        updated[1] = str(monomer_number)
        updated[2] = str(int(updated[2]) + offsets[monomer_number - 1]["atom"])
        updated[4:7] = product_fields[4:7]
        pre_atom_lines.append(base.format_record(updated, comment))

    post_atom_lines = []
    for global_id in ordered_global:
        fields, comment = product.atom_records()[global_id]
        updated = fields.copy()
        updated[0] = str(global_to_local[global_id])
        post_atom_lines.append(base.format_record(updated, comment))

    pre_topology = {section: [] for section in base.TOPOLOGY}
    for monomer_number, data in enumerate(monomers, 1):
        for section, (kind, _nrefs) in base.TOPOLOGY.items():
            for _record, old_type, refs in data.topology(section):
                global_refs = tuple(source_to_global[(monomer_number, ref)] for ref in refs)
                if all(ref in selected for ref in global_refs):
                    pre_topology[section].append((
                        old_type + offsets[monomer_number - 1][kind],
                        tuple(global_to_local[ref] for ref in global_refs),
                    ))
    post_topology = base.build_post_topology(product, global_to_local)

    pre_out = args.pre_output.expanduser().resolve()
    post_out = args.post_output.expanduser().resolve()
    map_out = args.atom_map.expanduser().resolve()
    report_out = args.report.expanduser().resolve()
    metadata_out = args.metadata.expanduser().resolve() if args.metadata else None
    output_paths = (pre_out, post_out, map_out, report_out)
    for path in output_paths + ((metadata_out,) if metadata_out else ()):
        path.parent.mkdir(parents=True, exist_ok=True)
    base.write_data(pre_out, "LAMMPS data: graph-derived pre-reaction template",
                    totals, pre_coeffs, pre_atom_lines, pre_topology)
    post_types = {kind: product.type_counts[kind] for kind in categories}
    base.write_data(post_out, "LAMMPS data: graph-derived post-reaction template",
                    post_types, base.copy_coefficients(product), post_atom_lines, post_topology)

    removed_local = sorted(tuple(sorted(global_to_local[a] for a in pair)) for pair in removed_global)
    added_local = sorted(tuple(sorted(global_to_local[a] for a in pair)) for pair in added_global)
    with map_out.open("w", encoding="utf-8") as handle:
        handle.write("template_id\tmonomer\tsource_atom_id\tproduct_atom_id\troles\n")
        for global_id in ordered_global:
            monomer_number, source_atom = global_to_source[global_id]
            roles = []
            if global_id in changed:
                roles.append("topology_change")
            if global_id in edges:
                roles.append("edge")
            if initiator_global and global_id in initiator_global:
                roles.append("initiator")
            handle.write(
                f"{global_to_local[global_id]}\t{monomer_number}\t{source_atom}\t"
                f"{global_id}\t{','.join(roles) or 'context'}\n"
            )

    candidate_text = [
        f"{global_to_source[a][0]}:{global_to_source[a][1]} -- "
        f"{global_to_source[b][0]}:{global_to_source[b][1]} "
        f"(template {global_to_local[a]}-{global_to_local[b]})"
        for a, b in added_global
    ]
    candidates = [
        {
            "source_refs": [
                f"{global_to_source[a][0]}:{global_to_source[a][1]}",
                f"{global_to_source[b][0]}:{global_to_source[b][1]}",
            ],
            "template_ids": [global_to_local[a], global_to_local[b]],
            "product_atom_ids": [a, b],
        }
        for a, b in added_global
    ]
    report_out.write_text(
        "# Graph-derived reaction template report\n\n"
        f"- Monomers: {len(monomers)}\n"
        f"- Template radius: {args.radius} bonds\n"
        f"- Selected atoms: {len(selected)}\n"
        f"- Edge atoms: {', '.join(map(str, sorted(global_to_local[a] for a in edges))) or 'none'}\n"
        f"- Removed template bonds: {removed_local}\n"
        f"- Added template bonds: {added_local}\n"
        f"- Initiator candidates: {candidate_text}\n"
        f"- Selected initiators: {args.initiator or 'not specified'}\n",
        encoding="utf-8",
    )
    if metadata_out:
        metadata_out.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "monomer_count": len(monomers),
                    "template_radius": args.radius,
                    "initiator_candidates": candidates,
                    "selected_source_refs": list(args.initiator or []),
                    "removed_template_bonds": [list(pair) for pair in removed_local],
                    "added_template_bonds": [list(pair) for pair in added_local],
                },
                indent=2,
            ) + "\n",
            encoding="utf-8",
        )
    print(f"Wrote {pre_out}")
    print(f"Wrote {post_out}")
    print(f"Wrote {map_out}")
    print(f"Wrote {report_out}")
    if metadata_out:
        print(f"Wrote {metadata_out}")
    print(f"Removed template bonds: {removed_local}")
    print(f"Added template bonds: {added_local}")
    print(f"Initiator candidates: {candidate_text}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
