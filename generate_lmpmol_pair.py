#!/usr/bin/env python3
"""Create synchronized LAMMPS bond/react molecule templates and their map file."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import data_to_bond_react_templates as converter


def infer_initiators(pre, post):
    """Use the only new bond, or stop when the graph is ambiguous."""
    added = sorted(converter.normalized_bonds(post) - converter.normalized_bonds(pre))
    if len(added) == 1:
        return added[0]
    raise ValueError(
        "Initiator selection is ambiguous. New template bonds are "
        f"{added}; select the two initiator source atoms in the GUI."
    )


def read_roles(path):
    rows = path.read_text(encoding="utf-8").splitlines()
    if not rows:
        raise ValueError(f"Empty template atom map: {path}")
    headers = rows[0].split("\t")
    required = {"template_id", "roles"}
    if not required <= set(headers):
        raise ValueError(f"Template atom map must contain columns {sorted(required)}")
    result = {}
    for line in rows[1:]:
        if not line.strip():
            continue
        values = dict(zip(headers, line.split("\t")))
        result[int(values["template_id"])] = set(values["roles"].split(","))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pre", type=Path, required=True)
    parser.add_argument("--post", type=Path, required=True)
    parser.add_argument("--atom-map", type=Path, required=True,
                        help="TSV metadata written by generate_reaction_data.py")
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--pre-name")
    parser.add_argument("--post-name")
    parser.add_argument("--map-name", default="reaction.map")
    parser.add_argument("--mapping-name", default="lmpmol_type_mapping.tsv")
    parser.add_argument(
        "--initiators", nargs=2, type=int, metavar=("ID1", "ID2"),
        help="pre-template initiator IDs; normally read from the atom-map TSV",
    )
    args = parser.parse_args(argv)

    pre_path = args.pre.expanduser().resolve()
    post_path = args.post.expanduser().resolve()
    atom_map_path = args.atom_map.expanduser().resolve()
    output_directory = args.output_directory.expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    pre = converter.read_data_file(pre_path)
    post = converter.read_data_file(post_path)
    family = converter.validate_force_field_family([pre, post])
    mapping, _representatives, labels, counts = converter.build_global_registry([pre, post])
    pre_graph = converter.reaction_graph(pre, mapping, labels)
    post_graph = converter.reaction_graph(post, mapping, labels)
    if pre_graph.atoms != post_graph.atoms:
        raise ValueError(
            "This synchronized pair requires identical pre/post atom ID sets; "
            f"pre-only={sorted(pre_graph.atoms - post_graph.atoms)}, "
            f"post-only={sorted(post_graph.atoms - pre_graph.atoms)}"
        )

    roles = read_roles(atom_map_path)
    mapped_initiators = sorted(atom_id for atom_id, values in roles.items()
                               if "initiator" in values)
    if args.initiators:
        initiators = tuple(args.initiators)
    elif len(mapped_initiators) == 2:
        initiators = tuple(mapped_initiators)
    else:
        initiators = infer_initiators(pre_graph, post_graph)
    added = sorted(converter.normalized_bonds(post_graph) - converter.normalized_bonds(pre_graph))
    removed = sorted(converter.normalized_bonds(pre_graph) - converter.normalized_bonds(post_graph))
    if initiators[0] == initiators[1] or not set(initiators) <= pre_graph.atoms:
        raise ValueError(
            f"Initiator pair {initiators} must contain two different pre-template atoms")
    edges = sorted(atom_id for atom_id, values in roles.items() if "edge" in values)
    converter.validate_map_ids(1, pre_graph, post_graph, initiators, edges, [], [])

    pre_output = output_directory / (args.pre_name or pre_path.with_suffix(".lmpmol").name)
    post_output = output_directory / (args.post_name or post_path.with_suffix(".lmpmol").name)
    map_output = output_directory / args.map_name
    mapping_output = output_directory / args.mapping_name
    converter.write_molecule_template(pre, pre_output, mapping)
    converter.write_molecule_template(post, post_output, mapping)
    map_output.write_text(
        converter.render_map_file(
            pre_graph, post_graph, initiators, edges, [], [], None),
        encoding="utf-8",
    )
    converter.write_mapping_report([pre, post], mapping_output, mapping, labels)

    print(f"Force-field family: {family}")
    print("Global type counts: " + ", ".join(f"{key}={value}" for key, value in counts.items()))
    print(f"Wrote {pre_output}")
    print(f"Wrote {post_output}")
    print(f"Wrote {map_output}")
    print(f"Initiator IDs: {initiators[0]}, {initiators[1]}")
    print(f"Edge IDs: {', '.join(map(str, edges)) if edges else 'none'}")
    print(f"Removed bonds: {removed}")
    print(f"Added bonds: {added}")
    print(f"Wrote {mapping_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
