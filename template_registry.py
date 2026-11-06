#!/usr/bin/env python3
"""Register a force-field-typed reaction template and detect equivalent templates."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import lammps_data_io as data_io


COEFF_FOR_TOPOLOGY = {
    "Bonds": "Bond Coeffs",
    "Angles": "Angle Coeffs",
    "Dihedrals": "Dihedral Coeffs",
    "Impropers": "Improper Coeffs",
}


def read_roles(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"Empty atom map: {path}")
    headers = lines[0].split("\t")
    result = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        values = dict(zip(headers, line.split("\t")))
        result[int(values["template_id"])] = sorted(values.get("roles", "").split(","))
    return result


def type_descriptors(data, section):
    result = {}
    for fields, comment in data.records(section):
        type_id = int(fields[0])
        result[type_id] = comment.strip() or "coeff:" + ",".join(fields[1:])
    return result


def normalized_interaction(section, type_label, colors):
    if section in {"Bonds", "Angles", "Dihedrals"}:
        reverse = tuple(reversed(colors))
        colors = min(tuple(colors), reverse)
    else:
        colors = tuple(sorted(colors))
    return (section, type_label, colors)


def canonical_references(section, references):
    references = tuple(references)
    if section in {"Bonds", "Angles", "Dihedrals"}:
        return min(references, tuple(reversed(references)))
    return tuple(sorted(references))


def topology_delta(pre, post):
    pre_atoms = pre.atom_records()
    post_atoms = post.atom_records()
    pre_atom_types = type_descriptors(pre, "Masses")
    post_atom_types = type_descriptors(post, "Masses")
    atom_changes = []
    for atom_id in sorted(set(pre_atoms) & set(post_atoms)):
        pre_fields = pre_atoms[atom_id][0]
        post_fields = post_atoms[atom_id][0]
        before_type = pre_atom_types.get(int(pre_fields[2]), f"atom_type:{pre_fields[2]}")
        after_type = post_atom_types.get(int(post_fields[2]), f"atom_type:{post_fields[2]}")
        before_charge = round(float(pre_fields[3]), 8)
        after_charge = round(float(post_fields[3]), 8)
        if before_type != after_type or before_charge != after_charge:
            atom_changes.append({
                "atom_id": atom_id,
                "before_type": before_type, "after_type": after_type,
                "before_charge": before_charge, "after_charge": after_charge,
            })
    interaction_changes = {}
    for section in data_io.TOPOLOGY:
        pre_labels = type_descriptors(pre, COEFF_FOR_TOPOLOGY[section])
        post_labels = type_descriptors(post, COEFF_FOR_TOPOLOGY[section])
        pre_set = {
            (canonical_references(section, refs), pre_labels.get(type_id, f"type:{type_id}"))
            for _record, type_id, refs in pre.topology(section)
        }
        post_set = {
            (canonical_references(section, refs), post_labels.get(type_id, f"type:{type_id}"))
            for _record, type_id, refs in post.topology(section)
        }
        interaction_changes[section] = {
            "removed": [{"atoms": list(refs), "type": label}
                        for refs, label in sorted(pre_set - post_set)],
            "added": [{"atoms": list(refs), "type": label}
                      for refs, label in sorted(post_set - pre_set)],
        }
    return {"atom_changes": atom_changes, "interactions": interaction_changes}


def write_report(path, identity):
    recommendation = identity["recommendation"]
    lines = [
        "# Reaction-template topology search",
        "",
        f"- Decision: **{recommendation['status'].replace('_', ' ')}**",
        f"- Reason: {recommendation['reason']}",
        f"- Typed template key: `{identity['typed_template_key']}`",
    ]
    if recommendation["matching_reactions"]:
        lines.append("- Existing matching reactions: " + ", ".join(
            recommendation["matching_reactions"]))
    lines.extend(["", "## Library topology changes", ""])
    lines.extend(f"- {item}" for item in identity.get("topology_changes", []))
    if not identity.get("topology_changes"):
        lines.append("- Manual reaction; see the typed topology differences below.")
    lines.extend(["", "## Atom type and charge changes", ""])
    for change in identity["typed_delta"]["atom_changes"]:
        lines.append(
            f"- Atom {change['atom_id']}: `{change['before_type']}` → "
            f"`{change['after_type']}`; charge {change['before_charge']:+.8f} → "
            f"{change['after_charge']:+.8f}")
    if not identity["typed_delta"]["atom_changes"]:
        lines.append("- None")
    for section, changes in identity["typed_delta"]["interactions"].items():
        lines.extend(["", f"## {section} changes", ""])
        for item in changes["removed"]:
            lines.append(f"- Remove {item['atoms']} type `{item['type']}`")
        for item in changes["added"]:
            lines.append(f"- Add {item['atoms']} type `{item['type']}`")
        if not changes["removed"] and not changes["added"]:
            lines.append("- None")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def typed_topology_signature(data, roles):
    atoms = data.atom_records()
    atom_types = type_descriptors(data, "Masses")
    adjacency = {atom_id: set() for atom_id in atoms}
    for _record, _bond_type, references in data.topology("Bonds"):
        atom1, atom2 = references
        adjacency[atom1].add(atom2)
        adjacency[atom2].add(atom1)

    colors = {}
    for atom_id, (fields, _comment) in atoms.items():
        atom_type = int(fields[2])
        charge = round(float(fields[3]), 8)
        role_text = ",".join(role for role in roles.get(atom_id, []) if role)
        colors[atom_id] = json.dumps(
            [atom_types.get(atom_type, f"atom_type:{atom_type}"), charge, role_text],
            separators=(",", ":"))
    for _iteration in range(6):
        colors = {
            atom_id: hashlib.sha256(json.dumps(
                [colors[atom_id], sorted(colors[neighbor] for neighbor in adjacency[atom_id])],
                separators=(",", ":")).encode("utf-8")).hexdigest()[:24]
            for atom_id in colors
        }

    interactions = []
    for section in data_io.TOPOLOGY:
        labels = type_descriptors(data, COEFF_FOR_TOPOLOGY[section])
        for _record, type_id, references in data.topology(section):
            interactions.append(normalized_interaction(
                section, labels.get(type_id, f"type:{type_id}"),
                [colors[atom_id] for atom_id in references]))
    payload = {
        "atom_colors": sorted(colors.values()),
        "interactions": sorted(interactions),
    }
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return {
        "hash": hashlib.sha256(serialized.encode("utf-8")).hexdigest(),
        "payload": payload,
    }


def load_json(path, default):
    if not path or not path.is_file():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--reaction-folder", required=True, type=Path)
    parser.add_argument("--pre", required=True, type=Path)
    parser.add_argument("--post", required=True, type=Path)
    parser.add_argument("--atom-map", required=True, type=Path)
    parser.add_argument("--selected-candidate", type=Path)
    args = parser.parse_args(argv)

    registry_path = args.registry.expanduser().resolve()
    reaction_folder = args.reaction_folder.expanduser().resolve()
    pre_path = args.pre.expanduser().resolve()
    post_path = args.post.expanduser().resolve()
    atom_map_path = args.atom_map.expanduser().resolve()
    candidate = load_json(
        args.selected_candidate.expanduser().resolve() if args.selected_candidate else None, {})
    registry = load_json(registry_path, {"schema_version": 1, "templates": []})
    registry["templates"] = [
        entry for entry in registry.get("templates", [])
        if entry.get("reaction_folder") != reaction_folder.name
    ]
    roles = read_roles(atom_map_path)
    pre_data = data_io.DataFile(pre_path)
    post_data = data_io.DataFile(post_path)
    pre_signature = typed_topology_signature(pre_data, roles)
    post_signature = typed_topology_signature(post_data, roles)
    typed_key = hashlib.sha256(json.dumps({
        "pre": pre_signature["hash"], "post": post_signature["hash"]
    }, sort_keys=True).encode("utf-8")).hexdigest()
    matches = [entry for entry in registry.get("templates", [])
               if entry.get("typed_template_key") == typed_key]
    recommendation = {
        "status": "reuse_existing_template" if matches else "new_template_required",
        "reason": (
            "The complete force-field-typed pre/post fix bond/react topology is "
            "identical to a registered template."
            if matches else
            "The complete force-field-typed pre/post fix bond/react topology has "
            "not been registered before."
        ),
        "matching_reactions": [entry.get("reaction_folder") for entry in matches],
        "authoritative": True,
    }
    identity = {
        "schema_version": 1,
        "reaction_folder": reaction_folder.name,
        "registered_at": datetime.now(timezone.utc).isoformat(),
        "rule_id": candidate.get("rule_id"),
        "preliminary_key": candidate.get("topology_identity", {}).get("preliminary_key"),
        "typed_template_key": typed_key,
        "pre_typed_hash": pre_signature["hash"],
        "post_typed_hash": post_signature["hash"],
        "topology_changes": candidate.get("action_summary", []),
        "typed_delta": topology_delta(pre_data, post_data),
        "recommendation": recommendation,
        "files": {
            "pre": str(pre_path), "post": str(post_path), "atom_map": str(atom_map_path)
        },
    }
    identity_path = reaction_folder / "template_identity.json"
    identity_path.write_text(json.dumps(identity, indent=2) + "\n", encoding="utf-8")
    report_path = reaction_folder / "topology_search_report.md"
    write_report(report_path, identity)
    registry.setdefault("templates", []).append(identity)
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    registry_path.write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    print(f"Typed template decision: {recommendation['status']}")
    print(recommendation["reason"])
    if recommendation["matching_reactions"]:
        print("Matching reactions: " + ", ".join(recommendation["matching_reactions"]))
    print(f"Wrote {identity_path}")
    print(f"Wrote {report_path}")
    print(f"Updated {registry_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
