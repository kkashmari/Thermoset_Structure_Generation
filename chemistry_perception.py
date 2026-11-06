#!/usr/bin/env python3
"""Use RDKit and an external SMARTS library to report reactive functional groups."""

from __future__ import annotations

import argparse
import hashlib
import html
import itertools
import json
from collections import deque
from pathlib import Path

import reaction_history

try:
    from rdkit import Chem, RDConfig, rdBase
except ImportError as exc:
    raise SystemExit(
        "RDKit is unavailable. Launch with the bundled .rdkit_env Python or install RDKit."
    ) from exc


def mol2_atom_ids(path: Path):
    ids = []
    in_atoms = False
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.strip() == "@<TRIPOS>ATOM"
            continue
        if in_atoms and line.strip():
            fields = line.split()
            if fields and fields[0].isdigit():
                ids.append(int(fields[0]))
    return ids


def read_molecule(path: Path):
    if path.suffix.lower() != ".mol2":
        raise ValueError(
            f"{path}: RDKit chemistry perception requires MOL2 input with elements and bond orders")
    molecule = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
    if molecule is None:
        raise ValueError(f"RDKit could not parse {path}")
    atom_ids = mol2_atom_ids(path)
    if len(atom_ids) != molecule.GetNumAtoms():
        raise ValueError(f"{path}: MOL2 atom table and RDKit atom count do not agree")
    for atom, source_id in zip(molecule.GetAtoms(), atom_ids):
        atom.SetIntProp("source_atom_id", source_id)
    failure = Chem.SanitizeMol(molecule, catchErrors=True)
    if failure != Chem.SanitizeFlags.SANITIZE_NONE:
        raise ValueError(f"{path}: RDKit sanitization failed at operation {failure}")
    return molecule


def load_library(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    groups = {}
    for entry in data.get("functional_groups", []):
        group_id = entry["id"]
        if group_id in groups:
            raise ValueError(f"Duplicate functional-group ID: {group_id}")
        query = Chem.MolFromSmarts(entry["smarts"])
        if query is None:
            raise ValueError(f"Invalid SMARTS for {group_id}: {entry['smarts']}")
        supported_selection_methods = {
            "minimum_heavy_atom_degree", "maximum_heavy_atom_degree",
            "minimum_total_hydrogens", "maximum_total_hydrogens",
        }
        for role, policy in entry.get("role_selection", {}).items():
            if role not in entry.get("reactive_query_atoms", {}):
                raise ValueError(
                    f"Functional group {group_id} selects unknown reactive role {role}")
            method = policy.get("method")
            if method not in supported_selection_methods:
                raise ValueError(
                    f"Functional group {group_id} has unsupported role-selection "
                    f"method {method!r}")
        groups[group_id] = (entry, query)
    for rule in data.get("reaction_rules", []):
        for participant in rule.get("participants", []):
            if participant["group"] not in groups:
                raise ValueError(
                    f"Reaction rule {rule['id']} references unknown group {participant['group']}")
    return data, groups


def _rdkit_hierarchy_entries():
    """Load RDKit's maintained functional-group hierarchy as detector entries."""
    hierarchy_path = Path(RDConfig.RDDataDir) / "Functional_Group_Hierarchy.txt"
    if not hierarchy_path.is_file():
        raise FileNotFoundError(
            f"RDKit functional-group hierarchy was not found: {hierarchy_path}")
    entries = []
    for line_number, line in enumerate(
            hierarchy_path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("//"):
            continue
        fields = [field.strip() for field in stripped.split("\t") if field.strip()]
        if len(fields) < 2:
            raise ValueError(
                f"Malformed RDKit hierarchy entry at {hierarchy_path}:{line_number}")
        hierarchy_id, smarts = fields[:2]
        label = fields[2] if len(fields) > 2 else hierarchy_id.rsplit(".", 1)[-1]
        entries.append({
            "id": f"rdkit:{hierarchy_id}",
            "name": label,
            "category": hierarchy_id.split(".", 1)[0],
            "hierarchy": hierarchy_id,
            "smarts": smarts,
            "source": "RDKit Functional Group Hierarchy",
        })
    return entries


def load_detection_catalog(path: Path | None):
    """Load broad detectors independently of reaction-capable library groups."""
    if path is None:
        return [], None
    data = json.loads(path.read_text(encoding="utf-8"))
    entries = []
    for source in data.get("sources", []):
        if not source.get("enabled", True):
            continue
        source_type = source.get("type")
        if source_type == "rdkit_functional_group_hierarchy":
            entries.extend(_rdkit_hierarchy_entries())
        else:
            raise ValueError(f"Unsupported functional-group catalog source: {source_type!r}")
    for entry in data.get("functional_groups", []):
        item = dict(entry)
        item.setdefault("source", path.name)
        item.setdefault("category", "Custom")
        entries.append(item)
    compiled = []
    seen = set()
    for entry in entries:
        detector_id = entry["id"]
        if detector_id in seen:
            raise ValueError(f"Duplicate detector ID in functional-group catalog: {detector_id}")
        seen.add(detector_id)
        query = Chem.MolFromSmarts(entry["smarts"])
        if query is None:
            raise ValueError(
                f"Invalid detector SMARTS for {detector_id}: {entry['smarts']}")
        compiled.append((entry, query))
    return compiled, data


def detect_catalog_groups(molecule, catalog_queries):
    detections = []
    for entry, query in catalog_queries:
        atom_matches = []
        for atom_indices in molecule.GetSubstructMatches(query, uniquify=True):
            atom_matches.append([
                molecule.GetAtomWithIdx(index).GetIntProp("source_atom_id")
                for index in atom_indices
            ])
        if atom_matches:
            detections.append({
                "id": entry["id"],
                "name": entry["name"],
                "category": entry.get("category", "Uncategorized"),
                "hierarchy": entry.get("hierarchy"),
                "source": entry.get("source", "external catalog"),
                "match_count": len(atom_matches),
                "atom_matches": atom_matches,
                "reaction_support": "detection_only",
            })
    return detections


def apply_role_selection(molecule, atom_ids, policy):
    """Filter role candidates using a library-configured connectivity criterion."""
    if len(atom_ids) < 2:
        return list(atom_ids), []
    atoms_by_source_id = {
        atom.GetIntProp("source_atom_id"): atom for atom in molecule.GetAtoms()
    }
    method = policy["method"]
    scored = []
    for atom_id in atom_ids:
        atom = atoms_by_source_id[atom_id]
        if method.endswith("heavy_atom_degree"):
            score = sum(neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors())
        elif method.endswith("total_hydrogens"):
            score = atom.GetTotalNumHs(includeNeighbors=True)
        else:  # Guarded by load_library validation.
            raise ValueError(f"Unsupported role-selection method {method!r}")
        scored.append((atom_id, score))
    target = (min if method.startswith("minimum") else max)(
        score for _atom_id, score in scored)
    selected = [atom_id for atom_id, score in scored if score == target]
    audit = [{"atom_id": atom_id, "score": score, "selected": score == target}
             for atom_id, score in scored]
    return selected, audit


def combined_molecule(molecules):
    combined = None
    ref_to_index = {}
    offset = 0
    for monomer_number in sorted(molecules):
        molecule = molecules[monomer_number]
        combined = Chem.Mol(molecule) if combined is None else Chem.CombineMols(combined, molecule)
        for atom in molecule.GetAtoms():
            source_id = atom.GetIntProp("source_atom_id")
            ref_to_index[f"{monomer_number}:{source_id}"] = offset + atom.GetIdx()
        offset += molecule.GetNumAtoms()
    return combined, ref_to_index


def neighborhood(molecule, centers, radius):
    distances = {center: 0 for center in centers}
    queue = deque(centers)
    while queue:
        atom_index = queue.popleft()
        if distances[atom_index] >= radius:
            continue
        for neighbor in molecule.GetAtomWithIdx(atom_index).GetNeighbors():
            neighbor_index = neighbor.GetIdx()
            if neighbor_index not in distances:
                distances[neighbor_index] = distances[atom_index] + 1
                queue.append(neighbor_index)
    return set(distances)


def fragment_signature(molecule, atom_indices, role_indices):
    marked = Chem.Mol(molecule)
    for atom in marked.GetAtoms():
        atom.SetAtomMapNum(0)
    for map_number, role in enumerate(sorted(role_indices), 1):
        marked.GetAtomWithIdx(role_indices[role]).SetAtomMapNum(map_number)
    return Chem.MolFragmentToSmiles(
        marked, atomsToUse=sorted(atom_indices), canonical=True,
        isomericSmiles=True, allBondsExplicit=True, allHsExplicit=True)


def candidate_topology_identity(molecules, selected, changes, rule_id, radius=3):
    pre, ref_to_index = combined_molecule(molecules)
    post = Chem.RWMol(pre)
    removed = []
    added = []
    for action in changes.get("remove", []):
        atom1, atom2 = (selected[role] for role in action["atoms"])
        index1, index2 = ref_to_index[atom1], ref_to_index[atom2]
        post.RemoveBond(index1, index2)
        removed.append([atom1, atom2])
    bond_types = {
        "1": Chem.BondType.SINGLE, "2": Chem.BondType.DOUBLE,
        "3": Chem.BondType.TRIPLE, "ar": Chem.BondType.AROMATIC,
        "am": Chem.BondType.SINGLE,
    }
    for action in changes.get("add", []):
        atom1, atom2 = (selected[role] for role in action["atoms"])
        order = str(action.get("order", "1"))
        if order not in bond_types:
            raise ValueError(f"Unsupported library bond order {order!r}")
        post.AddBond(ref_to_index[atom1], ref_to_index[atom2], bond_types[order])
        added.append([atom1, atom2, order])
    order_changes = []
    for action in changes.get("change_order", []):
        atom1, atom2 = (selected[role] for role in action["atoms"])
        order = str(action["order"])
        if order not in bond_types:
            raise ValueError(f"Unsupported library bond order {order!r}")
        bond = post.GetBondBetweenAtoms(ref_to_index[atom1], ref_to_index[atom2])
        if bond is None:
            raise ValueError(f"Cannot change order of missing predicted bond {atom1}-{atom2}")
        bond.SetBondType(bond_types[order])
        order_changes.append([atom1, atom2, order])
    formal_charge_changes = []
    for action in changes.get("formal_charge", []):
        atom_ref = selected[action["atoms"][0]]
        charge = int(action["charge"])
        post.GetAtomWithIdx(ref_to_index[atom_ref]).SetFormalCharge(charge)
        formal_charge_changes.append([atom_ref, charge])
    post = post.GetMol()
    failure = Chem.SanitizeMol(post, catchErrors=True)
    if failure != Chem.SanitizeFlags.SANITIZE_NONE:
        raise ValueError(f"Predicted product failed RDKit sanitization at operation {failure}")
    role_indices = {role: ref_to_index[atom_ref] for role, atom_ref in selected.items()}
    centers = set(role_indices.values())
    selected_atoms = neighborhood(pre, centers, radius) | neighborhood(post, centers, radius)
    pre_signature = fragment_signature(pre, selected_atoms, role_indices)
    post_signature = fragment_signature(post, selected_atoms, role_indices)
    payload = {
        "rule_id": rule_id,
        "radius": radius,
        "pre_signature": pre_signature,
        "post_signature": post_signature,
    }
    key = hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()[:20]
    return {
        "preliminary_key": key,
        "radius": radius,
        "pre_signature": pre_signature,
        "post_signature": post_signature,
        "removed_bonds": removed,
        "added_bonds": added,
        "bond_order_changes": order_changes,
        "formal_charge_changes": formal_charge_changes,
    }


def perceive(paths, library_path, topology_radius=3, catalog_path=None):
    library, group_queries = load_library(library_path)
    catalog_queries, _catalog = load_detection_catalog(catalog_path)
    input_histories = {
        number: reaction_history.load_atom_history(path)
        for number, path in enumerate(paths, 1)
    }
    monomers = []
    rdkit_molecules = {}
    all_matches = {}
    existing_bonds = set()
    existing_bond_orders = {}
    for monomer_number, path in enumerate(paths, 1):
        molecule = read_molecule(path)
        rdkit_molecules[monomer_number] = molecule
        atoms = []
        for atom in molecule.GetAtoms():
            atoms.append({
                "source_atom_id": atom.GetIntProp("source_atom_id"),
                "element": atom.GetSymbol(),
                "formal_charge": atom.GetFormalCharge(),
                "degree": atom.GetDegree(),
                "total_hydrogens": atom.GetTotalNumHs(includeNeighbors=True),
                "aromatic": atom.GetIsAromatic(),
                "hybridization": str(atom.GetHybridization()),
            })
        for bond in molecule.GetBonds():
            atom1 = bond.GetBeginAtom().GetIntProp("source_atom_id")
            atom2 = bond.GetEndAtom().GetIntProp("source_atom_id")
            pair = frozenset((f"{monomer_number}:{atom1}", f"{monomer_number}:{atom2}"))
            existing_bonds.add(pair)
            order_value = bond.GetBondTypeAsDouble()
            existing_bond_orders[pair] = (
                "ar" if bond.GetIsAromatic() else str(int(order_value))
                if order_value.is_integer() else str(order_value))
        matches = []
        catalog_detections = detect_catalog_groups(molecule, catalog_queries)
        for group_id, (entry, query) in group_queries.items():
            for match_number, atom_indices in enumerate(
                    molecule.GetSubstructMatches(query, uniquify=True), 1):
                atom_refs = [molecule.GetAtomWithIdx(index).GetIntProp("source_atom_id")
                             for index in atom_indices]
                reactive = {}
                selection_audit = {}
                for role, query_positions in entry.get("reactive_query_atoms", {}).items():
                    reactive[role] = [atom_refs[position] for position in query_positions]
                for role, specification in entry.get("neighbor_roles", {}).items():
                    neighbors = []
                    for anchor_id in reactive.get(specification["anchor_role"], []):
                        anchor_index = atom_refs.index(anchor_id)
                        anchor_atom = molecule.GetAtomWithIdx(atom_indices[anchor_index])
                        for neighbor in anchor_atom.GetNeighbors():
                            if neighbor.GetSymbol() == specification["element"]:
                                neighbors.append(neighbor.GetIntProp("source_atom_id"))
                    reactive[role] = sorted(set(neighbors))
                for role, policy in entry.get("role_selection", {}).items():
                    reactive[role], scored = apply_role_selection(
                        molecule, reactive.get(role, []), policy)
                    selection_audit[role] = {
                        "method": policy["method"],
                        "candidates": scored,
                        "selected_atom_ids": reactive[role],
                    }
                match = {
                    "group_id": group_id,
                    "group_name": entry["name"],
                    "match_number": match_number,
                    "atom_ids": atom_refs,
                    "reactive_atoms": reactive,
                    "role_selection": selection_audit,
                    "monomer": monomer_number,
                }
                matches.append(match)
                all_matches.setdefault(group_id, []).append(match)
        supported_detections = []
        for group_id, (entry, _query) in group_queries.items():
            group_matches = [
                match for match in matches if match["group_id"] == group_id
            ]
            if group_matches:
                supported_detections.append({
                    "id": f"reaction_library:{group_id}",
                    "name": entry["name"],
                    "category": "Reaction library",
                    "hierarchy": None,
                    "source": library_path.name,
                    "match_count": len(group_matches),
                    "atom_matches": [match["atom_ids"] for match in group_matches],
                    "reaction_support": "supported",
                })
        monomers.append({
            "number": monomer_number,
            "path": str(path),
            "name": path.name,
            "atom_count": molecule.GetNumAtoms(),
            "bond_count": molecule.GetNumBonds(),
            "formal_charge": sum(atom.GetFormalCharge() for atom in molecule.GetAtoms()),
            "matches": matches,
            "detected_functional_groups": catalog_detections + supported_detections,
            "atoms": atoms,
        })

    candidates = []
    for rule in library.get("reaction_rules", []):
        pools = [all_matches.get(item["group"], []) for item in rule.get("participants", [])]
        if not pools or any(not pool for pool in pools):
            continue
        for combination in itertools.product(*pools):
            monomer_numbers = [match["monomer"] for match in combination]
            if not rule.get("allow_intramolecular", False) and len(set(monomer_numbers)) < len(combination):
                continue
            sites = []
            role_options = {}
            valid = True
            for participant, match in zip(rule["participants"], combination):
                role_map = participant.get("role_map")
                if role_map:
                    mapped_roles = role_map
                    missing = [role for role in mapped_roles
                               if not match["reactive_atoms"].get(role)]
                    if missing:
                        valid = False
                        break
                    site_role = participant.get("site_role", next(iter(mapped_roles)))
                    site_atoms = match["reactive_atoms"][site_role]
                    sites.append({
                        "role": mapped_roles[site_role],
                        "group": match["group_name"],
                        "monomer": match["monomer"],
                        "match_number": match["match_number"],
                        "atom_options": [f"{match['monomer']}:{atom_id}" for atom_id in site_atoms],
                    })
                    available_items = (
                        (mapped_roles[group_role], match["reactive_atoms"][group_role])
                        for group_role in mapped_roles)
                else:
                    role = participant["role"]
                    role_atoms = match["reactive_atoms"].get(role, [])
                    if not role_atoms:
                        valid = False
                        break
                    sites.append({
                        "role": role,
                        "group": match["group_name"],
                        "monomer": match["monomer"],
                        "match_number": match["match_number"],
                        "atom_options": [f"{match['monomer']}:{atom_id}" for atom_id in role_atoms],
                    })
                    available_items = match["reactive_atoms"].items()
                for available_role, atom_ids in available_items:
                    refs = [f"{match['monomer']}:{atom_id}" for atom_id in atom_ids]
                    if available_role in role_options:
                        role_options[available_role] = sorted(
                            set(role_options[available_role]) & set(refs))
                    else:
                        role_options[available_role] = refs
            changes = rule.get("topology_changes", {})
            required_roles = set(rule.get("initiator_roles", []))
            for action in (changes.get("remove", []) + changes.get("add", [])
                           + changes.get("change_order", []) + changes.get("retype", [])
                           + changes.get("formal_charge", [])):
                required_roles.update(action["atoms"])
            if any(not role_options.get(role) for role in required_roles):
                valid = False
            if valid:
                ordered_roles = sorted(required_roles)
                option_sets = [role_options[role] for role in ordered_roles]
                for selected_values in itertools.product(*option_sets):
                    selected = dict(zip(ordered_roles, selected_values))
                    if len(set(selected.values())) != len(selected):
                        continue
                    commands = []
                    action_summary = []
                    action_valid = True
                    for action in changes.get("remove", []):
                        atom1, atom2 = (selected[role] for role in action["atoms"])
                        if frozenset((atom1, atom2)) not in existing_bonds:
                            action_valid = False
                            break
                        commands.append(f"remove {atom1} {atom2}")
                        action_summary.append(f"remove {atom1}–{atom2}")
                    if not action_valid:
                        continue
                    for action in changes.get("add", []):
                        atom1, atom2 = (selected[role] for role in action["atoms"])
                        if frozenset((atom1, atom2)) in existing_bonds:
                            action_valid = False
                            break
                        order = str(action.get("order", "1"))
                        commands.append(f"add {atom1} {atom2} {order}")
                        action_summary.append(f"add {atom1}–{atom2} ({order})")
                    if not action_valid:
                        continue
                    for action in changes.get("change_order", []):
                        atom1, atom2 = (selected[role] for role in action["atoms"])
                        pair = frozenset((atom1, atom2))
                        if pair not in existing_bonds:
                            action_valid = False
                            break
                        expected = str(action.get("from", ""))
                        if expected and existing_bond_orders.get(pair) != expected:
                            action_valid = False
                            break
                        order = str(action["order"])
                        commands.append(f"order {atom1} {atom2} {order}")
                        action_summary.append(
                            f"change {atom1}–{atom2} order "
                            f"{existing_bond_orders.get(pair)}→{order}")
                    if not action_valid:
                        continue
                    for action in changes.get("retype", []):
                        atom_ref = selected[action["atoms"][0]]
                        atom_type = str(action["atom_type"])
                        commands.append(f"type {atom_ref} {atom_type}")
                        action_summary.append(f"change {atom_ref} atom type to {atom_type}")
                    for action in changes.get("formal_charge", []):
                        atom_ref = selected[action["atoms"][0]]
                        action_summary.append(
                            f"set predicted formal charge of {atom_ref} to {int(action['charge']):+d}")
                    initiators = [selected[role] for role in rule.get("initiator_roles", [])]
                    identity = candidate_topology_identity(
                        rdkit_molecules, selected, changes, rule["id"], topology_radius)
                    candidate = {
                        "rule_id": rule["id"],
                        "rule_name": rule["name"],
                        "description": rule.get("description", ""),
                        "sites": sites,
                        "role_atoms": selected,
                        "initiator_source_refs": initiators,
                        "commands": commands,
                        "action_summary": action_summary,
                        "topology_identity": identity,
                    }
                    candidate["prior_reaction_participation"] = (
                        reaction_history.candidate_prior_participation(
                            candidate, input_histories))
                    candidate["history_priority"] = (
                        reaction_history.candidate_history_priority(candidate))
                    candidate["_generation_order"] = len(candidates)
                    candidates.append(candidate)
                    if len(candidates) >= 500:
                        break
            if len(candidates) >= 500:
                break
    candidates.sort(key=lambda candidate: (
        -candidate["history_priority"]["previously_reacted_initiator_count"],
        -candidate["history_priority"]["previously_reacted_role_atom_count"],
        -candidate["history_priority"]["prior_event_count"],
        candidate["_generation_order"],
    ))
    for rank, candidate in enumerate(candidates, 1):
        candidate.pop("_generation_order", None)
        candidate["history_priority"]["rank"] = rank
    return {
        "schema_version": 1,
        "rdkit_version": rdBase.rdkitVersion,
        "library": str(library_path),
        "functional_group_catalog": str(catalog_path) if catalog_path else None,
        "topology_radius": topology_radius,
        "monomers": monomers,
        "reaction_candidates": candidates,
    }


def render_html(result):
    sections = []
    for monomer in result["monomers"]:
        detected_rows = []
        for detection in monomer.get("detected_functional_groups", []):
            atom_text = "; ".join(
                ", ".join(map(str, atom_ids))
                for atom_ids in detection["atom_matches"])
            support = ("Reaction suggestions available"
                       if detection["reaction_support"] == "supported"
                       else "Detected only")
            detected_rows.append(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
                    html.escape(detection["name"]),
                    html.escape(detection["category"]),
                    detection["match_count"], html.escape(atom_text), support))
        rows = []
        for match in monomer["matches"]:
            reactive = "; ".join(
                f"{html.escape(role)}: " + ", ".join(
                    f"{match['monomer']}:{atom}" for atom in atom_ids)
                for role, atom_ids in match["reactive_atoms"].items())
            rows.append(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
                    html.escape(match["group_name"]), match["match_number"],
                    ", ".join(map(str, match["atom_ids"])), reactive or "—"))
        sections.append(
            f"<section><h2>Monomer {monomer['number']}: {html.escape(monomer['name'])}</h2>"
            f"<p>{monomer['atom_count']} atoms · {monomer['bond_count']} bonds · "
            f"formal charge {monomer['formal_charge']:+d}</p>"
            "<h3>Broad functional-group detection</h3>"
            "<p>These detections describe molecular features. They produce a reaction "
            "suggestion only when a validated reaction-library rule also matches.</p>"
            "<table><thead><tr><th>Detected group</th><th>Category</th>"
            "<th>Matches</th><th>Atoms in each match</th><th>Status</th></tr></thead><tbody>"
            + ("".join(detected_rows) if detected_rows else
               "<tr><td colspan='5'>No catalog matches</td></tr>")
            + "</tbody></table><h3>Reaction-supported functional groups</h3>"
            "<table><thead><tr><th>Functional group</th><th>Match</th>"
            "<th>Atoms in group (not necessarily reactive)</th>"
            "<th>Selected reactive roles</th></tr></thead><tbody>"
            + ("".join(rows) if rows else "<tr><td colspan='4'>No library matches</td></tr>")
            + "</tbody></table></section>")
    candidate_rows = []
    for index, candidate in enumerate(result["reaction_candidates"], 1):
        site_text = "<br>".join(
            f"{html.escape(role)}: {html.escape(atom_ref)}"
            for role, atom_ref in candidate["role_atoms"].items())
        changes = "<br>".join(html.escape(item) for item in candidate["action_summary"])
        prior = candidate.get("prior_reaction_participation", [])
        if prior:
            history_parts = []
            for participant in prior:
                events = "; ".join(
                    f"{event.get('reaction', 'earlier reaction')}"
                    + (f" as {', '.join(event.get('roles', []))}"
                       if event.get("roles") else "")
                    for event in participant["prior_events"])
                history_parts.append(
                    f"<strong>{html.escape(participant['atom_ref'])}</strong>: "
                    f"{html.escape(events)}")
            history_text = "<br>".join(history_parts)
        else:
            history_text = "None of these proposed reaction atoms has reacted before."
        priority = candidate.get("history_priority", {})
        priority_label = (
            "Prioritized: continues from a previously reacted atom"
            if priority.get("is_prioritized") else "No earlier reacted atom")
        candidate_rows.append(
            f"<tr><td>{index}</td><td>{html.escape(candidate['rule_name'])}<br>"
            f"<small>{html.escape(priority_label)}</small></td>"
            f"<td>{site_text}</td><td>{changes}<br><small>"
            f"{html.escape(candidate['description'])}</small></td>"
            f"<td>{history_text}</td></tr>")
    return """<!doctype html><html><head><meta charset="utf-8"><title>RDKit chemistry perception</title>
<style>body{font:15px system-ui;margin:28px;background:#f5f7fb;color:#172033}section{background:white;padding:18px;margin:16px 0;border-radius:10px}table{border-collapse:collapse;width:100%%}th,td{padding:9px;border-bottom:1px solid #d8dee9;text-align:left;vertical-align:top}th{background:#eaf1f8}.note{background:#fff7d6;padding:12px;border-left:4px solid #d39b00}</style></head><body>
<h1>RDKit chemistry-perception report</h1>
<p class="note">RDKit detected the functional groups and proposed atom-specific reaction sites and bond changes. Select and approve the complete change in the reaction editor. Template reuse is decided only after force-field typing by comparing the complete pre/post topology.</p>
<p>RDKit version: %s</p>%s
<section><h2>Compatible reaction candidates</h2><table><thead><tr><th>#</th><th>Rule</th><th>Selected role atoms</th><th>Proposed bond changes</th><th>Earlier reaction participation</th></tr></thead><tbody>%s</tbody></table></section>
</body></html>""" % (
        html.escape(result["rdkit_version"]), "".join(sections),
        "".join(candidate_rows) if candidate_rows else "<tr><td colspan='5'>No compatible library rule found</td></tr>")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--monomer", action="append", required=True, type=Path)
    parser.add_argument("--library", required=True, type=Path)
    parser.add_argument("--catalog", type=Path,
                        help="Broad functional-group detector catalog")
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--topology-radius", type=int, default=3)
    args = parser.parse_args(argv)
    paths = [path.expanduser().resolve() for path in args.monomer]
    library_path = args.library.expanduser().resolve()
    catalog_path = args.catalog.expanduser().resolve() if args.catalog else None
    if args.topology_radius < 0:
        raise ValueError("Topology radius must be non-negative")
    result = perceive(paths, library_path, args.topology_radius, catalog_path)
    json_path = args.json.expanduser().resolve()
    html_path = args.html.expanduser().resolve()
    json_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    html_path.write_text(render_html(result), encoding="utf-8")
    print(f"RDKit {result['rdkit_version']} sanitized {len(paths)} monomer(s)")
    for monomer in result["monomers"]:
        print(
            f"Monomer {monomer['number']}: "
            f"{len(monomer.get('detected_functional_groups', []))} detected group type(s); "
            f"{len(monomer['matches'])} reaction-supported match(es)")
    print(f"Compatible reaction candidates: {len(result['reaction_candidates'])}")
    print(f"Wrote {json_path}")
    print(f"Wrote {html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
