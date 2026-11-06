#!/usr/bin/env python3
"""Track atom participation across a sequential reaction series."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


HISTORY_FILENAME = "reaction_atom_history.json"


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


def load_atom_history(monomer_path: Path):
    """Return prior events keyed by atom ID when *monomer_path* is a tracked product."""
    monomer_path = monomer_path.expanduser().resolve()
    history_path = monomer_path.parent / HISTORY_FILENAME
    if not history_path.is_file():
        return {}
    data = json.loads(history_path.read_text(encoding="utf-8"))
    product_names = set(data.get("product_files", []))
    if monomer_path.name not in product_names:
        return {}
    input_ids = set(mol2_atom_ids(monomer_path))
    history = {}
    for atom_text, events in data.get("atoms", {}).items():
        atom_id = int(atom_text)
        if atom_id in input_ids and isinstance(events, list):
            history[atom_id] = events
    return history


def candidate_prior_participation(candidate, input_histories):
    """Describe candidate role atoms that have participated in earlier reactions."""
    roles_by_ref = {}
    for role, atom_ref in candidate.get("role_atoms", {}).items():
        roles_by_ref.setdefault(atom_ref, []).append(role)
    participants = []
    for atom_ref, current_roles in roles_by_ref.items():
        monomer_text, atom_text = atom_ref.split(":", 1)
        events = input_histories.get(int(monomer_text), {}).get(int(atom_text), [])
        if events:
            participants.append({
                "atom_ref": atom_ref,
                "current_roles": sorted(current_roles),
                "prior_events": events,
            })
    return participants


def candidate_history_priority(candidate):
    """Return a generic ranking score for continuity with earlier reactions.

    Previously reacted initiators rank ahead of other previously reacted role
    atoms because initiators define the principal reaction pair. No chemistry,
    material name, or atom ID is encoded here.
    """
    participants = candidate.get("prior_reaction_participation", [])
    prior_refs = {participant["atom_ref"] for participant in participants}
    initiators = set(candidate.get("initiator_source_refs", []))
    prior_event_count = sum(
        len(participant.get("prior_events", [])) for participant in participants)
    return {
        "previously_reacted_initiator_count": len(prior_refs & initiators),
        "previously_reacted_role_atom_count": len(prior_refs),
        "prior_event_count": prior_event_count,
        "is_prioritized": bool(prior_refs),
        "reason": (
            "Includes one or more atoms that participated in an earlier reaction."
            if prior_refs else
            "No proposed reaction atom participated in an earlier reaction."
        ),
    }


def write_product_history(monomer_paths, product_paths, changed_refs, reaction_name,
                          selected_candidate=None):
    """Carry input history to product IDs and append the current reaction event."""
    monomer_paths = [Path(path).expanduser().resolve() for path in monomer_paths]
    product_paths = [Path(path).expanduser().resolve() for path in product_paths]
    if not product_paths:
        raise ValueError("At least one product path is required for reaction history")

    product_history = {}
    offsets = {}
    offset = 0
    for monomer_number, path in enumerate(monomer_paths, 1):
        atom_ids = mol2_atom_ids(path)
        offsets[monomer_number] = offset
        prior = load_atom_history(path)
        for atom_id, events in prior.items():
            product_history[str(offset + atom_id)] = list(events)
        offset += len(atom_ids)

    roles_by_ref = {}
    if selected_candidate:
        for role, atom_ref in selected_candidate.get("role_atoms", {}).items():
            roles_by_ref.setdefault(atom_ref, []).append(role)

    for atom_ref in sorted(set(changed_refs)):
        monomer_text, atom_text = atom_ref.split(":", 1)
        monomer_number = int(monomer_text)
        product_atom_id = offsets[monomer_number] + int(atom_text)
        event = {
            "reaction": reaction_name,
            "source_ref": atom_ref,
            "roles": sorted(roles_by_ref.get(atom_ref, [])),
        }
        product_history.setdefault(str(product_atom_id), []).append(event)

    result = {
        "schema_version": 1,
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "reaction": reaction_name,
        "product_files": [path.name for path in product_paths],
        "atoms": product_history,
    }
    history_path = product_paths[0].parent / HISTORY_FILENAME
    history_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return history_path
