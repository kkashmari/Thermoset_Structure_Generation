#!/usr/bin/env python3
"""Generic parsing and writing helpers for LAMMPS full-style data files."""

from __future__ import annotations

import re
from pathlib import Path


SECTIONS = [
    "Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
    "Dihedral Coeffs", "Improper Coeffs", "BondBond Coeffs",
    "BondAngle Coeffs", "AngleAngle Coeffs", "AngleAngleTorsion Coeffs",
    "EndBondTorsion Coeffs", "MiddleBondTorsion Coeffs",
    "BondBond13 Coeffs", "AngleTorsion Coeffs", "Atoms", "Velocities",
    "Bonds", "Angles", "Dihedrals", "Impropers",
]

COEFF_SECTIONS = [
    "Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
    "Dihedral Coeffs", "Improper Coeffs", "BondBond Coeffs",
    "BondAngle Coeffs", "AngleAngle Coeffs", "AngleAngleTorsion Coeffs",
    "EndBondTorsion Coeffs", "MiddleBondTorsion Coeffs",
    "BondBond13 Coeffs", "AngleTorsion Coeffs",
]

SECTION_KIND = {
    "Masses": "atom", "Pair Coeffs": "atom", "Bond Coeffs": "bond",
    "Angle Coeffs": "angle", "BondBond Coeffs": "angle",
    "BondAngle Coeffs": "angle", "Dihedral Coeffs": "dihedral",
    "AngleAngleTorsion Coeffs": "dihedral", "EndBondTorsion Coeffs": "dihedral",
    "MiddleBondTorsion Coeffs": "dihedral", "BondBond13 Coeffs": "dihedral",
    "AngleTorsion Coeffs": "dihedral", "Improper Coeffs": "improper",
    "AngleAngle Coeffs": "improper",
}

TOPOLOGY = {
    "Bonds": ("bond", 2), "Angles": ("angle", 3),
    "Dihedrals": ("dihedral", 4), "Impropers": ("improper", 4),
}


def split_comment(line):
    if "#" in line:
        body, comment = line.split("#", 1)
        return body.strip(), comment.strip()
    return line.strip(), ""


def section_name(line):
    stripped = line.strip()
    for name in SECTIONS:
        if re.fullmatch(re.escape(name) + r"(?:\s*#.*)?", stripped):
            return name
    return None


class DataFile:
    """Read the sections and count headers needed by the template workflow."""

    def __init__(self, path):
        self.path = Path(path)
        self.sections = {}
        self.counts = {}
        self.type_counts = {}
        current = None
        for line in self.path.read_text(errors="replace").splitlines():
            found = section_name(line)
            if found:
                current = found
                self.sections.setdefault(found, [])
                continue
            if current is None:
                match = re.fullmatch(
                    r"\s*(\d+)\s+(atoms|bonds|angles|dihedrals|impropers)\s*", line)
                if match:
                    self.counts[match.group(2)] = int(match.group(1))
                match = re.fullmatch(
                    r"\s*(\d+)\s+(atom|bond|angle|dihedral|improper) types\s*", line)
                if match:
                    self.type_counts[match.group(2)] = int(match.group(1))
            else:
                self.sections[current].append(line)

    def records(self, section):
        result = []
        for line in self.sections.get(section, []):
            body, comment = split_comment(line)
            if not body or body.startswith("#"):
                continue
            fields = body.split()
            if fields and fields[0].lstrip("+-").isdigit():
                result.append((fields, comment))
        return result

    def atom_records(self):
        return {
            int(fields[0]): (fields, comment)
            for fields, comment in self.records("Atoms")
        }

    def topology(self, section):
        _kind, reference_count = TOPOLOGY[section]
        return [
            (int(fields[0]), int(fields[1]),
             tuple(map(int, fields[2:2 + reference_count])))
            for fields, _comment in self.records(section)
        ]


def format_record(fields, comment=""):
    line = " ".join(fields)
    return line + (f" # {comment}" if comment else "")


def copy_coefficients(data):
    return {
        section: [
            format_record(fields, comment)
            for fields, comment in data.records(section)
        ]
        for section in COEFF_SECTIONS
    }


def build_post_topology(post, global_to_local):
    selected = set(global_to_local)
    output = {}
    for section in TOPOLOGY:
        records = []
        for _old_id, type_id, references in post.topology(section):
            if all(reference in selected for reference in references):
                records.append((
                    type_id,
                    tuple(global_to_local[reference] for reference in references),
                ))
        output[section] = records
    return output


def coordinate_bounds(atom_lines, margin=5.0):
    if not atom_lines:
        raise ValueError("Cannot write a LAMMPS data file with no selected atoms")
    xyz = [
        tuple(map(float, split_comment(line)[0].split()[4:7]))
        for line in atom_lines
    ]
    return [
        (
            min(point[axis] for point in xyz) - margin,
            max(point[axis] for point in xyz) + margin,
        )
        for axis in range(3)
    ]


def write_data(path, title, type_counts, coefficients, atom_lines, topology):
    counts = {"atoms": len(atom_lines)}
    counts.update({section.lower(): len(topology[section]) for section in TOPOLOGY})
    bounds = coordinate_bounds(atom_lines)
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write(title + "\n\n")
        for key in ("atoms", "bonds", "angles", "dihedrals", "impropers"):
            handle.write(f"{counts[key]:8d} {key}\n")
        handle.write("\n")
        for kind in ("atom", "bond", "angle", "dihedral", "improper"):
            handle.write(f"{type_counts.get(kind, 0):8d} {kind} types\n")
        handle.write("\n")
        for (low, high), axis in zip(bounds, "xyz"):
            handle.write(f"{low:18.10f} {high:18.10f} {axis}lo {axis}hi\n")
        for section in COEFF_SECTIONS:
            lines = coefficients.get(section, [])
            if lines:
                handle.write(f"\n{section}\n\n")
                handle.write("\n".join(lines) + "\n")
        handle.write("\nAtoms # full\n\n")
        handle.write("\n".join(atom_lines) + "\n")
        for section in ("Bonds", "Angles", "Dihedrals", "Impropers"):
            handle.write(f"\n{section}\n\n")
            for record_id, (type_id, references) in enumerate(topology[section], 1):
                handle.write(
                    f"{record_id:7d} {type_id:5d} "
                    + " ".join(f"{reference:7d}" for reference in references)
                    + "\n"
                )
