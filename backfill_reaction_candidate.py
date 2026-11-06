#!/usr/bin/env python3
"""Recover the selected chemistry candidate for reaction folders made by older versions."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import chemistry_perception


def pairs_from_validation(values):
    result = set()
    for value in values:
        refs = re.findall(r"\d+:\d+", value)
        if len(refs) < 2:
            raise ValueError(f"Cannot read topology change: {value}")
        result.add(tuple(sorted(refs[:2])))
    return result


def pairs_from_commands(commands, command_name):
    result = set()
    for command in commands:
        fields = command.split()
        if fields and fields[0] == command_name:
            result.add(tuple(sorted(fields[1:3])))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reaction-folder", required=True, type=Path)
    parser.add_argument("--library", required=True, type=Path)
    args = parser.parse_args(argv)
    folder = args.reaction_folder.expanduser().resolve()
    old_perception = json.loads(
        (folder / "chemistry_perception_report.json").read_text(encoding="utf-8"))
    validation = json.loads(
        (folder / "reaction_validation_report.json").read_text(encoding="utf-8"))
    paths = [Path(item["path"]) for item in old_perception["monomers"]]
    current = chemistry_perception.perceive(
        paths, args.library.expanduser().resolve(), topology_radius=3)
    removed = pairs_from_validation(validation.get("removed_bonds", []))
    added = pairs_from_validation(validation.get("added_bonds", []))
    matches = [
        candidate for candidate in current["reaction_candidates"]
        if pairs_from_commands(candidate["commands"], "remove") == removed
        and pairs_from_commands(candidate["commands"], "add") == added
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one candidate matching the recorded topology change; found {len(matches)}")
    output = folder / "selected_reaction_candidate.json"
    output.write_text(json.dumps(matches[0], indent=2) + "\n", encoding="utf-8")
    print(f"Recovered {matches[0]['rule_name']}")
    print(f"Preliminary topology key: {matches[0]['topology_identity']['preliminary_key']}")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
