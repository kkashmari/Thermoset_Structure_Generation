# Automated thermoset reaction workflow

This folder contains the automated GUI and generic reaction-template tools. The
active template extractor uses its own material-independent LAMMPS data parser;
it does not import the older DETDA/EPON-specific template generator.

Launch the GUI by double-clicking `launch_reaction_workflow_gui.command`, or run:

```bash
/usr/bin/python3 reaction_workflow_gui.py
```

The launcher intentionally uses macOS `/usr/bin/python3`, because the currently
installed Homebrew Python does not include Tkinter.

Select any combination of stages. Selected stages run from top to bottom:

1. Build/edit the reaction with `react_monomers.py`, create the HTML graphs,
   minimize with Open Babel, validate the actual topology changes, and require
   explicit approval before continuing.
2. Run the existing force-field parameterization tool directly on the minimized
   `.mol2` file and write the selected LAMMPS `.data` file. Choose the `.frc`
   force-field file and charge method in the main workflow window; no second GUI
   must be opened or closed.
3. Generate reduced, genuine LAMMPS `read_data` pre/post reaction files.
4. Convert those reduced files to synchronized `.lmpmol` files and generate the
   matching LAMMPS reaction map using `data_to_bond_react_templates.py`.

Stage 4 writes these files into its selected output directory:

- one `.lmpmol` file named from each selected pre/post `.data` file;
- `reaction.map`;
- `lmpmol_type_mapping.tsv`.

Stage 3 detects bond changes by comparing the combined source-monomer graph with
the parameterized product graph. It extracts all atoms within the user-selected
graph radius and marks an atom as an edge only when it has a neighbor outside
that extracted graph. No functional-group, element, molecule-name, atom-count,
atom-ID, or force-field-type rule is used. If exactly one new bond exists, its
atoms can be used as initiators unambiguously. If multiple new bonds exist, the
GUI shows every graph-derived candidate and waits for you to select one. It does
not choose among multiple candidates automatically.

Stage 2 creates a new artifact directory beside the selected `.data` output on
every run. It retains `audit.json`, the generated NTA, typing scores, the typed
intermediate data, and the force-field-named final data. After verifying the
audit, the workflow copies that final file to the exact **Parameterized data**
path used by Stage 3. Exit status 3 means the parameterizer produced a file but
reported incomplete coefficient coverage; the workflow displays a prominent
warning and retains the audit for review.

During stage 1, enter `remove`, `add`, `type`, and `done` commands in the command
box below the live log. The GUI checks required inputs before each selected stage
and stops immediately if a stage fails.

After Stage 1, the GUI creates `reaction_validation_report.html` and
`reaction_validation_report.json`. The report compares the original monomer
connectivity with the edited product, lists every removed and added bond using
`monomer:atom` references, checks atom mapping and duplicate/self bonds, reports
per-atom coordination changes and charge conservation, and displays the
before/final interactive graphs side by side. The GUI will not start force-field
parameterization until you approve this report. Structural validation cannot
decide whether the proposed chemistry is scientifically intended; that remains
the purpose of the human approval step.
