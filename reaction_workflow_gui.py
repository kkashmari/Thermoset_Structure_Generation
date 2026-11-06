#!/usr/bin/env python3
"""workflow for reaction editing and REACTER file generation."""

from __future__ import annotations

import os
import json
import queue
import shutil
import subprocess
import sys
import tempfile
import threading
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk
from tkinter.scrolledtext import ScrolledText

import reaction_history


AUTOMATION_ROOT = Path(__file__).resolve().parent
ROOT = AUTOMATION_ROOT
GUI_PYTHON = sys.executable
WORKFLOW_PYTHON = os.environ.get(
    "REACTION_WORKFLOW_PYTHON",
    "/usr/local/bin/python3" if Path("/usr/local/bin/python3").is_file() else GUI_PYTHON,
)
RDKIT_PYTHON = AUTOMATION_ROOT / ".rdkit_env" / "bin" / "python"
CHEMISTRY_LIBRARY = AUTOMATION_ROOT / "chemistry_library.json"
FUNCTIONAL_GROUP_CATALOG = AUTOMATION_ROOT / "functional_groups_catalog.json"


class ReactionWorkflowGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Thermoset Reaction Workflow")
        self.geometry("1080x880")
        self.minsize(900, 720)
        self.configure(background="#eef3f8")
        self.events = queue.Queue()
        self.process: subprocess.Popen[str] | None = None
        self.worker: threading.Thread | None = None
        self.stop_requested = False
        self.approval_event = threading.Event()
        self.approval_result = False
        self.chemistry_event = threading.Event()
        self.chemistry_choice = None
        self.initiator_event = threading.Event()
        self.initiator_result = None
        self.session_root = None
        self.current_output_dir = None
        self.current_forcefield = None

        self.status_var = tk.StringVar(value="Ready")
        self.forcefield_var = tk.StringVar(value="")
        self.initiator_refs_var = tk.StringVar(value="")
        self.command_var = tk.StringVar()

        self._build_style()
        self._build_ui()
        self._set_default_monomers()
        self.after(100, self._drain_events)
        self.protocol("WM_DELETE_WINDOW", self._close)

    def _build_style(self):
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        style.configure("App.TFrame", background="#eef3f8")
        style.configure("Card.TLabelframe", background="#ffffff", padding=12,
                        borderwidth=1, relief="solid")
        style.configure("Card.TLabelframe.Label", background="#ffffff",
                        foreground="#17324d", font=("Helvetica", 11, "bold"))
        style.configure("Card.TLabel", background="#ffffff", foreground="#334e68")
        style.configure("Hint.TLabel", background="#ffffff", foreground="#6b7f93",
                        font=("Helvetica", 9))
        style.configure("Accent.TButton", background="#1479d1", foreground="white",
                        font=("Helvetica", 11, "bold"), padding=(18, 8))
        style.map("Accent.TButton", background=[("active", "#0d63ad"),
                                                ("disabled", "#9db4c8")])
        style.configure("Stop.TButton", foreground="#9b1c1c", padding=(12, 8))

    def _build_ui(self):
        header = tk.Frame(self, background="#102f4c", height=88)
        header.pack(fill="x")
        header.pack_propagate(False)
        tk.Label(header, text="Thermoset Reaction Workflow",
                 background="#102f4c", foreground="white",
                 font=("Helvetica", 21, "bold")).pack(anchor="w", padx=24, pady=(14, 2))
        tk.Label(header, text="RDKit perception → reaction editing → force field → REACTER templates",
                 background="#102f4c", foreground="#c8d8e8",
                 font=("Helvetica", 10)).pack(anchor="w", padx=24)

        outer = ttk.Frame(self, style="App.TFrame", padding=16)
        outer.pack(fill="both", expand=True)
        canvas = tk.Canvas(outer, background="#eef3f8", highlightthickness=0)
        scroll = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        self.content = ttk.Frame(canvas, style="App.TFrame")
        self.content.bind("<Configure>", lambda _e: canvas.configure(scrollregion=canvas.bbox("all")))
        window = canvas.create_window((0, 0), window=self.content, anchor="nw")
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(window, width=e.width))
        canvas.configure(yscrollcommand=scroll.set)
        canvas.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")

        inputs = ttk.LabelFrame(self.content, text="Reaction inputs", style="Card.TLabelframe")
        inputs.pack(fill="x", pady=5)
        inputs.columnconfigure(1, weight=1)
        ttk.Label(inputs, text="Monomer MOL2 or SMILES", style="Card.TLabel").grid(row=0, column=0, sticky="nw", padx=(0, 8))
        self.monomer_text = tk.Text(inputs, height=4, wrap="none", font=("Menlo", 10))
        self.monomer_text.grid(row=0, column=1, sticky="ew")
        button_column = ttk.Frame(inputs)
        button_column.grid(row=0, column=2, sticky="n", padx=(8, 0))
        ttk.Button(button_column, text="Add MOL2 files…", command=self._add_monomers).pack(fill="x")
        ttk.Button(button_column, text="Add SMILES…", command=self._add_smiles).pack(fill="x", pady=(5, 0))
        ttk.Button(button_column, text="Clear", command=lambda: self.monomer_text.delete("1.0", "end")).pack(fill="x", pady=5)
        self._path_row(inputs, 1, "Force-field .frc", self.forcefield_var,
                       save=False, extension=".frc")
        ttk.Label(inputs, text=(
            "Enter one MOL2 path or one SMILES string per line. Inputs may be mixed. "
            "SMILES are converted to explicit-hydrogen 3D MOL2 structures. "
            "All reports, intermediate files, parameterized monomers, products, templates, "
            "LMPMOL files, and reaction.map are created in a separate Reaction_### folder. "
            "The editable SMARTS library stays beside this GUI."),
            style="Hint.TLabel", wraplength=820).grid(
                row=2, column=1, sticky="w", pady=(6, 0))

        controls = ttk.Frame(self.content, style="App.TFrame")
        controls.pack(fill="x", pady=(10, 6))
        self.run_button = ttk.Button(controls, text="Run complete workflow", style="Accent.TButton", command=self._start)
        self.run_button.pack(side="left")
        self.stop_button = ttk.Button(controls, text="Stop", style="Stop.TButton", command=self._stop, state="disabled")
        self.stop_button.pack(side="left", padx=8)
        ttk.Button(controls, text="Open output folder", command=self._open_folder).pack(side="left")
        ttk.Label(controls, textvariable=self.status_var, style="App.TFrame").pack(side="right")

        terminal = ttk.LabelFrame(self.content, text="Live output and reaction commands", style="Card.TLabelframe")
        terminal.pack(fill="both", expand=True, pady=(5, 0))
        terminal.columnconfigure(0, weight=1)
        terminal.rowconfigure(0, weight=1)
        self.log = ScrolledText(terminal, height=15, wrap="word", background="#0d1b2a",
                                foreground="#e0eaf3", insertbackground="white",
                                font=("Menlo", 10))
        self.log.grid(row=0, column=0, columnspan=2, sticky="nsew")
        entry = ttk.Entry(terminal, textvariable=self.command_var)
        entry.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        entry.bind("<Return>", lambda _e: self._send_command())
        ttk.Button(terminal, text="Send command", command=self._send_command).grid(row=1, column=1, padx=(8, 0), pady=(8, 0))
        ttk.Label(terminal, text=("Examples: remove 1:17 1:18    add 1:17 2:23 1    "
                                  "order 1:13 1:14 2    done"),
                  style="Hint.TLabel").grid(row=2, column=0, columnspan=2, sticky="w", pady=(4, 0))

    def _path_row(self, parent, row, label, variable, save=False, extension=""):
        ttk.Label(parent, text=label, style="Card.TLabel").grid(row=row, column=0, sticky="w", padx=(0, 8), pady=4)
        ttk.Entry(parent, textvariable=variable).grid(row=row, column=1, sticky="ew", pady=4)
        command = (lambda: self._choose_save(variable, extension)) if save else (lambda: self._choose_file(variable, extension))
        ttk.Button(parent, text="Browse…", command=command).grid(row=row, column=2, padx=(8, 0), pady=4)

    def _set_default_monomers(self):
        pass

    def _add_monomers(self):
        paths = filedialog.askopenfilenames(
            initialdir=ROOT, filetypes=[("MOL2 structures", "*.mol2"), ("All files", "*")])
        if paths:
            existing = self.monomer_text.get("1.0", "end").strip()
            self.monomer_text.delete("1.0", "end")
            self.monomer_text.insert("1.0", "\n".join(([existing] if existing else []) + list(paths)))

    def _add_smiles(self):
        value = simpledialog.askstring(
            "Add SMILES",
            "Enter one or more SMILES strings, one per line:", parent=self)
        if not value:
            return
        new_values = [line.strip() for line in value.splitlines() if line.strip()]
        existing = self.monomer_text.get("1.0", "end").strip()
        self.monomer_text.delete("1.0", "end")
        self.monomer_text.insert(
            "1.0", "\n".join(([existing] if existing else []) + new_values))

    def _choose_file(self, variable, extension):
        path = filedialog.askopenfilename(initialdir=Path(variable.get()).parent if variable.get() else ROOT,
                                          filetypes=[(f"{extension} files", f"*{extension}"), ("All files", "*")])
        if path:
            variable.set(path)

    def _choose_save(self, variable, extension):
        path = filedialog.asksaveasfilename(initialdir=Path(variable.get()).parent if variable.get() else ROOT,
                                             initialfile=Path(variable.get()).name if variable.get() else "",
                                             defaultextension=extension,
                                             filetypes=[(f"{extension} files", f"*{extension}"), ("All files", "*")])
        if path:
            variable.set(path)

    def _choose_directory(self, variable):
        path = filedialog.askdirectory(initialdir=variable.get() or ROOT)
        if path:
            variable.set(path)

    def _start(self):
        molecular_inputs = [line.strip() for line in self.monomer_text.get("1.0", "end").splitlines() if line.strip()]
        if not self._validate_molecular_inputs(molecular_inputs):
            return
        mol2_paths = [Path(value).expanduser().resolve() for value in molecular_inputs
                      if self._existing_mol2(value)]
        monomer_parents = {path.parent for path in mol2_paths}
        if len(monomer_parents) > 1:
            messagebox.showerror(
                "Monomers must share one folder",
                "Place the initial monomer MOL2 files in the same folder. Numbered reaction-output folders will be created there.")
            return
        self.session_root = next(iter(monomer_parents), AUTOMATION_ROOT)
        forcefield = Path(self.forcefield_var.get()).expanduser().resolve()
        if not forcefield.is_file() or forcefield.suffix.lower() != ".frc":
            messagebox.showerror("Missing force field", "Choose an existing .frc force-field file.")
            return
        if not RDKIT_PYTHON.is_file():
            messagebox.showerror(
                "RDKit is not installed",
                f"The bundled RDKit environment is missing:\n{RDKIT_PYTHON}")
            return
        self.current_forcefield = forcefield
        output_dir = self._next_reaction_directory()
        self._launch_reaction(molecular_inputs, forcefield, output_dir, clear_log=True)

    @staticmethod
    def _existing_mol2(value):
        try:
            path = Path(value).expanduser()
            return path.is_file() and path.suffix.lower() == ".mol2"
        except OSError:
            return False

    def _validate_molecular_inputs(self, inputs, minimum=1):
        if len(inputs) < minimum:
            messagebox.showerror(
                "Missing molecular inputs",
                "Select the previous reacted-product MOL2 and at least one other monomer."
                if minimum > 1 else "Add the source monomer files used by the reaction.")
            return False
        missing_mol2 = [value for value in inputs
                        if value.lower().endswith(".mol2") and not self._existing_mol2(value)]
        if missing_mol2:
            messagebox.showerror(
                "Invalid molecular input",
                "These entries look like MOL2 paths but do not exist:\n\n"
                + "\n".join(missing_mol2))
            return False
        return True

    def _next_reaction_directory(self):
        if self.session_root is None:
            raise RuntimeError("Reaction-series root has not been initialized")
        index = 1
        while True:
            candidate = self.session_root / f"Reaction_{index:03d}"
            if not candidate.exists():
                candidate.mkdir(parents=True)
                return candidate
            index += 1

    def _launch_reaction(self, molecular_inputs, forcefield, output_dir, clear_log=False):
        self.current_output_dir = output_dir
        raw = output_dir / "reacted_product.mol2"
        minimized = output_dir / "reacted_product_minimized.mol2"
        self.stop_requested = False
        self.approval_event.clear()
        self.approval_result = False
        self.chemistry_event.clear()
        self.chemistry_choice = None
        self.initiator_event.clear()
        self.initiator_result = None
        if clear_log:
            self.log.delete("1.0", "end")
        else:
            self.log.insert("end", f"\n\n{'#' * 72}\nStarting {output_dir.name}\n{'#' * 72}\n")
            self.log.see("end")
        self.run_button.configure(state="disabled")
        self.stop_button.configure(state="normal")
        self.status_var.set("Running")
        config = {
            "stages": [1, 2, 3, 4], "molecular_inputs": molecular_inputs,
            "output_dir": str(output_dir),
            "raw": str(raw), "minimized": str(minimized),
            "forcefield": str(forcefield),
            "charge_method": "auto",
            "require_approval": True,
            "typed": str(output_dir / "reacted_product_parameterized.data"),
            "radius": 3, "initiators": [],
            "pre": str(output_dir / "pre_reaction.data"),
            "post": str(output_dir / "post_reaction.data"),
            "lmpmol_dir": str(output_dir),
            "template_registry": str(self.session_root / "reaction_template_registry.json"),
            "selected_candidate": str(output_dir / "selected_reaction_candidate.json"),
        }
        self.worker = threading.Thread(target=self._run_workflow, args=(config,), daemon=True)
        self.worker.start()

    def _start_followup_reaction(self, completed):
        previous_product = filedialog.askopenfilename(
            title="Select the previous reacted-product MOL2",
            initialdir=completed["output_dir"],
            filetypes=[("MOL2 structures", "*.mol2"), ("All files", "*")],
        )
        if not previous_product:
            return False
        needs_additional = messagebox.askyesno(
            "Next reaction",
            "Does the next reaction require one or more additional monomer molecules?\n\n"
            "Choose No for an intramolecular reaction such as ring closure.",
            parent=self,
        )
        additional = ()
        if needs_additional:
            use_files = messagebox.askyesnocancel(
                "Additional monomer format",
                "Choose Yes to select MOL2 files.\n"
                "Choose No to enter one or more SMILES strings.", parent=self)
            if use_files is None:
                return False
            if use_files:
                additional = filedialog.askopenfilenames(
                    title="Select the additional monomer MOL2 file(s)",
                    initialdir=self.session_root or ROOT,
                    filetypes=[("MOL2 structures", "*.mol2"), ("All files", "*")],
                )
            else:
                smiles = simpledialog.askstring(
                    "Additional monomer SMILES",
                    "Enter one or more SMILES strings, one per line:", parent=self)
                additional = tuple(
                    line.strip() for line in (smiles or "").splitlines() if line.strip())
            if not additional:
                return False
        molecular_inputs = [previous_product, *additional]
        if not self._validate_molecular_inputs(molecular_inputs, minimum=1):
            return False
        self.monomer_text.delete("1.0", "end")
        self.monomer_text.insert("1.0", "\n".join(molecular_inputs))
        output_dir = self._next_reaction_directory()
        self._launch_reaction(
            molecular_inputs, self.current_forcefield, output_dir, clear_log=False)
        return True

    def _run_workflow(self, config):
        try:
            input_dir = Path(config["output_dir"]) / "input_monomers"
            input_manifest = Path(config["output_dir"]) / "molecular_inputs.json"
            prepare_command = [
                str(RDKIT_PYTHON), "-u",
                str(AUTOMATION_ROOT / "prepare_molecular_inputs.py")]
            for molecular_input in config["molecular_inputs"]:
                prepare_command.extend(["--input", molecular_input])
            prepare_command.extend([
                "--output-dir", str(input_dir), "--manifest", str(input_manifest)])
            self._run_process(0, "Prepare MOL2 and SMILES inputs", prepare_command)
            input_data = json.loads(input_manifest.read_text(encoding="utf-8"))
            config["monomers"] = input_data["resolved_mol2"]
            self.events.put(("line", (
                "Resolved molecular inputs:\n  "
                + "\n  ".join(config["monomers"]) + "\n\n")))
            for stage in config["stages"]:
                if self.stop_requested:
                    raise RuntimeError("Workflow stopped by user")
                if stage == 1:
                    raw_parent = Path(config["raw"]).expanduser().resolve().parent
                    chemistry_html = raw_parent / "chemistry_perception_report.html"
                    chemistry_json = raw_parent / "chemistry_perception_report.json"
                    chemistry_command = [
                        str(RDKIT_PYTHON), "-u",
                        str(AUTOMATION_ROOT / "chemistry_perception.py")]
                    for monomer in config["monomers"]:
                        chemistry_command.extend(["--monomer", monomer])
                    chemistry_command.extend([
                        "--library", str(CHEMISTRY_LIBRARY),
                        "--catalog", str(FUNCTIONAL_GROUP_CATALOG),
                        "--json", str(chemistry_json), "--html", str(chemistry_html),
                        "--topology-radius", str(config["radius"]),
                    ])
                    self._run_process(0, "RDKit chemistry perception", chemistry_command)
                    chemistry_data = json.loads(chemistry_json.read_text(encoding="utf-8"))
                    self.events.put(("chemistry_choice", {
                        "report": str(chemistry_html),
                        "candidates": chemistry_data.get("reaction_candidates", []),
                        "monomers": chemistry_data.get("monomers", []),
                    }))
                    self.chemistry_event.wait()
                    if self.chemistry_choice is None:
                        raise RuntimeError("Reaction selection was cancelled")
                    scripted_input = None
                    if self.chemistry_choice["mode"] == "suggestion":
                        candidate = self.chemistry_choice["candidate"]
                        Path(config["selected_candidate"]).write_text(
                            json.dumps(candidate, indent=2) + "\n", encoding="utf-8")
                        scripted_input = [*candidate["commands"], "done"]
                        config["initiators"] = candidate["initiator_source_refs"]
                        self.events.put(("line", (
                            "Using approved chemistry-library suggestion:\n  "
                            + "\n  ".join(candidate["commands"]) + "\n\n")))
                    before_graph = raw_parent / "reaction_graph.html"
                    after_graph = raw_parent / "final_product_graph.html"
                    validation_html = raw_parent / "reaction_validation_report.html"
                    validation_json = raw_parent / "reaction_validation_report.json"
                    command = [WORKFLOW_PYTHON, "-u", str(ROOT / "react_monomers.py"), *config["monomers"],
                               "-o", config["raw"], "--minimized-output", config["minimized"],
                               "--graph-output", str(before_graph),
                               "--final-graph-output", str(after_graph)]
                    self._run_process(
                        stage, "Apply suggested reaction and minimize" if scripted_input
                        else "Edit reaction and minimize",
                        command, scripted_input=scripted_input)
                    validation_command = [
                        WORKFLOW_PYTHON, "-u", str(AUTOMATION_ROOT / "validate_reaction.py")]
                    for monomer in config["monomers"]:
                        validation_command.extend(["--monomer", monomer])
                    validation_command.extend([
                        "--product", config["raw"],
                        "--before-graph", str(before_graph),
                        "--after-graph", str(after_graph),
                        "--html", str(validation_html), "--json", str(validation_json),
                    ])
                    self._run_process(stage, "Validate edited reaction", validation_command)
                    validation_data = json.loads(validation_json.read_text(encoding="utf-8"))
                    if config["require_approval"]:
                        self.events.put(("approval", str(validation_html)))
                        self.approval_event.wait()
                        if not self.approval_result:
                            raise RuntimeError("Reaction was not approved; downstream stages were not run")
                        self.events.put(("line", "\n✓ Reaction approved by user\n\n"))
                    selected_candidate = None
                    selected_path = Path(config["selected_candidate"])
                    if selected_path.is_file():
                        selected_candidate = json.loads(selected_path.read_text(encoding="utf-8"))
                    history_path = reaction_history.write_product_history(
                        config["monomers"], [config["raw"], config["minimized"]],
                        validation_data.get("changed_atoms", []),
                        Path(config["output_dir"]).name, selected_candidate)
                    self.events.put(("line", f"Atom reaction history: {history_path}\n"))
                elif stage == 2:
                    source = Path(config["minimized"]).expanduser().resolve()
                    forcefield = Path(config["forcefield"]).expanduser().resolve()
                    target = Path(config["typed"]).expanduser().resolve()
                    if not source.is_file():
                        raise FileNotFoundError(f"Missing minimized MOL2 for parameterization: {source}")
                    if not forcefield.is_file() or forcefield.suffix.lower() != ".frc":
                        raise FileNotFoundError(f"Missing .frc force-field file: {forcefield}")
                    self._parameterize(stage, source, target, forcefield, config["charge_method"])
                    monomer_dir = Path(config["output_dir"]) / "parameterized_monomers"
                    monomer_dir.mkdir(parents=True, exist_ok=True)
                    config["data_monomers"] = []
                    for index, monomer in enumerate(config["monomers"], 1):
                        monomer_source = Path(monomer).expanduser().resolve()
                        monomer_target = monomer_dir / f"monomer_{index:02d}.data"
                        self._parameterize(
                            stage, monomer_source, monomer_target, forcefield,
                            config["charge_method"],
                            label=f"Parameterize monomer {index}")
                        config["data_monomers"].append(str(monomer_target))
                elif stage == 3:
                    for label, path in [
                            *[(f"monomer {index}", path)
                              for index, path in enumerate(config["data_monomers"], 1)],
                            ("typed product", config["typed"])]:
                        if not Path(path).is_file():
                            raise FileNotFoundError(f"Missing {label} file: {path}")
                    output_dir = Path(config["pre"]).parent
                    atom_map = output_dir / "reaction_template_atom_map.tsv"
                    metadata = output_dir / "reaction_template_metadata.json"
                    command = [WORKFLOW_PYTHON, "-u", str(AUTOMATION_ROOT / "generate_reaction_data.py")]
                    for monomer in config["data_monomers"]:
                        command.extend(["--monomer", monomer])
                    command.extend([
                        "--product", config["typed"], "--pre-output", config["pre"],
                        "--post-output", config["post"], "--atom-map", str(atom_map),
                        "--report", str(output_dir / "reaction_template_report.md"),
                        "--metadata", str(metadata),
                        "--radius", str(config["radius"]),
                    ])
                    if config["initiators"]:
                        command.extend(["--initiator", *config["initiators"]])
                    self._run_process(stage, "Generate reduced LAMMPS data", command)
                    if 4 in config["stages"]:
                        config["template_initiators"] = self._resolve_template_initiators(
                            config, atom_map, metadata)
                elif stage == 4:
                    for label, path in (("pre data", config["pre"]), ("post data", config["post"])):
                        if not Path(path).is_file():
                            raise FileNotFoundError(f"Missing {label}: {path}")
                    atom_map = Path(config["pre"]).parent / "reaction_template_atom_map.tsv"
                    metadata = Path(config["pre"]).parent / "reaction_template_metadata.json"
                    if "template_initiators" not in config:
                        config["template_initiators"] = self._resolve_template_initiators(
                            config, atom_map, metadata)
                    command = [WORKFLOW_PYTHON, "-u", str(AUTOMATION_ROOT / "generate_lmpmol_pair.py"),
                               "--pre", config["pre"], "--post", config["post"],
                               "--atom-map", str(atom_map),
                               "--output-directory", config["lmpmol_dir"]]
                    if config["template_initiators"]:
                        command.extend(["--initiators", *map(str, config["template_initiators"])])
                    self._run_process(stage, "Generate LMPMOL templates and reaction map", command)
                    register_command = [
                        WORKFLOW_PYTHON, "-u", str(AUTOMATION_ROOT / "template_registry.py"),
                        "--registry", config["template_registry"],
                        "--reaction-folder", config["output_dir"],
                        "--pre", config["pre"], "--post", config["post"],
                        "--atom-map", str(atom_map),
                    ]
                    if Path(config["selected_candidate"]).is_file():
                        register_command.extend([
                            "--selected-candidate", config["selected_candidate"]])
                    self._run_process(
                        stage,
                        "Authoritative fix bond/react topology check and registration",
                        register_command)
                    identity = json.loads(
                        (Path(config["output_dir"]) / "template_identity.json").read_text(
                            encoding="utf-8"))
                    config["template_recommendation"] = identity["recommendation"]
            self.events.put(("done", {
                "message": "Reaction workflow completed successfully.",
                "output_dir": config["output_dir"],
                "product": config["minimized"],
                "template_recommendation": config.get("template_recommendation", {}),
            }))
        except Exception as error:
            self.events.put(("error", str(error)))
        finally:
            self.process = None

    def _parameterize(self, stage, source, target, forcefield, charge_method,
                      label="Parameterize reacted product"):
        target.parent.mkdir(parents=True, exist_ok=True)
        artifact_dir = Path(tempfile.mkdtemp(
            prefix=f"{target.stem}_artifacts_", dir=str(target.parent)))
        self.events.put(("line", f"Parameterization artifacts: {artifact_dir}\n"))
        command = [
            WORKFLOW_PYTHON, "-u",
            str(ROOT / "Force_Field_Parametrization" / "parameterize.py"),
            str(source), "--ff", str(forcefield), "-o", str(artifact_dir),
            "--name", source.stem, "--frc-charge", charge_method,
        ]
        code = self._run_process(stage, label, command, accepted_codes=(0, 3))
        audit_path = artifact_dir / "audit.json"
        if not audit_path.is_file():
            raise FileNotFoundError(
                f"Parameterization did not create its audit file: {audit_path}")
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        generated = Path(audit.get("lammps_data", "")).expanduser()
        if not generated.is_file():
            raise FileNotFoundError(
                f"Parameterization did not create the reported LAMMPS data file: {generated}")
        shutil.copy2(generated, target)
        self.events.put(("line", f"Copied verified LAMMPS data to: {target}\n"))
        if code == 3:
            coverage = audit.get("total_coverage_percent", "unknown")
            self.events.put(("line", (
                "WARNING: the parameterizer reported incomplete coefficient coverage "
                f"({coverage}%). Review {audit_path}.\n")))

    def _resolve_template_initiators(self, config, atom_map, metadata):
        """Map a user choice to template IDs; never choose among multiple candidates."""
        if not atom_map.is_file():
            raise FileNotFoundError(f"Missing reaction template atom map: {atom_map}")
        rows = atom_map.read_text(encoding="utf-8").splitlines()
        if not rows:
            raise ValueError(f"Empty reaction template atom map: {atom_map}")
        headers = rows[0].split("\t")
        lookup = {}
        for line in rows[1:]:
            if not line.strip():
                continue
            values = dict(zip(headers, line.split("\t")))
            lookup[f"{values['monomer']}:{values['source_atom_id']}"] = int(
                values["template_id"])

        selected = config.get("initiators", [])
        if selected:
            try:
                return tuple(lookup[value] for value in selected)
            except KeyError as error:
                raise ValueError(
                    f"Initiator source atom is outside the generated template: {error.args[0]}") from error

        if not metadata.is_file():
            return None
        candidates = json.loads(metadata.read_text(encoding="utf-8")).get(
            "initiator_candidates", [])
        if len(candidates) <= 1:
            return tuple(candidates[0]["template_ids"]) if candidates else None

        self.initiator_event.clear()
        self.initiator_result = None
        self.events.put(("initiator_choice", candidates))
        self.initiator_event.wait()
        if not self.initiator_result:
            raise RuntimeError("No initiator pair was selected; Stage 4 was not run")
        config["initiators"] = self.initiator_result["source_refs"]
        return tuple(self.initiator_result["template_ids"])

    def _choose_chemistry_candidate(self, candidates, monomers):
        dialog = tk.Toplevel(self)
        dialog.title("Choose detected reaction")
        dialog.transient(self)
        dialog.geometry("820x560")
        frame = ttk.Frame(dialog, padding=18)
        frame.pack(fill="both", expand=True)
        ttk.Label(
            frame,
            text=("Choose a chemistry-library suggestion to apply automatically, or use "
                  "manual editing. Every suggestion below was matched by RDKit and will "
                  "still require validation. The final fix bond/react template decision "
                  "is made after force-field typing."),
            wraplength=760,
        ).pack(anchor="w", pady=(0, 10))
        detected_lines = []
        for monomer in monomers:
            names = sorted({
                item["name"]
                for item in monomer.get("detected_functional_groups", [])
            })
            supported = sorted({
                item["name"]
                for item in monomer.get("detected_functional_groups", [])
                if item.get("reaction_support") == "supported"
            })
            detected_lines.append(
                f"Monomer {monomer['number']}: detected "
                f"{', '.join(names) if names else 'no catalog groups'}; "
                f"reaction-supported {', '.join(supported) if supported else 'none'}")
        detection_summary = ScrolledText(
            frame, height=min(5, max(2, len(detected_lines) + 1)),
            wrap="word", font=("Menlo", 9))
        detection_summary.pack(fill="x", pady=(0, 10))
        detection_summary.insert("1.0", "\n".join(detected_lines))
        detection_summary.configure(state="disabled")
        listbox = tk.Listbox(frame, height=12, exportselection=False,
                             font=("Menlo", 10))
        listbox.pack(fill="both", expand=True)
        for index, candidate in enumerate(candidates, 1):
            roles = ", ".join(
                f"{role}={atom}" for role, atom in candidate["role_atoms"].items())
            prior = candidate.get("prior_reaction_participation", [])
            if prior:
                prior_refs = ", ".join(item["atom_ref"] for item in prior)
                history_prefix = f"[PRIORITY—REACTED BEFORE: {prior_refs}] "
            else:
                history_prefix = ""
            listbox.insert(
                "end", f"{index}. {history_prefix}{candidate['rule_name']} | {roles}")
        details = ScrolledText(frame, height=8, wrap="word", font=("Menlo", 10))
        details.pack(fill="x", pady=(10, 0))
        details.configure(state="disabled")

        def show_details(_event=None):
            selection = listbox.curselection()
            details.configure(state="normal")
            details.delete("1.0", "end")
            if selection:
                candidate = candidates[selection[0]]
                details.insert("end", candidate["description"] + "\n\n")
                details.insert("end", "Topology commands:\n  "
                               + "\n  ".join(candidate["commands"]))
                details.insert("end", "\n\nInitiators: "
                               + " — ".join(candidate["initiator_source_refs"]))
                prior = candidate.get("prior_reaction_participation", [])
                priority = candidate.get("history_priority", {})
                details.insert(
                    "end", "\n\nHistory priority: "
                    + ("yes" if priority.get("is_prioritized") else "no")
                    + f" (reacted initiators: "
                    f"{priority.get('previously_reacted_initiator_count', 0)}; "
                    f"reacted role atoms: "
                    f"{priority.get('previously_reacted_role_atom_count', 0)})")
                details.insert("end", "\n\nEarlier reaction participation:")
                if not prior:
                    details.insert("end", " none of the proposed atoms")
                for participant in prior:
                    details.insert("end", f"\n  {participant['atom_ref']}")
                    for event in participant["prior_events"]:
                        role_text = ", ".join(event.get("roles", [])) or "changed atom"
                        details.insert(
                            "end", f" — {event.get('reaction', 'earlier reaction')} "
                            f"as {role_text}")
            elif not candidates:
                details.insert(
                    "end", "No complete reaction rule matched. Choose Manual editing, or "
                    "add an appropriate rule to chemistry_library.json.")
            details.configure(state="disabled")

        listbox.bind("<<ListboxSelect>>", show_details)
        if candidates:
            listbox.selection_set(0)
        show_details()
        buttons = ttk.Frame(frame)
        buttons.pack(fill="x", pady=(12, 0))

        def use_selected():
            selection = listbox.curselection()
            if not selection:
                messagebox.showwarning("No suggestion", "Select one reaction suggestion.", parent=dialog)
                return
            self.chemistry_choice = {
                "mode": "suggestion", "candidate": candidates[selection[0]]}
            dialog.destroy()

        def manual():
            self.chemistry_choice = {"mode": "manual"}
            dialog.destroy()

        def cancel():
            self.chemistry_choice = None
            dialog.destroy()

        ttk.Button(buttons, text="Cancel", command=cancel).pack(side="right")
        ttk.Button(buttons, text="Manual editing", command=manual).pack(side="right", padx=8)
        use_button = ttk.Button(buttons, text="Use selected suggestion", command=use_selected)
        use_button.pack(side="right")
        if not candidates:
            use_button.configure(state="disabled")
        dialog.protocol("WM_DELETE_WINDOW", cancel)
        dialog.grab_set()
        dialog.wait_window()
        self.chemistry_event.set()

    def _choose_initiator(self, candidates):
        dialog = tk.Toplevel(self)
        dialog.title("Choose REACTER initiator pair")
        dialog.transient(self)
        dialog.resizable(False, False)
        frame = ttk.Frame(dialog, padding=18)
        frame.pack(fill="both", expand=True)
        ttk.Label(
            frame,
            text=("More than one new bond was detected. Choose which newly formed "
                  "bond is the REACTER initiator pair."),
            wraplength=560,
        ).pack(anchor="w", pady=(0, 12))
        choice = tk.IntVar(value=0)
        for index, candidate in enumerate(candidates):
            source = " — ".join(candidate["source_refs"])
            template = "–".join(map(str, candidate["template_ids"]))
            ttk.Radiobutton(
                frame,
                text=f"{source}   (template atoms {template})",
                variable=choice,
                value=index,
            ).pack(anchor="w", pady=4)
        buttons = ttk.Frame(frame)
        buttons.pack(fill="x", pady=(16, 0))

        def accept():
            self.initiator_result = candidates[choice.get()]
            self.initiator_refs_var.set(" ".join(self.initiator_result["source_refs"]))
            dialog.destroy()

        def cancel():
            self.initiator_result = None
            dialog.destroy()

        ttk.Button(buttons, text="Cancel", command=cancel).pack(side="right")
        ttk.Button(buttons, text="Use selected pair", command=accept).pack(
            side="right", padx=8)
        dialog.protocol("WM_DELETE_WINDOW", cancel)
        dialog.grab_set()
        dialog.wait_window()
        self.initiator_event.set()

    def _run_process(self, stage, name, command, accepted_codes=(0,), scripted_input=None):
        self.events.put(("stage", (stage, name)))
        self.events.put(("line", "$ " + " ".join(self._quote(value) for value in command) + "\n"))
        self.process = subprocess.Popen(command, cwd=ROOT, stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        text=True, bufsize=1)
        if scripted_input:
            assert self.process.stdin is not None
            self.process.stdin.write("\n".join(scripted_input) + "\n")
            self.process.stdin.flush()
            self.process.stdin.close()
        assert self.process.stdout is not None
        for line in self.process.stdout:
            self.events.put(("line", line))
        code = self.process.wait()
        self.process = None
        if code not in accepted_codes:
            raise RuntimeError(f"Stage {stage} ({name}) failed with exit code {code}")
        self.events.put(("line", f"\n✓ Stage {stage} complete: {name}\n\n"))
        return code

    @staticmethod
    def _quote(value):
        return f'"{value}"' if any(char.isspace() for char in str(value)) else str(value)

    def _send_command(self):
        value = self.command_var.get().strip()
        process = self.process
        if not value:
            return
        if (process is None or process.stdin is None or process.stdin.closed
                or process.poll() is not None):
            messagebox.showinfo("No interactive stage", "There is no running process waiting for input.")
            return
        try:
            process.stdin.write(value + "\n")
            process.stdin.flush()
            self.events.put(("line", f"> {value}\n"))
            self.command_var.set("")
        except (OSError, ValueError) as error:
            messagebox.showerror("Input failed", str(error))

    def _stop(self):
        self.stop_requested = True
        self.approval_event.set()
        self.chemistry_event.set()
        self.initiator_event.set()
        if self.process is not None and self.process.poll() is None:
            self.process.terminate()
        self.status_var.set("Stopping")

    def _open_folder(self):
        directory = self.current_output_dir or self.session_root or AUTOMATION_ROOT
        directory.mkdir(parents=True, exist_ok=True)
        if sys.platform == "darwin":
            subprocess.Popen(["open", str(directory)])
        elif os.name == "nt":
            os.startfile(directory)  # type: ignore[attr-defined]
        else:
            subprocess.Popen(["xdg-open", str(directory)])

    def _drain_events(self):
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == "line":
                    self.log.insert("end", value)
                    self.log.see("end")
                elif kind == "stage":
                    stage, name = value
                    self.status_var.set(f"Stage {stage}: {name}")
                    self.log.insert("end", f"\n{'=' * 72}\nSTAGE {stage}: {name}\n{'=' * 72}\n")
                    self.log.see("end")
                elif kind == "done":
                    self.status_var.set("Complete")
                    self.stop_button.configure(state="disabled")
                    recommendation = value.get("template_recommendation", {})
                    decision = recommendation.get(
                        "status", "not evaluated").replace("_", " ").upper()
                    matches = recommendation.get("matching_reactions", [])
                    decision_text = (
                        f"\n\nAuthoritative fix bond/react topology decision: {decision}\n"
                        + recommendation.get("reason", "")
                    )
                    if matches:
                        decision_text += "\nMatching template: " + ", ".join(matches)
                    another = messagebox.askyesno(
                        "Additional reaction?",
                        f"{value['message']}\n\nOutput folder:\n{value['output_dir']}\n\n"
                        f"{decision_text}\n\n"
                        "Do you want to generate another chemical reaction?\n\n"
                        "There is no fixed reaction count. Choose Yes to load the previous "
                        "product and continue while reactive groups remain.",
                    )
                    if another and self._start_followup_reaction(value):
                        continue
                    self.run_button.configure(state="normal")
                elif kind == "approval":
                    report = Path(value)
                    if sys.platform == "darwin":
                        subprocess.Popen(["open", str(report)])
                    elif os.name == "nt":
                        os.startfile(report)  # type: ignore[attr-defined]
                    else:
                        subprocess.Popen(["xdg-open", str(report)])
                    self.approval_result = messagebox.askyesno(
                        "Approve reaction?",
                        "The reaction validation report has opened.\n\n"
                        "Review the removed/added bonds, hydrogens, charges, and final structure.\n\n"
                        "Approve this reaction and continue to the selected downstream stages?",
                    )
                    self.approval_event.set()
                elif kind == "chemistry_choice":
                    report = Path(value["report"])
                    if sys.platform == "darwin":
                        subprocess.Popen(["open", str(report)])
                    elif os.name == "nt":
                        os.startfile(report)  # type: ignore[attr-defined]
                    else:
                        subprocess.Popen(["xdg-open", str(report)])
                    self._choose_chemistry_candidate(
                        value["candidates"], value.get("monomers", []))
                elif kind == "initiator_choice":
                    self._choose_initiator(value)
                elif kind == "error":
                    self.status_var.set("Stopped" if self.stop_requested else "Failed")
                    self.run_button.configure(state="normal")
                    self.stop_button.configure(state="disabled")
                    self.log.insert("end", f"\nERROR: {value}\n")
                    self.log.see("end")
                    if not self.stop_requested:
                        messagebox.showerror("Workflow failed", value)
        except queue.Empty:
            pass
        self.after(100, self._drain_events)

    def _close(self):
        if self.process is not None and self.process.poll() is None:
            if not messagebox.askyesno("Quit", "A workflow stage is running. Stop it and quit?"):
                return
            self.process.terminate()
        self.destroy()


if __name__ == "__main__":
    ReactionWorkflowGUI().mainloop()
