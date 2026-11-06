# RDKit-informed automated reaction workflow

Launch by double-clicking `launch_reaction_workflow_gui.command`.

The GUI has exactly two user inputs:

1. one or more monomer `.mol2` files;
2. one force-field `.frc` file.

The initial monomer files must be in the same directory. The GUI creates a new
`Reaction_001`, `Reaction_002`, ... directory there for every reaction, so no
reaction template overwrites an earlier one. Charge selection is automatic and
the reaction-template graph radius is 3 bonds.

## Automated sequence

1. The bundled RDKit environment reads and sanitizes every MOL2 file.
2. `functional_groups_catalog.json` loads RDKit's broad functional-group
   hierarchy and any user-added SMARTS detectors. Detection is reported even
   when no reaction mechanism is available. `chemistry_library.json` separately
   identifies reaction-supported groups and constructs complete, atom-specific
   reaction candidates. RDKit proposes the
   reacting atoms and bond changes; it does not display or make a template-reuse
   decision. For Reaction_002 and later, each suggestion identifies any proposed
   atom ID that participated in an earlier reaction and shows its earlier role.
3. `chemistry_perception_report.html` opens and the GUI asks whether to apply a
   detected suggestion or use manual editing.
4. For an approved suggestion, the GUI sends its library-defined bond removals,
   additions, hydrogen transfer, and initiators to the topology editor. In
   manual mode, the user supplies those commands.
5. Open Babel minimizes the complete reacted MOL2 product.
6. The structural validation report opens for approval.
7. The selected `.frc` file parameterizes every original monomer and the
   reacted product. Each parameterization keeps its own audit and artifacts.
8. The workflow derives reduced pre/post LAMMPS data templates from graph
   differences.
9. If more than one new bond exists, the GUI displays all graph-derived
   initiator candidates and waits for the user to select one.
10. The workflow creates synchronized `.lmpmol` files and `reaction.map`.
    It then performs the authoritative fix bond/react topology check. The check
    fingerprints the complete force-field-typed pre/post atom types, charges,
    bonds, angles, dihedrals, impropers, and reaction-role mapping. Only this
    signature decides whether a previous template can be reused or a new
    template is required.
11. The GUI asks whether another reaction is required. If yes, select the
    previous reacted-product MOL2 and one or more additional monomer MOL2 files.
    The same force field is reused, RDKit perception starts again, and all new
    results are written to the next numbered reaction directory.

RDKit never chooses a mechanism and never decides fix bond/react template
reuse. It provides element, bond-order, aromaticity,
ring, hybridization, hydrogen, charge, valence, and SMARTS-substructure
perception. The external library supplies named functional groups and compatible
reaction families and their topology actions. The user still selects the site
and approves the complete topology change.

## Files created in each `Reaction_###` directory

- `chemistry_perception_report.html` and `.json`
- `reaction_graph.html`
- `reacted_product.mol2`
- `reacted_product_minimized.mol2`
- `final_product_graph.html`
- `reaction_validation_report.html` and `.json`
- `reacted_product_parameterized.data`
- `parameterized_monomers/monomer_XX.data` and audit folders
- `pre_reaction.data` and `post_reaction.data`
- `reaction_template_atom_map.tsv`
- `reaction_template_metadata.json`
- `reaction_template_report.md`
- `pre_reaction.lmpmol` and `post_reaction.lmpmol`
- `reaction.map`
- `lmpmol_type_mapping.tsv`
- `selected_reaction_candidate.json` when a library suggestion was used
- `reaction_atom_history.json`, carrying each reacted product atom's reaction
  participation history into Reaction_002, Reaction_003, and later steps
- `template_identity.json`, containing the preliminary and typed topology keys
- `topology_search_report.md`, listing the reuse/new decision plus every atom
  type, charge, bond, angle, dihedral, and improper change

The shared monomer directory also contains `reaction_template_registry.json`.
It indexes every completed reaction folder and is used by later topology
searches. A template is reused only when the authoritative force-field-typed
pre/post topology signature is identical;
changed atom typing, charge, connectivity, bonded terms, or reaction-role
mapping produces a new fix bond/react topology.

## Follow-up reaction priority

For Reaction_002 and every later reaction, candidates containing atoms recorded
in `reaction_atom_history.json` are placed first. The ranking is generic:
previously reacted initiator atoms rank first, followed by other previously
reacted atoms used by the proposed topology change, followed by candidates with
no prior participation. Ties retain the chemistry library's normal order. The
GUI marks these entries with `PRIORITY—REACTED BEFORE` and preselects the first
ranked candidate; the user still approves the exact reaction. The ranking uses
only carried atom history and reaction roles—never material names or fixed atom
IDs.

## Chemistry library

`functional_groups_catalog.json` controls broad recognition. Its default
`rdkit_functional_group_hierarchy` source loads the hierarchy included with the
installed RDKit version. Additional detectors can be appended to its
`functional_groups` list with `id`, `name`, `category`, and `smarts`; Python and
the GUI do not need to be edited. These entries are detection-only and cannot
invent topology changes.

The GUI and HTML report explicitly separate **detected functional groups** from
**reaction-supported functional groups**. A molecule made from a new material
therefore needs no code update: if it contains already-supported chemistry, the
existing rules generate candidates. If RDKit recognizes a group for which no
validated rule exists, the group is still reported but no reaction is proposed.

`chemistry_library.json` is deliberately external to the Python source. Add or
revise SMARTS patterns and reaction-family compatibility there instead of adding
material names, atom IDs, or molecule-specific conditions to code.

The starter library recognizes epoxide, primary and secondary amine, hydroxyl,
thiol, isocyanate, cyanate ester, and acrylate-family alkene groups. Its starter
reaction rules define primary/secondary-amine–epoxide ring opening, including
the selected C–O cleavage, N–C formation, explicit N–H transfer to oxygen, and
the N/C initiator pair. The library also defines chemistry-informed cyanate
cyclotrimerization as three sequential states: first intermolecular link,
third-group addition, and intramolecular triazine-ring closure. Each state
appears in the same RDKit suggestion window as an epoxy reaction and is applied
only after the user selects one exact candidate.

These three entries describe topology states; they do not impose a three-reaction
limit. After every completed reaction the GUI asks whether another chemical
reaction is required. The sequence can repeat for additional cyanate groups,
additional monomers, other supported chemistry, or later network-growth events.

The cyanate rules are connectivity-based and contain no material names or atom
IDs. They generate bond-order and MOL2 atom-type changes together with each new
C–N bond. After an intermolecular stage, the follow-up dialog can load another
monomer; for the ring-closure stage choose `No` when asked whether an additional
monomer is required.

The epoxide entry uses the external `role_selection` policy
`minimum_heavy_atom_degree`. RDKit scores both ring carbons from their molecular
connectivity and retains only the least substituted carbon as the electrophile.
Symmetry-equivalent carbons remain tied and available. The policy contains no
material names or atom IDs and can be changed in `chemistry_library.json` for a
different reaction mechanism.

Manual editing supports changing an existing bond order without deleting and
recreating the bond. For example, `order 1:13 1:14 2` converts the selected
triple bond to a double bond. Equivalent commands are `order-bond
MONOMER:BOND_ID NEW_ORDER` and `order-number COMBINED_BOND_NUMBER NEW_ORDER`.

## Strict input behavior

Chemistry perception requires MOL2 input with meaningful elements and bond
orders. The GUI stops on RDKit parsing or sanitization failures. It does not
guess chemistry from a LAMMPS numeric atom type. Unsupported or ambiguous
REACTER choices also stop or request user selection rather than being inferred
from a material name.
