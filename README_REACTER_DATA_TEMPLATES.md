# Reduced LAMMPS data reaction structures

These are genuine LAMMPS `read_data` files, not lmpmol/native molecule files. Each includes box bounds, type counts, Masses, PCFF/Class-II coefficient sections, Atoms, Bonds, Angles, Dihedrals, and Impropers.

- Pre-reaction data: `detda_epon_pre_reaction.data`
- Post-reaction data: `detda_epon_post_reaction.data`
- Atom correspondence: `detda_epon_data_atom_map.tsv`

Both structures contain 20 atoms in identical order. The pre file combines the original DETDA and EPON coefficient namespaces by offsetting EPON types. The post file retains the types and coefficients from `reacted_product_minimized_pcff.data`.

Bond changes using reduced atom IDs: remove 8-9 and 14-15; add 8-15 and 9-14.

These reduced data files are suitable for `read_data`, inspection, or input to a data-to-REACTER-template conversion tool. LAMMPS `fix bond/react` itself ultimately loads native molecule templates through the `molecule` command, not data files.

Charge check: reduced pre=-0.45585000 e; reduced post=+0.45440000 e; complete reacted source=+0.75000000 e. Check the reacted charge assignment before production simulation.
