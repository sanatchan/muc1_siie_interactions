import glycosylator as gl
import numpy as np

INPUT_PDB = "/Users/sanaaatttt/Desktop/final_structural_proj/10rpts_muc_native.pdb"
OUTPUT_PDB = "/Users/sanaaatttt/Desktop/final_structural_proj/sial_10rpts_muc.pdb "

CORE1 = gl.Glycan.from_pdb("sialylated_core1.pdb")
CORE2 = gl.Glycan.from_pdb("sialylated_core2.pdb")

THR_LINK = gl.get_linkage("THR-glyco")
SER_LINK = gl.get_linkage("SER-glyco")


def find_atom(res, atom_id):
    return next(a for a in res.get_atoms() if a.id == atom_id)


def add_hydroxyl_hydrogen(protein, res, o_atom_id, h_atom_id):
    """Add the missing hydroxyl hydrogen (HG/HG1) off CB->O if not already present."""
    if h_atom_id in [a.id for a in res.get_atoms()]:
        return
    o_atom = find_atom(res, o_atom_id)
    cb = find_atom(res, "CB")
    direction = o_atom.coord - cb.coord
    direction = direction / np.linalg.norm(direction)
    h_coord = o_atom.coord + direction * 0.96
    protein.add_atoms(gl.Atom(h_atom_id, coord=h_coord, element="H"), residue=res)
    protein.add_bond(o_atom, find_atom(res, h_atom_id))


def glycosylate_residues(protein, residues, link, label):
    """Glycosylate residues in place, alternating CORE1/CORE2. Returns count attached."""
    attached = 0
    for i, res in enumerate(residues):
        glycan = CORE1 if i % 2 == 0 else CORE2
        glycan_label = "Core1" if i % 2 == 0 else "Core2"
        try:
            protein.glycosylate(glycan, link=link, residues=[res], inplace=True)
            attached += 1
            print(f"  \u2713 {label} {res.id[1]} \u2014 {glycan_label}")
        except Exception as e:
            print(f"  \u2717 {label} {res.id[1]} \u2014 {type(e).__name__}: {e!r}")
    return attached


# --- Load and prepare structure once ---
protein = gl.Protein.from_pdb(INPUT_PDB)
protein.reindex()
protein.infer_bonds()
protein.apply_standard_bonds()

all_res = list(protein.get_residues())

'''thr res'''

thr_res = [
    r for r in all_res
    if r.resname == "THR" and 1 <= r.id[1] <= 201 and r.id[1] % 20 == 1
]
print(f"Filtered THR residues: {len(thr_res)}")

#53-631 correspond to the range of extruding thr repeats spaced at 20 AA.

for res in thr_res:
    add_hydroxyl_hydrogen(protein, res, "OG1", "HG1")

thr_attached = glycosylate_residues(protein, thr_res, THR_LINK, "THR")
print(f"THR attached: {thr_attached}/{len(thr_res)}\n")

"""ser res"""

ser_res = [
    r for r in all_res
    if r.resname == "SER" and 10 <= r.id[1] <= 210 and r.id[1] % 20 == 10
]
print(f"Filtered SER residues: {len(ser_res)}")

for res in ser_res:
    add_hydroxyl_hydrogen(protein, res, "OG", "HG")

ser_attached = glycosylate_residues(protein, ser_res, SER_LINK, "SER")
print(f"SER attached: {ser_attached}/{len(ser_res)}\n")

print(f"Total glycans: {protein.count_glycans()}")

protein.to_pdb(OUTPUT_PDB)
print(f"Saved to {OUTPUT_PDB}")
