from Bio.PDB import MMCIFParser, PDBIO
import os

input_dir = "../inputs/my_structures/5UOI/my_cifs"
output_dir = "../inputs/my_structures/5UOI/my_pdbs"

parser = MMCIFParser(QUIET=True)
io = PDBIO()

for filename in os.listdir(input_dir):
    if filename.endswith(".cif"):
        structure_id = filename.rsplit(".", 1)[0]
        cif_path = os.path.join(input_dir, filename)
        pdb_path = os.path.join(output_dir, structure_id + ".pdb")
        print(f"Converting {filename} -> {structure_id}.pdb")
        structure = parser.get_structure(structure_id, cif_path)
        io.set_structure(structure)
        io.save(pdb_path)

