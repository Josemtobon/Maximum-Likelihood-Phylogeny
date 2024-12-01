from Bio import SeqIO

input_file = "grupo_08_alineamiento_recortado.fasta"
output_file = "grupo_08_alineamiento_recortado.nexus"

SeqIO.convert(input_file, "fasta", output_file, "nexus", molecule_type="DNA")
