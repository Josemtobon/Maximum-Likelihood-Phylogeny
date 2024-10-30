"""
Filter sequences of a genbank file based on gene type and specie
"""
from Bio import SeqIO
import pandas as pd


with open('sequences.gb', 'r') as gb:
    records = list(SeqIO.parse(gb, 'gb'))

with open('target_species.txt', 'r') as txt:
    target_species = [line.strip().split()[-1] for line in txt]

target_genera = ['Adelophryne', 'Eleutherodactylus', 'Diasporus', 'Phyzelaphryne', 'Ischnocnema']

seen_species = set()
records_to_write = []

def val_seq(seq):
    """
    Verifies if all characters in a sequence is valid
    """
    nucleotides = ['A', 'C', 'G', 'T']
    return all(nt in nucleotides for nt in seq)

for record in records:
    specie = None
    annotations = record.annotations

    if 'organism' in annotations and \
        len(record.seq) > 600 and \
        len(record.seq) < 700 and \
        val_seq(record.seq):
        organism = annotations['organism']
        organism_parts = organism.split()
        if organism_parts[1] in target_species and organism_parts[0] in target_genera:
            specie = annotations['organism']

    if specie and specie not in seen_species:
        for feature in record.features:
            if 'gene' in feature.qualifiers:
                gene = feature.qualifiers['gene'][0].lower()
                if gene == 'cox1' or gene == 'coi':
                    records_to_write.append(record)
                    seen_species.add(specie)
                    break


accession_df = pd.DataFrame(columns=["specie", "accession", "group"])
for record in records_to_write:
    specie = record.annotations["organism"].replace(" ", "_")
    accession = record.id
    if record.annotations["organism"].split()[0] == target_genera[-1]:
        group = "outgroup"
    else:
        group = "ingroup"

    accession_df.loc[len(accession_df)] = [specie, accession, group]


species_count = {}
for record in records_to_write:
    specie = record.annotations['organism']
    species_count[specie] =+ 1
    record.id = specie.replace(' ', '_')
    record.description = ''

for specie, count in species_count.items():
    print(f'{specie}: {count}')

SeqIO.write(records_to_write, 'sequences_to_align.fasta', 'fasta')
accession_df.to_csv("grupo_08_tabla_acceso.tsv", sep="\t")
