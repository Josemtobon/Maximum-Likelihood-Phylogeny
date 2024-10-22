"""
Filter sequences of a genbank file based on gene type and specie
"""
from Bio import SeqIO


with open('sequences.gb', 'r') as gb:
    records = list(SeqIO.parse(gb, 'gb'))

with open('target_species.txt', 'r') as txt:
    target_species = [line.strip().split()[-1] for line in txt]

target_genera = ['Adelophryne', 'Eleutherodactylus', 'Diasporus', 'Phyzelaphryne', 'Ischnocnema', 'Brachycephalus']

seen_species = set()
records_to_write = []

for record in records:
    specie = None
    annotations = record.annotations

    if 'organism' in annotations and len(record.seq) > 600 and len(record.seq) < 800:
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

species_count = {}
for record in records_to_write:
    specie = record.annotations['organism']
    species_count[specie] =+ 1
    record.id = specie.replace(' ', '_')
    record.description = ''

for specie, count in species_count.items():
    print(f'{specie}: {count}')

SeqIO.write(records_to_write, 'sequences_to_align.fasta', 'fasta')
