import re

def replace_nxv(protein_sequence):
    pattern = r'N(?=[^P]V)'
    changes = []
    modified_sequence = protein_sequence

    for match in re.finditer(pattern, protein_sequence):
        pos = match.start()
        original = protein_sequence[pos:pos+3]
        modified = 'J' + original[1:]
        changes.append((pos, original, modified))
        modified_sequence = modified_sequence[:pos] + 'J' + modified_sequence[pos+1:]
    
    return modified_sequence, changes

def process_fasta(filename, output_filename, changes_log_filename):
    total_changes = 0
    with open(filename, 'r') as file, \
         open(output_filename, 'w') as output_file, \
         open(changes_log_filename, 'w') as log_file:
        header = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    modified_sequence, changes = replace_nxv(''.join(sequence))
                    output_file.write(header + '\n')
                    output_file.write(modified_sequence + '\n')
                    if changes:
                        log_file.write(header + '\n')
                        log_file.write('Changes:\n')
                        for pos, original, modified in changes:
                            log_file.write(f"Position {pos}: {original} -> {modified}\n")
                            total_changes += 1
                header = line
                sequence = []
            else:
                sequence.append(line)
        if header:
            modified_sequence, changes = replace_nxv(''.join(sequence))
            output_file.write(header + '\n')
            output_file.write(modified_sequence + '\n')
            if changes:
                log_file.write(header + '\n')
                log_file.write('Changes:\n')
                for pos, original, modified in changes:
                    log_file.write(f"Position {pos}: {original} -> {modified}\n")
                    total_changes += 1

        log_file.write(f"\nTotal number of changes: {total_changes}\n")

fasta_file = 'proteome.fasta'
output_file = 'modified_proteome.fasta'
changes_log_file = 'changes_log.txt'

process_fasta(fasta_file, output_file, changes_log_file)
