import sys

if len(sys.argv) != 2:
    print("Usage: python fasta_length.py <fasta_file>")
    sys.exit(1)

fasta_file = sys.argv[1]

with open(fasta_file, 'r') as f:
    seq_id = None
    seq_len = 0

    for line in f:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            # Print previous sequence length
            if seq_id is not None:
                print(f"{seq_id}\t{seq_len}")

            seq_id = line[1:].split()[0]  # take first word as ID
            seq_len = 0
        else:
            seq_len += len(line)

    # Print last sequence
    if seq_id is not None:
        print(f"{seq_id}\t{seq_len}")
