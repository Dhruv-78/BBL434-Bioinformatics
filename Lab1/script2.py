import sys

if len(sys.argv) != 3:
    print("Usage: python kmer_count.py <fasta_file> <k>")
    sys.exit(1)

fasta_file = sys.argv[1]
k = int(sys.argv[2])

kmer_counts = {}
sequence = []

with open(fasta_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            # Process previous sequence
            seq = "".join(sequence)
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
            sequence = []
        else:
            sequence.append(line.upper())

    # Process last sequence
    seq = "".join(sequence)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1

# Print dictionary (optional)
for kmer, count in kmer_counts.items():
    print(f"{kmer}\t{count}")
