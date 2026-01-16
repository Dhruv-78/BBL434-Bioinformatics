import sys

if len(sys.argv) != 2:
    print("Usage: python count_fasta_records.py <mfa_file>")
    sys.exit(1)

mfa_file = sys.argv[1]
record_count = 0

with open(mfa_file, "r") as f:
    for line in f:
        if line.startswith(">"):
            record_count += 1

print(record_count)
