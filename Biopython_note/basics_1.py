### creating and manipulating DNA sequences, performing BLAST searches, analyzing sequence features like reverse complements, translating sequences to proteins, and working with Seq objects.

import Bio

######### 1 #########
# How many records are in the multi-FASTA file?

from Bio import SeqIO

# Open and count records in dna2.fasta
fasta_file = "dna2.fasta"
records = list(SeqIO.parse(fasta_file, "fasta"))

print("Number of records in file:", len(records))

######### 2 #########
# What is the length of the longest sequence in the file?

seq_lengths = [len(record.seq) for record in records]
longest = max(seq_lengths)
shortest = min(seq_lengths)
longest_ids = [record.id for record in records if len(record.seq) == longest]
shortest_ids = [record.id for record in records if len(record.seq) == shortest]

print("Longest sequence length:", longest, "IDs:", longest_ids)
print("Shortest sequence length:", shortest, "IDs:", shortest_ids)

######### 3 #########
# What is the length of the longest ORF appearing in reading frame 2 of any of the sequences?
start_codon = "ATG"
stop_codons = ["TAA", "TAG", "TGA"]
reading_frame = 3  # example: frame 2

longest_orf_length = 0
longest_orf_seq = ""
longest_orf_id = ""
longest_orf_start = 0

for record in records:
    seq = str(record.seq)
    seq = seq[reading_frame-1:]  # adjust for frame
    i = 0
    while i < len(seq)-2:
        codon = seq[i:i+3]
        if codon == start_codon:
            j = i
            while j < len(seq)-2:
                stop = seq[j:j+3]
                if stop in stop_codons:
                    orf_len = j + 3 - i
                    if orf_len > longest_orf_length:
                        longest_orf_length = orf_len
                        longest_orf_seq = seq[i:j+3]
                        longest_orf_id = record.id
                        longest_orf_start = i + 1  # 1-based index
                    break
                j += 3
            i = j
        else:
            i += 3

print("Longest ORF length in frame", reading_frame, ":", longest_orf_length)
print("Sequence ID:", longest_orf_id, "Start position:", longest_orf_start)

######### 4 #########
# What is the length of the longest ORF appearing in any sequence and in any forward reading frame?
longest_orf_length_all = 0
longest_orf_id_all = ""
longest_orf_frame_all = 0
longest_orf_start_all = 0

for record in records:
    seq = str(record.seq)
    for frame in range(3):
        sub_seq = seq[frame:]
        i = 0
        while i < len(sub_seq)-2:
            if sub_seq[i:i+3] == start_codon:
                j = i
                while j < len(sub_seq)-2:
                    if sub_seq[j:j+3] in stop_codons:
                        orf_len = j+3-i
                        if orf_len > longest_orf_length_all:
                            longest_orf_length_all = orf_len
                            longest_orf_id_all = record.id
                            longest_orf_frame_all = frame+1
                            longest_orf_start_all = i+1
                        break
                    j += 3
                i = j
            else:
                i += 3

print("Longest ORF in any frame:", longest_orf_length_all)
print("Sequence ID:", longest_orf_id_all, "Frame:", longest_orf_frame_all, "Start:", longest_orf_start_all)

######### 5 #########
# Find the most frequently occurring repeat of length 6 in all sequences. How many times does it occur in all?
repeat_length = 6  
repeat_counts = {}

for record in records:
    seq = str(record.seq)
    for i in range(len(seq) - repeat_length + 1):
        subseq = seq[i:i+repeat_length]
        repeat_counts[subseq] = repeat_counts.get(subseq, 0) + 1

# Most frequent repeat
max_count = max(repeat_counts.values())
most_freq_repeats = [seq for seq, count in repeat_counts.items() if count == max_count]

print("Most frequent repeat(s) of length", repeat_length, ":", most_freq_repeats)
print("Number of occurrences:", max_count)

######### 6 #########
# How many different 12-base sequences occur Max times?
repeat_length = 12
repeat_counts = {}

for record in records:
    seq = str(record.seq)
    for i in range(len(seq) - repeat_length + 1):
        subseq = seq[i:i+repeat_length]
        repeat_counts[subseq] = repeat_counts.get(subseq, 0) + 1

max_count = max(repeat_counts.values())
num_seqs_with_max = sum(1 for count in repeat_counts.values() if count == max_count)

print("Number of 12-base sequences with Max occurrences:", num_seqs_with_max)
