## DNA String, Sequences

Data Source: https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa

Functions will be used
```
# Define reverse complement function
def reverseComplement(s):
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    return ''.join(comp[b] for b in s[::-1])

# Define naive exact match (already provided)
def naive(p, t):
    occurrences = []
    for i in range(len(t)-len(p)+1):
        match = True
        for j in range(len(p)):
            if not t[i+j] == p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences

# Define function that checks both pattern and its reverse complement
def naive_with_rc(p, t):
    rc_p = reverseComplement(p)
    occurrences = naive(p, t)
    if rc_p != p:
        occurrences = sorted(list(set(occurrences + naive(rc_p, t))))
    return occurrences

# Read genome from FASTA file
def readGenome(filename):
    genome = ''
    with open(filename,'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
```

Questions
```
# Load lambda genome
genome = readGenome('lambda_virus.fa')

# Q1: how many times does "AGGT" or its reverse complement "ACCT" occur?
occ1 = naive_with_rc("AGGT", genome)
print("Q1: AGGT/ACCT occurrences =", len(occ1))

# Q2: how many times does "TTAA" or its reverse complement occur?
occ2 = naive_with_rc("TTAA", genome)
print("Q2: TTAA occurrences =", len(occ2))

# Q3: leftmost occurrence offset of "ACTAAGT" or its reverse complement
occ3 = naive_with_rc("ACTAAGT", genome)
print("Q3: leftmost ACTAAGT or RC offset =", min(occ3))

# Q4: leftmost occurrence offset of "AGTCGA" or its reverse complement
occ4 = naive_with_rc("AGTCGA", genome)
print("Q4: leftmost AGTCGA or RC offset =", min(occ4))
```

🧬 Question 5 & 6 — Approximate Matching (≤2 mismatches)

🔹 Step 1. Define naive_2mm function
```
def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        mismatches = 0
        for j in range(len(p)):
            if t[i + j] != p[j]:
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            occurrences.append(i)
    return occurrences
```
✅ Explanation
This version of the naïve algorithm allows up to 2 mismatches.
It loops through each position in the genome t, compares characters of pattern p, and counts mismatches.
If mismatches ≤ 2, that position is recorded as an approximate match.

🔹 Step 2. Search the lambda virus genome
```
# Q5: How many times does TTCAAGCC occur in the Lambda virus genome when allowing up to 2 mismatches?
pattern = 'TTCAAGCC'
matches = naive_2mm(pattern, genome)
print("Number of occurrences with ≤2 mismatches:", len(matches))

# Q6: What is the offset of the leftmost occurrence of AGGAGGTT
pattern = 'AGGAGGTT'
matches = naive_2mm(pattern, genome)
print("Number of approximate matches:", len(matches))
print("Leftmost occurrence (offset):", min(matches))
```

🧪 Question 7 — Detecting a Poor Sequencing Cycle in a FASTQ File
```
seqs, quals = readFastq('ERR037900_1.first1000.fastq')

import numpy as np
n_bases = len(quals[0])
avg_qual_per_cycle = np.zeros(n_bases)

for q in quals:
    for i in range(n_bases):
        avg_qual_per_cycle[i] += phred33ToQ(q[i])

avg_qual_per_cycle /= len(quals)
min_cycle = np.argmin(avg_qual_per_cycle)
min_value = avg_qual_per_cycle[min_cycle]
print(f"The lowest average quality is {min_value:.2f} at cycle {min_cycle}")
```
```
#PLOT:
seqs, quals = readFastq('ERR037900_1.first1000.fastq')

import numpy as np
n_bases = len(quals[0])
avg_qual_per_cycle = np.zeros(n_bases)

for q in quals:
    for i in range(n_bases):
        avg_qual_per_cycle[i] += phred33ToQ(q[i])

avg_qual_per_cycle /= len(quals)

import matplotlib.pyplot as plt
plt.plot(avg_qual_per_cycle)
plt.xlabel('Cycle (offset)')
plt.ylabel('Average quality score')
plt.title('Average read quality by sequencing cycle')
plt.show()
```
 
