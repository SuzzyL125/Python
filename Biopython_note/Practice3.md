## Edit distance between pattern and excerpt of human chromosome1
```
def editDistanceMatch(p, t):
    # Create distance matrix
    D = []
    for i in range(len(p)+1):
        D.append([0]*(len(t)+1))
    
    # Initialize first row and column
    for i in range(len(p)+1):
        D[i][0] = i
    for j in range(len(t)+1):
        D[0][j] = 0  # <-- allow matching anywhere in t

    # Fill in matrix
    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if p[i-1] == t[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    
    # The best match will be the smallest value in the last row
    return min(D[-1])
```
read data: https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/chr1.GRCh38.excerpt.fasta
```
def readGenome(filename):
    genome = ''
    with open(filename,'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
t = readGenome('chr1.GRCh38.excerpt.fasta')
```

q1: What is the edit distance of the best match between pattern GCTGATCGATCGTACG and the excerpt of human chromosome 1?
```
p = 'GCTGATCGATCGTACG'
print(editDistanceMatch(p,t))
# 3
```
q2: What is the edit distance of the best match between pattern GATTTACCAGATTGAG and the excerpt of human chromosome 1?
```
p = 'GATTTACCAGATTGAG'
print(editDistanceMatch(p,t))
# 2
```
Hints: Notice that this is for approximate match, use p to compare with string of t. In naive match of edit distance, it would be p compared with whole t.

## Overlap
```
def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match
def overlap_all_pairs(reads, k, map={}):
    def get_kmers(read, k):
        res = set()
        for i in range(0, len(read)-k+1):
            res.add(read[i:i+k])
        return res
    for read in reads:
        kmers = get_kmers(read, k)
        for kmer in kmers:
            if not kmer in map.keys():
                map[kmer] = set()
            map[kmer].add(read)
    pairs = []
    for head in reads:
        kmer = head[-k:]
        candidates = map[kmer]
        for tail in candidates:
            if (not head == tail and overlap(head, tail, k)):
                pairs.append((head, tail))
    return pairs
```
read data: https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR266411_1.for_asm.fastq
```
def readFastq(filename):
  sequences = []
  qualities = []
  with open(filename) as fh:
    while True:
      fh.readline()
      seq = fh.readline().rstrip()
      fh.readline()
      qual = fh.readline().rstrip()
      if len(seq) == 0:
        break
      sequences.append(seq)
      qualities.append(qual)
  return sequences, qualities
seqs, quals = readFastq('ERR266411_1.for_asm.fastq')
```
q3: How many distinct pairs of reads overlap
```
pairs = overlap_all_pairs(seqs, 30)
print(len(pairs))
# 904746
```
q4: how many reads have a suffix involved in an overlap?
```
print(len(set(pair[0] for pair in pairs)))
# 7161
```
