import random

######################################
###### Manipulating DNA strings ######
######################################
# Create a random seq
seq = ''.join([random.choice('ACGT') for _ in range(10)])
print(seq)
#GCCTATGTTA

seq[:3]
#GCC

seq[:-3]
#GCCTATG

seq[-3]
#T

def longestcommonprefix(s1,s2):
  i = 0
  while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
    i += 1
  return s1[:i]
longestcommonprefix('ACCATGT','ACCAGAC')
# ACCA

def match(s1,s2):
  if not len(s1) == len(s2):
    return False
  for i in range(len(s1)):
    if not s1[i] == s2[i]:
      return False
  return True
match('ACCATGT','ACCAGAC')
# False


def reverseComplement(s):
  complement = {'A':'T','T':'A','C':'G','G':'C'}
  t = ''
  for base in s:
    t = complement[base] + t
  return t
reverseComplement('ACGCTAGC')
# GCTAGCGT

######################################
###### Working with sequencing reads ######
######################################
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
seqs, quals = readFastq('1_control_psbA3_2019_minq7.fastq')
print(quals[:2])
# ["##$&$&/035881()'$0&*('-.=;685()$.%($'%%&#&)+..0,&+&%.-/+,%&()$3:0&@09BF=>CC8(78029F7=<=)+@+.6CCFFC@-8%2579<B8;88412134,,;:8./,#1#&(%((09;B=??48<=<@79*-:B540,8=B=444:<571-B5=ED2.56;110.5+,*)%%*", 
# ',&\'%/..0013618\'.\'(*\'#(%#$&&%%,-//&#\'$&&\'&%$(\'\'\'4+9;7<87756-1-+$$,3%%;)-%$%%$&)\'\'#067897:9$%)(<AC?,(**$&\'.6:<=394./41*,12((:?;3/9\'(-4<)=99D99BFAC9:;588+.&+&%()%-(59,,($,12+,,*6(-))&-$&0\'%)%&%\'""%(0((']

def phred33ToQ(qual):
    return ord(qual) - 33
phred33ToQ('#')
# 2

# Creating a Quality Score Histogram
def createHist(qualities):
    hist = [0]*53 #can customize this count
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist
h=createHist(quals)
print(h)
# [0, 1786, 29817, 71301, 70813, 58077, 48355, 41824, 
# 37500, 33766, 31333, 29527, 27735, 26588, 25171, 24274, 23521, 22669, 22459, 21796, 21316, 20863, 20304, 19680, 19114, 18263, 17350, 16411, 14863, 13370, 11701, 9984, 8404, 6856, 5179, 3869, 2740, 1871, 1175, 803, 458, 304, 189, 119, 69, 41, 24, 16, 14, 10, 3, 3, 0]

import matplotlib.pyplot as plt
plt.bar(range(len(h)),h)
plt.show()

# This function computes GC content at each base position across all reads:
def findGCByPos(reads):
    gc = [0] * 500
    totals = [0] * 500
    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
    return gc
gc=findGCByPos(seqs)
plt.plot(range(len(gc)),gc)
plt.show()

# gc[i]: counts the number of G or C bases at position i.
# totals[i]: counts the total number of bases at position i.

The GC ratio at each position = gc[i] / totals[i].


# Counting Nucleotide Frequencies
import collections
count = collections.Counter()
for seq in seqs:
    count.update(seq)
print(count)
# Counter({'T': 264026, 'A': 217776, 'G': 211076, 'C': 190800})


######################################
###### Matching ######################
######################################
# naïve exact string matching algorithm.
def naive(p,t):
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
t = 'AGCTTAGATAGC'
p = 'AGC'
naive(p,t)
# [0, 9]

# It searches for all occurrences of a pattern string p inside a longer text t. 
# For each possible starting position in t, it compares every character of p to the corresponding character in t. If all characters match, it records that position as a match. 
# Finally, it returns a list of all positions (offsets) where the pattern appears in the text.

# Approximate match (allow up to 2 mismatch
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
t = 'AGCTTAGATAGC'
p = 'AGC'
print(naive_2mm(p,t))
# [0, 5, 7, 9]

# reads a genome sequence from a FASTA file.
def readGenome(filename):
    genome = ''
    with open(filename,'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
genome = readGenome('phix.fa.1')

import random
# simulates sequencing reads from a genome
def generateReads(genome, numReads,readLen):
    reads = []
    for _ in range(numReads):
        start = random.randint(0,len(genome)-readLen)-1
        reads.append(genome[start:start+readLen])
    return reads
reads = generateReads(genome, 100, 100)
numMatched=0
for r in reads:
    matches=naive(r,genome)
    if len(matches)>0:
        numMatched += 1
print('%d / %d reads matched exactly!' % (numMatched,len(reads)))
# loops through all generated reads and uses the naive function to find where each read exactly matches the genome. 
# If a read matches at least once, it’s counted as a successful match. Finally, the program prints how many of the simulated reads had at least one perfect match in the genome (e.g., “85 / 100 reads matched exactly!”).

# use reverseComplement and naive match
# naive_with_rc(p, t):Finds matches of both p and its reverse complement.If they are identical (e.g., palindromic sequence), it avoids double-counting.Returns a sorted list of unique match offsets.
def naive_with_rc(p, t):
    rc_p = reverseComplement(p)
    occurrences = naive(p, t)

    # If reverse complement is different, search for it too
    if rc_p != p:
        rc_occurrences = naive(rc_p, t)
        # Combine results and remove duplicates
        occurrences = sorted(list(set(occurrences + rc_occurrences)))

    return occurrences

# Example usage
t = "ACGTACGTGACG"
p = "ACG"
print("Pattern matches:", naive_with_rc(p, t))
# Pattern matches: [0, 1, 4, 5, 9]


######################################
###### Boyer-Moore matching ######################
######################################
# Boyer-Moore preprocessing:http://d28rh4a8wq0iu5.cloudfront.net/ads1/code/bm_preproc.py
# preprocessing: bm_preproc.py
def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences

######################################
###### k-mer index ######################
######################################
import bisect

# For each position in the sequence: It extracts a substring of length k (a k-mer). Stores a tuple (k-mer, offset) in a list.
# Then it sorts all tuples alphabetically by k-mer, so binary search can be used
class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    # Takes a pattern p (e.g., 'TCTA'). Extracts its first k-mer, here 'TC'. Uses bisect_left to binary search the sorted list for where 'TC' first appears.
    # Then iterates forward to collect all matching positions (offsets) in the genome.
    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
def queryIndex(p,t,index):
  k = index.k
  offsets = []
  # For each possible match (i) returned by index.query(p): It compares the rest of the pattern (p[k:]) with the corresponding substring in t.
  # If the match is exact, adds i to the result list offsets.
  for i in index.query(p):
    if p[k:] == t[i+k:i+len(p)]:
      offsets.append(i)
  return offsets
  
t = 'GCTACGATCTAGAATCTA'
p = 'TCTA'
index = Index(t,2)
print(queryIndex(p,t,index))
# [7,14]

######################################
###### Approximate matching ######################
######################################
# naïve hamming
def naiveHamming(p,t,maxDistance):
    occurrences = []
    for i in xrange(len(t)-len(p)+1):
        nmm = 0
        match = True
        for j in xrange(len(p)):
            if not t[i+j] == p[j]:
                nmm += 1
                if nmm > maxDistance:
                    break
        if nmm <= maxDistance:
            occurrences.append(i)
    return occurrences

def approximate_match(p,t,n):
    segment_length = int(round(len(p)/(n+1)))
    all_matches =  set()
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length,len(p))
        p_bm = BoyerMoore(p[start:end],alpha='ACGT')
        matches = boyer_moore(p[start:end],p_bm,t)
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0,start):
                if not p[j] == t[m-start+j]:
                    mismatches = 1
                    if mismatches > n:
                        break
            for j in range(end,len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m-start)
    return list(all_matches)
t = 'CACTTAATTTG'
p = 'AACTTG'
print(approximate_match(p,t,2))
# [0,5]

######################################
###### Edit Distance ######################
######################################
def editDistance(x, y):
    # Create distance matrix: (len(x)+1) × (len(y)+1)
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
      
    # Initialize first row and column of matrix: D[i][0]-deletion/D[0][j]-insert
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
      
    # Fill in the rest of the matrix, distHor:Horizontal move-Insertion
    # D[0][0]=0
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]
x = "cat"
y = "but"
print(editDistance(x,y))
# 2

### The above code can be used when x and y of different length, but the result will target to let x and y have same length.
### The following code will be used if x and y of different length and check x meet any part of y and calculate it
### best approximate match of a short pattern p anywhere in a long text t
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
x = "cat"
y = "bat man"
print(editDistanceMatch(x,y))
# 1
print(editDistance(x,y))
# 5



######################################
###### Global Alignment ######################
######################################
alphabet = ['A','C','G','T']
# Match	A→A	0	No penalty
# Transition	A→G	2	Small penalty (purine↔purine)
# Transversion	A→C	4	Higher penalty (purine↔pyrimidine)
# Gap	A→-	8	Very high penalty (insertion/deletion)
score = [[0,4,2,4,8],[4,0,4,2,8],[2,4,0,4,8],[4,2,4,0,8],[8,8,8,8,8]]
def globalAlignment(x, y):
    # Create distance matrix: (len(x)+1) × (len(y)+1)
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
      
    # Initialize first row and column of matrix: D[i][0]-deletion/D[0][j]-insert
    for i in range(len(x)+1):
        D[i][0] = D[i-1][0] + score[alphabet.index(x[i-1])][-1]
    for i in range(len(y)+1):
        D[0][i] = D[0][i-1] + score[-1][alphabet.index(y[i-1])]
      
    # Fill in the rest of the matrix, distHor:Horizontal move-Insertion
    # D[0][0]=0
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + score[-1][alphabet.index(y[j-1])]
            distVer = D[i-1][j] + score[alphabet.index(x[i-1])][-1]
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]
x = "ACG"
y = "ACCT"
print(globalAlignment(x,y))
# 20


######################################
###### Overlap ######################
######################################
def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

# Limiation: It works for one pair only. In real genome assembly, we have thousands (or millions) of reads — we must find all overlapping pairs.
# If you called overlap(a, b, k) for every possible pair, that’s O(n²) comparisons — computationally impossible for large n.

###### overlap_all_pairs ######################
# Improving Efficiency Using k-mer Indexing, If k-mer "GTCCTA" doesn’t occur anywhere else in the dataset, we can skip checking overlaps between a and all other reads.
# finds all pairs of reads that overlap by at least k bases, without checking every pair explicitly.
def overlap_all_pairs(reads, k, map={}):
  
    # Extracts all k-mers (substrings of length k) from a read. e.g. get_kmers('AGTACCGT', 3) → {AGT, GTA, TAC, ACC, CCG, CGT}
    def get_kmers(read, k):
        res = set()
        for i in range(0, len(read)-k+1):
            res.add(read[i:i+k])
        return res

    # Build a k-mer index (hash map). {'CCG': {'AGTACCGT', 'CCGTTGA'}, 'CGT': {'AGTACCGT', 'CCGTTGA'}, ...}
    for read in reads:
        kmers = get_kmers(read, k)
        for kmer in kmers:
            if not kmer in map.keys():
                map[kmer] = set()
            map[kmer].add(read)
    pairs = []

    # Find overlap candidates efficiently. It compares only with candidates that share the same k-mer (suffix of head = prefix of tail).
    for head in reads:
        kmer = head[-k:]  # the last k bases of head (its suffix)
        candidates = map[kmer] # all reads that share that k-mer
        for tail in candidates:
            if (not head == tail and overlap(head, tail, k)):
                pairs.append((head, tail))
    return pairs
    # Returns a list of read pairs (head, tail) that overlap by ≥ k.
