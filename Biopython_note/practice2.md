---
Data source:http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/chr1.GRCh38.excerpt.fasta
```
t = readGenome('chr1.GRCh38.excerpt (1).fasta')
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
```

```
# ---- 1&2: How many alignments does the naive exact matching 
# algorithm try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG

def naive_with_counts(p, t):
    occurrences = []
    alignments = 0
    comparisons = 0
    for i in range(len(t) - len(p) + 1):
        alignments += 1
        match = True
        for j in range(len(p)):
            comparisons += 1
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences, alignments, comparisons

occurrences, alignments, comparisons = naive_with_counts(p, t)
print("Naive alignments:", alignments)
print("Naive comparisons:", comparisons)
print("Match positions:", occurrences)
```
---
- Naive alignments: 799954
- Naive comparisons: 984143
- Match positions: [56922]
---

```
# ---- 3: How many alignments does Boyer-Moore try when matching the string

p_bm = BoyerMoore(p) #from bm_preproc.py

def boyer_moore_with_counts(p, p_bm, t):
    occurrences = []
    alignments = 0
    i = 0
    while i < len(t) - len(p) + 1:
        alignments += 1
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
    return occurrences, alignments

occurrences_bm, alignments_bm = boyer_moore_with_counts(p, p_bm, t)
print("Boyer-Moore alignments:", alignments_bm)
print("Match positions (BM):", occurrences_bm)
```

---
- Boyer-Moore alignments: 127974
- Match positions (BM): [56922]
---

```
# ---- 4: How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, occur with up to 2 substitutions in the excerpt of 
# human chromosome 1?

# ---------- pigeonhole with counting index hits ----------
def pigeonhole_index_hit_counts(p, t, n, index):
    """Return:
       - total_index_hits: total number of index hits (sum over all partitions)
       - per_segment_counts: list of hit counts per partition
       - candidate_positions: set of candidate offsets produced by index (may include duplicates resolved by set)
    """
    seg_len = int(round(len(p) / (n + 1)))  # for len=24, n=2 -> seg_len=8
    total_index_hits = 0
    per_segment_counts = []
    candidate_positions = []  # collect raw candidate offsets (offset = index_hit - start)
    for i in range(n+1):
        start = i * seg_len
        # query k-mer beginning at start (we assume index.k == seg_len == 8)
        kmer = p[start:start + index.k]
        hits = index.query(kmer)
        per_segment_counts.append(len(hits))
        total_index_hits += len(hits)
        # convert each index hit position (where k-mer occurs in t) -> candidate alignment offset
        for hit_pos in hits:
            offset = hit_pos - start
            # keep offsets even if negative or out of range — they'll be filtered later if needed
            candidate_positions.append(offset)
    return total_index_hits, per_segment_counts, candidate_positions

# ---------- optional: verify approximate matches by full verification
# ----------
def verify_candidates(p, t, candidate_positions, n):
    matches = set()
    for offset in candidate_positions:
        if offset < 0 or offset + len(p) > len(t):
            continue
        # count mismatches
        mismatches = 0
        for j, ch in enumerate(p):
            if t[offset + j] != ch:
                mismatches += 1
                if mismatches > n:
                    break
        if mismatches <= n:
            matches.add(offset)
    return sorted(matches)



k = 8
index = Index(t, k)
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
total_hits, per_seg, candidates = pigeonhole_index_hit_counts(p, t, 2, index)
print("Per-segment index hit counts (segments 0..2):", per_seg)
print("Total index hits (sum over segments):", total_hits)
print("Total raw candidate positions (with duplicates):", len(candidates))

# Optional: verify and deduplicate to get actual approximate matches (<=2
# mismatches)
matches = verify_candidates(p, t, candidates, 2)
print("Verified approximate matches (<=2 mismatches), unique positions:", len(matches))
print("Positions:", matches)
```

---
- Per-segment index hit counts (segments 0..2): [13, 17, 60]
- Total index hits (sum over segments): 90
- Total raw candidate positions (with duplicates): 90
- Verified approximate matches (<=2 mismatches), unique positions: 19
- Positions: [56922, 84641, 147558, 160162, 160729, 191452, 262042, 273669, 364263, 421221, 429299, 465647, 551134, 635931, 657496, 681737, 717706, 724927, 747359]
---
