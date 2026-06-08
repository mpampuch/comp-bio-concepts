# Bioinformatics Sequence Modelling — Practice Problems (Julia)

>[!NOTE]
> First

A lecture-by-lecture problem set that builds up from probability fundamentals to
full Hidden Markov Models for gene finding. Each lecture has:

- **Concept check** — short questions to confirm you understood the lecture.
- **Math problems** — pencil-and-paper derivations / calculations.
- **Julia problems** — code to implement and test.

> **How to use this set**
>
> - Work the math problems first; they make the code obvious.
> - For the Julia problems, write the function _and_ a small test set
>   (`using Test`) before moving on. The expected answers in the math
>   problems double as your test oracles.
> - Each lecture's code builds on the previous one. By Lecture 10 you will
>   have a working gene-finding HMM assembled from pieces you wrote earlier.
> - Suggested project layout:
>   ```
>   SeqModel/
>     Project.toml          # add Test, (optionally) StatsBase, Plots
>     src/SeqModel.jl
>     test/runtests.jl
>   ```

A handful of problems are marked **(★ stretch)** — they go a little beyond the
lecture or require more design thinking.

---

## Lecture 1 — Intro: why probabilistic sequence models?

**Topics:** motivation (gene finding, ~20k genes), CpG dinucleotides and
methylation, CpG islands as a genomic landmark, why scores should be
_probabilities_ (a common language for comparing/combining models).

### Concept check

1. Why is "count how many times I saw this exact string inside vs. outside a
   CpG island" a poor strategy as the string gets long?
2. The lecture argues scores should be probabilities, not arbitrary confidence
   numbers. Give two concrete advantages of using probabilities.
3. In a uniformly random DNA sequence, what is the probability that a given
   adjacent pair of bases is exactly `CG`? Roughly how far apart would `CG`
   dinucleotides be on average?

### Math problems

1. Under a uniform i.i.d. model over `{A,C,G,T}`, compute the expected number of
   `CG` dinucleotides in a sequence of length `L`. (Answer for `L = 1000`.)
2. CpG islands are defined by `CG` being _denser_ than the background. If
   background spacing of `CG` is ~16 bp, and inside an island the spacing is
   ~4 bp, what is the ratio of `CG` densities (island : background)?

### Julia problems

1. **`gc_content(seq)`** — return the fraction of characters in `seq` that are
   `C` or `G`. Test on `"GATTACA"` (expect `1/7`) and `"GCGCGC"` (expect `1.0`).
2. **`count_dinucleotide(seq, dinuc)`** — count overlapping occurrences of a
   2-character pattern. Test: `count_dinucleotide("CGCGCG", "CG")` should be `3`.
3. **`cpg_count(seq)`** — count `CG` dinucleotides. Generate a random length-`L`
   DNA string and check empirically that the count is close to your Lecture-1
   math answer (average over many random strings).

```julia
# starter
const DNA = ('A','C','G','T')
random_dna(L) = String(rand(DNA, L))
```

---

## Lecture 2 — Probability review

**Topics:** sample space `Ω`, events, `P(A)`, joint `P(A,B)`, conditional
`P(A|B) = P(A,B)/P(B)`, multiplication rule `P(A,B)=P(A|B)P(B)`, repeated
application of the multiplication rule (chain of conditionals + one marginal),
independence `P(A,B)=P(A)P(B)`.

### Concept check

1. State the multiplication rule and explain the "zooming in / zooming out"
   intuition for conditional probability.
2. What does it mean, in words, for two events to be independent?

### Math problems (two-dice running example)

Sample space = rolling two fair dice (36 equally likely outcomes).
Let `A` = "die 1 is odd", `B` = "die 2 is even".

1. Compute `P(A)`, `P(B)`, `P(A,B)`, and `P(A|B)`. Confirm `A` and `B` are
   independent.
2. Write the joint `P(A,B,C,D)` decomposed fully via the multiplication rule
   into three conditionals and one marginal.
3. Give the general chain rule: for a string `x_1 x_2 … x_k`, write `P(x_1,…,x_k)`
   as a product of conditionals and one marginal (no Markov assumption yet).

### Julia problems

1. **`dice_space()`** — return all 36 `(d1,d2)` tuples.
2. **`prob(event, space)`** — given a predicate `event(outcome)::Bool` and a
   vector of equally-likely outcomes, return the naive probability (fraction of
   outcomes satisfying the predicate).
3. **`joint(a, b, space)`** and **`conditional(a, b, space)`** built on top of
   `prob`. Verify numerically: `P(A)=0.5`, `P(A,B)=0.25`, `P(A|B)=0.5`, and that
   `joint(a,b,space) ≈ prob(a,space)*prob(b,space)` (independence check).

```julia
# starter
isodd_d1((d1,d2))  = isodd(d1)
iseven_d2((d1,d2)) = iseven(d2)
```

---

## Lecture 3 — The Markov assumption

**Topics:** joint probabilities of strings are hard to estimate by raw counting
once `k` is large (exponential blow-up); the Markov assumption as a _locality_
assumption; conditional independence: `x_i ⫫ (everything before) | x_{i-1}`;
"the future is independent of the past given the present"; it is an _assumption_,
sometimes good (board game), sometimes bad (baseball season), sometimes
in-between (weather).

### Concept check

1. Rewrite `P(x_k | x_{k-1}, x_{k-2}, …, x_1)` under the (first-order) Markov
   assumption.
2. For each scenario, say whether the first-order Markov assumption is
   reasonable and why: (a) predicting tomorrow's baseball win from today's
   result only; (b) predicting the next Parcheesi board state from the current
   one; (c) predicting tomorrow's weather from today's only.

### Math problems

1. Starting from the full chain rule for `P(x_1,…,x_k)` (Lecture 2), apply the
   first-order Markov assumption and write the simplified factorization.
2. **Data-need blow-up.** Suppose you estimate `P(x)` for a length-`k` string by
   raw counting and you want at least ~10 expected observations of _each_
   possible length-`k` string. How large a dataset do you need as a function of
   `k` for a 4-letter alphabet? Evaluate for `k = 5, 10, 15`.

### Julia problems

1. **`naive_string_prob(query, examples)`** — estimate `P(query)` as
   `count(query in examples) / length(examples)` for a vector of equal-length
   strings. Demonstrate it returns `0.0` for a `query` never seen (the failure
   mode that motivates Markov).
2. **(★ stretch) `markov_factorization_terms(seq)`** — return the _list of
   conditional/marginal terms_ (as `(current, given)` pairs, with `given === nothing`
   for the marginal) that the first-order Markov assumption produces for `seq`.
   E.g. `"GATC"` → `[('G', nothing), ('A','G'), ('T','A'), ('C','T')]`. This is
   the scaffold you'll evaluate with a real Markov chain in Lecture 4.

---

## Lecture 4 — Markov chains (part 1): building & scoring

**Topics:** first-order Markov chain as the factorization
`P(x) = P(x_1)·∏_{i≥2} P(x_i | x_{i-1})`; estimating transition probabilities
from counts `P(b | a) = count(ab) / count(a·)`; the 4×4 transition table (16
dinucleotide conditionals); rows sum to 1, columns need not; interpretation
(outside-island `P(G|C)` is small); Markov chains as probabilistic automata;
worked example computing `P("GATC")`.

> **Notation reminder:** `P(G|C)` means "next base is `G` given previous base was
> `C`", i.e. the dinucleotide `CG` appears with `C` first. The letter on the
> _right_ of the bar occurs _first_ in the string.

### Math problems

Use this trained inside-island transition table `A[prev, curr]` (rows = previous
base, columns = current base; rows sum to 1):

|       | A    | C    | G    | T    |
| ----- | ---- | ---- | ---- | ---- |
| **A** | 0.18 | 0.27 | 0.43 | 0.12 |
| **C** | 0.17 | 0.37 | 0.27 | 0.19 |
| **G** | 0.16 | 0.34 | 0.38 | 0.12 |
| **T** | 0.08 | 0.36 | 0.38 | 0.18 |

Assume a uniform marginal `P(x_1) = 1/4`.

1. Verify each row sums to 1.
2. Compute `P("GATC")` by hand: `P(G)·P(A|G)·P(T|A)·P(C|T)`.
3. Compute `P("CGCG")`. Compare to `P("CCCC")` and interpret.

### Julia problems

1. **`transition_counts(seqs, alphabet=DNA)`** — return a `4×4` integer matrix of
   dinucleotide counts over a collection of training sequences. Index bases via a
   `Dict('A'=>1,'C'=>2,'G'=>3,'T'=>4)`.
2. **`transition_matrix(counts)`** — normalize counts into conditional
   probabilities (row-normalize). Add a test asserting every row sums to `1`.
3. **`seq_prob(seq, A; init=fill(0.25,4))`** — compute the joint probability of
   `seq` under transition matrix `A` and marginal `init`. Test against your
   hand-computed `P("GATC")` from the math section (use the table above).
4. **(★ stretch)** Train `transition_matrix` on a synthetic "island" generator
   (bias toward emitting `G` after `C`) and confirm the learned `A[C,G]` entry is
   noticeably higher than for a uniform-random generator.

---

## Lecture 5 — Markov chains (part 2): log-space & the CpG classifier

**Topics:** numerical **underflow** when multiplying many probabilities; fix by
working in **log space** (products → sums of logs); building two chains (inside
vs. outside) and classifying by the **log-likelihood ratio**
`S(x) = log( P_in(x) / P_out(x) ) = Σ_i [ log A_in(x_{i-1},x_i) − log A_out(x_{i-1},x_i) ]`;
the combined single table of per-dinucleotide log-ratios; decision rule
`S(x) > 0 ⇒ inside`; sanity checks (`"CGCGCG…"` scores high positive, a
`CG`-free string scores negative); train/test on different chromosomes; the two
overlapping histograms show the Markov assumption works well here.

### Math problems

1. Show algebraically that
   `log( P_in(x)/P_out(x) ) = Σ_{i=2}^{k} ( log A_in(x_{i-1},x_i) − log A_out(x_{i-1},x_i) )`
   (marginals cancel/ignored). Why does this let you store a _single_ table of
   log-ratios instead of two tables?
2. If `P_in(x)` is exactly twice `P_out(x)`, what is the ratio and the
   (natural / base-2) log-ratio? What if `P_out` is twice `P_in`?
3. **Underflow.** Roughly how many factors of `0.5` can a Float64 multiply before
   underflowing to `0`? (Recall `0.5^1100 == 0.0` in IEEE double.) Explain why
   the log-space sum does not have this problem.

### Julia problems

1. **`logspace_seq_logprob(seq, A; init)`** — same as `seq_prob` but returns a
   **log** probability by summing logs. Test that `exp(logspace_seq_logprob(...))`
   `≈ seq_prob(...)` for a short sequence.
2. Demonstrate underflow directly: show `prod(fill(0.5, 1100)) == 0.0` while
   `sum(fill(log(0.5), 1100))` is a large finite negative number.
3. **`logratio_table(A_in, A_out)`** — return the element-wise
   `log2.(A_in) .- log2.(A_out)` combined table.
4. **`cpg_score(seq, logratio)`** — sum the per-dinucleotide log-ratios along
   `seq`; positive ⇒ inside. Tests:
   - `cpg_score("CGCGCGCGCG", logratio) > 0`
   - `cpg_score("ATATATATAT", logratio) < 0`
5. **(★ stretch) evaluation.** Given labeled test sequences, compute the score
   for each, and report classification accuracy at threshold `0`. Optionally plot
   two overlapping histograms (inside vs. outside scores) like the lecture.

---

## Lecture 6 — Hidden Markov Models (part 1): definition

**Topics:** hidden state vs. observable emissions (baseball-outside-the-stadium
analogy); the trellis/lattice diagram; states `Q`, emissions `Σ`, paths;
vertical edges = emission probabilities `E[state, symbol]`, horizontal edges =
transition probabilities `A[state, state]`, plus initial vector `I`; the joint
`P(p_1..p_n, x_1..x_n) = I[p_1]·E[p_1,x_1]·∏_{i≥2} A[p_{i-1},p_i]·E[p_i,x_i]`;
the **occasionally dishonest casino** (fair/loaded coin) and filling in `A`, `E`,
`I` from the word problem.

### The occasionally dishonest casino (use throughout L6–L7)

- States `Q = {Fair, Loaded}`.
- Emissions `Σ = {H, T}`.
- Switch coins with prob `0.4` ⇒ stay with prob `0.6`.
- Fair coin: `P(H)=P(T)=0.5`. Loaded: `P(H)=0.8, P(T)=0.2`.
- Initial: `P(Fair)=P(Loaded)=0.5`.

### Math problems

1. Write out the full `A` (2×2), `E` (2×2), and `I` (length 2) from the word
   problem above. Confirm each row of `A` and `E` sums to 1.
2. Compute the joint probability of path `p = (Fair, Fair, Loaded)` with emissions
   `x = (H, H, T)`:
   `I[F]·E[F,H]·A[F,F]·E[F,H]·A[F,L]·E[L,T]`. Give the number.
3. How many numbers total define this HMM? Generalize: an HMM has
   `|Q|² + |Q||Σ| + |Q|` parameters — evaluate for this casino.

### Julia problems

1. Encode the casino HMM as Julia structures:
   ```julia
   struct HMM
       states::Vector{Symbol}
       symbols::Vector{Char}
       A::Matrix{Float64}   # |Q| x |Q|, A[i,j] = P(state j | state i)
       E::Matrix{Float64}   # |Q| x |Σ|, E[i,k] = P(symbol k | state i)
       I::Vector{Float64}   # |Q|, initial distribution
   end
   ```
   Build `casino::HMM`. Add tests that all rows of `A` and `E` sum to `1` and `I`
   sums to `1`.
2. **`joint_prob(hmm, path, emissions)`** — multiply `I·E·A·E·…`. Test against
   your hand answer from math problem 2.
3. **`joint_logprob(hmm, path, emissions)`** — the log-space version; check it
   agrees with `log(joint_prob(...))`.
4. **(★ stretch) `simulate(hmm, n)`** — sample a path of length `n` and its
   emissions from the HMM (use the initial vector, then transitions, emitting at
   each step). Useful for generating test data in Lecture 7.

---

## Lecture 7 — HMM decoding: the Viterbi algorithm

**Topics:** decoding = find `p* = argmax_p P(p | x) = argmax_p P(p, x)` (the
denominator `P(x)` is constant in `p`, so drop it); extending the emission string
by one base can ripple backward through the best path, so greedy doesn't work;
**optimal substructure** ⇒ dynamic programming; the Viterbi DP matrix
`V[state, i]` = max joint prob of any path ending in `state` at position `i`;
recurrence `V[s,i] = E[s, x_i] · max_{s'} ( V[s', i-1] · A[s', s] )`; first column
uses the initial vector `V[s,1] = I[s]·E[s,x_1]`; keep backpointers; **traceback**
from the best final cell; complexity `O(n·|Q|²)` time, traceback `O(n)`.

### Math problems (use the casino HMM, but with stickier coins for this part:

`stay = 0.9`, `switch = 0.1`, everything else as in L6)

1. For emission string `x = (H, H)`, fill in the `2×2` Viterbi matrix by hand
   (one column per emission, one row per state). Show the `max` you take in
   column 2 for each state and which previous state "wins".
2. Continue to `x = (H, T, T, H, H)` (you may do this in Julia, but write out the
   recurrence you are evaluating). What is `p*`?
3. **Ripple effect (conceptual).** The lecture shows that for `x = HTTHH` the best
   path is `FFFFF`, but for `x = HTTHHH` (one more `H`) the best path flips to all
   `L`. Explain in one or two sentences why decoding can't just append a state to
   the previous answer.
4. Why is each Viterbi cell `O(|Q|)` work, giving `O(n|Q|²)` overall?

### Julia problems

1. **`viterbi(hmm, emissions) -> (path, logprob)`** — implement in **log space**
   (work with `log.(A)`, `log.(E)`, `log.(I)`; sums instead of products; keep an
   integer backpointer matrix; traceback from the argmax of the last column).
   Return the most likely path as a `Vector{Symbol}` and its log-probability.
2. Tests:
   - On the sticky casino with `x = "HH"`, the decoded path should be
     `[:Loaded, :Loaded]` (loaded coin favors heads). Verify against your
     hand-filled matrix.
   - Reproduce the lecture's ripple example qualitatively: decode `"HTTHH"` and
     `"HTTHHH"` on an HMM tuned so the answers differ, and assert they are _not_
     simple extensions of each other.
3. **`viterbi_matrix(hmm, emissions)`** — return the full DP matrix (log values)
   and the backpointer matrix so you can inspect/print it like the lecture's
   trellis.
4. **Round-trip test.** Use `simulate(hmm, n)` from L6 to generate `(true_path,
emissions)`, run `viterbi`, and report the fraction of positions where the
   decoded state matches the true state.
5. **(★ stretch) complexity check.** Time `viterbi` for growing `n` and growing
   `|Q|`; confirm the empirical scaling is consistent with `O(n|Q|²)`.

---

## Lecture 8 — Applying HMMs to CpG-island finding

**Topics:** first (naive) design — states `{Inside, Outside}`, emissions
`{A,C,G,T}`; train `A` from the `I/O` label string and `E` from base-vs-label
counts; **the flaw**: this model never counts _dinucleotides_, so a run of `C`s
(no `CG`!) gets wrongly called an island; fix by **pushing the previous
nucleotide into the state space** — states become the Cartesian product
`{A,C,G,T} × {I,O}` (8 states); now transitions count dinucleotides; side effect:
emissions become **deterministic** (1 if the state's letter matches the emitted
symbol, else 0); interpreting the trained `8×8` transition heat-map (4 quadrants;
transitions into islands happen on `C`, out of islands on `G`).

### Math problems

1. **Why the naive model fails.** In the 2-state `{I,O}` model with emissions
   `{A,C,G,T}`, explain why the input `"CCCCCCCC"` can be decoded as _inside_ even
   though it contains no `CG`. Which information is the model missing?
2. **State-space size.** For the fixed model with states `{A,C,G,T}×{I,O}`, write
   the dimensions of `A` and `E`. How many entries of `E` are non-trivial
   (need training) vs. forced to 0/1?
3. Given the small labeled training example
   ```
   seq:   A C G C G T A C G C G ...
   label: O I I I I O O I I I I ...
   ```
   show how you would estimate `P(O | I)` (transition) and, in the _fixed_
   (8-state) model, `P(C_I | T_I)` (a dinucleotide transition). Write the
   count-ratio expression for each.

### Julia problems

1. **Naive model trainer.** Given parallel vectors `seq::String` and
   `labels::Vector{Char}` (`'I'`/`'O'`), build:
   - **`train_io_transitions(labels)`** → `2×2` matrix `P(label | prev label)`.
   - **`train_io_emissions(seq, labels)`** → `2×4` matrix `P(base | label)`.
     Assemble these into an `HMM` (reuse L6 struct) and run your L7 `viterbi`.
2. **Reproduce the bug.** Construct an HMM with the lecture's illustrative
   probabilities (inside emits more `C`/`G`; staying-inside prob `0.8`) and show
   `viterbi` mislabels a `CG`-free `C`-rich stretch as inside.
3. **Fixed (dinucleotide) model.** Implement the 8-state model:
   - state set `[:A_I,:C_I,:G_I,:T_I,:A_O,:C_O,:G_O,:T_O]`.
   - **`fixed_emissions()`** → `8×4` deterministic 0/1 matrix.
   - **`train_fixed_transitions(seq, labels)`** → `8×8` matrix that counts
     dinucleotides _with_ inside/outside labels.
     Decode the same `CG`-free input and assert it is now labeled _all outside_.
4. **(★ stretch) heat-map.** Visualize the trained `8×8` transition matrix
   (e.g. `heatmap` from `Plots`). Confirm the four-quadrant structure and that
   the into-island band sits on the `C_I` column and the out-of-island band on
   the `G_I` row, as described in the lecture.

---

## Lecture 9 — HMMs for gene finding (part 1): biology & signals

**Topics:** central dogma (DNA → mRNA → protein); the diner analogy
(menu = genome, orders = mRNAs, meals = proteins); **codons** (3 nt → 1 amino
acid; 64 codons, 20 amino acids ⇒ redundancy); the genetic code; **start codon**
`ATG`; **stop codons** `TAA`, `TAG`, `TGA`; **exons/introns**, splicing,
**donor** `GT` and **acceptor** `AG` motifs; UTRs; coding sequence length is a
multiple of 3; most `ATG`/`GT`/`AG` occurrences are _not_ real signals — context
matters.

### Math problems

1. Why must a codon be 3 nucleotides? Show that 1 or 2 nucleotides cannot encode
   20 amino acids but 3 can (and quantify the redundancy: 64 vs. 20).
2. A candidate coding sequence has length `L`. State the necessary (not
   sufficient) conditions involving: start codon, stop codon, and `L mod 3`.
3. The HBB example's coding sequence has length 444. Confirm it is a multiple of
   3 and say how many codons (including start and stop) it contains.

### Julia problems

1. **`codons(seq)`** — split a sequence into non-overlapping length-3 codons
   (assume `length(seq) % 3 == 0`; otherwise `@assert`). Test on a length-12
   string.
2. **`translate(seq, code)`** — given a `Dict{String,Char}` genetic code, return
   the amino-acid string. Provide at least the entries you need plus start/stop;
   `'*'` for stop. Test `translate("ATGAAATAA", code)` → `"MK*"` (M=start/Met,
   K=Lys, `*`=stop).
3. **`find_orfs(seq)`** — scan a single reading frame left-to-right and return all
   `(start, stop)` index pairs where an `ATG` is followed (in-frame) by the first
   stop codon. Test on a small constructed sequence with a known ORF.
4. **`is_donor(seq, i)` / `is_acceptor(seq, i)`** — return whether positions
   `i:i+1` equal `"GT"` / `"AG"`. Use them to list candidate splice sites and
   note (in a comment) that most candidates are _not_ real — motivating the HMM.
5. **(★ stretch) `gene_signal_report(seq)`** — given a putative gene region,
   check: starts with `ATG`, ends with a stop codon, interior introns begin `GT`
   and end `AG`, and exon-length-sum `% 3 == 0`. Return a struct/named tuple of
   booleans.

---

## Lecture 10 — HMMs for gene finding (part 2): building the model

**Topics:** building the gene model incrementally:

1. 3-state skeleton `N` (intergenic) → `E` (exon) ↔ `I` (intron), with self-loops
   and biologically-motivated _missing_ edges (no `N→I`, no `I→N`); start in `N`.
2. Split the exon state into **three codon-position states** `E0→E1→E2→E0`
   (different emission profiles per codon position; the 3rd "wobble" base).
3. Enforce **motifs** via deterministic-emission on-ramp/off-ramp states:
   start codon `ATG`, the three stop codons (a compact 5-node sub-graph), donor
   `GT`, acceptor `AG`.
4. Keep codons **in phase** across introns by **cloning the intron module 3×**
   (exit at position `k` ⇒ re-enter at `k+1`); this forces total exon length to
   be a multiple of 3.
5. Counting parameters: only non-deterministic emissions (N, E0, E1, E2, I) and
   only nodes with out-degree > 1 need training.
6. Reality check: real gene finding combines models like this with RNA
   sequencing + spliced alignment.

### Math problems

1. Draw (on paper) the 3-state `N/E/I` HMM. List every _allowed_ transition and
   explain why `N→I` and `I→N` are disallowed.
2. Explain why splitting `E` into `E0,E1,E2` as a directed cycle (no self-loops on
   the individual nodes) is what models "consecutive codons", and why entry must
   be at `E0` and exit must be after `E2`.
3. **Phase / multiple-of-three.** Argue that cloning the intron module three
   times (exit `E0`→re-enter `E1`, exit `E1`→re-enter `E2`, exit `E2`→re-enter
   `E0`) guarantees the concatenated exon length is a multiple of 3.
4. **Parameter count.** For the final model, count (a) the states with
   non-trivial emission tables, and (b) the states with out-degree > 1 (the only
   transitions you must train). Why are the `ATG`/`GT`/`AG` ramp states "free"?

### Julia problems

Build these as a small library so they can be decoded with your **L7 `viterbi`**.

1. **Skeleton model.** Define states `[:N, :E, :I]` and build the 3-state HMM with
   emissions over `{A,C,G,T}`. Encode the disallowed transitions as `0.0` in `A`
   and `I=[1,0,0]` (always start in `N`). Add a test asserting `A[:N,:I] == 0`
   and `A[:I,:N] == 0`.
2. **Codon-position model.** Replace `:E` with `[:E0,:E1,:E2]` wired as a cycle
   `E0→E1→E2→E0`, with `N→E0` and `E2→N`, and intron transitions allowed from/to
   each `Ek`. Provide a builder `build_codon_hmm(params)` returning an `HMM`.
3. **Motif ramp states.** Add deterministic-emission states for the start codon
   (`A→T→G`), a compact stop-codon sub-graph (the 5-node version that admits
   exactly `TAA`, `TAG`, `TGA`), and donor (`G→T`) / acceptor (`A→G`) ramps.
   Implement a helper **`deterministic_emission_row(base)`** returning a length-4
   row that is `1.0` at `base` and `0.0` elsewhere. Test that a path emitting
   anything other than the required base gets joint probability `0`.
4. **In-phase intron cloning.** Implement the three intron copies so that exiting
   at codon position `k` forces re-entry at `(k+1) mod 3`. Write a test: feed an
   emission string representing `exon(in-phase) — intron — exon` and assert the
   Viterbi path keeps codon phase (no `E0`-exit / `E2`-reentry that would skip a
   position).
5. **End-to-end gene parse.** Construct a synthetic genomic string:
   `intergenic … ATG [codons] GT [intron] AG [codons] STOP … intergenic`,
   train/seed reasonable emission profiles, run `viterbi`, and verify the decoded
   path: (a) starts and ends in `N`, (b) marks the `ATG` and stop codon as exonic,
   (c) marks the `GT…AG` intron as intronic, (d) total exonic length is a
   multiple of 3.
6. **(★ stretch) pooled intron training.** Show that pooling the training data
   across the three cloned intron modules (sharing one emission table) gives the
   same emission estimates as training a single intron state — i.e. cloning was
   only needed for the _transition_ structure (phase), not the emissions.

---

## Capstone — putting it all together

A single project that reuses everything above:

1. **Data.** Either (a) download a small annotated region (e.g. the HBB gene
   neighborhood) and parse FASTA + an exon/intron annotation, or (b) generate
   synthetic labeled data with your L6 `simulate` and the L10 gene model.
2. **Markov-chain CpG classifier (L4–L5).** Train inside/outside chains, build the
   log-ratio table, and evaluate the `S(x) > 0` classifier on held-out data;
   report accuracy and plot the two score histograms.
3. **CpG HMM (L8).** Train the 8-state dinucleotide-aware model; decode a test
   sequence with Viterbi; compare island calls to the Markov-chain classifier.
4. **Gene-finding HMM (L10).** Decode a synthetic (or real) region; report the
   predicted exon/intron/intergenic segmentation; check the gene-structure
   invariants (start/stop codons, donor/acceptor motifs, multiple-of-3 exon
   length).
5. **Write-up.** For each model, state where the Markov assumption helped and
   where it is a known simplification (e.g. introns' real motif diversity,
   higher-order context, phase memory).

### Suggested `test/runtests.jl` skeleton

```julia
using Test
using SeqModel

@testset "L2 probability" begin
    space = dice_space()
    @test prob(isodd_d1, space) ≈ 0.5
    @test joint(isodd_d1, iseven_d2, space) ≈ 0.25
    @test conditional(isodd_d1, iseven_d2, space) ≈ 0.5
end

@testset "L4 Markov chain" begin
    A = [0.18 0.27 0.43 0.12;
         0.17 0.37 0.27 0.19;
         0.16 0.34 0.38 0.12;
         0.08 0.36 0.38 0.18]
    @test all(≈(1.0), sum(A, dims=2))
    # P("GATC") = 0.25 * A[G,A] * A[A,T] * A[T,C]
    @test seq_prob("GATC", A) ≈ 0.25 * 0.16 * 0.12 * 0.36
end

@testset "L7 Viterbi (sticky casino)" begin
    hmm = sticky_casino()
    path, _ = viterbi(hmm, "HH")
    @test path == [:Loaded, :Loaded]
end
```

---

### Appendix: minimal genetic-code table for Lecture 9

You only need enough entries for the tests; here is a usable starter (extend as
needed):

```julia
const GENETIC_CODE = Dict(
    "ATG"=>'M',                              # start / Met
    "AAA"=>'K', "AAG"=>'K',                  # Lys
    "CCA"=>'P', "CCC"=>'P', "CCG"=>'P', "CCT"=>'P',  # Pro
    "GGT"=>'G', "GGC"=>'G', "GGA"=>'G', "GGG"=>'G',  # Gly
    "TAA"=>'*', "TAG"=>'*', "TGA"=>'*',      # stop
    # … add the rest of the 64 codons as you need them
)
```

> AI:
> If you'd like, I can scaffold the actual SeqModel Julia package (Project.toml, src/, test/) with function stubs and failing tests so you can start red-green-refactor immediately.
