# Practice problems for bioinformatics lectures
## Layout

```
SeqModel/
  Project.toml          # name + UUID, Random dep, Test as test target
  src/
    SeqModel.jl         # module: NotImplemented + @stub, includes, exports
    probability.jl      # L1–L3
    markov.jl           # L4–L5
    hmm.jl              # L6–L7  (defines the HMM struct)
    cpg_hmm.jl          # L8
    genes.jl            # L9–L10 (defines GENETIC_CODE)
  test/runtests.jl      # one @testset per lecture, asserting intended behavior
```

## How it works

- Every function is a documented **stub** that throws `NotImplemented` via a small `@stub "name"` macro, so each test surfaces a clear "implement `X`" message.
- `test/runtests.jl` encodes the *intended* behavior (with the exact expected numbers from the problem set, e.g. `seq_prob("GATC", A_TABLE) ≈ 0.25*0.16*0.12*0.36`, casino joint `≈ 0.006`, sticky-casino Viterbi on `"HH"` → `[:Loaded, :Loaded]`). Stub calls are inlined inside `@test` so each `NotImplemented` is recorded as a discrete failing test rather than aborting a whole testset.

## Verified red phase

```
Test Summary:        | Pass  Error  Total
SeqModel             |    4     48     52
```

The 4 passing are the non-stub sanity facts (e.g. `random_dna(50)` length, the IEEE underflow demonstration `0.5^1100 == 0.0`, the `A_TABLE` row-sum). The 48 errors are exactly the unimplemented stubs.

## Run it

```bash
cd /Users/markpampuch/SeqModel
julia --project=. --startup-file=no test/runtests.jl
```

> Note: `julia --project=. -e 'using Pkg; Pkg.test()'` tripped on a Revise precompile directive coming from your global `startup.jl`; running the test file directly (above) sidesteps that. If you want `] test` to work, we can pin/skip Revise or adjust your startup file.

## Suggested workflow
Implement top-down — Lecture 1 first (`gc_content`, `count_dinucleotide`, `cpg_count`), re-run, watch L1 go green, then proceed. The later integration tests (L8 fixed-model "C-run is not an island", L10 skeleton-HMM invariants) will pass once their building-block functions plus `viterbi` are done, giving you a natural end-to-end checkpoint.

One design note: `build_skeleton_hmm`/`build_codon_hmm` take a `params` NamedTuple whose schema I documented in the test (e.g. `p_stay_N`, `p_E_to_I`, `emit_N`); adjust that schema to taste when you implement them — just keep the asserted invariants (`states == [:N,:E,:I]`, `I == [1,0,0]`, `A[N,I] == A[I,N] == 0`).

