# Bioinformatics Sequence Modelling — Practice Problems (Julia)

## Layout

```
SeqModel/
  Project.toml          # name + UUID, Random dep, Test as test target
  notebooks/            # Interactive Pluto.jl notebooks (Lectures 1-10)
    00_overview.jl
    01_dna_and_probability.jl
    02_markov_chains_and_cpg.jl
    03_hmms_and_viterbi.jl
    04_cpg_island_hmms.jl
    05_gene_finding_hmms.jl
  src/
    SeqModel.jl         # module: NotImplemented + @stub, includes, exports
    probability.jl      # L1–L3
    markov.jl           # L4–L5
    hmm.jl              # L6–L7  (defines the HMM struct)
    cpg_hmm.jl          # L8
    genes.jl            # L9–L10 (defines GENETIC_CODE)
  test/runtests.jl      # testsets for module functions & notebook execution
```

## How it works — Interactive Red-Green Refactoring in Pluto

- Every function starts as a documented **stub** throwing `NotImplemented` via `@stub "name"`.
- In the **Pluto.jl notebooks** ([`notebooks/`](file:///Users/markpampuch/Downloads/tmp/pHMMs/comp-bio-concepts/HMMs/notebooks)), each exercise includes:
  - **Concept Check & Math Problems**: Conceptual background and mathematical formulations.
  - **Function Stub Cell**: Starter code in the 🔴 **Red Phase**.
  - **Live Answer Check Cell**: Reactive test box rendering `Markdown.Admonition` boxes (`still_missing()` / `keep_working()`).
- As you write a working function implementation in the cell, the test box immediately turns into a 🟢 **Green Phase** success box (`correct()`).

## Run it

### Interactive Pluto Notebooks
Start Pluto in Julia and open any notebook in `notebooks/`:
```julia
using Pluto
Pluto.run()
```

### Command Line Test Suite
Run the test runner to verify both module stubs and notebook execution:
```bash
julia --project=. --startup-file=no test/runtests.jl
```


