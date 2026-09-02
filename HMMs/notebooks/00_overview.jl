### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> chapter = 0
#> section = 0.0
#> title = "00: HMMs Curriculum Overview"
#> description = "Interactive Red-Green Refactoring with Pluto.jl Notebooks"

using Markdown
using InteractiveUtils

# ╔═╡ 81654bca-22ad-498f-8813-8e16600597f2
md"""
# **Bioinformatics Sequence Modelling — Interactive Notebooks**

Welcome! This interactive course refactors bioinformatics sequence modelling into a hands-on **Red-Green Refactoring** learning experience using Pluto.jl notebooks, inspired by the MIT 18.S191 Computational Thinking notebooks.

### **How Red-Green Refactoring Works Here**
- 🔴 **Red Phase (Not Implemented)**: When you first open a notebook, every function starts as an unimplemented stub (`@stub "function_name"`). The live answer check cell directly below it displays a **yellow or red alert box** detailing what's missing and what expected behavior to implement.
- 🟢 **Green Phase (Correct)**: Edit the function code cell right above the check box. As soon as you write a correct implementation, Pluto reactively re-evaluates the check cell and turns it into a **green alert box** ("Got it! 🎉").

---

### **Course Notebooks Map**
1. **[01_dna_and_probability.jl](file:///Users/markpampuch/Downloads/tmp/pHMMs/comp-bio-concepts/HMMs/notebooks/01_dna_and_probability.jl)** — *Lectures 1–3*: DNA Alphabet, CpG Counts, Naive Probability & The Markov Assumption.
2. **[02_markov_chains_and_cpg.jl](file:///Users/markpampuch/Downloads/tmp/pHMMs/comp-bio-concepts/HMMs/notebooks/02_markov_chains_and_cpg.jl)** — *Lectures 4–5*: Markov Chains, Log-Space Scoring & The CpG Log-Likelihood-Ratio Classifier.
3. **[03_hmms_and_viterbi.jl](file:///Users/markpampuch/Downloads/tmp/pHMMs/comp-bio-concepts/HMMs/notebooks/03_hmms_and_viterbi.jl)** — *Lectures 6–7*: Hidden Markov Models, Dishonest Casino & Viterbi Trellis Decoding.
4. **[04_cpg_island_hmms.jl](file:///Users/markpampuch/Downloads/tmp/pHMMs/comp-bio-concepts/HMMs/notebooks/04_cpg_island_hmms.jl)** — *Lecture 8*: Naive 2-State HMM vs. Dinucleotide-Aware Fixed 8-State CpG Model.
5. **[05_gene_finding_hmms.jl](file:///Users/markpampuch/Downloads/tmp/pHMMs/comp-bio-concepts/HMMs/notebooks/05_gene_finding_hmms.jl)** — *Lectures 9–10*: Codons, ORFs, Splice Motifs & Multi-State Gene-Finding HMM Construction.

---

### **Running in Pluto**
To launch any notebook interactively:
```julia
using Pluto
Pluto.run()
```
Or run the full test suite in terminal:
```bash
julia --project=. test/runtests.jl
```
"""

# ╔═╡ Cell order:
# ╟─81654bca-22ad-498f-8813-8e16600597f2
