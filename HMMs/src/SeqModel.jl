"""
    SeqModel

Scaffold for the bioinformatics sequence-modelling problem set
(`bioinformatics-sequence-modelling-problems-to-solve-in-Julia.md`).

Every function below is a **stub** that throws `NotImplemented`. The test suite
(`test/runtests.jl`) is written against the intended behaviour, so it currently
**fails on purpose** (the "red" of red-green-refactor). Implement one function,
re-run `] test`, watch it go green, then move on.

Source is split per lecture group:
  - `probability.jl`  : Lectures 1-3
  - `markov.jl`       : Lectures 4-5
  - `hmm.jl`          : Lectures 6-7 (defines the `HMM` type)
  - `cpg_hmm.jl`      : Lecture 8
  - `genes.jl`        : Lectures 9-10
"""
module SeqModel

using Random

"""Sentinel error thrown by every unimplemented stub."""
struct NotImplemented <: Exception
    what::String
end
Base.showerror(io::IO, e::NotImplemented) =
    print(io, "NotImplemented: implement `", e.what, "` (see problem set)")

"""`@stub "name"` — body for an unimplemented function; throws `NotImplemented`."""
macro stub(name)
    return :(throw(NotImplemented($(esc(name)))))
end

# Order matters: hmm.jl defines `HMM`, used by cpg_hmm.jl and genes.jl.
include("probability.jl")
include("markov.jl")
include("hmm.jl")
include("cpg_hmm.jl")
include("genes.jl")

# ---- exports --------------------------------------------------------------
# Lectures 1-3
export DNA, random_dna,
       gc_content, count_dinucleotide, cpg_count,
       dice_space, prob, joint, conditional, isodd_d1, iseven_d2,
       naive_string_prob, markov_factorization_terms
# Lectures 4-5
export transition_counts, transition_matrix, seq_prob,
       logspace_seq_logprob, logratio_table, cpg_score
# Lectures 6-7
export HMM, casino, sticky_casino,
       joint_prob, joint_logprob, simulate,
       viterbi, viterbi_matrix
# Lecture 8
export train_io_transitions, train_io_emissions, build_io_hmm,
       fixed_states, fixed_emissions, train_fixed_transitions, build_fixed_hmm
# Lectures 9-10
export GENETIC_CODE, codons, translate, find_orfs,
       is_donor, is_acceptor, gene_signal_report,
       deterministic_emission_row, build_skeleton_hmm, build_codon_hmm

end # module
