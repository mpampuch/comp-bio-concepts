# =============================================================================
# Lectures 4-5 : Markov chains, log-space scoring, the CpG classifier
# =============================================================================

"""
    transition_counts(seqs; alphabet=DNA) -> Matrix{Int}

`4×4` matrix of dinucleotide counts over a collection of training sequences.
`C[prev, curr]` counts occurrences of `prev` immediately followed by `curr`
(rows/cols indexed by `BASE_INDEX`).
"""
function transition_counts(seqs; alphabet = DNA)
    @stub "transition_counts"
end

"""
    transition_matrix(counts) -> Matrix{Float64}

Row-normalize a count matrix into conditional probabilities `P(curr | prev)`.
Every row must sum to 1.
"""
function transition_matrix(counts::AbstractMatrix)
    @stub "transition_matrix"
end

"""
    seq_prob(seq, A; init=fill(0.25, 4)) -> Float64

Joint probability of `seq` under a first-order Markov chain:
`init[x_1] * prod(A[x_{i-1}, x_i])`. Indices via `BASE_INDEX`.
"""
function seq_prob(seq::AbstractString, A::AbstractMatrix; init = fill(0.25, 4))
    @stub "seq_prob"
end

# ---- Lecture 5 : log space + log-likelihood-ratio classifier ----------------

"""
    logspace_seq_logprob(seq, A; init=fill(0.25, 4)) -> Float64

Same as [`seq_prob`](@ref) but returns a **log** probability by summing logs
(avoids underflow). `exp(logspace_seq_logprob(...)) ≈ seq_prob(...)`.
"""
function logspace_seq_logprob(seq::AbstractString, A::AbstractMatrix; init = fill(0.25, 4))
    @stub "logspace_seq_logprob"
end

"""
    logratio_table(A_in, A_out) -> Matrix{Float64}

Combined per-dinucleotide log-ratio table `log2.(A_in) .- log2.(A_out)`.
"""
function logratio_table(A_in::AbstractMatrix, A_out::AbstractMatrix)
    @stub "logratio_table"
end

"""
    cpg_score(seq, logratio) -> Float64

Sum the per-dinucleotide log-ratios along `seq`. Positive ⇒ looks like it came
from inside a CpG island; negative ⇒ outside.
"""
function cpg_score(seq::AbstractString, logratio::AbstractMatrix)
    @stub "cpg_score"
end
