# =============================================================================
# Lectures 1-3 : DNA basics, naive probability, the Markov assumption
# =============================================================================

"""DNA alphabet, fixed order A,C,G,T (this order indexes all matrices)."""
const DNA = ('A', 'C', 'G', 'T')

"""Map a base character to its 1-based index in `DNA` (A=>1, C=>2, G=>3, T=>4)."""
const BASE_INDEX = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4)

"""
    random_dna(L) -> String

Return a uniformly random DNA string of length `L`.
"""
random_dna(L::Integer) = String(rand(DNA, L))

# ---- Lecture 1 -------------------------------------------------------------

"""
    gc_content(seq) -> Float64

Fraction of characters in `seq` that are `C` or `G`.
`gc_content("GATTACA") == 1/7`, `gc_content("GCGCGC") == 1.0`.
"""
function gc_content(seq::AbstractString)
    @stub "gc_content"
end

"""
    count_dinucleotide(seq, dinuc) -> Int

Number of **overlapping** occurrences of the 2-character pattern `dinuc` in `seq`.
`count_dinucleotide("CGCGCG", "CG") == 3`.
"""
function count_dinucleotide(seq::AbstractString, dinuc::AbstractString)
    @stub "count_dinucleotide"
end

"""
    cpg_count(seq) -> Int

Number of `CG` dinucleotides in `seq` (overlapping). Convenience wrapper.
"""
function cpg_count(seq::AbstractString)
    @stub "cpg_count"
end

# ---- Lecture 2 : naive probability over a finite, equally-likely space ------

"""
    dice_space() -> Vector{Tuple{Int,Int}}

All 36 `(d1, d2)` outcomes of rolling two fair dice.
"""
function dice_space()
    @stub "dice_space"
end

"`isodd_d1((d1,d2))` — event A: first die is odd."
isodd_d1((d1, d2)) = isodd(d1)
"`iseven_d2((d1,d2))` — event B: second die is even."
iseven_d2((d1, d2)) = iseven(d2)

"""
    prob(event, space) -> Float64

Naive probability: fraction of equally-likely `space` outcomes for which the
predicate `event(outcome)::Bool` holds.
"""
function prob(event, space)
    @stub "prob"
end

"""
    joint(a, b, space) -> Float64

`P(a, b)` — fraction of outcomes where both predicates hold.
"""
function joint(a, b, space)
    @stub "joint"
end

"""
    conditional(a, b, space) -> Float64

`P(a | b) = P(a, b) / P(b)`.
"""
function conditional(a, b, space)
    @stub "conditional"
end

# ---- Lecture 3 : why raw counting fails; Markov factorization ---------------

"""
    naive_string_prob(query, examples) -> Float64

Estimate `P(query)` as `count(query in examples) / length(examples)` over a
vector of (equal-length) strings. Returns `0.0` when `query` is unseen — the
failure mode that motivates the Markov assumption.
"""
function naive_string_prob(query::AbstractString, examples::AbstractVector{<:AbstractString})
    @stub "naive_string_prob"
end

"""
    markov_factorization_terms(seq) -> Vector{Tuple{Char,Union{Char,Nothing}}}

Terms produced by the first-order Markov factorization of `P(seq)`, as
`(current, given)` pairs (`given === nothing` marks the initial marginal).
`markov_factorization_terms("GATC")` →
`[('G', nothing), ('A','G'), ('T','A'), ('C','T')]`.
"""
function markov_factorization_terms(seq::AbstractString)
    @stub "markov_factorization_terms"
end
