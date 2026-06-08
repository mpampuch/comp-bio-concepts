# =============================================================================
# Lecture 8 : HMMs for CpG-island finding (naive 2-state + fixed 8-state)
# =============================================================================

# ---- Naive {Inside, Outside} model (the one with the dinucleotide blind spot)

"""
    train_io_transitions(labels) -> Matrix{Float64}

`2×2` transition matrix `P(label | prev label)` estimated from a label string of
`'I'`/`'O'` characters. Row/col order: `[:Inside, :Outside]`.
"""
function train_io_transitions(labels::AbstractString)
    @stub "train_io_transitions"
end

"""
    train_io_emissions(seq, labels) -> Matrix{Float64}

`2×4` emission matrix `P(base | label)` from parallel `seq` / `labels`.
Row order `[:Inside, :Outside]`, column order `DNA`.
"""
function train_io_emissions(seq::AbstractString, labels::AbstractString)
    @stub "train_io_emissions"
end

"""
    build_io_hmm(seq, labels) -> HMM

Assemble the naive 2-state CpG HMM (`[:Inside, :Outside]`, emissions `DNA`)
from trained transitions/emissions; uniform initial distribution.
"""
function build_io_hmm(seq::AbstractString, labels::AbstractString)
    @stub "build_io_hmm"
end

# ---- Fixed model: previous nucleotide pushed into the state space -----------

"""
    fixed_states() -> Vector{Symbol}

The 8 states of the dinucleotide-aware model:
`[:A_I,:C_I,:G_I,:T_I, :A_O,:C_O,:G_O,:T_O]` (`_I` inside, `_O` outside island).
"""
function fixed_states()
    @stub "fixed_states"
end

"""
    fixed_emissions() -> Matrix{Float64}

`8×4` **deterministic** emission matrix: row for state `b_X` is `1.0` in the
column for base `b` and `0.0` elsewhere (emission is determined by the state).
"""
function fixed_emissions()
    @stub "fixed_emissions"
end

"""
    train_fixed_transitions(seq, labels) -> Matrix{Float64}

`8×8` transition matrix for the fixed model. Counts labelled dinucleotides:
`P(curr_base&label | prev_base&label)`. Rows/cols ordered as [`fixed_states`](@ref).
"""
function train_fixed_transitions(seq::AbstractString, labels::AbstractString)
    @stub "train_fixed_transitions"
end

"""
    build_fixed_hmm(seq, labels) -> HMM

Assemble the 8-state dinucleotide-aware CpG HMM (deterministic emissions +
trained transitions).
"""
function build_fixed_hmm(seq::AbstractString, labels::AbstractString)
    @stub "build_fixed_hmm"
end
