# =============================================================================
# Lectures 6-7 : Hidden Markov Models + Viterbi decoding
# =============================================================================

"""
    HMM

A discrete Hidden Markov Model.

Fields:
  - `states::Vector{Symbol}`  — the state set `Q`.
  - `symbols::Vector{Char}`   — the emission alphabet `Σ`.
  - `A::Matrix{Float64}`      — `|Q|×|Q|`, `A[i,j] = P(state j | state i)`.
  - `E::Matrix{Float64}`      — `|Q|×|Σ|`, `E[i,k] = P(symbol k | state i)`.
  - `I::Vector{Float64}`      — `|Q|`, initial distribution over states.
"""
struct HMM
    states::Vector{Symbol}
    symbols::Vector{Char}
    A::Matrix{Float64}
    E::Matrix{Float64}
    I::Vector{Float64}
end

"1-based index of `s` in `hmm.states`."
state_index(hmm::HMM, s::Symbol) = findfirst(==(s), hmm.states)
"1-based index of emission `c` in `hmm.symbols`."
symbol_index(hmm::HMM, c::Char) = findfirst(==(c), hmm.symbols)

"""
    casino() -> HMM

The occasionally-dishonest casino from Lecture 6:
states `[:Fair, :Loaded]`, symbols `['H','T']`, switch prob `0.4` (stay `0.6`),
fair `P(H)=P(T)=0.5`, loaded `P(H)=0.8, P(T)=0.2`, initial `0.5/0.5`.
"""
function casino()
    @stub "casino"
end

"""
    sticky_casino() -> HMM

Same as [`casino`](@ref) but stickier coins (stay `0.9`, switch `0.1`) — the
variant used for the hand-worked Viterbi examples in Lecture 7.
"""
function sticky_casino()
    @stub "sticky_casino"
end

"""
    joint_prob(hmm, path, emissions) -> Float64

`I[p_1]·E[p_1,x_1]·prod(A[p_{i-1},p_i]·E[p_i,x_i])` for `path` (states) and
`emissions` (chars) of equal length.
"""
function joint_prob(hmm::HMM, path::AbstractVector{Symbol}, emissions::AbstractString)
    @stub "joint_prob"
end

"""
    joint_logprob(hmm, path, emissions) -> Float64

Log-space version of [`joint_prob`](@ref); `≈ log(joint_prob(...))`.
"""
function joint_logprob(hmm::HMM, path::AbstractVector{Symbol}, emissions::AbstractString)
    @stub "joint_logprob"
end

"""
    simulate(hmm, n; rng=Random.default_rng()) -> (path, emissions)

Sample a length-`n` state path and emission string from `hmm`
(`path::Vector{Symbol}`, `emissions::String`).
"""
function simulate(hmm::HMM, n::Integer; rng = Random.default_rng())
    @stub "simulate"
end

# ---- Lecture 7 : Viterbi ----------------------------------------------------

"""
    viterbi(hmm, emissions) -> (path::Vector{Symbol}, logprob::Float64)

Most likely state path for `emissions` via the Viterbi DP (implement in **log
space**). `logprob` is the joint log-probability of `(path, emissions)`.
"""
function viterbi(hmm::HMM, emissions::AbstractString)
    @stub "viterbi"
end

"""
    viterbi_matrix(hmm, emissions) -> (V::Matrix{Float64}, back::Matrix{Int})

Full Viterbi DP table `V` (log values, `|Q|×n`) and integer backpointer matrix
`back`, for inspection/printing of the trellis.
"""
function viterbi_matrix(hmm::HMM, emissions::AbstractString)
    @stub "viterbi_matrix"
end
