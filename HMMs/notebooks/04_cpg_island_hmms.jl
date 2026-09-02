### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> title = "04: CpG-Island HMMs (Naive vs Fixed 8-State)"
#> description = "Lecture 8"
#> chapter = 4
#> section = 4.0

using Markdown
using InteractiveUtils

# ╔═╡ 43b11222-9ce4-46d8-8dc1-ebdc92b97839
md"""
# **Notebook 04: CpG-Island HMMs**
`Bioinformatics Sequence Modelling — Lecture 8`

This notebook contrasts the naive 2-state CpG HMM with the dinucleotide-aware 8-state model that incorporates the previous nucleotide into the state space.
"""

# ╔═╡ d615d65b-047f-42ec-aa8a-ef3ab780cbd4
begin
	using Random

	struct NotImplemented <: Exception
		what::String
	end
	Base.showerror(io::IO, e::NotImplemented) = print(io, "NotImplemented: implement `", e.what, "`")
	macro stub(name)
		return :(throw(NotImplemented($(esc(name)))))
	end

	still_missing(text=md"Replace `@stub` with your implementation.") = 
		Markdown.MD(Markdown.Admonition("warning", "Not Implemented Yet (Red Phase)", [text]))

	keep_working(text=md"The answer is not quite right.") = 
		Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]))

	yays = [md"Fantastic!", md"Splendid!", md"Great! 🎉", md"Well done!", md"Keep it up!", md"Good job!", md"Awesome!", md"You got the right answer!"]
	correct(text=rand(yays)) = 
		Markdown.MD(Markdown.Admonition("correct", "Got it! (Green Phase)", [text]))

	const DNA = ('A', 'C', 'G', 'T')
	const BASE_INDEX = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4)
	random_dna(L::Integer) = String(rand(DNA, L))
end

# ╔═╡ 4aec8970-3ea6-479e-be7e-b4d8037ec4e3
struct HMM
    states::Vector{Symbol}
    symbols::Vector{Char}
    A::Matrix{Float64}
    E::Matrix{Float64}
    I::Vector{Float64}
end
state_index(hmm::HMM, s::Symbol) = findfirst(==(s), hmm.states)
symbol_index(hmm::HMM, c::Char) = findfirst(==(c), hmm.symbols)

function joint_logprob(hmm::HMM, path::AbstractVector{Symbol}, emissions::AbstractString)
    n = length(emissions)
    n == 0 && return 0.0
    p1 = state_index(hmm, path[1])
    x1 = symbol_index(hmm, emissions[1])
    lp = log(hmm.I[p1]) + log(hmm.E[p1, x1])
    for i in 2:n
        pi_prev = state_index(hmm, path[i-1])
        pi_curr = state_index(hmm, path[i])
        xi = symbol_index(hmm, emissions[i])
        lp += log(hmm.A[pi_prev, pi_curr]) + log(hmm.E[pi_curr, xi])
    end
    return lp
end

function viterbi(hmm::HMM, emissions::AbstractString)
    n = length(emissions)
    n == 0 && return (Symbol[], 0.0)
    n_states = length(hmm.states)
    V = fill(-Inf, n_states, n)
    back = zeros(Int, n_states, n)
    x1 = symbol_index(hmm, emissions[1])
    for s in 1:n_states
        if hmm.I[s] > 0 && hmm.E[s, x1] > 0
            V[s, 1] = log(hmm.I[s]) + log(hmm.E[s, x1])
        end
    end
    for t in 2:n
        xt = symbol_index(hmm, emissions[t])
        for s in 1:n_states
            e_prob = hmm.E[s, xt]
            e_prob <= 0 && continue
            log_e = log(e_prob)
            best_val = -Inf
            best_prev = 0
            for prev in 1:n_states
                V[prev, t-1] == -Inf && continue
                a_prob = hmm.A[prev, s]
                a_prob <= 0 && continue
                val = V[prev, t-1] + log(a_prob) + log_e
                if val > best_val
                    best_val = val
                    best_prev = prev
                end
            end
            V[s, t] = best_val
            back[s, t] = best_prev
        end
    end
    last_state = argmax(V[:, n])
    best_logprob = V[last_state, n]
    best_path = Vector{Symbol}(undef, n)
    curr = last_state
    for t in n:-1:1
        best_path[t] = hmm.states[curr]
        curr = back[curr, t]
    end
    return (best_path, best_logprob)
end

# ╔═╡ dcc58ee1-6024-4ef4-a6f2-911ba86f45f0
md"""
## **Part 1** — *Naive 2-State HMM* (`[:Inside, :Outside]`)

The naive 2-state model treats single nucleotide emissions independently, creating a dinucleotide blind spot (it misclassifies C-runs as CpG islands!).

- `train_io_transitions(labels)`: estimate $2 \times 2$ matrix $P(label \mid prev\_label)$ over `[:Inside, :Outside]`.
- `train_io_emissions(seq, labels)`: estimate $2 \times 4$ emission matrix $P(base \mid label)$ over `DNA`.
- `build_io_hmm(seq, labels)`: assemble HMM with uniform initial distribution.
"""

# ╔═╡ 6f1bfb4d-f9fd-446b-98d3-3850f802268a
function train_io_transitions(labels::AbstractString)
	@stub "train_io_transitions"
end

function train_io_emissions(seq::AbstractString, labels::AbstractString)
	@stub "train_io_emissions"
end

function build_io_hmm(seq::AbstractString, labels::AbstractString)
	@stub "build_io_hmm"
end

# ╔═╡ 1f872826-6517-4b70-9b72-2692e9e959e2
let
	res_t = try train_io_transitions("IIOOIO") catch e; e end
	res_e = try train_io_emissions("ACGTAC", "IIOOIO") catch e; e end
	res_h = try build_io_hmm("ACGTAC", "IIOOIO") catch e; e end
	if res_t isa Exception || res_e isa Exception || res_h isa Exception || res_t isa Missing
		still_missing(md"Implement naive 2-state functions: `train_io_transitions`, `train_io_emissions`, `build_io_hmm`.")
	elseif size(res_t) == (2,2) && all(≈(1.0), sum(res_t, dims=2)) && size(res_e) == (2,4) && res_h isa HMM
		correct()
	else
		keep_working(md"Check naive 2-state shapes: transition 2x2, emission 2x4, HMM initial distribution uniform.")
	end
end

# ╔═╡ 31bbfa04-d36d-4232-a376-d534e3151144
md"""
## **Part 2** — *Fixed 8-State Dinucleotide-Aware Model*

The 8 states are `[:A_I,:C_I,:G_I,:T_I, :A_O,:C_O,:G_O,:T_O]`.
Emissions are **deterministic**: state `b_X` emits base `b` with probability $1.0$.

- `fixed_states()` $\to$ 8 state symbols.
- `fixed_emissions()` $\to 8 \times 4$ deterministic matrix.
- `train_fixed_transitions(seq, labels)` $\to 8 \times 8$ transition matrix.
- `build_fixed_hmm(seq, labels)` $\to$ full 8-state HMM.
"""

# ╔═╡ 91f6904e-993d-42cf-a45f-23dd792b8f8c
function fixed_states()
	@stub "fixed_states"
end

function fixed_emissions()
	@stub "fixed_emissions"
end

function train_fixed_transitions(seq, labels)
	@stub "train_fixed_transitions"
end

function build_fixed_hmm(seq, labels)
	@stub "build_fixed_hmm"
end

# ╔═╡ cbb5cc8b-fedc-4385-a41e-faabe5b06de9
let
	st = try fixed_states() catch e; e end
	em = try fixed_emissions() catch e; e end
	tr = try train_fixed_transitions("ACGTACGT", "IIIIOOOO") catch e; e end
	hm = try build_fixed_hmm("ACGTACGT", "IIIIOOOO") catch e; e end
	if st isa Exception || em isa Exception || tr isa Exception || hm isa Exception || st isa Missing
		still_missing(md"Implement fixed 8-state functions.")
	elseif length(st) == 8 && size(em) == (8,4) && size(tr) == (8,8) && hm isa HMM
		# Integration test check: C-run test
		seq = "AAAACGCGCGCGAAAA"
		labels = "OOOOIIIIIIIIOOOO"
		trained_hmm = try build_fixed_hmm(seq, labels) catch e; e end
		if trained_hmm isa HMM
			path_try = try viterbi(trained_hmm, "AAAACCCCCCCCAAAA") catch e; e end
			if path_try isa Tuple && all(s -> endswith(String(s), "_O"), path_try[1])
				correct(md"Fantastic! The 8-state model correctly avoids the C-run false positive!")
			else
				keep_working(md"The fixed 8-state model should decode a pure C-run as outside island (`_O` states).")
			end
		else
			still_missing(md"Implement `build_fixed_hmm`.")
		end
	else
		keep_working(md"Check 8-state specifications: 8 states, 8x4 deterministic emissions, 8x8 transition matrix.")
	end
end

# ╔═╡ Cell order:
# ╟─43b11222-9ce4-46d8-8dc1-ebdc92b97839
# ╠═d615d65b-047f-42ec-aa8a-ef3ab780cbd4
# ╠═4aec8970-3ea6-479e-be7e-b4d8037ec4e3
# ╟─dcc58ee1-6024-4ef4-a6f2-911ba86f45f0
# ╠═6f1bfb4d-f9fd-446b-98d3-3850f802268a
# ╠═1f872826-6517-4b70-9b72-2692e9e959e2
# ╟─31bbfa04-d36d-4232-a376-d534e3151144
# ╠═91f6904e-993d-42cf-a45f-23dd792b8f8c
# ╠═cbb5cc8b-fedc-4385-a41e-faabe5b06de9
