### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> title = "03: Hidden Markov Models & Viterbi Decoding"
#> description = "Lectures 6-7"
#> chapter = 3
#> section = 3.0

using Markdown
using InteractiveUtils

# ╔═╡ 2de1eb4a-8ad3-4597-9896-a3362829e118
md"""
# **Notebook 03: Hidden Markov Models & Viterbi Decoding**
`Bioinformatics Sequence Modelling — Lectures 6 to 7`

This notebook introduces discrete Hidden Markov Models (`HMM`), joint path-emission scoring, simulation, and exact log-space Viterbi trellis decoding.
"""

# ╔═╡ b3541a27-0169-44f9-ac09-e82be527837d
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

# ╔═╡ 3a826430-1c9f-4190-bc86-2f43f07fa394
md"""
## **Lecture 6** — *Hidden Markov Model Definition & Casino Example*

An `HMM` contains:
- `states::Vector{Symbol}` ($Q$)
- `symbols::Vector{Char}` ($\Sigma$)
- `A::Matrix{Float64}` ($|Q| \times |Q|$ transition matrix)
- `E::Matrix{Float64}` ($|Q| \times |\Sigma|$ emission matrix)
- `I::Vector{Float64}` ($|Q|$ initial state distribution)
"""

# ╔═╡ 90b9ac9d-1852-4c5e-8b18-758f3d797dd3
struct HMM
    states::Vector{Symbol}
    symbols::Vector{Char}
    A::Matrix{Float64}
    E::Matrix{Float64}
    I::Vector{Float64}
end

state_index(hmm::HMM, s::Symbol) = findfirst(==(s), hmm.states)
symbol_index(hmm::HMM, c::Char) = findfirst(==(c), hmm.symbols)

# ╔═╡ 3684e23c-d406-46f4-a710-578da2f0800e
md"""
### **Exercise 6.1** — `casino()`
👉 Construct the dishonest casino HMM:
- States: `[:Fair, :Loaded]`
- Symbols: `['H', 'T']`
- Switch prob: 0.4 (stay 0.6)
- Fair emissions: $P(H)=P(T)=0.5$
- Loaded emissions: $P(H)=0.8, P(T)=0.2$
- Initial: $0.5 / 0.5$
"""

# ╔═╡ 8aafc432-c7e1-40c3-b1a5-10e5cec8d47f
function casino()
	@stub "casino"
end

# ╔═╡ 9c7f7a80-1e17-4617-84dd-e32f4a98bd2f
let
	res = try casino() catch e; e end
	if res isa Exception || res isa Missing
		still_missing(md"Implement `casino()`.")
	elseif res isa HMM && res.states == [:Fair, :Loaded] && all(≈(1.0), sum(res.A, dims=2)) && all(≈(1.0), sum(res.E, dims=2)) && sum(res.I) ≈ 1.0
		correct()
	else
		keep_working(md"Verify casino parameters: states [:Fair, :Loaded], row sums of A and E equal 1.0, sum(I) == 1.0.")
	end
end

# ╔═╡ 98e5f56e-f5db-4b97-8334-05e460edebe2
md"""
### **Exercise 6.2** — `joint_prob` & `joint_logprob`
👉 Joint probability of path $p$ and emissions $x$:
$$P(p, x) = I[p_1] E[p_1, x_1] \prod_{i=2}^n A[p_{i-1}, p_i] E[p_i, x_i]$$
"""

# ╔═╡ 2cdb3928-8f7c-4d47-9595-90755c1fefa2
function joint_prob(hmm::HMM, path::AbstractVector{Symbol}, emissions::AbstractString)
	@stub "joint_prob"
end

function joint_logprob(hmm::HMM, path::AbstractVector{Symbol}, emissions::AbstractString)
	@stub "joint_logprob"
end

# ╔═╡ ba28654f-bffc-431c-a6e9-d98aeff3a555
let
	c = try casino() catch ex; ex end
	if c isa Exception
		still_missing(md"Implement `casino()` first.")
	else
		p = [:Fair, :Fair, :Loaded]
		e = "HHT"
		res = try joint_prob(c, p, e) catch ex; ex end
		res_log = try joint_logprob(c, p, e) catch ex; ex end
		if res isa Exception || res_log isa Exception || res isa Missing
			still_missing(md"Implement `joint_prob` and `joint_logprob`.")
		elseif res ≈ 0.006 && res_log ≈ log(0.006)
			correct()
		else
			keep_working(md"`joint_prob(casino(), [:Fair,:Fair,:Loaded], 'HHT')` should be 0.006 and logprob ≈ log(0.006).")
		end
	end
end

# ╔═╡ 4d664ab3-9d8e-4545-b66b-2f3d0a69dd89
md"""
### **Exercise 6.3** — `simulate(hmm, n; rng)`
👉 Sample a length-$n$ state path and emission string from `hmm`.
"""

# ╔═╡ 5a0655a9-19c5-4a6e-ba00-eaee5d24316b
function simulate(hmm::HMM, n::Integer; rng = Random.default_rng())
	@stub "simulate"
end

# ╔═╡ efa1fadc-7fdf-4120-b63a-ae015ef26807
let
	c = try casino() catch ex; ex end
	if c isa Exception
		still_missing(md"Implement `casino()` first.")
	else
		res = try simulate(c, 10) catch e; e end
		if res isa Exception || res isa Missing
			still_missing(md"Implement `simulate(hmm, n)`.")
		elseif res isa Tuple && length(res[1]) == 10 && length(res[2]) == 10
			correct()
		else
			keep_working(md"`simulate` should return `(path, emissions)` of length n.")
		end
	end
end

# ╔═╡ 276612ae-c464-48f9-b283-c0319f2adf0e
md"""
## **Lecture 7** — *Viterbi Trellis Decoding*

### **Exercise 7.1** — `sticky_casino()`
👉 Same as `casino()` but stickier coins: stay probability $0.9$, switch probability $0.1$.
"""

# ╔═╡ 3378d4e8-7d98-44f3-b656-76908ef21de8
function sticky_casino()
	@stub "sticky_casino"
end

# ╔═╡ 70c526c4-0e2b-4f76-a3d4-d3016ae7de16
let
	res = try sticky_casino() catch e; e end
	if res isa Exception || res isa Missing
		still_missing(md"Implement `sticky_casino()`.")
	elseif res isa HMM && res.A[1,1] == 0.9 && res.A[1,2] == 0.1
		correct()
	else
		keep_working(md"`sticky_casino()` stay prob should be 0.9.")
	end
end

# ╔═╡ 5fe7de74-9ccc-46bb-ac1e-e4efcf18061d
md"""
### **Exercise 7.2** — `viterbi(hmm, emissions)`
👉 Find the most likely state path via Viterbi DP in **log space**. Return `(path, logprob)`.
"""

# ╔═╡ 92910a9b-aa6d-45b2-bd16-0b328f027cee
function viterbi(hmm::HMM, emissions::AbstractString)
	@stub "viterbi"
end

# ╔═╡ 7dc90a27-0e8c-4fcc-881b-98b6ad6ebec2
let
	sc = try sticky_casino() catch e; e end
	if sc isa Exception
		still_missing(md"Implement `sticky_casino()` first.")
	else
		res = try viterbi(sc, "HH") catch e; e end
		if res isa Exception || res isa Missing
			still_missing(md"Implement `viterbi(hmm, emissions)`.")
		elseif res isa Tuple && res[1] == [:Loaded, :Loaded]
			correct()
		else
			keep_working(md"`viterbi(sticky_casino(), 'HH')` should decode to `[:Loaded, :Loaded]`.")
		end
	end
end

# ╔═╡ 5c60f829-032d-438e-945d-82c4956aadb9
md"""
### **Exercise 7.3** — `viterbi_matrix(hmm, emissions)`
👉 Return full Viterbi DP table `V` ($|Q| \times n$) and backpointer matrix `back`.
"""

# ╔═╡ e77e40b4-056d-4411-8c40-09c3ec696efb
function viterbi_matrix(hmm::HMM, emissions::AbstractString)
	@stub "viterbi_matrix"
end

# ╔═╡ f18f8c23-b0a7-4e5b-b01b-062f3ad45601
let
	sc = try sticky_casino() catch e; e end
	if sc isa Exception
		still_missing(md"Implement `sticky_casino()` first.")
	else
		res = try viterbi_matrix(sc, "HTTHH") catch e; e end
		if res isa Exception || res isa Missing
			still_missing(md"Implement `viterbi_matrix`.")
		elseif res isa Tuple && size(res[1]) == (2, 5) && size(res[2]) == (2, 5)
			correct()
		else
			keep_working(md"`viterbi_matrix` should return (V, back) matrices of size (2, 5) for length 5 emissions.")
		end
	end
end

# ╔═╡ Cell order:
# ╟─2de1eb4a-8ad3-4597-9896-a3362829e118
# ╠═b3541a27-0169-44f9-ac09-e82be527837d
# ╟─3a826430-1c9f-4190-bc86-2f43f07fa394
# ╠═90b9ac9d-1852-4c5e-8b18-758f3d797dd3
# ╟─3684e23c-d406-46f4-a710-578da2f0800e
# ╠═8aafc432-c7e1-40c3-b1a5-10e5cec8d47f
# ╠═9c7f7a80-1e17-4617-84dd-e32f4a98bd2f
# ╟─98e5f56e-f5db-4b97-8334-05e460edebe2
# ╠═2cdb3928-8f7c-4d47-9595-90755c1fefa2
# ╠═ba28654f-bffc-431c-a6e9-d98aeff3a555
# ╟─4d664ab3-9d8e-4545-b66b-2f3d0a69dd89
# ╠═5a0655a9-19c5-4a6e-ba00-eaee5d24316b
# ╠═efa1fadc-7fdf-4120-b63a-ae015ef26807
# ╟─276612ae-c464-48f9-b283-c0319f2adf0e
# ╠═3378d4e8-7d98-44f3-b656-76908ef21de8
# ╠═70c526c4-0e2b-4f76-a3d4-d3016ae7de16
# ╟─5fe7de74-9ccc-46bb-ac1e-e4efcf18061d
# ╠═92910a9b-aa6d-45b2-bd16-0b328f027cee
# ╠═7dc90a27-0e8c-4fcc-881b-98b6ad6ebec2
# ╟─5c60f829-032d-438e-945d-82c4956aadb9
# ╠═e77e40b4-056d-4411-8c40-09c3ec696efb
# ╠═f18f8c23-b0a7-4e5b-b01b-062f3ad45601
