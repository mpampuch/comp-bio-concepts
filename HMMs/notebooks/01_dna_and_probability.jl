### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> chapter = 1
#> section = 1.0
#> title = "01: DNA Basics, Naive Probability & The Markov Assumption"
#> description = "Lectures 1-3"

using Markdown
using InteractiveUtils

# ╔═╡ c4876862-4cea-4ef0-9b94-3e7de74f672c
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

# ╔═╡ de6ade27-d38a-4dbc-9f1c-502835c301cb
md"""
# **Notebook 01: DNA Basics, Naive Probability & The Markov Assumption**
`Bioinformatics Sequence Modelling — Lectures 1 to 3`

This notebook features live answer checks. Function cells start as `@stub` (🔴 Red Phase). Edit the functions to turn the answer check boxes into green alerts (🟢 Green Phase)!
"""

# ╔═╡ c0b07b58-d822-496d-8e17-3385d32bb38d
md"""
## **Lecture 1** — *Intro: Why probabilistic sequence models?*

**Topics:** motivation (gene finding, ~20k genes), CpG dinucleotides and methylation, CpG islands as a genomic landmark, why scores should be *probabilities*.

### **Concept Check**
1. *Why is "count how many times I saw this exact string inside vs. outside a CpG island" a poor strategy as the string gets long?*
   - **Answer:** The number of possible strings of length $L$ is $4^L$, which grows exponentially. For $L=100$, $4^{100} \approx 1.6 \times 10^{60}$, far exceeding any genomic dataset size. Most strings will be observed $0$ times.
2. *The lecture argues scores should be probabilities. Give two concrete advantages.*
   - **Answer:** (a) Probabilities provide a universal, normalized scale $[0, 1]$ allowing comparison of different length sequences and models. (b) They combine naturally via standard probability laws (e.g. multiplication rule, Bayes rule).
3. *In a uniformly random DNA sequence, what is the probability that a given adjacent pair of bases is exactly `CG`? Spacing?*
   - **Answer:** $P(\text{C}) \times P(\text{G}) = \frac{1}{4} \times \frac{1}{4} = \frac{1}{16} = 0.0625$. Average spacing is $16$ bp.

### **Math Problems**
1. *Under uniform i.i.d., expected `CG` count in $L=1000$ bp?*
   - There are $L - 1 = 999$ adjacent positions. Expected count $= 999 \times \frac{1}{16} \approx 62.44$.
2. *Island spacing ~4 bp vs background ~16 bp. Density ratio?*
   - Island density $= 1/4 = 0.25$, background $= 1/16 = 0.0625$. Ratio $= 0.25 / 0.0625 = 4.0$.
"""

# ╔═╡ 610c5082-2882-42e2-8022-beeb56e1f33b
md"""
### **Exercise 1.1** — `gc_content(seq)`
👉 Return the fraction of characters in `seq` that are `'C'` or `'G'`.
"""

# ╔═╡ 7834b6a6-907e-4bbd-bffd-de8de1f2f986
function gc_content(seq::AbstractString)
	@show split(seq, dlm="G")
	# @stub "gc_content"
	println("T")
end

# ╔═╡ ad44f845-249a-4992-a29c-507c4529e36e
let
	res = try gc_content("GATTACA") catch e; e end
	if res isa NotImplemented || res isa Missing || res isa Exception
		still_missing(md"Implement `gc_content(seq)` to compute the fraction of 'C' and 'G' characters.")
	elseif !(res isa Real)
		keep_working(md"`gc_content` should return a number (`Float64`).")
	elseif gc_content("GATTACA") ≈ 1/7 && gc_content("GCGCGC") ≈ 1.0 && gc_content("ATAT") == 0.0
		correct()
	else
		keep_working(md"`gc_content('GATTACA')` should be ≈ 1/7 (0.1428) and `gc_content('GCGCGC')` should be 1.0.")
	end
end

# ╔═╡ ef215f8e-eb6a-473c-a3b0-3930ef581de1
md"""
### **Exercise 1.2** — `count_dinucleotide(seq, dinuc)`
👉 Count **overlapping** occurrences of 2-character pattern `dinuc` in `seq`.
"""

# ╔═╡ 4ebc1a1d-7ffa-49bc-863c-6670f5e17128
function count_dinucleotide(seq::AbstractString, dinuc::AbstractString)
	@stub "count_dinucleotide"
end

# ╔═╡ 8ec94807-9238-4261-8bda-d2fbd79d263c
let
	res = try count_dinucleotide("CGCGCG", "CG") catch e; e end
	if res isa NotImplemented || res isa Missing || res isa Exception
		still_missing(md"Implement `count_dinucleotide(seq, dinuc)` for overlapping 2-character patterns.")
	elseif !(res isa Integer)
		keep_working(md"`count_dinucleotide` should return an `Integer`.")
	elseif count_dinucleotide("CGCGCG", "CG") == 3 && count_dinucleotide("AAAA", "CG") == 0 && count_dinucleotide("CGCG", "CGC") == 0
		correct()
	else
		keep_working(md"`count_dinucleotide('CGCGCG', 'CG')` should be 3 and `count_dinucleotide('AAAA', 'CG')` should be 0.")
	end
end

# ╔═╡ c4445cbc-5053-4355-a990-c7cb6fc67e84
md"""
### **Exercise 1.3** — `cpg_count(seq)`
👉 Convenience wrapper counting `"CG"` dinucleotides.
"""

# ╔═╡ faf74d30-7c0f-472c-9ec8-a2e99de4d912
function cpg_count(seq::AbstractString)
	@stub "cpg_count"
end

# ╔═╡ 93e76b03-914f-455d-ac41-1c83fb0712dc
let
	res = try cpg_count("CGATCG") catch e; e end
	if res isa NotImplemented || res isa Missing || res isa Exception
		still_missing(md"Implement `cpg_count(seq)`.")
	elseif !(res isa Integer)
		keep_working(md"`cpg_count` should return an `Integer`.")
	elseif cpg_count("CGATCG") == 2 && cpg_count("GATTACA") == 0
		correct()
	else
		keep_working(md"`cpg_count('CGATCG')` should return 2.")
	end
end

# ╔═╡ 8009e2c4-7e88-4c8f-86f3-0f4298029302
md"""
## **Lecture 2** — *Probability review*

**Topics:** Sample space $\Omega$, events, $P(A)$, joint $P(A,B)$, conditional $P(A|B) = \frac{P(A,B)}{P(B)}$, multiplication rule, independence.

### **Math Problems (Two-Dice Running Example)**
Sample space $\Omega = \{(d_1, d_2) \mid d_1, d_2 \in \{1..6\}\}$ ($36$ equally-likely outcomes).
Let $A$ = "die 1 is odd", $B$ = "die 2 is even".
- $P(A) = 18/36 = 0.5$, $P(B) = 18/36 = 0.5$.
- $P(A,B) = 9/36 = 0.25$.
- $P(A|B) = \frac{P(A,B)}{P(B)} = \frac{0.25}{0.5} = 0.5$.
- Note $P(A,B) = P(A)P(B) = 0.25$, so $A$ and $B$ are independent.
"""

# ╔═╡ 2433c197-58cc-4acd-bf48-bf9c38752b53
isodd_d1((d1, d2)) = isodd(d1)
iseven_d2((d1, d2)) = iseven(d2)

# ╔═╡ 4410d99a-8888-4f89-8281-dca0bbcd27dc
md"""
### **Exercise 2.1** — `dice_space()`
👉 Return all 36 `(d1, d2)` tuples of rolling two fair dice.
"""

# ╔═╡ c72280db-ff0f-4e8e-9692-189e5a0edcdf
function dice_space()
	@stub "dice_space"
end

# ╔═╡ 86a5a341-9d33-4489-8939-ebf32551786e
let
	res = try dice_space() catch e; e end
	if res isa NotImplemented || res isa Missing || res isa Exception
		still_missing(md"Implement `dice_space()` to return a vector of all 36 `(d1, d2)` tuples.")
	elseif !(res isa AbstractVector) || length(res) != 36
		keep_working(md"`dice_space()` must return a vector of 36 elements.")
	elseif all(t -> t isa Tuple{Int,Int} && 1 <= t[1] <= 6 && 1 <= t[2] <= 6, res) && length(unique(res)) == 36
		correct()
	else
		keep_working(md"Check that all 36 unique pairs from (1,1) to (6,6) are generated.")
	end
end

# ╔═╡ b0ef6a89-429f-474c-a0d9-18fbc8af9a33
md"""
### **Exercise 2.2** — `prob(event, space)`
👉 Naive probability: fraction of `space` outcomes for which `event(outcome)` is `true`.
"""

# ╔═╡ 8e96e7aa-f96f-405e-81bc-dfdd8bf4164e
function prob(event, space)
	@stub "prob"
end

# ╔═╡ d4c3105f-2aac-4b1c-9888-c480d2787e7f
let
	ds = try dice_space() catch e; e end
	if ds isa Exception
		still_missing(md"Implement `dice_space()` first.")
	else
		res = try prob(isodd_d1, ds) catch e; e end
		if res isa NotImplemented || res isa Missing || res isa Exception
			still_missing(md"Implement `prob(event, space)`.")
		elseif !(res isa Real)
			keep_working(md"`prob` should return a Float64.")
		elseif prob(isodd_d1, ds) ≈ 0.5 && prob(iseven_d2, ds) ≈ 0.5
			correct()
		else
			keep_working(md"`prob(isodd_d1, dice_space())` should be 0.5.")
		end
	end
end

# ╔═╡ 7cdaa812-46ae-48c5-af99-2adef792a935
md"""
### **Exercise 2.3** — `joint(a, b, space)` & `conditional(a, b, space)`
👉 `joint(a,b,space)` computes $P(a,b)$ and `conditional(a,b,space)` computes $P(a|b) = P(a,b)/P(b)$.
"""

# ╔═╡ f7f82c38-6af5-4efa-b5af-7d592a7a40f2
function joint(a, b, space)
	@stub "joint"
end

function conditional(a, b, space)
	@stub "conditional"
end

# ╔═╡ 6cc43c06-4b36-4ab9-90ef-9756f959f4d3
let
	ds = try dice_space() catch e; e end
	if ds isa Exception
		still_missing(md"Implement `dice_space()` first.")
	else
		res_j = try joint(isodd_d1, iseven_d2, ds) catch e; e end
		res_c = try conditional(isodd_d1, iseven_d2, ds) catch e; e end
		if res_j isa Exception || res_c isa Exception || res_j isa Missing || res_c isa Missing
			still_missing(md"Implement both `joint` and `conditional`.")
		elseif res_j ≈ 0.25 && res_c ≈ 0.5
			correct()
		else
			keep_working(md"`joint` should return 0.25 and `conditional` should return 0.5 for the dice events.")
		end
	end
end

# ╔═╡ fc82378a-10e7-4f40-8c3e-210cd2d67f14
md"""
## **Lecture 3** — *The Markov assumption*

**Topics:** exponential data explosion of raw string counting; Markov assumption as a locality assumption ($x_i \perp \text{past} \mid x_{i-1}$); "the future is independent of the past given the present".
"""

# ╔═╡ 4ee4d8e2-ae06-448b-bf01-ce62c4f4113f
md"""
### **Exercise 3.1** — `naive_string_prob(query, examples)`
👉 Estimate $P(\text{query})$ as `count(query in examples) / length(examples)`. Returns `0.0` when unseen.
"""

# ╔═╡ 5801aba5-37c7-478d-92e7-a1c9db3fa619
function naive_string_prob(query::AbstractString, examples::AbstractVector{<:AbstractString})
	@stub "naive_string_prob"
end

# ╔═╡ 8c6943dc-d180-4f7b-b2f1-d29708ade3e3
let
	ex = ["ACGT", "ACGT", "TTTT", "GGGG"]
	res1 = try naive_string_prob("ACGT", ex) catch e; e end
	res2 = try naive_string_prob("CCCC", ex) catch e; e end
	if res1 isa Exception || res2 isa Exception || res1 isa Missing
		still_missing(md"Implement `naive_string_prob(query, examples)`.")
	elseif res1 ≈ 0.5 && res2 == 0.0
		correct()
	else
		keep_working(md"`naive_string_prob('ACGT', ex)` should be 0.5 and for unseen 'CCCC' it should be 0.0.")
	end
end

# ╔═╡ 2ad35a91-e77c-419d-a26c-e6179b90fcdd
md"""
### **Exercise 3.2** — `markov_factorization_terms(seq)`
👉 Return the list of `(current, given)` pairs for first-order Markov factorization (`given === nothing` for initial marginal).
E.g. `"GATC"` $\to$ `[('G', nothing), ('A','G'), ('T','A'), ('C','T')]`.
"""

# ╔═╡ 59a7f298-e939-4374-b130-06703c454e56
function markov_factorization_terms(seq::AbstractString)
	@stub "markov_factorization_terms"
end

# ╔═╡ 60e8ca1d-e44a-4948-b7c4-8e3d8e1784c9
let
	res = try markov_factorization_terms("GATC") catch e; e end
	if res isa Exception || res isa Missing
		still_missing(md"Implement `markov_factorization_terms(seq)`.")
	elseif res == [('G', nothing), ('A', 'G'), ('T', 'A'), ('C', 'T')]
		correct()
	else
		keep_working(md"`markov_factorization_terms('GATC')` should be `[('G', nothing), ('A','G'), ('T','A'), ('C','T')]`.")
	end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "9e64cf17f9522d20edabe6f2b4ec85252943fcae"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"
"""

# ╔═╡ Cell order:
# ╟─de6ade27-d38a-4dbc-9f1c-502835c301cb
# ╠═c4876862-4cea-4ef0-9b94-3e7de74f672c
# ╟─c0b07b58-d822-496d-8e17-3385d32bb38d
# ╟─610c5082-2882-42e2-8022-beeb56e1f33b
# ╠═7834b6a6-907e-4bbd-bffd-de8de1f2f986
# ╠═ad44f845-249a-4992-a29c-507c4529e36e
# ╟─ef215f8e-eb6a-473c-a3b0-3930ef581de1
# ╠═4ebc1a1d-7ffa-49bc-863c-6670f5e17128
# ╠═8ec94807-9238-4261-8bda-d2fbd79d263c
# ╟─c4445cbc-5053-4355-a990-c7cb6fc67e84
# ╠═faf74d30-7c0f-472c-9ec8-a2e99de4d912
# ╠═93e76b03-914f-455d-ac41-1c83fb0712dc
# ╟─8009e2c4-7e88-4c8f-86f3-0f4298029302
# ╠═2433c197-58cc-4acd-bf48-bf9c38752b53
# ╟─4410d99a-8888-4f89-8281-dca0bbcd27dc
# ╠═c72280db-ff0f-4e8e-9692-189e5a0edcdf
# ╠═86a5a341-9d33-4489-8939-ebf32551786e
# ╟─b0ef6a89-429f-474c-a0d9-18fbc8af9a33
# ╠═8e96e7aa-f96f-405e-81bc-dfdd8bf4164e
# ╠═d4c3105f-2aac-4b1c-9888-c480d2787e7f
# ╟─7cdaa812-46ae-48c5-af99-2adef792a935
# ╠═f7f82c38-6af5-4efa-b5af-7d592a7a40f2
# ╠═6cc43c06-4b36-4ab9-90ef-9756f959f4d3
# ╟─fc82378a-10e7-4f40-8c3e-210cd2d67f14
# ╟─4ee4d8e2-ae06-448b-bf01-ce62c4f4113f
# ╠═5801aba5-37c7-478d-92e7-a1c9db3fa619
# ╠═8c6943dc-d180-4f7b-b2f1-d29708ade3e3
# ╟─2ad35a91-e77c-419d-a26c-e6179b90fcdd
# ╠═59a7f298-e939-4374-b130-06703c454e56
# ╠═60e8ca1d-e44a-4948-b7c4-8e3d8e1784c9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
