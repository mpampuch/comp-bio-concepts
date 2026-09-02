### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> title = "02: Markov Chains, Log-Space Scoring & CpG Classifier"
#> description = "Lectures 4-5"
#> chapter = 2
#> section = 2.0

using Markdown
using InteractiveUtils

# ╔═╡ bd490a62-8931-457a-bf7b-85a3ec2ecc21
md"""
# **Notebook 02: Markov Chains, Log-Space Scoring & CpG Classifier**
`Bioinformatics Sequence Modelling — Lectures 4 to 5`

This notebook covers building, training, and scoring Markov chains in linear and log space, building up to the CpG island log-likelihood-ratio classifier.
"""

# ╔═╡ e7d77d4d-bbf2-4888-bb15-39ce3920db7d
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

# ╔═╡ 1bb053e8-b44a-467a-b441-43398c3fcc8d
md"""
## **Lecture 4** — *Markov chains (part 1): building & scoring*

Transition matrix $A[prev, curr] = P(curr \mid prev)$. Rows sum to 1.
Indexing order: $A \to 1, C \to 2, G \to 3, T \to 4$ (`BASE_INDEX`).

Given trained transition table $A_{TABLE}$:
$$A_{TABLE} = \begin{bmatrix}
0.18 & 0.27 & 0.43 & 0.12 \\
0.17 & 0.37 & 0.27 & 0.19 \\
0.16 & 0.34 & 0.38 & 0.12 \\
0.08 & 0.36 & 0.38 & 0.18
\end{bmatrix}$$
"""

# ╔═╡ 35f44fa8-31cd-423d-9dae-2b2c11ff66b2
const A_TABLE = [0.18 0.27 0.43 0.12;
                 0.17 0.37 0.27 0.19;
                 0.16 0.34 0.38 0.12;
                 0.08 0.36 0.38 0.18]

# ╔═╡ 1ef95154-a5d4-4832-81cf-924cf9483c58
md"""
### **Exercise 4.1** — `transition_counts(seqs; alphabet=DNA)`
👉 Return a $4 \times 4$ integer count matrix of dinucleotide occurrences across sequences.
"""

# ╔═╡ b0eae793-8f87-4457-a7c0-97dd327e9012
function transition_counts(seqs; alphabet = DNA)
	@stub "transition_counts"
end

# ╔═╡ 1994a2c3-6feb-40a2-9c08-9adc530ab5e0
let
	res = try transition_counts(["ACGTACGT"]) catch e; e end
	if res isa Exception || res isa Missing
		still_missing(md"Implement `transition_counts(seqs)`.")
	elseif !(res isa AbstractMatrix) || size(res) != (4, 4)
		keep_working(md"`transition_counts` must return a 4x4 Matrix.")
	else
		correct()
	end
end

# ╔═╡ 6f9ab85e-3855-4a3f-9ffb-18676654cdbb
md"""
### **Exercise 4.2** — `transition_matrix(counts)`
👉 Row-normalize count matrix into transition probabilities $P(curr \mid prev)$. Every row must sum to $1.0$.
"""

# ╔═╡ 6b724724-395b-4dd8-b01c-01a075cf9d82
function transition_matrix(counts::AbstractMatrix)
	@stub "transition_matrix"
end

# ╔═╡ 8cb17070-cc17-4ab1-8222-9188493a2fb8
let
	cnts = try transition_counts(["ACGTACGT"]) catch e; e end
	if cnts isa Exception
		still_missing(md"Implement `transition_counts` first.")
	else
		res = try transition_matrix(cnts) catch e; e end
		if res isa Exception || res isa Missing
			still_missing(md"Implement `transition_matrix(counts)`.")
		elseif size(res) == (4,4) && all(≈(1.0), sum(res, dims=2))
			correct()
		else
			keep_working(md"Check that `transition_matrix` row-normalizes counts so each row sums to 1.0.")
		end
	end
end

# ╔═╡ b644595c-1808-4f6c-9332-4f9074660223
md"""
### **Exercise 4.3** — `seq_prob(seq, A; init=fill(0.25, 4))`
👉 Joint probability $P(seq) = init[x_1] \cdot \prod_{i=2}^n A[x_{i-1}, x_i]$.
"""

# ╔═╡ ca54b3e1-618d-4c07-adbd-e7d31d83416c
function seq_prob(seq::AbstractString, A::AbstractMatrix; init = fill(0.25, 4))
	@stub "seq_prob"
end

# ╔═╡ f4d6adef-a3dd-4814-84c9-44b8489899e5
let
	res = try seq_prob("GATC", A_TABLE) catch e; e end
	expected = 0.25 * 0.16 * 0.12 * 0.36
	if res isa Exception || res isa Missing
		still_missing(md"Implement `seq_prob(seq, A)`.")
	elseif res ≈ expected
		correct()
	else
		keep_working(md"`seq_prob('GATC', A_TABLE)` should be 0.25 * 0.16 * 0.12 * 0.36.")
	end
end

# ╔═╡ 9bfae159-9ae0-4ee8-9627-f8bb2b2955b7
md"""
## **Lecture 5** — *Log space & CpG log-likelihood-ratio classifier*

In floating point math, $0.5^{1100} = 0.0$ (IEEE 754 underflow). Log-space converts products into sums:
$$\log P(seq) = \log init[x_1] + \sum_{i=2}^n \log A[x_{i-1}, x_i]$$

The CpG classifier computes:
$$\text{Score}(seq) = \sum_{i=2}^n \log_2 \frac{A_{in}[x_{i-1}, x_i]}{A_{out}[x_{i-1}, x_i]}$$
"""

# ╔═╡ b7110275-2c3e-4a36-9119-46e53d7c0636
const A_IN = [0.20 0.30 0.30 0.20;
              0.10 0.10 0.70 0.10;
              0.20 0.30 0.30 0.20;
              0.20 0.30 0.30 0.20]
const A_OUT = [0.30 0.20 0.20 0.30;
               0.30 0.30 0.10 0.30;
               0.30 0.20 0.20 0.30;
               0.30 0.20 0.20 0.30]

# ╔═╡ 36e95d93-c816-4faf-98c4-a08326c1948d
md"""
### **Exercise 5.1** — `logspace_seq_logprob(seq, A; init=fill(0.25, 4))`
👉 Compute joint log probability by summing logs. Verify $\exp(\text{logprob}) \approx \text{seq\_prob}$.
"""

# ╔═╡ 69cb0e0f-3480-4b8c-9ea0-63994cc92e75
function logspace_seq_logprob(seq::AbstractString, A::AbstractMatrix; init = fill(0.25, 4))
	@stub "logspace_seq_logprob"
end

# ╔═╡ 2a2c8b61-5d44-47d1-a9ab-c95c2f8e4a67
let
	res = try logspace_seq_logprob("GATC", A_TABLE) catch e; e end
	sp = try seq_prob("GATC", A_TABLE) catch e; e end
	if res isa Exception || sp isa Exception || res isa Missing
		still_missing(md"Implement `logspace_seq_logprob`.")
	elseif exp(res) ≈ sp
		correct()
	else
		keep_working(md"`exp(logspace_seq_logprob(...))` should equal `seq_prob(...)`.")
	end
end

# ╔═╡ e32582be-ba26-4f83-b6c4-36f92b8d728c
md"""
### **Exercise 5.2** — `logratio_table(A_in, A_out)`
👉 Matrix of per-dinucleotide log-ratios $\log_2(A_{in}) - \log_2(A_{out})$.
"""

# ╔═╡ b3b7e6e6-159b-4136-8ef6-9bb916c4ee25
function logratio_table(A_in::AbstractMatrix, A_out::AbstractMatrix)
	@stub "logratio_table"
end

# ╔═╡ 44037e5a-f85b-4d50-a94a-a4b2289771c6
let
	res = try logratio_table(A_IN, A_OUT) catch e; e end
	if res isa Exception || res isa Missing
		still_missing(md"Implement `logratio_table`.")
	elseif size(res) == (4,4) && res[2,3] ≈ log2(0.70 / 0.10)
		correct()
	else
		keep_working(md"`logratio_table` should be 4x4 with entries log2(A_in) - log2(A_out).")
	end
end

# ╔═╡ f5e5bbd8-9a2e-422a-a794-9be98777e8d4
md"""
### **Exercise 5.3** — `cpg_score(seq, logratio)`
👉 Sum log-ratios along `seq`. Positive $\implies$ CpG island; negative $\implies$ background.
"""

# ╔═╡ 7856211c-4cac-4845-81cd-f154539958e6
function cpg_score(seq::AbstractString, logratio::AbstractMatrix)
	@stub "cpg_score"
end

# ╔═╡ 64778101-fca0-480a-be43-21ff03f84d9f
let
	lr = try logratio_table(A_IN, A_OUT) catch e; e end
	if lr isa Exception
		still_missing(md"Implement `logratio_table` first.")
	else
		s_cpg = try cpg_score("CGCGCGCGCG", lr) catch e; e end
		s_at = try cpg_score("ATATATATAT", lr) catch e; e end
		if s_cpg isa Exception || s_at isa Exception || s_cpg isa Missing
			still_missing(md"Implement `cpg_score`.")
		elseif s_cpg > 0 && s_at < 0
			correct()
		else
			keep_working(md"`cpg_score('CGCGCGCGCG', lr)` should be > 0 and `cpg_score('ATATATATAT', lr)` should be < 0.")
		end
	end
end

# ╔═╡ Cell order:
# ╟─bd490a62-8931-457a-bf7b-85a3ec2ecc21
# ╠═e7d77d4d-bbf2-4888-bb15-39ce3920db7d
# ╟─1bb053e8-b44a-467a-b441-43398c3fcc8d
# ╠═35f44fa8-31cd-423d-9dae-2b2c11ff66b2
# ╟─1ef95154-a5d4-4832-81cf-924cf9483c58
# ╠═b0eae793-8f87-4457-a7c0-97dd327e9012
# ╠═1994a2c3-6feb-40a2-9c08-9adc530ab5e0
# ╟─6f9ab85e-3855-4a3f-9ffb-18676654cdbb
# ╠═6b724724-395b-4dd8-b01c-01a075cf9d82
# ╠═8cb17070-cc17-4ab1-8222-9188493a2fb8
# ╟─b644595c-1808-4f6c-9332-4f9074660223
# ╠═ca54b3e1-618d-4c07-adbd-e7d31d83416c
# ╠═f4d6adef-a3dd-4814-84c9-44b8489899e5
# ╟─9bfae159-9ae0-4ee8-9627-f8bb2b2955b7
# ╠═b7110275-2c3e-4a36-9119-46e53d7c0636
# ╟─36e95d93-c816-4faf-98c4-a08326c1948d
# ╠═69cb0e0f-3480-4b8c-9ea0-63994cc92e75
# ╠═2a2c8b61-5d44-47d1-a9ab-c95c2f8e4a67
# ╟─e32582be-ba26-4f83-b6c4-36f92b8d728c
# ╠═b3b7e6e6-159b-4136-8ef6-9bb916c4ee25
# ╠═44037e5a-f85b-4d50-a94a-a4b2289771c6
# ╟─f5e5bbd8-9a2e-422a-a794-9be98777e8d4
# ╠═7856211c-4cac-4845-81cd-f154539958e6
# ╠═64778101-fca0-480a-be43-21ff03f84d9f
