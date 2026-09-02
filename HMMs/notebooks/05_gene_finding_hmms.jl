### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> title = "05: Gene Signals & Building the Gene-Finding HMM"
#> description = "Lectures 9-10"
#> chapter = 5
#> section = 5.0

using Markdown
using InteractiveUtils

# ╔═╡ 6023c619-5e4f-42a8-b188-e036f9aaf480
md"""
# **Notebook 05: Gene Signals & Building the Gene-Finding HMM**
`Bioinformatics Sequence Modelling — Lectures 9 to 10`

This final notebook covers reading frame signal analysis (codons, ORFs, translation, splice motifs) and constructs full multi-state HMM architectures for gene prediction.
"""

# ╔═╡ f2929f38-831e-4c3f-8f2b-6101401bbb69
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

# ╔═╡ 382aee3b-7ea0-47f9-abd5-4a4501e7175f
struct HMM
    states::Vector{Symbol}
    symbols::Vector{Char}
    A::Matrix{Float64}
    E::Matrix{Float64}
    I::Vector{Float64}
end
state_index(hmm::HMM, s::Symbol) = findfirst(==(s), hmm.states)
symbol_index(hmm::HMM, c::Char) = findfirst(==(c), hmm.symbols)

# ╔═╡ 345dee65-d9b6-4080-9c17-5e60751ee700
md"""
## **Lecture 9** — *Codons, Translation, ORFs & Splice Motifs*
"""

# ╔═╡ d3bf3dfa-42d9-4571-a938-ec67d3e17b2a
const GENETIC_CODE = Dict(
    "ATG" => 'M', "AAA" => 'K', "AAG" => 'K',
    "CCA" => 'P', "CCC" => 'P', "CCG" => 'P', "CCT" => 'P',
    "GGT" => 'G', "GGC" => 'G', "GGA" => 'G', "GGG" => 'G',
    "TAA" => '*', "TAG" => '*', "TGA" => '*'
)
const STOP_CODONS = ("TAA", "TAG", "TGA")
const START_CODON = "ATG"

# ╔═╡ e8af9030-36ff-4cc1-9ac1-cbf92fe96ece
md"""
### **Exercise 9.1 & 9.2** — `codons(seq)` & `translate(seq)`
👉 `codons(seq)` splits into 3-bp chunks. `translate(seq)` converts codons to amino acid characters (`'*'` for stop).
"""

# ╔═╡ a734470d-e658-44e8-ba40-957b31a7d710
function codons(seq::AbstractString)
	@stub "codons"
end

function translate(seq::AbstractString, code = GENETIC_CODE)
	@stub "translate"
end

# ╔═╡ 7efd51b5-1802-4ea2-a40c-2c44cca329f3
let
	res_c = try codons("ATGAAATAA") catch e; e end
	res_t = try translate("ATGAAATAA") catch e; e end
	if res_c isa Exception || res_t isa Exception || res_c isa Missing
		still_missing(md"Implement `codons` and `translate`.")
	elseif res_c == ["ATG", "AAA", "TAA"] && res_t == "MK*"
		correct()
	else
		keep_working(md"`codons('ATGAAATAA')` should be `['ATG', 'AAA', 'TAA']` and `translate` should return `'MK*'`. ")
	end
end

# ╔═╡ a04d1f6d-1dcd-431f-857b-370ee8775fb9
md"""
### **Exercise 9.3** — `is_donor(seq, i)` & `is_acceptor(seq, i)`
👉 Check for 5' donor `GT` motif and 3' acceptor `AG` motif.
"""

# ╔═╡ 44ba3953-877d-44fe-b46f-81fe179fd723
function is_donor(seq::AbstractString, i::Integer)
	@stub "is_donor"
end

function is_acceptor(seq::AbstractString, i::Integer)
	@stub "is_acceptor"
end

# ╔═╡ 5ac549f6-1076-4d29-8194-45444661ea2d
let
	d = try is_donor("AGTC", 2) catch e; e end
	a = try is_acceptor("CAGT", 2) catch e; e end
	if d isa Exception || a isa Exception || d isa Missing
		still_missing(md"Implement `is_donor` and `is_acceptor`.")
	elseif d == true && a == true && !is_donor("AAAA", 1)
		correct()
	else
		keep_working(md"`is_donor('AGTC', 2)` (checking pos 2:3 == 'GT') should be true.")
	end
end

# ╔═╡ e4ddf188-d916-400f-ab34-b405d9eead2d
md"""
### **Exercise 9.4 & 9.5** — `find_orfs(seq)` & `gene_signal_report(seq)`
👉 `find_orfs` returns `(start, stop)` 1-based index pairs for frame-1 `ATG` to first stop codon.
👉 `gene_signal_report` checks `(; starts_atg, ends_stop, multiple_of_three)`.
"""

# ╔═╡ 8578da38-da35-4b59-a26d-f1eb0131120d
function find_orfs(seq::AbstractString)
	@stub "find_orfs"
end

function gene_signal_report(seq::AbstractString)
	@stub "gene_signal_report"
end

# ╔═╡ 03076ec9-224f-4eba-a07d-af92139a898e
let
	orfs = try find_orfs("ATGAAATAA") catch e; e end
	rep = try gene_signal_report("ATGAAATAA") catch e; e end
	if orfs isa Exception || rep isa Exception || orfs isa Missing
		still_missing(md"Implement `find_orfs` and `gene_signal_report`.")
	elseif (1, 9) in orfs && rep.starts_atg && rep.ends_stop && rep.multiple_of_three
		correct()
	else
		keep_working(md"`find_orfs('ATGAAATAA')` should include (1, 9) and `gene_signal_report` should return true for all fields.")
	end
end

# ╔═╡ ac505cb8-da96-45cd-ba12-e1abe4d78021
md"""
## **Lecture 10** — *Gene-Finding HMM Construction*

- `deterministic_emission_row(base)`: length-4 vector with $1.0$ at `BASE_INDEX[base]` and $0.0$ elsewhere.
- `build_skeleton_hmm(params)`: 3-state skeleton `[:N, :E, :I]` (intergenic, exon, intron) with disallowed transitions $N \to I, I \to N$ set to $0.0$ and $I=[1,0,0]$.
- `build_codon_hmm(params)`: elaborated gene HMM with 3 codon phase states `[:E0, :E1, :E2]`.
"""

# ╔═╡ bf53de9c-3503-4f3b-bd2a-afbeb047d06a
function deterministic_emission_row(base::Char)
	@stub "deterministic_emission_row"
end

function build_skeleton_hmm(params)
	@stub "build_skeleton_hmm"
end

function build_codon_hmm(params)
	@stub "build_codon_hmm"
end

# ╔═╡ f679be14-101e-488c-bdaf-cab498d7cd72
let
	row_g = try deterministic_emission_row('G') catch e; e end
	params = (p_stay_N = 0.99, p_enter_gene = 0.01,
	          p_stay_E = 0.90, p_E_to_I = 0.05, p_E_to_N = 0.05,
	          p_stay_I = 0.90, p_I_to_E = 0.10,
	          emit_N = fill(0.25, 4), emit_E = fill(0.25, 4), emit_I = fill(0.25, 4))
	skel = try build_skeleton_hmm(params) catch e; e end
	cod = try build_codon_hmm(params) catch e; e end

	if row_g isa Exception || skel isa Exception || cod isa Exception || row_g isa Missing
		still_missing(md"Implement `deterministic_emission_row`, `build_skeleton_hmm`, and `build_codon_hmm`.")
	elseif row_g == [0.0, 0.0, 1.0, 0.0] && skel isa HMM && cod isa HMM
		ni = findfirst(==(:N), skel.states); ii = findfirst(==(:I), skel.states)
		if skel.states == [:N, :E, :I] && skel.I == [1.0, 0.0, 0.0] && skel.A[ni, ii] == 0.0 && skel.A[ii, ni] == 0.0 && :E0 in cod.states
			correct()
		else
			keep_working(md"Check skeleton invariants: states [:N,:E,:I], start intergenic I=[1,0,0], A[N,I] == A[I,N] == 0.0.")
		end
	else
		keep_working(md"Check deterministic emission row and HMM builds.")
	end
end

# ╔═╡ Cell order:
# ╟─6023c619-5e4f-42a8-b188-e036f9aaf480
# ╠═f2929f38-831e-4c3f-8f2b-6101401bbb69
# ╠═382aee3b-7ea0-47f9-abd5-4a4501e7175f
# ╟─345dee65-d9b6-4080-9c17-5e60751ee700
# ╠═d3bf3dfa-42d9-4571-a938-ec67d3e17b2a
# ╟─e8af9030-36ff-4cc1-9ac1-cbf92fe96ece
# ╠═a734470d-e658-44e8-ba40-957b31a7d710
# ╠═7efd51b5-1802-4ea2-a40c-2c44cca329f3
# ╟─a04d1f6d-1dcd-431f-857b-370ee8775fb9
# ╠═44ba3953-877d-44fe-b46f-81fe179fd723
# ╠═5ac549f6-1076-4d29-8194-45444661ea2d
# ╟─e4ddf188-d916-400f-ab34-b405d9eead2d
# ╠═8578da38-da35-4b59-a26d-f1eb0131120d
# ╠═03076ec9-224f-4eba-a07d-af92139a898e
# ╟─ac505cb8-da96-45cd-ba12-e1abe4d78021
# ╠═bf53de9c-3503-4f3b-bd2a-afbeb047d06a
# ╠═f679be14-101e-488c-bdaf-cab498d7cd72
