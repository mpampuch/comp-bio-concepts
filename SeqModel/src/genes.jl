# =============================================================================
# Lectures 9-10 : gene biology signals + gene-finding HMM construction
# =============================================================================

"""
Minimal genetic code (extend to all 64 codons as needed). `'*'` marks a stop.
"""
const GENETIC_CODE = Dict(
    "ATG" => 'M',                                    # start / Met
    "AAA" => 'K', "AAG" => 'K',                      # Lys
    "CCA" => 'P', "CCC" => 'P', "CCG" => 'P', "CCT" => 'P',  # Pro
    "GGT" => 'G', "GGC" => 'G', "GGA" => 'G', "GGG" => 'G',  # Gly
    "TAA" => '*', "TAG" => '*', "TGA" => '*',        # stop
)

"""The three stop codons (DNA alphabet)."""
const STOP_CODONS = ("TAA", "TAG", "TGA")
"""The start codon (DNA alphabet)."""
const START_CODON = "ATG"

# ---- Lecture 9 : codons / translation / ORFs / splice motifs ----------------

"""
    codons(seq) -> Vector{String}

Split `seq` into non-overlapping length-3 codons. `@assert length(seq) % 3 == 0`.
"""
function codons(seq::AbstractString)
    @stub "codons"
end

"""
    translate(seq, code=GENETIC_CODE) -> String

Translate `seq` codon-by-codon to amino acids (`'*'` for stop).
`translate("ATGAAATAA") == "MK*"`.
"""
function translate(seq::AbstractString, code = GENETIC_CODE)
    @stub "translate"
end

"""
    find_orfs(seq) -> Vector{Tuple{Int,Int}}

Scan reading frame 1 left-to-right; return `(start, stop)` 1-based index pairs,
each an `ATG` followed in-frame by the first stop codon (inclusive of the stop).
"""
function find_orfs(seq::AbstractString)
    @stub "find_orfs"
end

"""
    is_donor(seq, i) -> Bool

True iff `seq[i:i+1] == "GT"` (the intron 5' donor motif).
"""
function is_donor(seq::AbstractString, i::Integer)
    @stub "is_donor"
end

"""
    is_acceptor(seq, i) -> Bool

True iff `seq[i:i+1] == "AG"` (the intron 3' acceptor motif).
"""
function is_acceptor(seq::AbstractString, i::Integer)
    @stub "is_acceptor"
end

"""
    gene_signal_report(seq) -> NamedTuple

Check gene-structure hallmarks of a putative single-exon coding region:
`(; starts_atg, ends_stop, multiple_of_three)`.
(★ stretch: extend to multi-exon with donor/acceptor checks.)
"""
function gene_signal_report(seq::AbstractString)
    @stub "gene_signal_report"
end

# ---- Lecture 10 : building the gene-finding HMM -----------------------------

"""
    deterministic_emission_row(base) -> Vector{Float64}

Length-4 emission row that is `1.0` at `base`'s `BASE_INDEX` and `0.0` elsewhere
(used for motif on-ramp/off-ramp states: start codon, stop codons, donor `GT`,
acceptor `AG`).
"""
function deterministic_emission_row(base::Char)
    @stub "deterministic_emission_row"
end

"""
    build_skeleton_hmm(params) -> HMM

The 3-state skeleton `[:N, :E, :I]` (intergenic / exon / intron) with the
biologically-disallowed transitions (`N→I`, `I→N`) fixed to `0.0` and `I=[1,0,0]`
(always start intergenic). `params` supplies the trainable probabilities
(e.g. a NamedTuple of self-loop / switch probabilities and per-state emissions).
"""
function build_skeleton_hmm(params)
    @stub "build_skeleton_hmm"
end

"""
    build_codon_hmm(params) -> HMM

The elaborated gene model: exon split into codon-position states
`[:E0,:E1,:E2]` wired as a cycle, plus deterministic motif ramp states (start
codon, stop-codon sub-graph, donor/acceptor) and the 3×-cloned intron module
that keeps codon phase (forcing multiple-of-three exon length).
"""
function build_codon_hmm(params)
    @stub "build_codon_hmm"
end
