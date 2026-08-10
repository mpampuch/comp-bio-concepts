using Test
using SeqModel

# These tests describe the *intended* behaviour of every stub in `src/`.
# With the scaffold's unimplemented stubs they all fail/error on purpose
# (the "red" of red-green-refactor). Implement a function, re-run `] test`,
# and watch its testset turn green.
#
# Stub calls are inlined inside `@test` so a `NotImplemented` throw is recorded
# as a single failing test rather than aborting an entire `@testset`.

# --- shared fixtures (plain data, never throws) ------------------------------

# Lecture-4 trained inside-island transition table, rows = prev, cols = curr,
# base order A,C,G,T. Rows sum to 1.
const A_TABLE = [0.18 0.27 0.43 0.12;
                 0.17 0.37 0.27 0.19;
                 0.16 0.34 0.38 0.12;
                 0.08 0.36 0.38 0.18]

# Lecture-5 inside/outside chains chosen so CpG-rich scores > 0 and AT-rich < 0.
const A_IN = [0.20 0.30 0.30 0.20;
              0.10 0.10 0.70 0.10;
              0.20 0.30 0.30 0.20;
              0.20 0.30 0.30 0.20]
const A_OUT = [0.30 0.20 0.20 0.30;
               0.30 0.30 0.10 0.30;
               0.30 0.20 0.20 0.30;
               0.30 0.20 0.20 0.30]

@testset "SeqModel" begin

    # =====================================================================
    @testset "L1 DNA basics" begin
        @test gc_content("GATTACA") ≈ 1 / 7
        @test gc_content("GCGCGC") ≈ 1.0
        @test count_dinucleotide("CGCGCG", "CG") == 3
        @test count_dinucleotide("AAAA", "CG") == 0
        @test cpg_count("CGATCG") == 2
        @test length(random_dna(50)) == 50         # not a stub: sanity
    end

    # =====================================================================
    @testset "L2 probability" begin
        @test length(dice_space()) == 36
        @test prob(isodd_d1, dice_space()) ≈ 0.5
        @test prob(iseven_d2, dice_space()) ≈ 0.5
        @test joint(isodd_d1, iseven_d2, dice_space()) ≈ 0.25
        @test conditional(isodd_d1, iseven_d2, dice_space()) ≈ 0.5
        # independence: P(A,B) ≈ P(A)·P(B)
        @test joint(isodd_d1, iseven_d2, dice_space()) ≈
              prob(isodd_d1, dice_space()) * prob(iseven_d2, dice_space())
    end

    # =====================================================================
    @testset "L3 Markov assumption" begin
        examples = ["ACGT", "ACGT", "TTTT", "GGGG"]
        @test naive_string_prob("ACGT", examples) ≈ 0.5
        @test naive_string_prob("CCCC", examples) == 0.0   # unseen -> 0
        @test markov_factorization_terms("GATC") ==
              [('G', nothing), ('A', 'G'), ('T', 'A'), ('C', 'T')]
    end

    # =====================================================================
    @testset "L4 Markov chains" begin
        # counts -> row-normalized probabilities
        @test all(≈(1.0), sum(transition_matrix(transition_counts(["ACGTACGT"])), dims = 2))
        # rows of the given table sum to 1
        @test all(≈(1.0), sum(A_TABLE, dims = 2))
        # P("GATC") = 0.25 * A[G,A] * A[A,T] * A[T,C]
        @test seq_prob("GATC", A_TABLE) ≈ 0.25 * 0.16 * 0.12 * 0.36
    end

    # =====================================================================
    @testset "L5 log space + CpG score" begin
        # log-space round-trips to the linear probability
        @test exp(logspace_seq_logprob("GATC", A_TABLE)) ≈ seq_prob("GATC", A_TABLE)
        # underflow motivation (pure Julia facts, always pass)
        @test prod(fill(0.5, 1100)) == 0.0
        @test sum(fill(log(0.5), 1100)) < 0 && isfinite(sum(fill(log(0.5), 1100)))
        # combined log-ratio table shape
        @test size(logratio_table(A_IN, A_OUT)) == (4, 4)
        # decision rule: CpG-rich > 0, AT-rich < 0
        @test cpg_score("CGCGCGCGCG", logratio_table(A_IN, A_OUT)) > 0
        @test cpg_score("ATATATATAT", logratio_table(A_IN, A_OUT)) < 0
    end

    # =====================================================================
    @testset "L6 HMM definition (casino)" begin
        @test casino().states == [:Fair, :Loaded]
        @test all(≈(1.0), sum(casino().A, dims = 2))
        @test all(≈(1.0), sum(casino().E, dims = 2))
        @test sum(casino().I) ≈ 1.0
        # joint of (Fair,Fair,Loaded) emitting "HHT"
        @test joint_prob(casino(), [:Fair, :Fair, :Loaded], "HHT") ≈ 0.006
        @test joint_logprob(casino(), [:Fair, :Fair, :Loaded], "HHT") ≈ log(0.006)
        # simulate returns equal-length path/emissions
        @test length(simulate(casino(), 10)[1]) == 10
        @test length(simulate(casino(), 10)[2]) == 10
    end

    # =====================================================================
    @testset "L7 Viterbi decoding" begin
        # sticky casino: "HH" decodes to two loaded coins (hand-worked)
        @test viterbi(sticky_casino(), "HH")[1] == [:Loaded, :Loaded]
        # logprob agrees with joint_logprob of the returned path
        let (p, lp) = viterbi(sticky_casino(), "HH")
            @test lp ≈ joint_logprob(sticky_casino(), p, "HH")
        end
        # DP matrix has |Q| rows and n columns
        @test size(viterbi_matrix(sticky_casino(), "HTTHH")[1]) == (2, 5)
    end

    # =====================================================================
    @testset "L8 CpG-island HMMs" begin
        # naive 2-state training shapes
        @test size(train_io_transitions("IIOOIO")) == (2, 2)
        @test all(≈(1.0), sum(train_io_transitions("IIOOIO"), dims = 2))
        @test size(train_io_emissions("ACGTAC", "IIOOIO")) == (2, 4)

        # fixed (dinucleotide) model
        @test length(fixed_states()) == 8
        @test size(fixed_emissions()) == (8, 4)
        @test all(≈(1.0), sum(fixed_emissions(), dims = 2))  # deterministic rows
        @test size(train_fixed_transitions("ACGTACGT", "IIIIOOOO")) == (8, 8)

        # integration: the fixed model must NOT call a CpG-free C-run an island
        # (this is the bug the naive model has; see Lecture 8).
        seq    = "AAAACGCGCGCGAAAA"
        labels = "OOOOIIIIIIIIOOOO"
        hmm = build_fixed_hmm(seq, labels)
        # a C-run has no CG dinucleotide -> should decode all outside (_O states)
        path = viterbi(hmm, "AAAACCCCCCCCAAAA")[1]
        @test all(s -> endswith(String(s), "_O"), path)
    end

    # =====================================================================
    @testset "L9 gene signals" begin
        @test codons("ATGAAATAA") == ["ATG", "AAA", "TAA"]
        @test translate("ATGAAATAA") == "MK*"
        @test is_donor("AGTC", 2)        # positions 2:3 == "GT"
        @test is_acceptor("CAGT", 2)     # positions 2:3 == "AG"
        @test !is_donor("AAAA", 1)
        # ORF: ATG ... first in-frame stop, inclusive
        @test (1, 9) in find_orfs("ATGAAATAA")
        let r = gene_signal_report("ATGAAATAA")
            @test r.starts_atg
            @test r.ends_stop
            @test r.multiple_of_three
        end
    end

    # =====================================================================
    @testset "L10 gene-finding HMM construction" begin
        @test deterministic_emission_row('G') == [0.0, 0.0, 1.0, 0.0]
        @test deterministic_emission_row('A') == [1.0, 0.0, 0.0, 0.0]

        # skeleton model invariants (params schema documented in build_skeleton_hmm)
        params = (p_stay_N = 0.99, p_enter_gene = 0.01,
                  p_stay_E = 0.90, p_E_to_I = 0.05, p_E_to_N = 0.05,
                  p_stay_I = 0.90, p_I_to_E = 0.10,
                  emit_N = fill(0.25, 4), emit_E = fill(0.25, 4), emit_I = fill(0.25, 4))
        skel = build_skeleton_hmm(params)
        @test skel.states == [:N, :E, :I]
        @test skel.I == [1.0, 0.0, 0.0]                       # always start intergenic
        ni = findfirst(==(:N), skel.states); ii = findfirst(==(:I), skel.states)
        @test skel.A[ni, ii] == 0.0                            # no N -> I
        @test skel.A[ii, ni] == 0.0                            # no I -> N
        @test all(≈(1.0), sum(skel.A, dims = 2))

        # codon-position model: returns an HMM that still starts intergenic
        cod = build_codon_hmm(params)
        @test cod isa HMM
        @test :E0 in cod.states && :E1 in cod.states && :E2 in cod.states
        @test cod.states[argmax(cod.I)] == :N
    end

    # =====================================================================
    @testset "Pluto Notebook Execution" begin
        notebooks_dir = joinpath(@__DIR__, "..", "notebooks")
        nb_files = filter(f -> endswith(f, ".jl"), readdir(notebooks_dir, join = true))
        @test length(nb_files) == 6
        for nb in sort(nb_files)
            @testset "Notebook: $(basename(nb))" begin
                @test_nowarn include(nb)
            end
        end
    end

end

