context('test extract_random_seqs_from_genome')

test_that("No filtering Ns", {
    genome_fa <- system.file('seqs', 'sbj_aa.fa', package = "metablastr")
    output_fa <- tempfile()
    
    set.seed(123)
    extract_random_seqs_from_genome(subject_genome = genome_fa,
                                    interval_width = 10,
                                    size = 20,
                                    file_name = output_fa)
    
    rand_seq <- Biostrings::readAAStringSet(output_fa)

    expect_equal(length(rand_seq), 20)
    expect_true(all(lengths(rand_seq) == 10))
})

test_that("Filtering sequences with count N > n_max", {
    
    genome_fa <- system.file('seqs', 'sbj_aa.fa', package = "metablastr")
    output_fa <- tempfile()
    
    set.seed(123)
    extract_random_seqs_from_genome(subject_genome = genome_fa,
                                    interval_width = 10,
                                    size = 20,
                                    file_name = output_fa,
                                    n_max = 1)
    
    rand_seq <- Biostrings::readAAStringSet(output_fa)
    nn <- Biostrings::vcountPattern('N', rand_seq)
    
    expect_equal(length(rand_seq), sum(nn <= 1))
    expect_true(all(lengths(rand_seq) == 10))
})
