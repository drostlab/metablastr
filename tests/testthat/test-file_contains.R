context("Test: file_contains_*() ")

test_that("file_contains_aa() works properly with input aa fasta file ...", {
  expect_null(file_contains_aa(system.file('seqs/qry_aa.fa', package = 'metablastr'), "query"))
  
})

test_that("file_contains_aa() throws error when input file contains dna rather than aa ...", {
  expect_error(file_contains_aa(system.file('seqs/qry_nn.fa', package = 'metablastr'), "query"))
  
})

test_that("file_contains_dna() works properly with input dna fasta file ...", {
  expect_null(file_contains_dna(system.file('seqs/qry_nn.fa', package = 'metablastr'), "query"))
  
})

test_that("file_contains_dna() throws error when input file contains aa rather than dna ...", {
  expect_error(file_contains_dna(system.file('seqs/qry_aa.fa', package = 'metablastr'), "query")
)
  
})

