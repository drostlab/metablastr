context("Test: outformat2num()")

test_that("The outformat2num() function works properly..",{
    expect_equal(outformat2num("pair"), 0)
    expect_equal(outformat2num("qa.ident"), 1)
    expect_equal(outformat2num("qa.nonident"), 2)
    expect_equal(outformat2num("fq.ident"), 3)
    expect_equal(outformat2num("fq.nonident"), 4)
    expect_equal(outformat2num("xml"), 5)
    expect_equal(outformat2num("tab"), 6)
    expect_equal(outformat2num("tab.comment"), 7)
    expect_equal(outformat2num("ASN.1.text"), 8)
    expect_equal(outformat2num("ASN.1.binary"), 9)
    expect_equal(outformat2num("csv"), 10)
    expect_equal(outformat2num("ASN.1"), 11)
    expect_equal(outformat2num("json.seq.aln"), 12)
    expect_equal(outformat2num("json.blast.multi"), 13)
    expect_equal(outformat2num("xml2.blast.multi"), 14)
    expect_equal(outformat2num("json.blast.single"), 15)
    expect_equal(outformat2num("xml2.blast.single"), 16)
    expect_equal(outformat2num("report"), 18)
    
    
})
