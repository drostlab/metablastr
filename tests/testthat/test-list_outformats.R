context("Test: list_outformats()")

test_that("The list_outformats() function works properly..",{
    expect_equal(list_outformats(), c(
        "pair",
        "qa.ident",
        "qa.nonident",
        "fq.ident",
        "fq.nonident",
        "xml",
        "tab",
        "tab.comment",
        "ASN.1.text",
        "ASN.1.binary",
        "csv",
        "ASN.1",
        "json.seq.aln",
        "json.blast.multi",
        "xml2.blast.multi",
        "json.blast.single",
        "xml2.blast.single",
        "SAM",
        "report"
    ))
})
