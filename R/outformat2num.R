# helper function to map human readable output format to numeric specification
# used as input by the BLAST command line tool
outformat2num <- function(out.format) {
    if (is_outformat(out.format)) {
        return(switch(out.format,
               pair = 0,
               qa.ident = 1,
               qa.nonident = 2,
               fq.ident = 3,
               fq.nonident = 4,
               xml = 5,
               tab = 6,
               tab.comment = 7,
               ASN.1.text = 8,
               ASN.1.binary = 9,
               csv = 10,
               ASN.1 = 11,
               json.seq.aln = 12,
               json.blast.multi = 13,
               xml2.blast.multi = 14,
               json.blast.single = 15,
               xml2.blast.single = 16,
               report = 18
               ))
        
    }
}
