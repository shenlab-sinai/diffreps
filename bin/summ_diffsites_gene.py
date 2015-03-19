#!/usr/bin/env python

if __name__ == "__main__":

    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Summarize diffReps sites for \
                                                  each gene.",
                                     prog="summ_diffsites_gene.py")

    parser.add_argument("-a", "--abs", help="Use absolute log2FC", 
                        action="store_true")
    parser.add_argument("-f", "--filter", help="Filter based on location",
                        choices=["Genebody", "Promoter", "Promoter1k", 
                                 "Promoter3k", "ProximalPromoter"])
    parser.add_argument("diffreps", help="Annotated diffReps file", type=str)

    args = parser.parse_args()

    try:
        diff_f = open(args.diffreps)
    except IOError:
        print "Open diffReps file: {0} error.".format(diff_f)
        sys.exit()

    # header = ["Chrom", "Start", "End", "Length", "Treatment.cnt", "Control.cnt", 
    #           "Treatment.avg", "Control.avg", "Treatment.enr", "Control.enr", 
    #           "Event", "logFC", "pval", "padj", "winSta", "winEnd", "winFC", 
    #           "winP", "winQ", "GName", "TName", "Strand", "TSS", "TES", "Feature", 
    #           "D2TSS"]

    header = diff_f.readline().strip().split("\t")

    gene_tbl = {}
    for rec in diff_f:
        tokens = rec.strip().split("\t")
        site_tbl = dict(zip(header, tokens))
        feature = site_tbl["Feature"]
        gene = site_tbl["GName"]
        if gene == "" or args.filter != None and args.filter not in feature:
            continue
        if "log2FC" in header:
            log2fc = float(site_tbl["log2FC"])
        else:
            log2fc = float(site_tbl["logFC"])
        if args.abs:
            log2fc = abs(log2fc)

        if gene in gene_tbl:
            gene_tbl[gene].append(log2fc)
        else:
            gene_tbl[gene] = [log2fc]

    print "Gene\tAvgLog2FC"
    for k in gene_tbl.viewkeys():
        print "{0}\t{1}".format(k, sum(gene_tbl[k]) / len(gene_tbl[k]))

