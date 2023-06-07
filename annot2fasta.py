def process(annotation, out):
    annot_dict = {}
    with open(annotation) as annot:
        for l in annot:
            s = l.split("\t")
            guide = s[0]
            gene = s[1].split("_")[1]
            if gene not in annot_dict:
                annot_dict[gene] = [guide]
            else:
                annot_dict[gene].append(guide)
    with open(out, 'w') as o:
        for g in annot_dict:
            i = 1
            for guide in annot_dict[g]:
                o.write(">%s_%s\n%s\n" % (g, i, guide) )
                i = i+1
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("parse MOFA library guide anntotation and convert to CB2 fasta")
    parser.add_argument("annotation", help="")
    parser.add_argument("out", help="")
    args = parser.parse_args()
    process(**vars(args))
