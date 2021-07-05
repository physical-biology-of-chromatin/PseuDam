#!/usr/local/bin/python
import os
import gffutils
import argparse

def validate_file(f):
    if not os.path.exists(f):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(f))
    return f


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="create transcript to genes file from a gtf file."
    )
    parser.add_argument(
        "-g", "--gtf", dest="gtf", required=True, type=validate_file,
        help="gtf file", metavar="FILE"
    )
    args = parser.parse_args()

    db = gffutils.create_db(
        args.gtf,
        dbfn=":memory:",
        force=True,
        merge_strategy="merge",
        disable_infer_transcripts=True,
        disable_infer_genes=True
    )
    with open("t2g.txt", "w") as t2g:
        for gene in db.all_features():
            for transcript in db.children(
              gene, featuretype='transcript', order_by='start'):
                t2g.write(
                    str(transcript["transcript_id"][0]) +
                    "\t" +
                    str(gene["gene_id"][0]) +
                    "\n"
                )
