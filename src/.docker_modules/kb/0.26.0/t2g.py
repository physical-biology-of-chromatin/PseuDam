#!/usr/local/bin/python
import os
import re
import gzip
import argparse


def validate_file(f):
    if not os.path.exists(f):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(f))
    return f


def t2g_line(transcript, gene):
    return str(transcript) + "\t" + str(gene) + "\n"


def build_gene_re():
    return re.compile(".*gene_id\s+\"(\S+)\";.*")


def build_transcript_re():
    return re.compile(".*transcript_id\s+\"(\S+)\";.*")


def get_gene(line, gene_re):
    return gene_re.match(line)


def get_transcript(line, transcript_re):
    return transcript_re.match(line)


def gtf_line(line, transcript_re, gene_re):
    transcript_id = get_transcript(line, transcript_re)
    gene_id = get_gene(line, gene_re)
    return {'transcript_id': transcript_id, 'gene_id': gene_id}


def write_t2g_line(t2g, line, transcript_re, gene_re):
    results = gtf_line(line, transcript_re, gene_re)
    if results['transcript_id']:
        t2g.write(
            t2g_line(
                results['transcript_id'].group(1),
                results['gene_id'].group(1)
            )
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="create transcript to genes file from a gtf file."
    )
    parser.add_argument(
        "-g", "--gtf", dest="gtf", required=True, type=validate_file,
        help="gtf file", metavar="FILE"
    )
    args = parser.parse_args()
    gene_re = build_gene_re()
    transcript_re = build_transcript_re()

    with gzip.open(args.gtf, "rb") as gtf:
        with open("t2g_dup.txt", "w") as t2g:
            for line in gtf:
                write_t2g_line(t2g, str(line), transcript_re, gene_re)