#!/usr/bin/env python3

import argparse
import Bio.SeqIO
import sys
# custom lib
import pylib


def get_args():
	ap = argparse.ArgumentParser(description = "filter KO mapping by bin fasta;"
		" only records associated with headers in fasta will be in the output")
	ap.add_argument("input", type = str, nargs = "?", default = "-",
		help = "2-column tsv input as genes to KO mappings; compatible with "
			"blastkoala output and kofamscan mapping output format "
			"(default: <stdin>)")
	ap.add_argument("-f", "--bin-fasta", type = str, required = True,
		metavar = "fasta",
		help = "input bin fasta (required)")
	ap.add_argument("-o", "--output", type = str, default = "-",
		metavar = "tsv",
		help = "output filtered ko mapping (default: <stdout>)")
	# parse and refine args
	args = ap.parse_args()
	if args.input == "-":
		args.input = sys.stdin
	if args.output == "-":
		args.output = sys.stdout
	return args


def load_fasta_headers(file) -> set:
	ret = set()
	for seq in Bio.SeqIO.parse(file, format = "fasta"):
		ret.add(seq.name)
	return ret


def parse_contig_name_from_gene_name(gene_name, delimiter = "_"):
	*s, idx = gene_name.split(delimiter)
	return delimiter.join(s)


def filter_ko_mapping_by_contig_name(infile, outfile, contig_name_set: set):
	with pylib.file_util.get_fp(infile, "r") as ifp:
		with pylib.file_util.get_fp(outfile, "w") as ofp:
			for line in ifp:
				gene_name = line.split("\t")[0]
				contig_name = parse_contig_name_from_gene_name(gene_name)
				if contig_name in contig_name_set:
					ofp.write(line)
	return


def main():
	args = get_args()
	# load bin fasta info
	fasta_headers = load_fasta_headers(args.bin_fasta)
	# filter and save
	filter_ko_mapping_by_contig_name(args.input, args.output,
		contig_name_set = fasta_headers)
	return


if __name__ == "__main__":
	main()
