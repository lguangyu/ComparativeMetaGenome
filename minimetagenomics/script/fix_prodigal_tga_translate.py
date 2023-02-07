#!/usr/bin/env python3

import argparse
import sys
# bio format i/o
import Bio
import Bio.SeqIO
# custom lib
import pylib


def get_args():
	ap = argparse.ArgumentParser(description = "prodigal (v2.6.3) treats codon "
		"TGA as W in its translation table 11; while in other parts of the "
		"world it is stop (*); this script 'hard-fix' this issue by replacing "
		"all stop (*) in the middle of faa sequences to 'W'")
	ap.add_argument("input", type = str, nargs = "?", default = "-",
		metavar = "faa",
		help = "input fasta (default: <stdin>)")
	ap.add_argument("-o", "--output", type = str, default = "-",
		metavar = "faa",
		help = "output fasta (default: <stdout>)")
	# parse and refine args
	args = ap.parse_args()
	if args.input == "-":
		args.input = sys.stdin
	if args.output == "-":
		args.output = sys.stdout
	return args


def main():
	args = get_args()
	with pylib.file_util.get_fp(args.input, "r") as ifp:
		with pylib.file_util.get_fp(args.output, "w") as ofp:
			for faa in Bio.SeqIO.parse(ifp, format = "fasta"):
				seq = str(faa.seq)
				seq = seq[:-1].replace("*", "W") + seq[-1]
				faa.seq = Bio.Seq.Seq(seq)
				Bio.SeqIO.write(faa, ofp, format = "fasta")
	return


if __name__ == "__main__":
	main()
