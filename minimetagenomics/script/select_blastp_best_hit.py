#!/usr/bin/env python3

import argparse
#import collections
import sys
# custom lib
import pylib
from pylib.format_io.field_named_list import FieldNamedList


def get_args():
	ap = argparse.ArgumentParser(description = "find the best hit for each cds "
		"in blastp output")
	ap.add_argument("input", type = str, nargs = "?", default = "-",
		metavar = "tsv",
		help = "input blastp output (default: <stdin>)")
	ap.add_argument("-o", "--output", type = str, default = "-",
		metavar = "txt",
		help = "output filtered blastp output (default: <stdout>)")
	# parse and refine args
	args = ap.parse_args()
	if args.input == "-":
		args.input = sys.stdin
	if args.output == "-":
		args.output = sys.stdout
	return args


class BlastpRecord(FieldNamedList):
	query_id	= FieldNamedList.add_field(0)
	hit_accs	= FieldNamedList.add_field(1)
	evalue		= FieldNamedList.add_field(2, float)
	score		= FieldNamedList.add_field(3, int)
	bitscore	= FieldNamedList.add_field(4, float)
	length		= FieldNamedList.add_field(5, int)
	mismatch	= FieldNamedList.add_field(6, int)
	pident		= FieldNamedList.add_field(7, float)

	@property
	def is_qualified(self):
		# fitering method
		return True

	def is_better_than(self, other):
		return self.evalue < other.evalue


def load_blastp(f) -> list:
	ret = list()
	with pylib.file_util.get_fp(f, "r") as fp:
		for line in fp:
			line = line.rstrip("\r\n")
			if line:
				ret.append(BlastpRecord.from_str(line))
	return ret


def select_blastp_best_hit(blastp_recs) -> dict:
	ret = dict()
	for r in blastp_recs:
		if r.is_qualified:
			if r.query_id not in ret:
				ret[r.query_id] = r
			elif r.is_better_than(ret[r.query_id]):
				ret[r.query_id] = r
	return ret


def save_best_hit(f, best_hits: dict):
	with pylib.file_util.get_fp(f, "w") as fp:
		keys = sorted(best_hits.keys())
		for k in keys:
			print(best_hits[k].to_str(), file = fp)
	return


def main():
	args = get_args()
	blastp_recs	= load_blastp(args.input)
	best_hits	= select_blastp_best_hit(blastp_recs)
	save_best_hit(args.output, best_hits)
	return


if __name__ == "__main__":
	main()
