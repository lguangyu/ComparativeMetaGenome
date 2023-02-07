#!/usr/bin/env python3

import argparse
import collections
import sys
# custom lib
import pylib


def get_args():
	ap = argparse.ArgumentParser(description = "combine gene prediction from "
		"multile GFF files based on input priority")
	ap.add_argument("inputs", type = str, nargs = "+",
		metavar = "gff [gff ...]",
		help = "input GFF tables; when CDS regions from multiple GFF's overlap "
			"with each other, the selection will be decided according to this "
			"input list as priority; all other tying CDS's will be discarded")
	ap.add_argument("-l", "--min-length", type = int, default = "108",
		metavar = "int",
		help = "minimum CDS filter length; this filtering will be performed "
			"preceeding the overlapping check (default: 108, i.e. 36 aa)")
	ap.add_argument("-o", "--output", type = str, default = "-",
		metavar = "gff",
		help = "output GFF table (default: <stdout>)")
	# parse and refine args
	args = ap.parse_args()
	if args.output == "-":
		args.output = sys.stdout
	return args


class GffRecordWithPriority(pylib.format_io.gff.GffRecord):
	def __init__(self, *ka, priority: int, **kw):
		super().__init__(*ka, **kw)
		self.priority = priority
		return


def add_gff_cds_inplace(gff, dest_dict: collections.defaultdict, priority: int,
		min_length: int):
	"""
	parse input gff and append each record <rec> into list
	dest_dict[rec.sequence].
	"""
	priority = int(priority)
	with pylib.file_util.get_fp(gff, "r") as fp:
		for line in fp:
			line = line.rstrip("\r\n")
			if (not line) or line.startswith("#"):
				continue
			record = GffRecordWithPriority.from_str(line, priority = priority)
			# check min length
			if record.end_pos - record.start_pos + 1 >= min_length:
				dest_dict[record.sequence].append(record)
	return


def load_gffs_by_sequence(gff_list: list, min_length: int) -> dict:
	ret = collections.defaultdict(list)
	for i, gff in enumerate(gff_list):
		add_gff_cds_inplace(gff, dest_dict = ret, priority = i,
			min_length = min_length)
	return ret


class Position(object):
	def __init__(self, pos, is_end: bool, *ka, cds, **kw):
		super().__init__(*ka, **kw)
		self.pos	= pos
		self.is_end	= is_end
		self.cds	= cds
		return

	@property
	def is_start(self):
		return not self.is_end

	def in_prior_to(self, other) -> bool:
		return self.cds.priority < other.cds.priority


def resolve_conflict_cds(cds_list) -> list:
	# this list stored all start/end positions of all cds records
	pos_list = list()
	for cds in cds_list:
		pos_list.append(Position(cds.start_pos, is_end = False, cds = cds))
		pos_list.append(Position(cds.end_pos, is_end = True, cds = cds))
	# now we sort this list inplace by the position
	pos_list.sort(key = lambda x: x.pos)
	# travel through the sorted list, find wanted cds records
	ret = list()
	querying_start_pos = None
	for pos in pos_list:
		if pos.is_start:
			if querying_start_pos is None:
				querying_start_pos = pos
			elif pos.in_prior_to(querying_start_pos):
				# replace the querying_start_pos for higher priority
				# this will result in discarding the low-priority cds
				querying_start_pos = pos
		else:
			if querying_start_pos is None:
				# reached here because the mate start pos is discarded some time
				# earlier
				continue
			elif pos.cds is querying_start_pos.cds:
				# found a paired mate start, time to take that cds
				ret.append(pos.cds)
				querying_start_pos = None # important! reset this for next round
			# else, indicate the mate start pos is discarded some time earlier
	return ret


def main():
	args = get_args()
	# load gff cds
	cds = load_gffs_by_sequence(args.inputs, args.min_length)
	resolved_cds = {k: resolve_conflict_cds(v) for k, v in cds.items()}
	# save output
	with pylib.file_util.get_fp(args.output, "w") as fp:
		for seq in sorted(resolved_cds.keys()):
			for i in resolved_cds[seq]:
				print(i.to_str(), file = fp)
	return


if __name__ == "__main__":
	main()
