#!/usr/bin/env python3

import argparse
import collections
# bio format i/o
import Bio
import Bio.SeqIO
# custom lib
import pylib


def get_args():
	ap = argparse.ArgumentParser(description = "extract GFF-defined CDS from "
		"nucleotide fasta")
	gp = ap.add_argument_group("input")
	gp.add_argument("-g", "--gff", type = str, required = True,
		metavar = "gff",
		help = "input GFF table (required)")
	gp.add_argument("-f", "--fasta", type = str, required = True,
		metavar = "fna",
		help = "input fasta file, assumed to be nucleotides (required)")
	gp = ap.add_argument_group("output options")
	gp.add_argument("-n", "--cds-nucleotide", type = str,
		metavar = "fna",
		help = "fasta output as cds in nucleotides (optional)")
	gp.add_argument("-a", "--cds-amino-acid", type = str,
		metavar = "faa",
		help = "fasta output as cds in amino acids (optional)")
	# parse and refine args
	args = ap.parse_args()
	return args


class NameIdAssigner(collections.Counter):
	"""
	a name-counter that can provide identical & serial ids to each str called
	through; for example:

	>>> nia = NameIdAssigner()
	>>> nia.update("a") # called the 1st time with name "a"
	1
	>>> nia("a") # called the 2nd time with "a"
	a_2

	>>> nia("b") # called the 1st time with "b"
	b_1
	>>> nia.update("b") # called the 2nd time with name "b"
	2
	"""
	def update(self, s: str) -> None:
		"""
		update the counter with key a; unlike collections.Counter.update(), this
		function also returns the values of counter after update
		"""
		if (s is not None) and (not isinstance(s, str)):
			raise TypeError("s must be str, not '%s'" % type(s).__name__)
		super().update([s])
		return self[s]

	def __call__(self, s: str, fmt = "%u"):
		return (s + "_" + fmt) % self.update(s)


def extract_gff_cds(gff_in, fna_in) -> (list, list):
	"""
	extract cds sequences (as nucleotides) from <fna_in> based on cds records in
	<gff_in>; returns a list of cds (as nucleotides);
	"""
	# load fna as dict
	with pylib.file_util.get_fp(fna_in, "r") as fp:
		fna_ref = Bio.SeqIO.to_dict(Bio.SeqIO.parse(fp, format = "fasta"))
	# 
	seq_name_resolver = NameIdAssigner()
	cds_nuc_list = list()
	with pylib.file_util.get_fp(gff_in, "r") as fp:
		for line in fp:
			line = line.rstrip("\n\r")
			if (not line) or line.startswith("#"):
				continue
			gff = pylib.format_io.gff.GffRecord.from_str(line)
			if gff.feature.lower() == "cds":
				# get cds nucleotide sequence
				cds_nuc = fna_ref[gff.sequence][gff.start_pos - 1 : gff.end_pos]
				if gff.strand == "-":
					cds_nuc = cds_nuc.reverse_complement()
				cds_nuc = cds_nuc[len(cds_nuc) % 3:]
				# add seq information
				cds_nuc.id = seq_name_resolver(gff.sequence, fmt = "%03u")
				cds_nuc.description = ""
				cds_nuc_list.append(cds_nuc)
	return cds_nuc_list


def save_cds_fna(fasta, cds_list):
	if fasta: # also skip if fasta == ""
		with pylib.file_util.get_fp(fasta, "w") as fp:
			for cds in cds_list:
				Bio.SeqIO.write(cds, fp, format = "fasta")
	return


def save_cds_faa(fasta, cds_list, *, table = 11):
	if fasta: # also skip if fasta == ""
		with pylib.file_util.get_fp(fasta, "w") as fp:
			for cds in cds_list:
				if len(cds) % 3:
					print(cds)
				cds_trans = cds.translate(table = table)
				cds_trans.id = cds.id
				cds_trans.description = cds.description
				Bio.SeqIO.write(cds_trans, fp, format = "fasta")
	return


def main():
	args = get_args()
	if (args.cds_nucleotide is args.cds_amino_acid is None):
		exit("terminating: no output specified")
	# load gff cds
	cds_nuc_list = extract_gff_cds(args.gff, args.fasta)
	# save 
	save_cds_fna(args.cds_nucleotide, cds_nuc_list)
	save_cds_faa(args.cds_amino_acid, cds_nuc_list)
	return


if __name__ == "__main__":
	main()
