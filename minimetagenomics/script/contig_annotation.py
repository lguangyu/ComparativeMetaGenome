#!/usr/bin/env python3

import argparse
import collections
import sys
# custom lib
import pylib


def get_args():
	ap = argparse.ArgumentParser(description = "combine diamond blastp and "
		"blastdbcmd results to assign phylogeny to each contigs (LCA)")
	ap.add_argument("-c", "--seq-to-contig", type = str, required = True,
		metavar = "tsv",
		help = "seq id to contig map (required)")
	ap.add_argument("-a", "--seq-to-accs", type = str, required = True,
		metavar = "tsv",
		help = "seq id to accession map (required)")
	ap.add_argument("-t", "--accs-to-taxid", type = str, required = True,
		metavar = "tsv",
		help = "accession to tax id map (required)")
	ap.add_argument("-d", "--nodes-dump", type = str, required = True,
		metavar = "nodes.dmp",
		help = "nodes.dmp file in ncbi taxonomy db (required)")
	ap.add_argument("-n", "--names-dump", type = str, required = True,
		metavar = "names.dmp",
		help = "names.dmp file in ncbi taxonomy db (required)")
	ap.add_argument("--full-path", action = "store_true",
		help = "if set, all ranks in taxonomic tree will be reported; by "
			"default, only kingdom, phylumn, class, order, family, genus, "
			"species and strains will be reported (default: off)")
	ap.add_argument("--output-tax-path", type = str, required = True,
		metavar = "tsv",
		help = "output contig-sequence taxonomic path (required)")
	ap.add_argument("--output-lca", type = str, required = True,
		metavar = "tsv",
		help = "output contig LCA annotation (required)")
	# parse and refine args
	args = ap.parse_args()
	return args


def load_map_file(f, delimiter = "\t"):
	ret = dict()
	with pylib.file_util.get_fp(f, "r") as fp:
		for line in fp:
			line = line.rstrip("\r\n")
			if not line:
				continue
			k, v = line.split(delimiter)
			ret[k] = v
	return ret


def get_seq_to_taxid_map(seq_to_accs, accs_to_taxid):
	"""
	return a dict of str->int, with keys as seq names, and values are tax ids
	of each seq
	"""
	s2a	= load_map_file(seq_to_accs)
	a2t	= load_map_file(accs_to_taxid)
	s2t = {k : int(a2t.get(v, 0)) for k, v in s2a.items()} # use 0 as missing
	return s2t


def get_contig_seqs_dict(seq_to_contig):
	"""
	return a dict of str->set, with keys as contig names, and values are sets
	of seq names in each contig
	"""
	c2s = collections.defaultdict(set)
	s2c = load_map_file(seq_to_contig)
	for k, v in s2c.items():
		c2s[v].add(k)
	return c2s


def assign_tax_path(contig_seqs: dict, seq_to_taxid: dict, tax_tree, *,
		full_path = None):
	ret = collections.defaultdict(dict)
	for ctg, seqs in contig_seqs.items():
		for seq in seqs:
			# filter if tax id is 0 (guess to be error/missing in database)
			taxid = seq_to_taxid[seq]
			if (taxid != 0) and tax_tree.has_node(taxid):
				nd = tax_tree.get_node(taxid)
				ret[ctg][seq] = nd.tax_path if full_path else nd.simp_tax_path
	return ret


def lca(seq_tax_path: dict, proceed_frac = 0.7) -> dict:
	"""
	run non-strict LCA on each contig; if one tax has more than 70% votes at any
	tax-level, it will be selected and proceed to the next level;
	"""
	if (proceed_frac <= 0.5) or (proceed_frac > 1.0):
		raise ValueError("invalid proceed_frac, expected 0.5 < (value) <= 1.0")
	paths = list(seq_tax_path.values())
	count_thres = proceed_frac * len(paths)
	ret = "unclassified"
	for taxons in zip(*paths):
		counts = collections.Counter(taxons)
		for t, c in counts.items():
			if c >= count_thres:
				# assign as the current taxon
				# and continue to next level (outer loop)
				ret = t
				break
		else:
			# reaches here when inner loop finishes without break
			# i.e. no any taxon passed the threshold
			break
	return ret


def assign_contigs_lca(contig_seq_tax: dict, proceed_frac = 0.7) -> dict:
	ret = {k : lca(v) for k, v in contig_seq_tax.items()}
	return ret


def result_to_str_recursive(res, *, delimiter = "\t", prefix = "") -> str:
	if isinstance(res, dict):
		ret = ""
		for k in sorted(res.keys()):
			ret += result_to_str_recursive(res[k], delimiter = delimiter,
				prefix = k + delimiter)
	elif isinstance(res, list) or isinstance(res, set):
		ret = prefix + (";").join(res) + "\n"
	else:
		ret = prefix + str(res) + "\n"
	return ret


def save_output(f, res: dict, *, delimiter = "\t"):
	with pylib.file_util.get_fp(f, "w") as fp:
		fp.write(result_to_str_recursive(res, delimiter = delimiter))
	return


def main():
	args = get_args()
	# load ncbi tax tree
	ncbi_tax_tree	= pylib.format_io.ncbi_taxdump.NcbiTaxTree.parse_dump_files(
		nodes_file = args.nodes_dump, names_file = args.names_dump)
	# load mappings
	seq_to_taxid	= get_seq_to_taxid_map(args.seq_to_accs, args.accs_to_taxid)
	contig_seqs		= get_contig_seqs_dict(args.seq_to_contig)
	# assign tax
	contig_seq_tax	= assign_tax_path(contig_seqs, seq_to_taxid, ncbi_tax_tree,
		full_path = args.full_path)
	contig_lca		= assign_contigs_lca(contig_seq_tax, proceed_frac = 0.5)
	# save output
	save_output(args.output_tax_path, contig_seq_tax)
	save_output(args.output_lca, contig_lca)
	return


if __name__ == "__main__":
	main()
