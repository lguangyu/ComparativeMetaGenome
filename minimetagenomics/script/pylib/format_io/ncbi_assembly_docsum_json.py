#!/usr/bin/env python3

import collections
import json
from ..import file_util


class NcbiAssemblyDocsumJsonRecord(dict):
	_quality_order_ = {
		"Complete Genome": 0,
		"Contig": 1,
		"Scaffold": 2,
	}

	@property
	def accession(self):
		return self.get("lastmajorreleaseaccession", None) or\
			self["assemblyaccession"]
	@property
	def assembly_name(self):
		return self["assemblyname"]
	@property
	def n50(self):
		return int(self["contign50"])
	@property
	def coverage(self):
		return float(self["coverage"])
	@property
	def taxid(self):
		return int(self["taxid"])
	@property
	def organism(self):
		return self["organism"]
	@property
	def speciesname(self):
		return self["speciesname"]

	def is_better_than(self, other):
		if not isinstance(other, NcbiAssemblyDocsumJsonRecord):
			raise TypeError("can't campare between '%s' and '%s'"\
				% (type(self).__name__, type(other).__name__))
		return (self.n50 > other.n50) or (self.coverage > other.coverage)

	def to_tabluar(self, delimiter = "\t"):
		return delimiter.join([
			self.accession,
			self.assembly_name,
			str(self.taxid),
			self.speciesname,
			self.organism,
			self["assemblystatus"],
			str(self.n50),
			self["ftppath_genbank"],
			self["ftppath_refseq"],
		])

	def __lt__(self, other):
		return self.is_better_than(other)
		

class NcbiAssemblyDocsumJson(dict):
	@classmethod
	def from_files(cls, *fnames, **kw):
		new = cls(**kw)
		new.setdefault("uids", list())
		for f in fnames:
			with file_util.get_fp(f, "r") as fp:
				raw = json.load(fp)
			new["uids"].extend(raw["result"]["uids"])
			for i in raw["result"]["uids"]:
				new[i] = NcbiAssemblyDocsumJsonRecord(raw["result"][i])
		return new

	def group_by_species(self) -> dict:
		ret = collections.defaultdict(list)
		for i in self["uids"]:
			rec = self[i]
			ret[rec.speciesname].append(rec)
		return ret

	def pick_best_record_per_species(self) -> list:
		asm_by_species = self.group_by_species()
		return [sorted(v)[0] for v in asm_by_species.values()]

