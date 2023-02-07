#!/usr/bin/env python3

from .field_named_list import FieldNamedList
from ..import file_util


class NcbiTaxdumpRecordBase(FieldNamedList):
	@classmethod
	def from_str(cls, *ka, delimiter = "\t|\t", **kw):
		return super().from_str(*ka, delimiter = delimiter, **kw)

	def to_str(self, delimiter = "\t|\t"):
		return super().to_str(delimiter = delimiter)


class NcbiTaxdumpNodeRecord(NcbiTaxdumpRecordBase):
	tax_id			= NcbiTaxdumpRecordBase.add_field(0, int) # tax_id
	parent_tax_id	= NcbiTaxdumpRecordBase.add_field(1, int) # parent tax_id
	rank			= NcbiTaxdumpRecordBase.add_field(2) # rank
	embl_code		= NcbiTaxdumpRecordBase.add_field(3) # embl code
	div_id			= NcbiTaxdumpRecordBase.add_field(4) # division id
	inh_div_flag	= NcbiTaxdumpRecordBase.add_field(5) # inherited div flag (1 or 0)
	gnt_code_id		= NcbiTaxdumpRecordBase.add_field(6, int) # genetic code id
	inh_gc_flag		= NcbiTaxdumpRecordBase.add_field(7) # inherited GC flag (1 or 0)
	mt_gnt_code_id	= NcbiTaxdumpRecordBase.add_field(8, int) # mitochondrial genetic code id
	inh_mgc_flag	= NcbiTaxdumpRecordBase.add_field(9) # inherited MGC flag (1 or 0)
	gbk_hid_flag	= NcbiTaxdumpRecordBase.add_field(10) # GenBank hidden flag (1 or 0)
	hid_subtrt_flag	= NcbiTaxdumpRecordBase.add_field(11) # hidden subtree root flag (1 or 0)
	comments		= NcbiTaxdumpRecordBase.add_field(12) # comments


class NcbiTaxdumpNameRecord(NcbiTaxdumpRecordBase):
	tax_id			= FieldNamedList.add_field(0, int) # tax_id
	name_txt		= FieldNamedList.add_field(1) # name_txt
	unique_name		= FieldNamedList.add_field(2) # unique name
	name_class		= FieldNamedList.add_field(3) # name class


class NcbiTaxdumpNameMixin(object):
	"""
	this mixin class is used specifically with class NcbiTaxdumpNode to handle
	the taxonomy (specifically, scientific name) features; therefore must be
	used as a mixin class with the later one;
	"""
	def __init__(self, *ka, sci_name = None, **kw):
		super().__init__(*ka, **kw)
		self.sci_name = sci_name
		return

	@property
	def tax_path(self) -> list:
		if self.is_root:
			ret = list()
		else:
			ret = self.parent.tax_path
			ret.append(self.sci_name)
		return ret

	SIMP_RANKS = ("superkingdom", "kingdom", "phylum", "class", "order",
		"family", "genus", "species", "strain")

	@property
	def simp_tax_path(self) -> list:
		"""
		similar to .tax_path, but only extract ranks in:
		kingdom, phylum, class, order, family, genus, species, strain
		"""
		if self.is_root:
			ret = list()
		else:
			ret = self.parent.simp_tax_path
			if self.rank in self.SIMP_RANKS:
				ret.append(self.sci_name)
		return ret


class NcbiTaxdumpNode(NcbiTaxdumpNodeRecord, NcbiTaxdumpNameMixin):
	def __init__(self, *ka, parent = None, **kw):
		super().__init__(*ka, **kw)
		self.children = dict()
		self.set_parent(parent)
		return

	def set_parent(self, parent = None):
		if parent is None:
			self.parent = parent
		elif isinstance(parent, NcbiTaxdumpNode):
			self.parent = parent
			self.parent.children[self.tax_id] = self
		else:
			raise TypeError("parent must be NcbiTaxdumpNode, not '%s'"\
				% type(parent).__name__)
		return

	@property
	def is_root(self):
		# root as parent is none
		return self.parent is None


class NcbiTaxTree(object):
	def __init__(self, *ka, **kw):
		super().__init__(*ka, **kw)
		self.nodes_dict = dict()
		self.root = None
		return

	def has_node(self, tax_id):
		return tax_id in self.nodes_dict

	def get_node(self, tax_id):
		return self.nodes_dict[tax_id]

	def _add_node(self, node):
		# add to tree's node dict
		self.nodes_dict[node.tax_id] = node
		return

	def _all_nodes_set_parent(self):
		for node in self.nodes_dict.values():
			# add parent-child relationship
			# current known tax_id == parent_tax_id is the root
			if node.tax_id == node.parent_tax_id:
				node.set_parent(None)
				self.root = node
			else:
				node.set_parent(self.get_node(node.parent_tax_id))
		return

	def _parse_nodes_dump(self, nodes_file):
		with file_util.get_fp(nodes_file, "r") as fp:
			for line in fp:
				line = line.rstrip("\r\n")
				if line.endswith("\t|"):
					line = line[:-2]
				# parse as node
				node = NcbiTaxdumpNode.from_str(line)
				self._add_node(node)
		self._all_nodes_set_parent()
		return

	def _parse_names_dump(self, names_file):
		with file_util.get_fp(names_file, "r") as fp:
			for line in fp:
				line = line.rstrip("\r\n")
				if line.endswith("\t|"):
					line = line[:-2]
				# parse name and add to node
				name = NcbiTaxdumpNameRecord.from_str(line)
				if name.name_class == "scientific name":
					self.get_node(name.tax_id).sci_name\
						= (name.unique_name or name.name_txt)
		return

	@classmethod
	def parse_dump_files(cls, nodes_file, names_file, *ka, **kw):
		new = cls(*ka, **kw)
		new._parse_nodes_dump(nodes_file)
		new._parse_names_dump(names_file)
		return new
