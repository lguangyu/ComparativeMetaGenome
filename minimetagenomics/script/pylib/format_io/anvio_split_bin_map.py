#!/usr/bin/env python3

import re
from .field_named_list import FieldNamedList


class AnvioSplitBinMap(FieldNamedList):
	split	= FieldNamedList.add_field(0)
	bin		= FieldNamedList.add_field(1)

	@property
	def contig(self):
		"""
		parse contig from split id
		"""
		ret = re.sub(r"_split_\d+$", "", self.split)
		return ret
