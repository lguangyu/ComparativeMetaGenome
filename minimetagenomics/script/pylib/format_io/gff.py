#!/usr/bin/env python3

from .field_named_list import FieldNamedList


class GffRecord(FieldNamedList):
	sequence	= FieldNamedList.add_field(0)
	source		= FieldNamedList.add_field(1)
	feature		= FieldNamedList.add_field(2)
	start_pos	= FieldNamedList.add_field(3, int)
	end_pos		= FieldNamedList.add_field(4, int)
	score		= FieldNamedList.add_field(5, float)
	strand		= FieldNamedList.add_field(6)
	phase		= FieldNamedList.add_field(7, int)
	attributes	= FieldNamedList.add_field(8)

	@property
	def on_forward_strand(self):
		return self.strand == "+"

	@property
	def on_reverse_strand(self):
		return not self.on_forward_strand
