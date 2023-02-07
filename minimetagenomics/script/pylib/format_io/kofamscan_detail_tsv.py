#!/usr/bin/env python3

from .field_named_list import FieldNamedList


class KofamscanDetailTsvRecord(FieldNamedList):
	over_thres	= FieldNamedList.add_field(0)
	gene_name	= FieldNamedList.add_field(1)
	ko			= FieldNamedList.add_field(2)
	threshold	= FieldNamedList.add_field(3, float)
	score		= FieldNamedList.add_field(4, float)
	e_val		= FieldNamedList.add_field(5, float)
	ko_def		= FieldNamedList.add_field(6)

	@property
	def is_over_threshold(self):
		return self.over_thres == "*"
