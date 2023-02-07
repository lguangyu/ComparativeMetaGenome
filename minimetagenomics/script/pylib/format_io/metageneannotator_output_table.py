#!/usr/bin/env python3

from .field_named_list import FieldNamedList


class MetaGeneAnnotatorOutputTableRecord(FieldNamedList):
	gene_id		= FieldNamedList.add_field(0)
	start_pos	= FieldNamedList.add_field(1, int)
	end_pos		= FieldNamedList.add_field(2, int)
	strand		= FieldNamedList.add_field(3)
	frame		= FieldNamedList.add_field(4, int)
	complete	= FieldNamedList.add_field(5)
	score		= FieldNamedList.add_field(6, float)
	model		= FieldNamedList.add_field(7)
	rbs_start	= FieldNamedList.add_field(8)
	rbs_end		= FieldNamedList.add_field(9)
	rbs_score	= FieldNamedList.add_field(10)

	def compile_gff_attributes(self, sequence: str = ""):
		gff_feat = [
			"ID=" + (sequence + "_" if sequence else "") + self.gene_id,
			"complete=" + self.complete,
			"model=" + self.model,
			"rbs_start=" + self.rbs_start,
			"rbs_end=" + self.rbs_end,
			"rbs_score=" + self.rbs_score,
		]
		return (";").join(gff_feat)

	def to_gff(self, sequence: str = "", feat_type: str = ""):
		gff_field = [
			sequence,
			"MetaGeneAnnotator",
			feat_type,
			str(self.start_pos),
			str(self.end_pos),
			str(self.score),
			self.strand,
			str(self.frame),
			self.compile_gff_attributes(sequence),
		]
		return ("\t").join(gff_field)
