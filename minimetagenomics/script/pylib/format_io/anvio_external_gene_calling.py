#!/usr/bin/env python3

from .field_named_list import FieldNamedList
from . import gff


class AnvioExternalGeneCallingRecord(FieldNamedList):
	gid			= FieldNamedList.add_field(0, int)
	contig		= FieldNamedList.add_field(1)
	start		= FieldNamedList.add_field(2, int)
	stop		= FieldNamedList.add_field(3, int)
	direction	= FieldNamedList.add_field(4)
	partial		= FieldNamedList.add_field(5, int)
	call_type	= FieldNamedList.add_field(6, int)
	source		= FieldNamedList.add_field(7)
	version		= FieldNamedList.add_field(8)
	aa_sequence	= FieldNamedList.add_field(9)

	@classmethod
	def parse_from_gff_record(cls, gff_record: gff.GffRecord, *,
			gid: int, aa_sequence: str, call_type = 1,
			source: str = "", version: str = "", **kw):
		if not isinstance(gff_record, gff.GffRecord):
			raise TypeError("gff_record must be GffRecord, not '%s'"\
				% type(gff_record).__name__)
		new = cls([""] * 10, **kw)
		new.gid			= gid
		new.contig		= gff_record.sequence
		new.start		= gff_record.start_pos - 1
		new.stop		= gff_record.end_pos
		new.direction	= ("f" if gff_record.on_forward_strand else "r")
		new.partial		= (0 if gff_record.phase == 0 else 1)
		new.call_type	= call_type # for bacteria, it's almost 1
		new.source		= source
		new.version		= version
		new.aa_sequence	= aa_sequence
		return new
