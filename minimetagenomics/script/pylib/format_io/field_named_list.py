#!/usr/bin/env python3

from . import base


class FieldNamedList(list, base.FormatIOBase):
	def __init__(self, *ka, **kw):
		super().__init__(*ka, **kw)
		return

	@classmethod
	def from_str(cls, s: str, *ka, delimiter = "\t", **kw):
		new = cls(s.split(delimiter), *ka, **kw)
		return new

	def to_str(self, delimiter = "\t"):
		return delimiter.join(self)

	@staticmethod
	def add_field(index, cast = None):
		def fget(self):
			return self[index] if cast is None else cast(self[index])
		def fset(self, value):
			self[index] = str(value)
			return
		return property(fget = fget, fset = fset)
