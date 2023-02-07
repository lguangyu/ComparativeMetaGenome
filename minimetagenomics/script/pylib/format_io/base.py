#!/usr/bin/env python3

import abc


class FormatIOBase(abc.ABC):
	@classmethod
	@abc.abstractmethod
	def from_str(cls, s: str, *ka, **kw):
		pass

	@abc.abstractmethod
	def to_str(self, *ka, **kw):
		pass
