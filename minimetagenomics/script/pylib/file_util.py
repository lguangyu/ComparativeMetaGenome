#!/usr/bin/env python3

import io
import sys


def get_fp(f, *ka, factory = open, **kw) -> io.IOBase:
	"""
	wrapper function (originally over open()) for safely getting file handle
	when opened file handled used instead of file name

	ARGUMENTS
	---------
	f: [io.IOBase or str]
		opened file handle or str;
	factory: [callable]
		wrapped function to open a file handle;
	*ka, **kw:
		other arguments passed to factory

	RETURN
	------
	io.IOBase:
		if argument f is a file handle, directly return f; if f is a str, return
		factory(f, *ka, **kw);

	EXCEPTION
	---------
	raise TypeError if f is neither io.IOBase nor str;
	"""
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = factory(f, *ka, **kw)
	else:
		raise TypeError("f must be file handle or str, not '%s'"\
			% type(f).__name__)
	return ret
