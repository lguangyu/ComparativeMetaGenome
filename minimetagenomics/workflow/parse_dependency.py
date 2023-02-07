#!/usr/bin/env python3

import re
import sys

dep_ids = set([i for i in sys.argv[1:] if re.match("^[0-9]+$", i)])

if dep_ids:
	sys.stdout.write("--dependency=aftercorr:" + ((",").join(dep_ids)))
