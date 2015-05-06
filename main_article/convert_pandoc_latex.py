#!/usr/bin/env python2.7

import io
import sys
from uni2tex import *

def usage():
	print("")
	print("usage:\n\tpython " + sys.argv[0] + " <unicode_latex_file> <new_latex_file>")
	print("")

if (len(sys.argv) != 3):
	usage()
	quit()

fhandle = io.open(sys.argv[1], mode = 'r')
ohandle = io.open(sys.argv[2], mode = 'w')
for line in fhandle:
	ohandle.write(uni2tex(line))
fhandle.close()
ohandle.close()

