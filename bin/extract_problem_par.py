#!/usr/bin/env python

import re
import sys

print(sys.argv)
if len(sys.argv) != 2:
    print("USAGE: extract_problem.par LOGFILE")
    sys.exit(1)

f = open(sys.argv[1])
data = f.readlines()
f.close()

namelist = re.compile('\&[A-Z]')
change = re.compile('\*\s[A-Z]')
slash = re.compile('0:\s{1,4}/$')

for line in data:
    line = line.strip()
    if re.search(namelist, line):
        print('$%s' % (line.split('&')[-1]))
    if re.search(change, line):
        l = line.split('* ')[-1].strip(',')
        l = re.sub(' *"', '"', l)
        l = re.sub(' *, *', ', ', l)
        l = re.sub(' *= *', ' = ', l)
        print('  %s' % l)
    if re.search(slash, line.strip()):
        print('/\n')
