#!/usr/bin/env python

"""USAGE:

OPTIONS:
    -h, -?, --help: obtain help

TODO:
"""
import sys, getopt
from pygsl.sf import debye_1
# checks options
def checkopt():
    flags={}
    sopt='h?'
    lopt=['help']
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], sopt, lopt)
    except getopt.GetoptError:
        sys.stderr.write(__doc__)
    for opt, optarg in opts:
        if opt in ('-h', '-?', '--help'):
            sys.stdout.write(__doc__)
            sys.exit()
    return opts, args, flags

opts, args, flags=checkopt()

if __name__=='__main__':
    # Import Psyco if available
    try:
        import psyco
        psyco.full()
    except ImportError:
        pass

    print "# 320/x [K] debye(x) []"
    for T in range(1,10000):
        print "%s\t%s" % (0.1*T, debye_1(3200./T)[0])
    # ...your code here...
