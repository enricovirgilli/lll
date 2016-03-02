#!/usr/bin/env python

"""USAGE:

OPTIONS:
    -h, -?, --help: obtain help

TODO:
"""
import sys, getopt

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
    # ...your code here...
    INFILE=file('spibgsmoothed.dat')
    OUFILE=file('spibg.dat', 'w')
    data=[line for line in INFILE]
    data=[(float(x),float(y)) for x,y in [line.split() for line in data]]

    for i in range(20):
        OUFILE.write("1.,\n")
    keV=20
    for e, bg in data:
        if e==keV:
            OUFILE.write("%s,\n" % bg)
            keV+=1
