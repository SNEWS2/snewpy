# -*- coding: utf-8 -*-
"""
snewpy.scripts.to_snowglobes
============================
"""

from argparse import ArgumentParser

def main(options=None):
    p = ArgumentParser(description='Convert to SNOwGLoBES format.')
    p.add_argument('infile', nargs=1,
                   help='Supernova model input file.')

    if options is None:
        args = p.parse_args()
    else:
        args = p.parse_args(options)
