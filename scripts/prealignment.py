#!/usr/bin/env python

import argparse
import pandas
from Bio.SeqIO import parse
import os

parser = argparse.ArgumentParser()

parser.add_argument( '-i', '--hmmsearch',
                     type=str,
                     dest='hmmsearch',
                     help='hmmsearch output' )

parser.add_argument( '-p', '--proteins',
                     type=str,
                     dest='proteins',
                     help='FASTA file of protein sequences' )

parser.add_argument( '-o', '--outdir',
                     type=str,
                     dest='outdir',
                     help='output directory' )

args = parser.parse_args()

# borrowed from CheckV
def parse_hmmsearch( path ) :
    with open(path) as f:
        names = [
            "qname",
            "qacc",
            "tname",
            "tacc",
            "eval",
            "score",
            "bias",
            "beval",
            "bscore",
            "bbias",
        ]
        formats = [str, str, str, str, float, float, float, float, float, float]
        for line in f:
            if not line.startswith("#"):
                values = line.split()
                yield dict([(names[i], formats[i](values[i])) for i in range(10)])

hits = pandas.DataFrame( [ hits for hits in parse_hmmsearch( args.hmmsearch ) ] )

hits['genome'] = [ '_'.join(p.split('_')[:-1]) for p in hits['qname'] ]

print( hits['tname'].value_counts() )

orthologs = []
for model in set( hits['tname'] ) :
    for genome in set( hits['genome'] ) :
        mhits = hits[ ( hits['tname']  == model  ) &
                      ( hits['genome'] == genome ) ]
        if len( mhits ) == 0 : continue
        tophit = mhits.sort_values( 'bscore', ascending=False ).iloc[0]
        orthologs.append( tophit.name )

ortholog_hits = hits.loc[ orthologs ]
ortholog_names = list( ortholog_hits['qname'] )

# create output directory if it doesn't exist
if not os.path.exists( args.outdir ) :
    os.mkdir( args.outdir )

handles = {}
for n,rec in enumerate( parse( open( args.proteins ), 'fasta' ) ) :
    if not rec.id in ortholog_names : continue
    h = ortholog_hits[ ortholog_hits['qname'] == rec.id ].iloc[0]
    if not h['tname'] in handles :
        handles[ h['tname'] ] = open( 'proteins/' + 
                                      h['tname'] +
                                      '.faa', 'w' )
    handles[ h['tname'] ].write( rec.format( 'fasta' ) )

# close file handles
for key in handles.keys() :
    handles[key].close()
