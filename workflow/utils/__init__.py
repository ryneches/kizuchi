import os
from Bio.SeqIO import parse
from itertools import chain, islice, product
from collections import defaultdict

# itertools.batched won't be available until 3.12.
def batched( iterable, size=10 ) :
  iterator = iter(iterable)
  for first in iterator :
    yield chain( [first], islice( iterator, size - 1 ) )

# resource use projection functions : estimation of 
# empirical constants can be found in docs/Benchmarks.ipynb

def clustalo_memory( wildcards, attempt ) :
  a = 3e-5
  b = 0.4
  c = 500
  fastapath = os.path.join( 'clusters',
                            wildcards.ani,
                            wildcards.cid,
                            'genes',
                            wildcards.hmm + '.fna' )

  max_gene_length = max( [ len(rec)
                           for rec
                           in parse( open( fastapath ), 'fasta' ) ] )
  
  x = max_gene_length

  return a * x**2 + b * x + c * attempt

def trimal_memory( wildcards, attempt ) :
  a = 0
  b = 0.05
  c = 100
  fastapath = os.path.join( 'clusters',
                            wildcards.ani,
                            wildcards.cid,
                            'alignments',
                            wildcards.hmm + '_aln.fna' )
  
  n_sequences = len( [ rec
                       for rec
                       in parse( open( fastapath ), 'fasta' ) ] )
  
  x = n_sequences

  return a * x**2 + b * x + c * attempt

def fasttree_memory( wildcards, attempt ) :
  a = 0
  b = 1.8e-4
  c = 70
  fastapath = os.path.join( 'clusters',
                             wildcards.ani,
                             wildcards.cid,
                             'filtered_trimmed_alignments',
                             wildcards.hmm + '_aln.fna' )
  
  sequence_lengths = [ len(rec)
                       for rec
                       in parse( open( fastapath ), 'fasta' ) ]
  
  if len( sequence_lengths ) == 0 : return 0

  x = max( sequence_lengths ) * len( sequence_lengths )

  return a * x**2 + b * x + c * attempt

# borrowed from CheckV
def parse_hmmsearch( path ) :
  with open(path) as f:
    names = [ 'tname', 'qacc', 'qname', 'tacc',   'eval', 
              'score', 'bias', 'beval', 'bscore', 'bbias' ]
    formats = [ str,   str,   str,   str,   float, 
                float, float, float, float, float ]
    for line in f:
      if not line.startswith( '#' ) :
        values = line.split()
        yield dict( [ ( names[i], formats[i](values[i]) ) 
                      for i in range(10) ] )

def bipartition_compatibility( A, B ) :
    '''
    Test for compatibility of two bipartitions.
    
    [1] Salichos, Leonidas, Alexandros Stamatakis, and
        Antonis Rokas. "Novel information theory-based 
        measures for quantifying incongruence among 
        phylogenetic trees." Molecular biology and 
        evolution 31, no. 5 (2014): 1261-1271.
    [2] Salichos, Leonidas, and Antonis Rokas. "Inferring
        ancient divergences requires genes with strong 
        phylogenetic signals." Nature 497, no. 7449 
        (2013): 327-331.
    '''
    A0,A1 = A
    B0,B1 = B
    return ( not bool( A0 & B0 ) ) \
         | ( not bool( A1 & B0 ) ) \
         | ( not bool( A0 & B1 ) ) \
         | ( not bool( A1 & B1 ) )

def bipartition_compatibility_ratio( T1, T2, links ) :
    '''
    Compute the ratio of compatible bipartitions of connected
    leaf nodes between two gene trees. Connected leaf names 
    are mapped into a common namespace, so only bipartitions
    of the connected leaf nodes are included in the ratio.
    '''
    L1,L2 = zip(*links)
    
    B1 = frozenset( frozenset((a,b)) for a,b in 
                    [ [ frozenset( L1.index(i) for i in a if i in L1 ), 
                        frozenset( L1.index(i) for i in b if i in L1 ) ]
                      for a,b in T1.bipartitions() ]
                    if a and b )
    
    B2 = frozenset( frozenset((a,b)) for a,b in 
                    [ [ frozenset( L2.index(i) for i in a if i in L2 ), 
                        frozenset( L2.index(i) for i in b if i in L2 ) ]
                      for a,b in T2.bipartitions() ]
                    if a and b )
    
    S = [ bipartition_compatibility( b1, b2 ) for b1,b2 in product( B1, B2 ) ]
    
    return { 'ratio'   : sum(S)/len(S),
             'jaccard' : len( B1 & B2 ) / len( B1 | B2 ) }

def find_links( genes1, genes2 ) :
    '''
    For two sets of genes called from the same set of genomes, find the linkage
    relationships connecting the genes.
    '''
    t1_genomes = { name : name.rsplit('_', 3)[0] for name in genes1 }
    t2_genomes = { name : name.rsplit('_', 3)[0] for name in genes2 }
    
    t1_genes = defaultdict(list)
    for name in genes1 :
        genome = name.rsplit( '_', 3 )[0]
        t1_genes[genome].append( name )
    
    t2_genes = defaultdict(list)
    for name in genes2 :
        genome = name.rsplit( '_', 3 )[0]
        t2_genes[genome].append( name )
    
    links = []
    for genome in set( t1_genes.keys() ) & set( t2_genes.keys() ) :
        
        # if both genes are single copy a genome, it is linked
        # if either gene has duplications, each copy of each
        # gene is linked only if they are on the same contig
        
        if len( t1_genes[ genome] ) == len( t2_genes[ genome ] ) == 1 :
            
            links.append( ( t1_genes[genome][0], t2_genes[genome][0] ) )
            
        else :
            for gene1 in t1_genes[ genome ] :
                for gene2 in t2_genes[ genome ] :
                    contig1 = gene1.rsplit('_', 2)[0]
                    contig2 = gene2.rsplit('_', 2)[0]
                    if contig1 == contig2 :
                        links.append( ( gene1, gene2 ) )
    
    return links
