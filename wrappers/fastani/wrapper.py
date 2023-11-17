__author__    = 'Russell Neches'
__copyright__ = 'Copyright 2023, Russell Neches'
__email__     = 'russell@vort.org'
__license__   = 'MIT'


from snakemake.shell import shell
from collections import defaultdict
import tempfile
import os
import pysam

LOG = snakemake.log_fmt_shell( stdout=True, stderr=True )
log = snakemake.params.prelog
extra = snakemake.params.get( 'extra', '' )

with tempfile.TemporaryDirectory() as tmpdir :
    
    open( log, 'w' ).write( 'created temp dir : {t}\n'.format( t=tmpdir ) )
    
    reflist = os.path.join( tmpdir, 'references.txt' )
    qrylist = os.path.join( tmpdir, 'queries.txt'   )
    
    R = snakemake.input.fastas[:len(snakemake.input.fastas)//2]
    Q = snakemake.input.fastas[len(snakemake.input.fastas)//2:]
    
    open( log, 'a' ).write( 'references : {r}\n'.format( r=str(','.join(R)) ) )
    open( log, 'a' ).write( 'queries    : {q}\n'.format( q=str(','.join(Q)) ) )
    
    for contigs,batches,glist,group in zip( [ defaultdict(list), defaultdict(list) ],
                                            [ R, Q ],
                                            [ reflist, qrylist ],
                                            [ 'reference', 'query' ] ) :
        
        open( log, 'a' ).write( 'writing group {g}...\n'.format( g=group ) )
                               
        for batch in batches :
            with pysam.FastaFile( batch ) as fasta :
                open( log, 'a' ).write( '   batch : {b}\n'.format( b=batch ) )
                for ref in fasta.references :
                    genome = ref.rsplit('_',2)[0]
                    with open( os.path.join( tmpdir, genome + '.fna' ), 'a' ) as f :
                        n = len( contigs[ genome ] )
                        contigs[ genome ].append( ref )
                        f.write( '>{contig}\n{data}\n'.format( contig=genome + '_contig_' + str(n),
                                                               data=fasta.fetch( ref ) ) )
        
        open( log, 'a' ).write( '{n} {g} genomes\n'.format( n=str(len(contigs)),
                                                            g=group ) )
        open( log, 'a' ).write( '{n} {g} contigs\n'.format( n=str( sum( [ len(c) for c in contigs.values() ] ) ),
                                                            g=group ) )
        with open( glist, 'w' ) as f :
            f.write( '\n'.join( [ os.path.join( tmpdir, g + '.fna' ) for g in contigs.keys() ] ) )

        open( log, 'a' ).write( 'wrote {l}\n'.format( l=glist ) )
    
    open( log, 'a' ).write( 'pre-ANI staging complete, starting fastANI\n\n\n' )
    
    shell(
        'fastANI '
        '{extra} '
        '-t {snakemake.threads} '
        '--rl {reflist} '
        '--ql {qrylist} '
        '-o {snakemake.output[0]} '
        '{LOG} '
    )
    
    open( log, 'a' ).write( '\n\nfastANI run complete, cleaned up temp directory.\n' )
