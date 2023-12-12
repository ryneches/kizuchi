'''Snakemake wrapper for MAFFT.'''

__author__    = 'Russell Neches'
__copyright__ = 'Copyright 2023, Russell Neches'
__email__     = 'russell@vort.org'
__license__   = 'MIT'


from snakemake.shell import shell

# Placeholder for optional parameters
extra = snakemake.params.get( 'extra', '' )
# Formats the log redrection string
log = snakemake.log_fmt_shell( stdout=False, stderr=True )

# if there are fewer than two records in the input file, write an
# empty output file and exit.
if open( snakemake.input[0] ).read().count( '>' ) < 3 :
    
    with open( snakemake.log[0], 'w' ) as f :
        f.write( 'Too few records to perform an alignment.' ) 
    with open( snakemake.output[0], 'w' ) as f :
        f.write( '' )
    
else :
    
    # Executed shell command
    shell(
        'mafft {extra} '
        '{snakemake.input[0]} '
        '> {snakemake.output[0]} '
        '{log}'
    )
