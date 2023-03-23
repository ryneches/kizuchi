"""Snakemake wrapper for clustal omega."""

__author__ = "Michael Hall"
__copyright__ = "Copyright 2019, Michael Hall"
__email__ = "mbhall88@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Placeholder for optional parameters
extra = snakemake.params.get( 'extra', '' )
# Formats the log redrection string
log = snakemake.log_fmt_shell( stdout=True, stderr=True )

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
        'clustalo {extra}'
        ' --threads={snakemake.threads}'
        ' --in {snakemake.input[0]}'
        ' --out {snakemake.output[0]}'
        ' --log={snakemake.log[0]} --verbose'
    )
