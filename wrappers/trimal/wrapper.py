__author__ = 'Russell Neches'
__copyright__ = 'Copyright 2022, Russell Neches'
__email__ = 'russel@vort.org'
__license__ = 'MIT'


from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get( 'extra', '' )

# if there are fewer than two records in the input file, write an
# empty output file and exit.
if open( snakemake.input[0] ).read().count( '>' ) < 3 :
    
   # with open( snakemake.log, 'w' ) as f :
   #     f.write( 'Too few records to perform an alignment.' )
    
    with open( snakemake.output[0], 'w' ) as f :
        f.write( '' )

else :

    shell(
        'trimal '
        '-in {snakemake.input} '
        '-out {snakemake.output} '
        '{extra}'
    )
