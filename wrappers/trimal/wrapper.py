__author__    = 'Russell Neches'
__copyright__ = 'Copyright 2022, Russell Neches'
__email__     = 'russell@vort.org'
__license__   = 'MIT'


from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell( stdout=True, stderr=True )
extra = snakemake.params.get( 'extra', '' )

empty_output_warning = 'WARNING: Output alignment has not been generated. It is empty.'

# if there are fewer than two records in the input file, write an
# empty output file and exit.
if open( snakemake.input[0] ).read().count( '>' ) < 3 :
    
   with open( snakemake.log[0], 'w' ) as f :
      f.write( 'Too few records to perform an alignment.' )
   with open( snakemake.output[0], 'w' ) as f :
      f.write( '' )

else :

    shell(
        'trimal '
        ' -in {snakemake.input}'
        ' -out {snakemake.output}'
        ' {extra}'
        ' {log}'
    )

    if empty_output_warning in open( snakemake.log[0] ).read() :
       with open( snakemake.log[0], 'w' ) as f :
          f.write( '\n Snakemake wrapper : writing empty file anyway.' )
       with open( snakemake.output[0], 'w' ) as f :
          f.write( '' )
