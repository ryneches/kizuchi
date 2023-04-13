__author__ = "Nikos Tsardakas Renhuldt"
__copyright__ = "Copyright 2021, Nikos Tsardakas Renhuldt"
__email__ = "nikos.tsardakas_renhuldt@tbiokem.lth.se"
__license__ = "MIT"


from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell( stdout=True, stderr=True )
extra = snakemake.params.get( 'extra', '' )

# if there are fewer than two records in the input file, write an
# empty output file and exit.
if open( snakemake.input[0] ).read().count( '>' ) < 3 :
    
   with open( snakemake.log[0], 'w' ) as f :
      f.write( 'Too few records to perform an alignment.' )
   with open( snakemake.output[0], 'w' ) as f :
      f.write( '' )

else :

    shell(
        'fasttree'
        ' {extra}'
        ' -out {snakemake.output.tree}'
        ' -log {snakemake.log[0]}'
        ' {snakemake.input.alignment}'
    )
