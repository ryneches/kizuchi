__author__    = "Russell Neches"
__copyright__ = "Copyright 2022, Russell Neches"
__email__     = "russell@vort.org"
__license__   = "MIT"


from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell( stdout=False, stderr=True )
extra = snakemake.params.get( 'extra', '' )

gff = snakemake.output.get( 'gff' )
if not gff :
    gff = '/dev/null'

fna = snakemake.output.get( 'fna' )
if fna :
    fna = '-d ' + fna
else :
    fna = ''

shell(
    "prodigal-gv "
    "-m "
    "-p meta "
    "-f gff "
    "-i {snakemake.input.fna} "
    "-o {gff} "
    "-a {snakemake.output.faa} "
    "{fna} "
    "{extra} "
    "2> {snakemake.log} "
)
