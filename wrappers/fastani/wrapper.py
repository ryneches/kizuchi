__author__    = 'Russell Neches'
__copyright__ = 'Copyright 2023, Russell Neches'
__email__     = 'russell@vort.org'
__license__   = 'MIT'


from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell( stdout=True, stderr=True )
extra = snakemake.params.get( 'extra', '' )

shell(
    'fastANI '
    '{extra} '
    '-t {snakemake.threads} '
    '--rl {snakemake.input.reference} '
    '--ql {snakemake.input.query} '
    '-o {snakemake.output[0]} '
    '{log} '
)
