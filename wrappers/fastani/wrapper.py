__author__ = 'Russell Neches'
__copyright__ = 'Copyright 2023, Russell Neches'
__email__ = 'russell@vort.org'
__license__ = 'MIT'


from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get('extra', '')

shell(
    'fastANI '
    '{extra} '
    '-r {snakemake.input.reference} '
    '-q {snakemake.input.query} '
    '-o {snakemake.params.name} '
    '2> {snakemake.log} '
)
