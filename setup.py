from setuptools import setup

setup(
    name='cobra-meta',
    version='1.3.0',
    py_modules=['cobra'],
    install_requires=[
        'biopython',
        'pysam'
    ],
    author='LINXING CHEN',
    description='COBRA (Contig Overlap Based Re-Assembly) is a bioinformatics tool to get higher quality viral genomes assembled from metagenomes of short paired-end reads. COBRA was written in Python. COBRA has so far only been tested on assembled contigs from metaSPAdes, IDBA_UD, and MEGAHIT.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/linxingchen/cobra',
    entry_points={
        'console_scripts': [
            'cobra-meta=cobra:main',
        ],
    },
)
