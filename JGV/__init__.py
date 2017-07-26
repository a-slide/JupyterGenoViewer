# -*- coding: utf-8 -*-

# Define self package variable
__version__ = '1.0a1'
__all__ = ["JVG", "JGV_helper_fun", "JGV_Alignment", "JGV_Annotation", "JGV_Level", "JGV_Reference"]
__doc__='JGV is an embed genomic viewer for Jupyter notebook written in python3'
__author__= 'Adrien Leger'
__email__ = 'aleg@ebi.ac.uk'
__url__ = "https://github.com/a-slide/JupyterGenoViewer"
__licence__ = 'GPLv3'
__classifiers__ = [
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',]
__install_requires__ = ['numpy>=1.11.1', 'pandas>=0.18.1', 'matplotlib>=1.5.1', 'pysam>= 0.9.0', 'notebook>=4.0.0']
__package_data__ =  ['data/yeast.bam', 'data/yeast.fa.gz', "data/yeast.gtf.gz", "data/yeast.gff3.gz", "JGV_Test_Notebook.ipynb"]

# Collect info in a dictionnary for setup.py
setup_dict = {}

setup_dict = {
    "name": __name__,
    "version": __version__,
    "description": __doc__,
    "url": __url__,
    "author": __author__,
    "author_email": __email__,
    "license": __licence__,
    "classifiers": __classifiers__,
    "install_requires": __install_requires__,
    "packages": [__name__],
    "package_dir": {__name__: __name__},
    "package_data": {__name__: __package_data__},
    }
