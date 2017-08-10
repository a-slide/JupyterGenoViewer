# -*- coding: utf-8 -*-

# Define self package variable
__version__ = '1.0a2'
__all__ = ["JVG", "JGV_helper_fun", "JGV_Alignment", "JGV_Annotation", "JGV_Level", "JGV_Reference"]
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
__install_requires__ = ['numpy>=1.11.1', 'pandas>=0.18.1', 'matplotlib>=1.5.1', 'pysam>= 0.9.0', 'notebook>=4.0.0', 'pycl>=1.0.3']
__dependency_links__ = ['https://github.com/a-slide/pycl/archive/1.0.3.tar.gz#egg=pycl-1.0.3']
__package_data__ =  ['data/yeast.bam', 'data/yeast.fa.gz', "data/yeast.gtf.gz", "data/yeast.gff3.gz", "JGV_Test_Notebook.ipynb"]
__python_requires__='>=3'

__description__='JGV is an embed genomic viewer for Jupyter notebook written in python3'
__long_description__="""JGV is a lightweight genomic viewer, taking advantage of maplotlib python library to generate annotation and
sequencing coverage plots. The genomic interval plotting method is higly customizable and allow users to analyse their results in a jupyter
notebook directly. The package can parse a variety of standard annotation file (bed, gff3, gtf...) and compute the sequencing coverage
from SAM or BAM files as well as from bed coverage files"""

# Collect info in a dictionnary for setup.py
setup_dict = {
    "name": __name__,
    "version": __version__,
    "description": __description__,
    "long_description": __long_description__,
    "url": __url__,
    "author": __author__,
    "author_email": __email__,
    "license": __licence__,
    "classifiers": __classifiers__,
    "install_requires": __install_requires__,
    "packages": [__name__],
    "package_dir": {__name__: __name__},
    "package_data": {__name__: __package_data__},
    "dependency_links": __dependency_links__,
    "python_requires": __python_requires__,
    }
