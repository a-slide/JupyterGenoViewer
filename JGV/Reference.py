# -*- coding: utf-8 -*-

"""
  Reference.py
  JGV is a Python3 package for an embed genomic viewer in Jupyter notebook. Do not import the package in a
  non-interactive environment

  Copyright 2016 Adrien Leger <aleg@ebi.ac.ul>
  [Github](https://github.com/a-slide)

  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
  License as published by the Free Software Foundation; either version 3 of the License, or(at your option) any later
  version

  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
  (http://www.gnu.org/licenses/gpl-3.0.html).

  You should have received a copy of the GNU General Public License along with this program; if not, write to the Free
  Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
"""

# Strandard library imports
from collections import OrderedDict, Counter, namedtuple
from os import access, R_OK
import gzip
import csv

# Third party import
import pandas as pd

# Local lib import
from JGV_helper_fun import extensions, file_basename, dir_path
from JGV_helper_fun import jprint as print

#~~~~~~~CLASS~~~~~~~#
class Reference(object):
    """
    Parse a fasta reference file or a fasta index and save the list of reference sequences ids and their lengths
    """
    version = "0.0.1"

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#
    def __init__ (self, fp, name=None, verbose=False, ref_list=[], output_index=False):
        """
         * fp
            A fasta file containing the reference sequences OR an tab separated index file containing at least 2 columns
            with the refid and the length in bases (like a .fa.fai file generated by samtools faidx, or with the
            output_index option of this function)
            The fasta option will take more time as the file has to be parsed to get the refid and length of sequences.
            Both fasta and infex file can be gziped
        *  name
            Name of the data file that will be used as track name for plotting. If not given, will be deduced from fp
            file name
        * verbose
            If True, will print more information during initialisation and calls of all the object methods
            [ DEFAULT: False ]
        * ref_list
            list of reference sequence id to select from the data file, by default all [ DEFAULT: [] ]
        * output_index
            If True will write a simple A 2 column index tsv file containing the Reference sequence ids and their
            lengths [ DEFAULT: False ]
        """
        # Verify that the file is readable
        assert access(fp, R_OK), "{} is not readable".format(fp)

        #Save self variable
        self.fp = fp
        self.name = name if name else file_basename(fp)
        self.verbose = verbose
        self.ext = extensions(fp)[0]

        # If the file is in fasta format
        if self.ext in ["fa", "fasta"]:
            if self.verbose: print ("Parsing fasta file")

            # File handling for both uncompressed or compressed fasta file
            if fp.endswith(".gz"):
                open_fun, open_mode = gzip.open, "rt"
            else:
                open_fun, open_mode = open, "r"

            # Parse fasta file refid and count the length of each sequence if in the ref_list
            with open_fun(fp, open_mode) as f:
                d = OrderedDict()
                last_ref = None
                for l in f:
                    if l.startswith(">"):
                        refid = l[1:].split()[0].strip()
                        if not ref_list or refid in ref_list:
                            d[refid] = 0
                            last_ref = refid
                        else:
                            last_ref = None
                    elif last_ref:
                        d[last_ref] += len(l.strip())

            # Transform the counter in a Dataframe and sort by length
            self.d = pd.Series(d, name="length", dtype = "int64")
            self.d.sort_values(inplace=True, ascending=False)

            # Write the index in a file for quicker loading next time
            if output_index:
                index_file = "{}/{}.tsv".format(dir_path(fp), file_basename(fp))
                if self.verbose: print ("Write a fasta index file: {}".format(index_file))
                self.d.to_csv(index_file, sep="\t")

        # In the case the file is not in fasta format, try to parse it as a 2 columns tabulated file with refid and length for each sequence
        else:
            if self.verbose: print ("Assume the file is a fasta index")
            self.d = pd.read_csv(fp, sep="\t", squeeze=True, comment="#",  usecols=[0,1], index_col=0, header=None)
            if ref_list: self.d = self.d[(self.d.index.isin(ref_list))]
            self.d.name= "length"
            self.d.sort_values(inplace=True, ascending=False)

        if self.verbose:
            print ("\tFound {} reference sequences".format(self.refid_count))

    def __str__(self):
        """readable description of the object"""
        msg = "{} instance\n".format(self.__class__.__name__)
        # list all values in object dict in alphabetical order
        for k,v in OrderedDict(sorted(self.__dict__.items(), key=lambda t: t[0])).items():
            if k == "d":
                for refid, length in v.items():
                    msg+="\t\t{}\tlength: {}\n".format(refid, length)
            else:
                msg+="\t{}\t{}\n".format(k, v)
        return (msg)

    def __repr__ (self):
        return ("{}-{} / Refid count {}".format(self.__class__.__name__, self.name, self.refid_count))

    #~~~~~~~PROPERTY METHODS~~~~~~~#
    @property
    def refid_list(self):
        """List of all the sequence ids found in the annotation file"""
        return list(self.d.index)

    @property
    def refid_count(self):
        """Number of unique reference sequence ids found"""
        return len(self.d)

    #~~~~~~~PUBLIC METHODS~~~~~~~#
    def get_refid_len (self, refid):
        """ Return the length of a given refid, If the reference is not found return None"""
        if refid not in self.d:
            if self.verbose: print("The reference sequence {} was not found in the reference list".format(refid))
            return None
        else:
            return self.d[refid]
