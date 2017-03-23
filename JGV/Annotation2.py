# -*- coding: utf-8 -*-

"""
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
from collections import OrderedDict
from os import access, R_OK

# Third party import
import pandas as pd

# Local lib import
from JGV_helper_fun import *
from JGV_helper_fun import jprint as print

#~~~~~~~CLASS~~~~~~~#

class Annotation(object):

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, fp, name=None, verbose=False):
        """
         * fp
            An URL to a standard genomic file containing features annotations among the following format:
            gff3 file (http://www.ensembl.org/info/website/upload/gff3.html)
            gtf file  (http://www.ensembl.org/info/website/upload/gff.html)
            bed http://www.ensembl.org/info/website/upload/bed.html
            Valid URL schemes include http, ftp, s3, and file.
            The file can eventually be compressed in ‘gzip’, ‘bz2’, ‘zip’ or ‘xz’
        *  name
            Name of the data file that will be used as track name for plotting. If not given, will be deduced from fp
            file name
        * verbose
            If True, will print more information during initialisation and calls of all the object methods.
        """

        # Verify that the file is readable
        if not access(fp, R_OK):
            raise IOError ("{} is not readable".format(fp))

        #Save self variable
        self.fp = fp
        self.name = name if name else file_basename(fp)
        self.verbose = verbose
        self.ext = extensions(fp)[0]

        if self.ext == "gtf":
            if self.verbose: print("Use GTF parser to parse annotations in", self.name)
            self.df = self._gtf_parser(fp)
        elif self.ext == "gff3":
            if self.verbose: print("Use GFF3 parser to parse annotations in ", self.name)
            self.df = self._gff3_parser(fp)
        elif self.ext == "bed":
            if self.verbose: print("Use BED parser to parse annotations in ", self.name)
            self.df = self._bed_parser(fp)
        else:
            msg = "The file is not in gtf/gff3/bed format (.gff3/.gtg/.bed +-.gz). Please provide a correctly formated file"
            raise ValueError(msg)

        # Sort the dataframe
        if self.verbose: print("Sort annotation features by coordinates".format(self.ext))
        self.df.sort_values(by=["seqid","start","end"], inplace=True)

        # Verbose informations
        if self.verbose:
            print("Annotation features count:{}".format(self.feature_count))
            print("Seqid count:{}".format(self.seqid_count))
            print("Feature type count:{}".format(self.type_count))

    def __str__(self):
        """readable description of the object"""
        msg = "{} instance\n".format(self.__class__.__name__)
        # list all values in object dict in alphabetical order
        for k,v in OrderedDict(sorted(self.__dict__.items(), key=lambda t: t[0])).items():
            if k != "df":
                msg+="\t{}\t{}\n".format(k, v)
        return (msg)


    def __repr__ (self):
        return ("{} {} {} features".format(self.__class__.__name__, self.name, self.feature_count))

    #~~~~~~~PROPERTY METHODS~~~~~~~#

    @property
    def feature_count(self):
        """List of all the sequence ids found in the annotation file"""
        return len(self.df)

    @property
    def seqid_count(self):
        """Number of unique seqid found"""
        return self.df["seqid"].nunique()

    @property
    def seqid_list(self):
        """List of unique seqid found"""
        return self.df["seqid"].unique()

    @property
    def seqid_count_uniq(self):
        """List of all the sequence ids found in the annotation file"""
        return pd.DataFrame(self.df.groupby("seqid").size().sort_values(ascending=False), columns=["count"])

    @property
    def type_count(self):
        """Number of unique seqid found"""
        return self.df["type"].nunique()

    @property
    def type_list(self):
        """List of unique seqid found"""
        return self.df["type"].unique()

    @property
    def type_count_uniq(self):
        """List of all the sequence ids found in the annotation file"""
        return pd.DataFrame(self.df.groupby("type").size().sort_values(ascending=False), columns=["count"])

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _bed_parser(self, fp):
        # Import the file in a dataframe
        col_names = ["seqid","start","end","ID","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]
        df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#")
        # Cast the start and end field in integer
        df[['start', 'end']] = df[['start', 'end']].astype(int)
        # Extract the ID field = first field of attribute
        df['type'] = "unknown"
        # Drop unecessary fields and reorganise order
        df = df[["seqid","start","end","strand","ID", "type"]]
        return df

    def _gff3_parser(self, fp):
        # Import the file in a dataframe
        col_names = ["seqid","source","type","start","end","score","strand","frame","attribute"]
        df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#")
        # Cast the start and end field in integer
        df[['start', 'end']] = df[['start', 'end']].astype(int)
        # Extract the ID field = first field of attribute
        df['ID'] = df["attribute"].str.split(';').str[0].str[3:]
        # Drop unecessary fields and reorganise order
        df = df[["seqid","start","end","strand","ID","type"]]
        return df

    def _gtf_parser(self, fp):
        # Import the file in a dataframe
        col_names = ["seqid","source","type","start","end","score","strand","frame","attribute"]
        df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#")
        # Cast the start and end field in integer
        df[['start', 'end']] = df[['start', 'end']].astype(int)
        # Extract the ID field = first field of attribute=
        df['ID'] = df["attribute"].str.split('\"').str[1]
        # Drop unecessary fields and reorganise order
        df = df[["seqid","start","end","strand","ID","type"]]
        return df

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def interval_features (self, seqid, start=None, end=None, feature_type=None):
        """
        Parse the annotation file for the given seqid and interval and return a dataframe containing all the features
        found for each original line file. Features are identified by their ID field for gff3 files, by the entire
        attribute field for the bad files and by the firt element in the attribute field for the gtf files
        * seqid
            Name of the sequence from the initial fasta file to display
        * start
            Start of the window to display. If not given will start from 0 [ DEFAULT: None ]
        * end
            End of the window to display. If not given will start from end of the sequence [ DEFAULT: None ]
        * feature_types
            Name of a valid feature type ( "exon"|"transcript"|"gene"|"CDS"...) or list of names of feature type for
            which a row will be returned. The option is not available for bed files. If not given, all features type
            found in the interval will be returned [ DEFAULT: None ]
        """
        # Verify that the sequence is in the seqid list
        if seqid not in self.seqid_list:
            if self.verbose: print("Seqid ({}) not found in the list for the annotation {}".format(seqid, self.name))
            return pd.DataFrame()

        # Select the seqid
        df = self.df[(self.df["seqid"] == seqid)]

        # Select start, end and feature types if required
        if feature_type and self.ext != "bed":
            if type(feature_type) == list:
                df = df[(df["type"].isin(feature_type))]
            elif type(feature_type) == str:
                df = df[(df["type"] == feature_type)]
        if start:
            df = df[(df["end"] > start)]
        if end:
            df = df[(df["start"] < end)]

        return df.copy().reset_index(drop=True)
