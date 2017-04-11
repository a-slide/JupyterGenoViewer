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
    """
    Parse data from a file containing genomic annotation in GFF3, GTF or BED format.
    Can return the list of annotations for a given interval
    """

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, fp, name=None, verbose=False, ref_list=[]):
        """
         * fp
            An URL to a standard genomic file containing features annotations among the following format:
              gff3: http://www.ensembl.org/info/website/upload/gff3.html
              gtf:  http://www.ensembl.org/info/website/upload/gff.html
              bed:  http://www.ensembl.org/info/website/upload/bed.html
            Valid URL schemes include http, ftp, s3, and file.
            The file can eventually be compressed in ‘gzip’, ‘bz2’, ‘zip’ or ‘xz’
        *  name
            Name of the data file that will be used as track name for plotting. If not given, will be deduced from fp
            file name  [ DEFAULT: None ]
        * verbose
            If True, will print more information during initialisation and calls of all the object methods
            [ DEFAULT: False ]
        * ref_list
            list of reference sequence id to select from the data file, by default all [ DEFAULT: [] ]
        """

        # Verify that the file is readable
        assert access(fp, R_OK), "{} is not readable".format(fp)

        #Save self variable
        self.fp = fp
        self.name = name if name else file_basename(fp)
        self.verbose = verbose
        self.ext = extensions(fp)[0]

        if self.ext == "gtf":
            if self.verbose: print("Use GTF parser to parse annotations in", self.name)
            self.df = self._gtf_parser(fp, ref_list)
        elif self.ext == "gff3":
            if self.verbose: print("Use GFF3 parser to parse annotations in ", self.name)
            self.df = self._gff3_parser(fp, ref_list)
        elif self.ext == "bed":
            if self.verbose: print("Use BED parser to parse annotations in ", self.name)
            self.df = self._bed_parser(fp, ref_list)
        else:
            msg = "The file is not in gtf/gff3/bed format (.gff3/.gtg/.bed +-.gz). Please provide a correctly formated file"
            raise ValueError(msg)

        # Sort the dataframe
        #if self.verbose: print("\tSort annotation features by coordinates".format(self.ext))
        #self.df.sort_values(by=["refid","start","end"], inplace=True)

        if self.verbose:
            print ("\tFound {} features in {} reference sequences".format( self.feature_count, self.refid_count))

    def __str__(self):
        """readable description of the object"""
        msg = "{} instance\n".format(self.__class__.__name__)
        # list all values in object dict in alphabetical order
        for k,v in OrderedDict(sorted(self.__dict__.items(), key=lambda t: t[0])).items():
            if k != "df":
                msg+="\t{}\t{}\n".format(k, v)
        return (msg)

    def __repr__ (self):
        return ("{}-{} / Feature count {}".format(self.__class__.__name__, self.name, self.feature_count))

    #~~~~~~~PROPERTY METHODS~~~~~~~#

    @property
    def feature_count(self):
        """Number of features collected"""
        return len(self.df)

    @property
    def refid_count(self):
        """Number of unique reference sequence ids found"""
        return self.df["refid"].nunique()

    @property
    def type_count(self):
        """Number of unique feature type found"""
        return self.df["type"].nunique()

    @property
    def refid_list(self):
        """List of unique reference sequence ids found"""
        return self.df["refid"].unique().tolist()

    @property
    def type_list(self):
        """List of unique feature type found"""
        return self.df["type"].unique().tolist()

    @property
    def refid_count_uniq(self):
        """List of unique reference sequence ids with count of associated features"""
        return pd.DataFrame(self.df.groupby("refid").size().sort_values(ascending=False), columns=["count"])

    @property
    def type_count_uniq(self):
        """List of unique feature types with count of associated features"""
        return pd.DataFrame(self.df.groupby("type").size().sort_values(ascending=False), columns=["count"])

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _bed_parser(self, fp, ref_list=[]):
        """Parse a bed formated file"""
        # Import the file in a dataframe
        col_names = ["refid","start","end","ID","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]
        df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#")
        # Select only features from the reference list
        if ref_list:
            df = df[(df["refid"].isin(ref_list))]
        # Cast the start and end field in integer
        df[['start', 'end']] = df[['start', 'end']].astype(int)
        # Extract the ID field = first field of attribute
        df['type'] = "unknown"
        # Drop unecessary fields and reorganise order
        df = df[["refid","start","end","strand","ID", "type"]]
        return df

    def _gff3_parser(self, fp, ref_list=[]):
        """Parse a gff3 formated file"""
        # Import the file in a dataframe
        col_names = ["refid","source","type","start","end","score","strand","frame","attribute"]
        df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#")
        # Select only features from the reference list
        if ref_list:
            df = df[(df["refid"].isin(ref_list))]
        # Cast the start and end field in integer
        df[['start', 'end']] = df[['start', 'end']].astype(int)
        # Extract the ID field = first field of attribute
        df['ID'] = df["attribute"].str.split(';').str[0].str[3:]
        # Drop unecessary fields and reorganise order
        df = df[["refid","start","end","strand","ID","type"]]
        return df

    def _gtf_parser(self, fp, ref_list=[]):
        """Parse a gtf formated file"""
        # Import the file in a dataframe
        col_names = ["refid","source","type","start","end","score","strand","frame","attribute"]
        df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#")
        # Select only features from the reference list
        if ref_list:
            df = df[(df["refid"].isin(ref_list))]
        # Cast the start and end field in integer
        df[['start', 'end']] = df[['start', 'end']].astype(int)
        # Extract the ID field = first field of attribute=
        df['ID'] = df["attribute"].str.split('\"').str[1]
        # Drop unecessary fields and reorganise order
        df = df[["refid","start","end","strand","ID","type"]]
        return df

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def interval_features (self, refid, start, end, feature_types=None, max_features_per_type=None):
        """
        Parse the annotation file for the given refid and interval and return a dataframe containing all the features
        found for each original line. Features are identified by their ID field for gff3 files, by the entire
        attribute field for the bed files and by the first element in the attribute field for the gtf files
        * refid
            Name of the sequence from the original fasta file to display
        * start
            Start of the window to display. The coordinate is not verified, if outside of the range it will
            return an empty dataframe
        * end
            End of the window to display. The coordinate is not verified, if outside of the range it will
            return an empty dataframe
        * feature_types
            Name of a valid feature type ( "exon"|"transcript"|"gene"|"CDS"...) or list of names of feature type for
            which a row will be returned. The option is not available for bed files. If not given, all features type
            found in the interval will be returned [ DEFAULT: None ]
        * max_features_per_type
            Maximal total number of features for a particular feature type. If more are found, a random sampling will
            be performed. If None, all the features will be returned [ DEFAULT: None ]
        """
        # Verifications and auto adjustment of coordinates
        if not refid in self.refid_list:
            if self.verbose: print ("The reference {} is not in the list of references with alignment".format(refid))
            return pd.DataFrame(columns=["refid","start","end","strand","ID","type"])
        if feature_types and self.ext == "bed":
            if self.verbose: print ("Incompatible options. Bed files do not allow to identify feature_type")
            return pd.DataFrame(columns=["refid","start","end","strand","ID","type"])

        # Select the refid and coordinates
        df = self.df[(self.df["refid"] == refid)&(self.df["end"] > start)&(self.df["start"] < end)]
        if df.empty:
            if self.verbose: print("No feature found in the requested interval")
            return pd.DataFrame(columns=["refid","start","end","strand","ID","type"])

        # Cast str to list
        if type(feature_types) == str: feature_types = [feature_types]

        # Filter_df by type and max number per type
        select_list = []
        for type_name, type_df in df.groupby("type"):
            # Filter out if not in the list
            if not feature_types or type_name in feature_types:
                sdf = df[(df["type"] == type_name)]
                if sdf.empty:
                    if self.verbose: print("No feature of type {} found in the requested interval".format(type_name))
                elif max_features_per_type and len(sdf)>max_features_per_type:
                    select_list.append(sdf.sample(max_features_per_type))
                else:
                    select_list.append(sdf)
        # Merge the selected features in a single df
        if select_list:
            df = pd.concat(select_list)
        else:
            if self.verbose: print("No feature found in the requested interval")
            return pd.DataFrame(columns=["refid","start","end","strand","ID","type"])

        # Return a sorted copy of the  df
        df.sort_values(by=["refid","start","end"], inplace=True)
        return df.copy().reset_index(drop=True)
