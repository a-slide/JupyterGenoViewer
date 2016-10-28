# -*- coding: utf-8 -*-

"""
  JGV_Annotation.py
  JGV is a Python3 package for an embed genomic viewer in Jupyter notebook
  Do not try to import the package in a non-interactive environment

  Copyright 2016 Adrien Leger <aleg@ebi.ac.ul>
  [Github](https://github.com/a-slide)
  
  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or(at your option) any later version
  
  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
  (http://www.gnu.org/licenses/gpl-3.0.html).
  
  You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
"""

# Strandard library imports
from collections import OrderedDict, namedtuple, Counter
from os import access, R_OK
import csv

# Third party import
import pysam
import pandas as pd

# Local lib import
from JGV_helper_fun import *

#~~~~~~~CLASS~~~~~~~#

class Annotation(object):

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#
    
    def __init__ (self, fp, name=None, verbose=False):
        """
         * fp
            A standard gff3 file (http://www.ensembl.org/info/website/upload/gff3.html) or gtf
            (http://www.ensembl.org/info/website/upload/gff.html)containing features annotations. Could be uncompressed or archived in gz
            format. Ideally the  file would be already indexed with tabix bgzip. If not the program will sort the features and index the 
            file (can take time)
        *  name
            Name of the data file that will be used as track name for plotting.
        * verbose
            If True, will print more information during initialisation and calls of all the object methods.
        """
        #Save self variable
        self.name= name if name else file_basename(fp)
        self.verbose=verbose
        
        # Verify that the file is readable
        if not access(fp, R_OK):
            raise IOError ("{} is not readable".format(fp))
        
        # Define file format and attributes field parser 
        if extensions(fp)[0] == "gtf":
            self.format = "gtf"
            self.get_ID = self._get_gtf_ID
        elif extensions(fp)[0] == "gff3":
            self.format = "gff3"
            self.get_ID = self._get_gff3_ID
        else:
            raise ValueError ("The file is not in gtf or gff3 format (.gff3.gz, .gff3, .gtf.gz, .gtf). Please provide a correctly formated file")
        
        # Save the file path list
        self.fp = fp
        
        # If not indexed, sort and compress and index the original file with tabix 
        if not access(fp+".tbi", R_OK):
           
            # Import in a panda dataframe, remove empty rows, and convert coordinates in integer
            if self.verbose:  print("Indexing file with tabix\n\tImport annotation file and clean data")
            df = pd.read_csv(fp, names=["seqname","source","feature","start","stop","score","strand","frame","attribute"], sep="\t")
            df.dropna(inplace=True)
            df[['start', 'stop']] = df[['start', 'stop']].astype(int)
            
            # Sort the dataframe
            if self.verbose: print("\tSort lines by coordinates")
            df.sort_values(by=["seqname","start","stop"], inplace=True)
            
            # Remove the extension, name the output file and write in file
            if self.verbose: print("\tWrite a new sorted annotation file")
            temp_file = "{}/{}_sorted.{}".format(dir_path(fp), file_basename(fp),extensions(fp)[0])
            df.to_csv(temp_file, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
            
            # Compress and index the sorted file with tabix
            if self.verbose: print("\tCompress and index with tabix")
            self.fp = pysam.tabix_index(temp_file, preset="gff", force=True)
            
    def __str__(self):
        """readable description of the object"""
        msg = "{} instance\n".format(self.__class__.__name__) 
        msg+= "\tParameters list\n"
        # list all values in object dict in alphabetical order
        sorted_d = OrderedDict(sorted(self.__dict__.items(), key=lambda t: t[0]))
        for k,v in sorted_d.items():
            msg+="\t{}\t{}\n".format(k, v)
        return (msg)
        
    #~~~~~~~PROPERTY METHODS~~~~~~~#
        
    @property
    def seqid_count(self):
        """List of all the sequence ids found in the annotation file"""
        with pysam.TabixFile(self.fp, parser=pysam.asGTF()) as tbf:
            c = Counter()
            for line in tbf.fetch():
                c[line.contig]+=1
        df = pd.DataFrame.from_dict(c, orient='index', dtype=int)
        df.columns = ['count']
        df.sort_values(by="count", inplace=True, ascending=False)
        return df

    @property
    def feature_type_count(self):
        """List of all the feature types found in the annotation file"""
        with pysam.TabixFile(self.fp, parser=pysam.asGTF()) as tbf:
            c = Counter()
            for line in tbf.fetch():
                c[line.feature]+=1
        df = pd.DataFrame.from_dict(c, orient='index', dtype=int)
        df.columns = ['count']
        df.sort_values(by="count", inplace=True, ascending=False)
        return df
    
    #~~~~~~~PRIVATE METHODS~~~~~~~#
    
    def _get_gtf_ID (self, line):
        """
        Parse a gtf line and extract the feature ID corresponding to the feature type of the line, if possible.
        If not found, then an empty string will be returned.
        """
        if line.feature == "exon":
            id_field = "exon_id"
        elif line.feature == "CDS":
            id_field = "cdsid"
        elif line.feature == "transcript":
            id_field = "transcript_id"
        elif line.feature == "gene":
            id_field = "gene_id"
        else:
            return ""
        
        for a in line.attributes.strip().split(";"):           
            if a:
                l = a.strip().split(" ")
                if len(l) == 2:
                    if l[0] ==  id_field:
                        return l[1][1:-1]
        return ""

    def _get_gff3_ID (self, line):
        """
        Parse a gff3 line and extract the feature ID.
        If not found, then an empty string will be returned.
        """
        for a in line.attributes.strip().split(";"):
            if a:
                l = a.strip().split("=")
                if len(l) == 2:
                    if l[0] == "ID":
                        return l[1]
        return ""

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def get_features (self, seqid, start=None, end=None, feature_types=[]):
        """
        Parse the annotation file for the given seqid and interval and return a feature level dictionary containing a list of NamedTuples
        for each original line in the gff or gtf file. Each features are identified by their ID for gff3 file. For gtf file the ID is given
        only for exon, cds, transcript or genes if found.
        * seqid
            Name of the sequence from the initial fasta file to display
        * start
            Start of the window to display. If not given will start from 0 [ DEFAULT: None ]
        * end
            End of the window to display. If not given will start from end of the sequence [ DEFAULT: None ]
        * feature_types
            List of features types for which a track will be displayed if at least 1 feature of this type was found in the requested
            interval ( "exon"|"transcript"|"gene"|"CDS"...). If not given, all features type found in the interval will be displayed
            [ DEFAULT: [] ]
        """
        # Init dict to collect data and create a custom named tuple pattern
        feature_dict = OrderedDict()
        feature = namedtuple('feature', ["ID","start","end","strand","frame"])
 
        # Iterate over the indexed file containing the features
        with pysam.TabixFile(self.fp, parser=pysam.asGTF()) as f:
            
            for line in f.fetch(seqid, start, end, parser=pysam.asGTF()):
                if not feature_types or line.feature in feature_types:
                    
                    if not line.feature in feature_dict:
                        feature_dict[line.feature] = []
                    
                    feature_dict[line.feature].append(feature(
                            self.get_ID (line),
                            line.start,
                            line.end,
                            line.strand,
                            line.frame))
                    
            if self.verbose:
                for k,v in feature_dict.items():
                    print ("{}: {}".format(k, len(v)))
            
            return feature_dict
            
