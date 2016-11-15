# -*- coding: utf-8 -*-

"""
  JGV_Alignment.py
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
from collections import OrderedDict
from os import access, R_OK

# Third party import
import pysam
import pandas as pd

# Local lib import
from JGV_helper_fun import *
from JGV_helper_fun import jprint as print

#~~~~~~~CLASS~~~~~~~#

class Alignment(object):

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#
    
    def __init__ (self, fp, name=None, verbose=False):
        """
         * fp
             A standard bam file already SORTED BY COORDINATES (http://samtools.sourceforge.net/SAM1.pdf) containing aligned reads.
             Ideally the file would also be already indexed with samtools index. If not the program will index the file (can take time)
        *  name
            Name of the data file that will be used as track name for plotting. If not given, will be deduced from fp file name
        * verbose
            If True, will print more information during initialisation and calls of all the object methods.
        """
        #Save self variable
        self.name= name if name else file_basename(fp)
        self.verbose=verbose
        
        # Verify that the file is readable
        if not access(fp, R_OK):
            raise IOError ("{} is not readable".format(fp))

        # Verify that the file is in bam format
        if extensions(fp)[0] != "bam":
            raise ValueError ("The file is not in BAM format (.bam). Please provide a correctly formated file")
        
        with pysam.AlignmentFile(fp) as bam:
            # If not bam index available
            if not bam.has_index():
                # Verify if file is sorted
                if bam.header["HD"]['SO'] == "coordinate":
                    if self.verbose: print ("Indexing bam file with samtools index")
                    pysam.index(fp, catch_stdout=False)
                else:
                    raise UserWarning ("The bam file needs to be sorted by coordinates to be processed")
        
         # The index should be available at this stage
        with pysam.AlignmentFile(fp) as bam:
            assert bam.has_index(), "Bam file is not indexed ?"
            
            # save the list of the sequences found in the file
            self.seqid_list = bam.references
            self.n_seq = len(self.seqid_list)


        # Save the file path list
        self.fp = fp
            
    def __str__(self):
        """readable description of the object"""
        msg = "{} instance\n".format(self.__class__.__name__) 
        msg+= "\tParameters list\n"
        # list all values in object dict in alphabetical order
        sorted_d = OrderedDict(sorted(self.__dict__.items(), key=lambda t: t[0]))
        for k,v in sorted_d.items():
            msg+="\t{}\t{}\n".format(k, v)
        return (msg)
        
    def __repr__ (self):
        return ("{}-{}".format(self.__class__.__name__, self.name))
        
        
    #~~~~~~~PUBLIC METHODS~~~~~~~#
    
    def seqid_read_count (self):
        """
        Return a dataframe containing 3 colums ("length", "read_count", "read_density") for all seq_id reference sequence in the BAM file.
        The read count is computationally intensive and can take time for large BAM files
        """       
        with pysam.AlignmentFile(self.fp, "rb") as bam:
            df = pd.DataFrame(columns=["length", "read_count", "rpkm"], index=bam.references)
            total_read = 0
            for seqid, length in zip(bam.references, bam.lengths):
                read_count = bam.count(seqid)
                df.loc[seqid]["length"] = length
                df.loc[seqid]["read_count"] = read_count
                total_read += read_count
            
        for seqid, val in df.iterrows():
            df.loc[seqid]["rpkm"] = val.read_count / (val.length/1000) / (total_read/1000000)
        
        return df

    def interval_coverage (self, seqid, start=None, end=None, n_step=500, mode="auto"):
        """
        Parse the alignment file for a given seqid and interval. Compute the coverage at read level in mode "read_count" or at base level
        in mode "base_coverage". In both case the interval is splited in a number of windows equal to n_step, for which the coverage in computed.
        The method return a tuple of 2 synchronised lists containing the starting positions of the windows and the coverage.
        * seqid
            Name of the sequence from the initial fasta file to display
        * start
            Start of the window to display. If not given will start from 0 [ DEFAULT: None ]
        * end
            End of the window to display. If not given will start from end of the sequence [ DEFAULT: None ]
       * n_step
            Number of alignment count bins to divide the displayed window. Low number will result in low resolution but will be faster. The
            default value is a decent compromise between speed and precision [ DEFAULT: 500 ]
        * mode
            Mode to parse and represent the alignment track ("auto"|"read_count"|"base_coverage"). The "read_count" is fast for large
            regions, but do not account for spliced reads. In "auto" mode it is automatically selected for intervals larger than 5000000.
            The "base_coverage" mode compute the median base coverage per step. It is faster for low coverage and small intervals.
            In "auto" mode it is automatically selected for intervals larger smaller than 5000000 [ DEFAULT: "auto" ]
        """
        
        # Init lists to collect data
        pos_list = []
        cov_list = []
        
        #Auto define start and stop if not given
        if not start:
            start = 0
            if self.verbose: print ("Autodefine start position: {}".format(start))
        if not end:
            end = self._len_seqid(seqid)
            if self.verbose: print ("Autodefine end position: {}".format(end))
        
        # Autocorect the number of steps if needed to be always superior or equal to the requested interval
        len_range = end-start
        if n_step > len_range:
            if self.verbose: warnings.warn("Short interval = autocorect n_step to {}".format(n_step))
            n_step = len_range
            
        step_size = (end-start) / float(n_step)

        # Read_count method, quicker for long intervals
        if mode == "read_count" or (mode == "auto" and len_range > 5000000):
            if self.verbose: print("Select read_count mode") 
            
            start_pos = start
            with pysam.AlignmentFile(self.fp) as bam:
                for i in range(n_step):

                    pos_list.append(round(start_pos))
                    cov_list.append(bam.count(seqid, start_pos, start_pos+step_size))

                    start_pos+=step_size

        # Coverage method, quicker for short intervals
        elif mode == "base_coverage" or (mode == "auto" and len_range <= 5000000):
            if self.verbose: print("Select base_coverage mode")
            with pysam.AlignmentFile(self.fp) as bam:
                mat = bam.count_coverage(seqid, start, end)
                
                start_pos = 0
                while start_pos <= len_range-step_size:
                    
                    s = round(start_pos)
                    e = round(start_pos+step_size)

                    pos_list.append(round(start+start_pos))
                    cov_list.append(sum(mat[0][s:e])+sum(mat[1][s:e])+sum(mat[2][s:e])+sum(mat[3][s:e]))

                    start_pos+=step_size
        
        # Return a tuple of the 2 lists
        return (pos_list, cov_list)
        

