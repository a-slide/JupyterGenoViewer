# -*- coding: utf-8 -*-

"""
  Alignment.py
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
from collections import OrderedDict, Counter
from os import access, R_OK
from gzip import open as gopen

# Third party import
import pysam
import pandas as pd

# Local lib import
from JGV_helper_fun import *
from JGV_helper_fun import jprint as print

#~~~~~~~CLASS~~~~~~~#

class Alignment(object):

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, fp, name=None, verbose=False, min_coverage=5, ref_list=[], output_bed=False):
        """
         * fp
             A standard BAM or SAM (http://samtools.sourceforge.net/SAM1.pdf) containing aligned reads.
             The files do not need to be sorted or indexed
        *  name
            Name of the data file that will be used as track name for plotting. If not given, will be deduced from fp file name
        * verbose
            If True, will print more information during initialisation and calls of all the object methods.
        * min_coverage
            [ DEFAULT: 5 ]
        * ref_list
            [ DEFAULT: [] ]
        * output_coverage_file
            [ DEFAULT: None ]
        """
        # Verify that the file is readable
        assert access(fp, R_OK), "{} is not readable".format(fp)

        #Save self variable
        self.fp = fp
        self.name = name if name else file_basename(fp)
        self.verbose = verbose
        self.ext = extensions(fp)[0]
        self.nbases = 0

        if self.ext in ["bam","sam"]:
            if self.verbose: print("Compute coverage from bam/sam file ", self.fp)
            self.d = self._bam_parser(fp, min_coverage, ref_list)
        elif self.ext == "bed":
            if self.verbose: print("Extract coverage from bed file", self.fp)
            self.d = self._bed_parser(fp, min_coverage, ref_list)
        else:
            msg = "The file is not in SAM/BAM/BED format. Please provide a correctly formated file"
            raise ValueError(msg)

        if output_bed:
            #assert access(output_file, W_OK), "{} is not writable".format(fp)
            outfp = "{}/{}.bed.gz".format(dir_path(fp), file_basename(fp))
            if self.verbose: print("Write coverage data in file ", outfp)
            self._write_coverage_file (outfp)
            self.outfp=outfp

    def __str__(self):
        """readable description of the object"""
        msg = "{} instance\n".format(self.__class__.__name__)
        # list all values in object dict in alphabetical order
        for k,v in OrderedDict(sorted(self.__dict__.items(), key=lambda t: t[0])).items():
            if k == "d":
                for refid, refval in v.items():
                    msg+="\t{}\tlength:{}\tnbases:{}\n".format(refid, refval["length"], refval["nbases"])
            else:
                msg+="\t{}\t{}\n".format(k, v)
        return (msg)

    def __repr__ (self):
        return ("{}-{}".format(self.__class__.__name__, self.name))

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _bam_parser(self, fp, min_coverage=5, ref_list= []):

        d = OrderedDict()

        # Parse Bam/Sam file
        with pysam.AlignmentFile(fp) as bam:

            # Save information for References based on the bam header
            for refid, length in zip(bam.references, bam.lengths):
                if not ref_list or refid in ref_list:
                    d[refid] = {"length":length, "nbases":0, "+":Counter(), "-":Counter()}

            # Tally the coverage for each base
            for line in bam:
                refid = line.reference_name
                if refid in d:
                    strand = "-" if line.is_reverse else "+"
                    for position in line.get_reference_positions():
                        d[refid][strand][position] += 1

        # Remove base with coverage below threshold and transform in Pandas Series
        for refid in d.keys():
            for strand in ["+","-"]:
                s = OrderedDict()
                for position, coverage in d[refid][strand].items():
                    if coverage >= min_coverage:
                        s[position] = coverage
                        self.nbases += coverage
                        d[refid]["nbases"] += coverage
                d[refid][strand] = pd.Series(s, dtype="int64")
                d[refid][strand].sort_index(inplace=True)

        return d

    def _bed_parser (self, fp, min_coverage=5, ref_list = []):

        d = OrderedDict()

        try:
            fin = gopen(fp, "rt") if fp[-2:].lower() == "gz" else open (file, "r")

            for line in fin:
                if line.startswith("#"):
                    sl = line[1:-1].split("\t")
                    refid = sl[0]
                    if not ref_list or refid in ref_list:
                        length = int(sl[1])
                        d[refid] = {"length":length, "nbases":0, "+":OrderedDict(), "-":OrderedDict()}
                else:
                    sl = line[0:-1].split("\t")
                    refid = sl[0]
                    if refid in d:
                        position = int(sl[1])
                        coverage = int(sl[4])
                        strand = sl[5]
                        d[refid][strand][position] = coverage
                        d[refid]["nbases"] += coverage
                        self.nbases += coverage

            for refid in d.keys():
                for strand in ["+","-"]:
                    d[refid][strand] = pd.Series(d[refid][strand], dtype="int64")

        # close the file properly
        finally:
            return d
            try:
                f.close()
            except:
                pass

    def _write_coverage_file (self, outfp, buf_size=8192 ):
        with gopen (outfp, "wt") as out:

            # Write header containing chromosome information
            for refid, refval in self.d.items():
                out.write("#{}\t{}\n".format(refid, refval["length"]))

            # Bufferized writing of lines
            i = 0
            str_buf = ""
            for refid, refval in self.d.items():
                for strand in ["+","-"]:
                    for position, coverage in refval[strand].items():
                        i+=1
                        str_buf += "{0}\t{1}\t{1}\tpos{2}\t{3}\t{4}\n".format(refid, position, i, coverage, strand)
                        if i%buf_size == 0:
                            out.write(str_buf)
                            str_buf = ""

                    out.write(str_buf)
                    str_buf = ""

    #~~~~~~~PROPERTY METHODS~~~~~~~#

    @property
    def refid_len(self):
        """List of all the sequence ids found in the annotation file"""
        s = pd.Series(name="length")
        for refid, refval in self.d.items():
            s.loc[refid] = refval["length"]
        return s.sort_values(ascending=False)

    @property
    def refid_nbases(self):
        """List of all the sequence ids found in the annotation file"""
        s = pd.Series(name="nbases")
        for refid, refval in self.d.items():
            s.loc[refid] = refval["nbases"]
        return s.sort_values(ascending=False)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def refid_count (self):
        """
        Return a dataframe containing 3 colums ("reads", "rpkm") for all seq_id reference
        sequence in the BAM file. The read count is computationally intensive and can take time for large BAM files
        """

        df = pd.DataFrame(columns=["nbases", "bpkm"])

        for refid, refval in self.d.items():
            df.loc[refid] = [
                refval["nbases"],
                refval["nbases"]/(refval["length"]/1000)/(self.nbases/1000000)]

        return df

    def interval_coverage (self, refid, start=None, end=None, bins=500, bin_repr_fun = "max"):
        """
        Parse the alignment file for a given seqid and interval. Compute the coverage at read level in mode "read_count"
        or at base level in mode "base_coverage". In both case the interval is splited in a number of windows equal to
        n_step, for which the coverage in computed. The method return a dataframe containing the starting positions of
        the windows and the coverage.
        * seqid
            Name of the sequence from the original fasta file to display
        * start
            Start of the window to display. If not given will start from 0 [ DEFAULT: None ]
        * end
            End of the window to display. If not given will start from end of the sequence [ DEFAULT: None ]
       * n_step
            Number of alignment count bins to divide the displayed window. Low number will result in low resolution but
        """
        # Verifications and auto adjustment of coordinates
        assert refid in self.d, "The reference {} is not in the list of references with alignment".format(refid)
        if not start or start < 0:
            start = 0
            if self.verbose: print ("Auto adjust start position: {}".format(start))
        if not end or end > self.d[refid]["length"]:
            end = self.d[refid]["length"]
            if self.verbose: print ("Auto adjust end position: {}".format(end))
        assert start < end, "Start coordinate must be smaller than end"
        if bins > end-start:
            bins = end-start
            if self.verbose: print ("Auto adjust the number of bins to match the interval: {}".format(bins))
        step = (end-start)/bins
        if self.verbose: print ("Define size of each bin: {}".format(step))

        # Select positions windows and get maximun
        df = pd.DataFrame(columns=["+", "-"], dtype=int)

        for i in np.arange (start, end, step):
            winstart = int(i)
            winend = int(i+step)

            for strand in ["+","-"]:
                l = self.d[refid][strand][(self.d[refid][strand].index >= winstart)&(self.d[refid][strand].index < winend)]
                if l.empty:
                    df.loc[int(i), strand] =  0
                elif bin_repr_fun == "max":
                    df.loc[int(i), strand] = l.max()
                elif bin_repr_fun == "sum":
                    df.loc[int(i), strand] = l.sum()
                elif bin_repr_fun == "mean":
                    df.loc[int(i), strand] = l.sum()/step

        return df
