# -*- coding: utf-8 -*-

# Standard library imports
from collections import OrderedDict, Counter
import gzip

# Third party import
import pysam
import pandas as pd
import numpy as np
from pycl.pycl import extensions_list, file_basename, dir_path, jprint, is_readable_file

#~~~~~~~CLASS~~~~~~~#
class Alignment(object):
    """
    Parse data and compute the base resolution coverage from a file containing aligned reads in BAM, SAM or BED format.
    Can return the coverage for a given interval
    """
    
    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#
    def __init__ (self, fp, name=None, min_coverage=5, refid_list=[], output_bed=False, verbose=False, **kwargs):
        """
         * fp
             A standard BAM or SAM (http://samtools.sourceforge.net/SAM1.pdf) containing aligned reads and a standard
             header. The files do not need to be sorted or indexed.
             One can also use a 6 fields bed (chrom, chromStart, chromEnd, name, score, strand, where score is the
             coverage value (Much faster than from a Bam/Sam file, can be gzipped). http://www.ensembl.org/info/website/upload/bed.html
        *  name
            Name of the data file that will be used as track name for plotting. If not given, will be deduced from fp
            file name  [ DEFAULT: None ]
        * min_coverage
            Minimal coverage to compute the data. If less, the coverage will be considered null. Not used for
            if fp is a bed coverage file [ DEFAULT: 5 ]
        * refid_list
            list of reference sequence id to select from the data file, by default all, Not used for if fp is a bed
            coverage file [ DEFAULT: [] ]
        * output_bed
            If True will be write a 6 columns compressed bed file containing the coverage values for + and - strand
            excluding positions with coverage lesser than min_coverage.the option will apply only is the input file is
            BAM or SAM. [ DEFAULT: False ]
        """
        # Verify that the file is readable
        is_readable_file(fp)

        #Save self variable
        self.fp = fp
        self.name = name if name else file_basename(fp)
        self.ext = extensions_list(fp)[0]
        self.nbases = 0

        if self.ext in ["bam","sam"]:
            if verbose: jprint ("Compute coverage from bam/sam file ", self.fp)
            self.d = self._bam_parser(fp, min_coverage, refid_list)
            if output_bed:
                outfp = "{}/{}.bed.gz".format(dir_path(fp), file_basename(fp))
                if verbose: jprint ("Write coverage data in file ", outfp)
                self._write_coverage_file (outfp)
                self.outfp=outfp

        elif self.ext == "bed":
            if verbose: jprint ("Extract coverage from bed file", self.fp)
            self.d = self._bed_parser(fp, min_coverage, refid_list)
        
        else:
            msg = "The file is not in SAM/BAM/BED format. Please provide a correctly formated file"
            raise ValueError(msg)

        if verbose:
            jprint ("\tTotal base coverage {} in {} reference sequences".format( self.nbases, self.refid_count))

    def __str__(self):
        msg = "{} instance\n".format(self.__class__.__name__)
        # list all values in object dict in alphabetical order
        for k,v in OrderedDict(sorted(self.__dict__.items(), key=lambda t: t[0])).items():
            if k == "d":
                for refid, refval in v.items():
                    msg+="\t{}\tnbases:{}\n".format(refid, refval["nbases"])
            else:
                msg+="\t{}\t{}\n".format(k, v)
        return (msg)

    def __repr__ (self):
        return ("{}: {} - Base coverage {}".format(self.__class__.__name__, self.name, self.name, self.nbases))

    #~~~~~~~PROPERTY METHODS~~~~~~~#
    @property
    def refid_list(self):
        """List of unique reference sequence ids found
        """
        return list(self.d.keys())

    @property
    def refid_nbases(self):
        """List of unique reference sequence ids found associated with their base coverage
        """
        s = pd.Series(name="nbases")
        for refid, refval in self.d.items():
            s.loc[refid] = refval["nbases"]
        return s.sort_values(ascending=False)

    @property
    def refid_count (self):
        """Number of unique reference sequence ids found
        """
        return len(self.d)

    #~~~~~~~PRIVATE METHODS~~~~~~~#
    def _bam_parser(self, fp, min_coverage=5, refid_list= [], verbose=False, **kwargs):
        """Parse a sam or bam formated file
        """
        d = OrderedDict()
        with pysam.AlignmentFile(fp) as bam:
            # Compute the genomic coverage for each reads
            if verbose: jprint ("\tTally coverage for each base")
            for line in bam:
                refid = line.reference_name
                # If not refid filter or if the refid is in the autozized list
                if not refid_list or refid in refid_list:
                    # Create a new entry if not in the dict
                    if not refid in d:
                        d[refid] = {"nbases":0, "+":Counter(), "-":Counter()}
                    # Save coverage
                    strand = "-" if line.is_reverse else "+"
                    for position in line.get_reference_positions():
                        d[refid][strand][position] += 1
        
        d = self._clean_d (d=d, min_coverage=min_coverage, verbose=verbose)
        return d
        
    def _bed_parser (self, fp, min_coverage=5, refid_list= [], verbose=False, **kwargs):
        """Extract data from a coverage bad file
        """
        d = OrderedDict()

        # File handling for both uncompressed or compressed fasta file
        if fp.endswith(".gz"):
            open_fun, open_mode = gzip.open, "rt"
        else:
            open_fun, open_mode = open, "r"
        # Parse fasta file refid
        with open_fun(fp, open_mode) as fin:
            if verbose: jprint ("\tExtract base coverage data")
            for line in fin:
                sl = line[0:-1].split("\t")
                refid = sl[0]
                if not refid_list or refid in refid_list:
                    # Create a new entry if not in the dict
                    if not refid in d:
                        d[refid] = {"nbases":0, "+":Counter(), "-":Counter()}
                    position = int(sl[1])
                    coverage = int(sl[4])
                    strand = sl[5]
                    d[refid][strand][position] = coverage
                    d[refid]["nbases"] += coverage
                    self.nbases += coverage
                    
        d = self._clean_d (d=d, min_coverage=min_coverage, verbose=verbose)
        return d
    
    def _clean_d (self, d, min_coverage=5, verbose=False, **kwargs):
        """ Remove base with coverage below threshold and transform in Pandas Series
        """
        if verbose: jprint ("\tFilter and sort the coverage results by position")
        for refid in d.keys():
            for strand in ["+","-"]:
                s = OrderedDict()
                for position, coverage in d[refid][strand].items():
                    if coverage >= min_coverage:
                        s[position] = coverage
                        self.nbases += coverage
                        d[refid]["nbases"] += coverage
                d[refid][strand] = pd.Series(s)
                d[refid][strand].sort_index(inplace=True)
        return d
    
    def _write_coverage_file (self, outfp, buf_size=8192, verbose=False, **kwargs):
        """Bufferized writer for the coverage bed file
        """
        with gzip.open (outfp, "wt") as out:
            # Write header containing chromosome information\t
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
            # Empty the rest of the buffer
            out.write(str_buf)
            str_buf = ""

    #~~~~~~~PUBLIC METHODS~~~~~~~#
    def interval_coverage (self, refid, start, end, bins=500, bin_repr_fun = "max", verbose=False, **kwargs):
        """
        Parse the alignment file for a given refid and interval. The interval is splited in a number of windows equal to
        bins, for which the coverage in computed. The method return a dataframe containing the starting positions of
        the windows and the coverage for the + and - strands. If the refid or the coordinates are invalid a zero filled
        dataframe will be returned.
        * refid
            Name of the sequence from the original fasta file to display
        * start
            Start of the window to display. The coordinate is not verified, if outside of the range it will
            return empty bins
        * end
            End of the window to display. The coordinate is not verified, if outside of the range it will
            return empty bins
       * bins
            Number of alignment count bins to divide the displayed window. Low number will result in low resolution
            high value could result in a long ploting time. The value is automatically adjusted if lower than base
            resolution, ie if the requested interval is lower than the number of bins [ DEFAULT: 500 ]
        * bin_repr_fun
            Function to represent each bin ("max", "mean" and "sum") [ DEFAULT: "max" ]
        """
        if verbose: jprint ("Compute coverage from the windows: {}:{}-{}".format(refid, start, end))
        df = pd.DataFrame(columns=["+", "-"], dtype=int)

        # Adjust number of bins and calculate step
        if bins > end-start:
            bins = end-start
            if verbose: jprint ("\tAuto adjust the number of bins to match the interval: {}".format(bins))
        step = (end-start)/bins
        if verbose: jprint ("\tDefine size of each bin: {}".format(step))

        # If refid is not in the self refid-list
        if not refid in self.refid_list:
            if verbose: jprint ("\tThe reference {} is not in the list of references with alignment".format(refid))
            for i in np.arange (start, end, step):
                for strand in ["+","-"]:
                    df.loc[int(i), strand] =  0
            return df

        # Select positions windows and get maximun
        if verbose: jprint ("\tCompute coverage...")
        for i in np.arange (start, end, step):
            winstart = int(i)
            winend = int(i+step)
            for strand in ["+","-"]:
                l = self.d[refid][strand][(self.d[refid][strand].index >= winstart)&(self.d[refid][strand].index < winend)]
                if l.empty:
                    df.loc[winstart, strand] =  0
                elif bin_repr_fun == "max":
                    df.loc[winstart, strand] = l.max()
                elif bin_repr_fun == "sum":
                    df.loc[winstart, strand] = l.sum()
                elif bin_repr_fun == "mean":
                    df.loc[winstart, strand] = l.sum()/step
        if verbose:
            if df["+"].sum() + df["-"].sum() == 0:
                jprint ("\tNull coverage for both strands in the requested interval")
            elif df["+"].sum() == 0:
                jprint ("\tNull coverage for the positive strand in the requested interval")
            elif df["-"].sum() == 0:
                jprint ("\tNull coverage for the negative strand in the requested interval")
        return df
