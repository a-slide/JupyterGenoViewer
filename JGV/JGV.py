# -*- coding: utf-8 -*-

"""
  JGV.py
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

# Standard library imports
from collections import OrderedDict
import warnings
from sys import exit as sysexit

# Third party imports
try:
    current= "Numpy 1.11.1"
    import numpy as np
    current = "Matplotlib 1.5.1"
    from matplotlib.patches import FancyArrowPatch as Arrow
    from matplotlib.gridspec import GridSpec
    import pylab as pl
    current = "Pandas 0.18.1"
    import pandas as pd
    current = "Pysam 0.9.0"
    import pysam
    current = "Jupyter 4.1.0"
    cfg = get_ipython()
    from IPython.core.display import display
except (NameError, ImportError):
    print ("The third party package {}+ is required by JVG. Please verify your dependencies".format(current))
    sysexit()

# Local lib import
try:
    current = "JGV_helper_fun"
    from JGV_helper_fun import extensions, file_basename, dir_path, color_palette
    from JGV_helper_fun import jprint as print
    current = "Reference"
    from Reference import Reference
    current = "Annotation"
    from Annotation import Annotation
    current = "Alignment"
    from Alignment import Alignment
    current = "Level"
    from Level import Level
except ImportError:
    print ("Can not import the local packages {}. Please verify JVG source code directory".format(current))
    sysexit()

#~~~~~~~CLASS~~~~~~~#
class JGV(object):
    version = "0.0.1"

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#
    def __init__ (self, fp, name=None, verbose=False, ref_list=[], output_index=False):
        """
         * fp
            A fasta file containing the reference sequences OR an tab separated index file containing at least 2 columns
            with the refid and the length in bases (like a .fa.fai file generated by samtools faidx).
            The fasta option will take more time as the file has to be parsed to get the refid and length of sequences.
            A 2 column index tsv file will be automatically generated for latter usage as an index file.
            Both fasta and infex file can be gziped
        *  name
            Name of the data file that will be used as track name for plotting. If not given, will be deduced from fp
            file name
        * verbose
            If True, will print more information during initialisation and calls of all the object methods.
        * ref_list
            list of reference sequence id to select from the data file, by default all [ DEFAULT: [] ]
        * output_index
            If True will write a simple A 2 column index tsv file containing the Reference sequence ids and their
            lengths [ DEFAULT: False ]
        """
        #Save self variable
        self.verbose = verbose
        self.ref_list = ref_list

        # Store the reference genome informations
        if self.verbose: print("Add reference genome file", bold=True)
        self.reference = Reference(
            fp=fp,
            name=name,
            verbose=self.verbose,
            ref_list=self.ref_list,
            output_index=output_index)

        # List to store annotations and alignment tracks
        self.annotations = []
        self.alignments = []

    def __str__(self):
        """readable description of the object"""
        msg = "{} instance\n".format(self.__class__.__name__)
        msg+= "\tParameters list\n"
        # list all values in object dict in alphabetical order
        sorted_d = OrderedDict(sorted(self.__dict__.items(), key=lambda t: t[0]))
        for k,v in sorted_d.items():
            msg+="\t{}\t{}\n".format(k, v)
        return (msg)

    #~~~~~~~PUBLIC SET METHODS~~~~~~~#
    def add_annotation(self, fp, name=None):
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
        """
        if self.verbose: print("Add annotation file", bold=True)
        a = Annotation(
            fp=fp,
            name=name,
            verbose=self.verbose,
            ref_list=self.ref_list)

        if self.verbose:
            for refid in self.reference.refid_list:
                if refid not in  a.refid_list:
                    warnings.warn("No annotation found for {}".format(refid))

        self.annotations.append(a)

    def add_alignment(self, fp, name=None, min_coverage=5, output_bed=False):
        """
         * fp
             A standard BAM or SAM (http://samtools.sourceforge.net/SAM1.pdf) containing aligned reads and a standard
             header. The files do not need to be sorted or indexed.
             One can also use a 6 fields bed (chrom, chromStart, chromEnd, name, score, strand) file with a hastaged
             commented header listing the reference sequences id and length, similar to the format generated by the
             output_bed option (Much faster than from a Bam/Sam file, can be gzipped).
             http://www.ensembl.org/info/website/upload/bed.html
        *  name
            Name of the data file that will be used as track name for plotting. If not given, will be deduced from fp
            file name  [ DEFAULT: None ]
        * min_coverage
            Minimal coverage to compute the data. If less, the coverage will be considered null. Not used for
            if fp is a bed coverage file [ DEFAULT: 5 ]
        * output_bed
            If True will be write a 6 columns compressed bed file containing the coverage values for + and - strand
            excluding positions with coverage lesser than min_coverage.the option will apply only is the input file is
            BAM or SAM. The file starts with a header consisting of a list of the ID of the reference sequences and
            their length [ DEFAULT: False ]. Example:
              #chr20	64444167
              #chr21	46709983
              chr20	276516	276516	pos1	5	+
              chr20	276517	276517	pos2	5	+
        """
        if self.verbose: print("Add alignment file", bold=True)
        a = Alignment(
            fp=fp,
            name=name,
            verbose=self.verbose,
            min_coverage=min_coverage,
            ref_list=self.ref_list,
            output_bed=output_bed)

        if self.verbose:
            for refid in self.reference.refid_list:
                if refid not in  a.refid_list:
                    warnings.warn("No coverage found for {}".format(refid))

        self.alignments.append(a)

    def annotation_summary (self):
        """Display table summarizing annotation file information"""
        if not self.annotations:
            warnings.warn("No annotation track loaded")
            return None

        count_df = pd.DataFrame(columns=["Feature count", "Refid count", "Feature type count"])
        rcu_df = pd.DataFrame()
        tcu_df = pd.DataFrame()

        for a in self.annotations:
            count_df.loc[a.name] = [a.feature_count, a.refid_count, a.type_count]
            rcu = a.refid_count_uniq
            rcu.columns=[a.name]
            rcu_df = pd.merge(left=rcu_df, right=rcu, how='outer', right_index=True, left_index=True)
            tcu = a.type_count_uniq
            tcu.columns=[a.name]
            tcu_df = pd.merge(left=tcu_df, right=tcu, how='outer', right_index=True, left_index=True)

        print("Counts per Annotation file", bold=True)
        display(count_df)
        print("Counts per Reference sequence", bold=True)
        rcu_df.sort_index()
        display(rcu_df)
        print("Counts per feature types", bold=True)
        tcu_df.sort_index()
        display(tcu_df)

    def alignment_summary (self):
        """Display table summarizing annotation file information"""
        if not self.alignments:
            warnings.warn("No alignment track loaded")
            return None

        count_df = pd.DataFrame(columns=["Refid count", "Base coverage"])
        rbc_df = pd.DataFrame()

        for a in self.alignments:
            count_df.loc[a.name] = [a.refid_count, a.nbases]
            rbc = pd.DataFrame(a.refid_nbases)
            rbc.columns=[a.name]
            rbc_df = pd.merge(left=rbc_df, right=rbc, how='outer', right_index=True, left_index=True)

        print("Counts per Alignment file", bold=True)
        display(count_df)
        print("Counts per Reference sequence", bold=True)
        rbc_df.sort_index()
        display(rbc_df)

    def refid_coverage_plot (self,
        norm_depth = True,
        norm_len =  True,
        plot_style="ggplot",
        figwidth = 20,
        figheight = 5,
        log = False,
        ref_list = [],
         **kwargs):
        """
        * norm_len
            If True, for each refid, the base counts are normalised by the length of refid in bases
        * norm_depth [ DEFAULT: True ]
            If True, for each track, the base counts are normalised by the overall number of bases mapped
        * plot_style [ DEFAULT: True ]
            Default plot style for pyplot ('grayscale'|'bmh'|'ggplot'|'dark_background'|'classic'|'fivethirtyeight'...)
            [ DEFAULT: "ggplot" ]
        * figwidth
             Width of the ploting area in inches [ DEFAULT: 20 ]
        * figheight
             height of the ploting area in inches [ DEFAULT: 5 ]
        * log
            if True the yscale will be log10 else it will be linear [ DEFAULT: True ]
        * ref_list
            list of reference sequence id to display, by default all. The list is also used to reorder the reference in
            the dataframe and in the graph [ DEFAULT: [] ]
        * kwargs
            Additional parameters for plot appearance derived from pylab basic plot arguments such as: color, alpha,
            fontsize...
        """
        if not self.alignments:
            warnings.warn("No alignment track loaded")
            return None

        df = pd.DataFrame()
        for a in self.alignments:
            refid_df = pd.DataFrame(a.refid_nbases)
            refid_df.columns=[a.name]
            df = pd.merge(left=df, right=refid_df, how='outer', right_index=True, left_index=True)

        # Filter no listed refid is requested + reorder index
        if ref_list:
            df = df[(df.index.isin(ref_list))]
            df = df.reindex(ref_list)
        else:
            df.sort_index(inplace=True)

        # Normalize by lenth and depth is requested
        if norm_depth:
            for name, col in df.iteritems():
                df[name] = col/col.sum()*1000
        if norm_len:
            for refid in df.index:
                ref_len = self.reference.get_refid_len(refid)
                df.loc[refid] = df.loc[refid]/ref_len*1000000

        # Prepare default plotting options
        fontsize = kwargs["fontsize"] if "fontsize" in kwargs else 12
        color =  kwargs["color"] if "color" in kwargs else None
        alpha = kwargs["alpha"] if "alpha" in kwargs else 1
        pl.style.use(plot_style)

        # Plot data
        fig, ax = pl.subplots(figsize=(figwidth, figheight))
        df.plot.bar(ax=ax, alpha=alpha, color = color)

        # Tweak plot
        if log: ax.set_yscale("log")
        a = ax.set_title("Distribution of bases per reference sequence", fontsize=fontsize+2, y=1.05)
        a = ax.tick_params(axis='both', which='major', labelsize=fontsize, bottom=True, top=False, left=True, right=False)
        a = ax.legend(bbox_to_anchor=(1, 1), loc=2,frameon=False, fontsize=fontsize)

        return df

    def interval_plot (self,
        refid,
        start=None,
        end=None,
        plot_style="ggplot",
        figwidth = 30,
        alignment_track_height=5,
        annotation_track_height=2,
        alignment_bins = 500,
        alignment_bin_repr_fun = "max",
        alignment_log=True,
        alignment_color=("dodgerblue", "darkorange"),
        alignment_alpha=0.5,
        feature_types=[],
        max_features_per_type=500,
        annotation_offset=None,
        annotation_label=False,
        annotation_color="grey",
        **kwargs):
        """
        * refid
            Name of the sequence from the original fasta file to display
        * start
            Start of the window to display. If not given, will be set to 0 [ DEFAULT: None ]
        * end
            End of the window to display. If not given, will be set to the length of refid [ DEFAULT: None ]
        * plot_style [ DEFAULT: True ]
            Default plot style for pyplot ('grayscale'|'bmh'|'ggplot'|'dark_background'|'classic'|'fivethirtyeight'...)
            [ DEFAULT: "ggplot" ]
        * figwidth
             Width of the ploting area in inches [ DEFAULT: 20 ]
        * alignment_track_height
            Height of individual aligment tracks [DEFAULT : 5 ]
        * annotation_track_height
            Height of individual annotation tracks for each feature types [DEFAULT : 2 ]
        * alignment_bins
            Number of alignment count bins to divide the displayed window. Low number will result in low resolution
            high value could result in a long ploting time. The value is automatically adjusted if lower than base
            resolution, ie if the requested interval is lower than the number of bins [ DEFAULT: 500 ]
        * alignment_bin_repr_fun
            Function to represent each bin ("max", "mean" and "sum") [ DEFAULT: "max" ]
        * alignment_log
            if True the yscale will be log10 else it will be linear [ DEFAULT: True ]
        * alignment_color
            Tuple of 2 color for the alignment + and - tracks [DEFAULT : ("dodgerblue", "darkorange") ]
        * alignment_alpha
            Transparency of the alignment coverage area between 0 and 1 [ DEFAULT: 0.5 ]
        * feature_types
            Name of a valid feature type ( "exon"|"transcript"|"gene"|"CDS"...) or list of names of feature type for
            which a row will be returned. The option is not available for bed files. If not given, all features type
            found in the interval will be returned [ DEFAULT: None ]
        * max_features_per_type
            Maximal total number of features for a particular feature type. If more are found, a random sampling will
            be performed. If None, all the features will be returned [ DEFAULT: 500 ]
        * annotation_offset
            Minimal distance between 2 contigous annotation features on the same level. If not given, will be
            automatically set to 1/400 of the windows to display [DEFAULT : None ]
        * annotation_label
            If True, labels of features will be plotted. To be avoid when expecting many features [DEFAULT : False ]
        * annotation_color
            Color of the annotation arrows [DEFAULT : "grey" ]
        * kwargs
        """
        # Verify that the sequence is in the refid list and that at least one alignment or annotation file was loaded
        if refid not in self.reference.refid_list:
            warnings.warn("Requested reference sequence not found: {}".format(refid))
            return None
        if not self.alignments and not self.annotations:
            warnings.warn("No annotation and alignment track loaded")
            return None

        # Auto define start and stop and overlapping annotation offset if not given
        if not start:
            start = 0
            if self.verbose: print ("Autodefine start position: {}".format(start))
        if not end:
            end = self.reference.get_refid_len(refid)-1
            if self.verbose: print ("Autodefine end position: {}".format(end))
        if start >= end:
            raise ValueError ("Invalid coordinates (start: {}, end :{}) start has to be greater than end")
        if not annotation_offset:
            annotation_offset = int((end-start)/400)
            if self.verbose:print ("Estimated overlap offset: {}".format(annotation_offset))

        figheight = 0

        # Extract alignment coverage data and compute the coverage tracks height
        alignments_dict = OrderedDict()
        if self.alignments:
            if self.verbose: print ("Extract alignment data", bold=True)
            for a in self.alignments:
                alignments_dict[a.name] = a.interval_coverage(
                    refid=refid, start=start, end=end, bins=alignment_bins, bin_repr_fun=alignment_bin_repr_fun)
                # +1 for space bewtween tracks
                figheight += alignment_track_height+1

        # Extract feature annotation data and compute the feature tracks height
        annot_tracks_heigth = 0
        annotation_dict = OrderedDict()
        if self.annotations:
            if self.verbose: print ("Extract annotation data", bold=True)
            for a in self.annotations:
                annotation_dict[a.name] = a.interval_features(
                    refid=refid, start=start, end=end, feature_types=feature_types,
                    max_features_per_type=max_features_per_type)
                # Take empty df into account for ploting
                n = 1 if annotation_dict[a.name].empty else annotation_dict[a.name].type.nunique()
                # +1 for space bewtween tracks
                figheight += n*annotation_track_height+1

        # Create a pylot figure object with an empty grid
        fig = pl.figure (figsize= (figwidth, figheight))
        grid = GridSpec (nrows=figheight, ncols=1, hspace=0.5)
        pl.style.use (plot_style)

        # Curent height marker
        h = 0
        if self.alignments:
            for track_name, track_df in alignments_dict.items():
                if self.verbose: print ("\tAlignment track name: {}".format(track_name))

                # Prepare the subplot grid
                ax = pl.subplot(grid[h:h+alignment_track_height])
                h+=alignment_track_height
                ax.set_xlim((start, end))
                ax.ticklabel_format(useOffset=False, style='plain')
                if alignment_log: ax.set_yscale("log")
                ax.yaxis.set_tick_params(left=True, right=False, labelleft=True, labelright=False)
                ax.xaxis.set_tick_params(bottom=False, top=False, labelbottom=False, labeltop=False)
                ax.set_ylabel(track_name)

                # Plot the positive strand
                if track_df["+"].sum () == 0:
                    ax.text(0.5, 0.6,'No Coverage on the + strand', ha='center', va='center', transform=ax.transAxes)
                else:
                    ax.fill_between(x=track_df.index, y1=0, y2=list(track_df["+"]),
                        alpha=alignment_alpha, color=alignment_color[0], label="Positive strand")
                # Plot the negative strand
                if track_df["-"].sum () == 0:
                    ax.text(0.5, 0.4,'No Coverage on the - strand', ha='center', va='center', transform=ax.transAxes)
                else:
                    ax.fill_between(x=track_df.index, y1=0, y2=list(track_df["-"]),
                        alpha=alignment_alpha, color=alignment_color[1], label="Negative strand")
                # If elements were added to the ax
                if ax.collections: ax.legend(bbox_to_anchor=(1, 1), loc=2,frameon=False)

            # Add x labels if last element
            ax.xaxis.set_tick_params(bottom=True, labelbottom=True)

        if self.annotations:
            for track_name, track_df in annotation_dict.items():
                h+=1
                if self.verbose: print ("\tAlignment track name: {}".format(track_name))

                # No feature case
                if track_df.empty:
                    ax = pl.subplot(grid[h:h+annotation_track_height])
                    ax.set_xlim((start, end))
                    ax.text(0.5, 0.5,'No feature found', ha='center', va='center', transform=ax.transAxes)
                    ax.yaxis.set_tick_params(left=False, right=False, labelleft=False, labelright=False)
                    ax.xaxis.set_tick_params(bottom=True, labelbottom=True)
                    ax.ticklabel_format(useOffset=False, style='plain')
                    ax.grid(axis="y", b=False)
                    ax.set_title (track_name)
                    h+=annotation_track_height

                # General case
                else:
                    first=True
                    for feature_type, feature_df in track_df.groupby("type"):

                        # Prepare the ploting area
                        ax = pl.subplot(grid[h:h+annotation_track_height])
                        h+=annotation_track_height
                        ax.set_xlim((start, end))
                        ax.yaxis.set_tick_params(left=False, right=False, labelleft=False, labelright=False)
                        ax.xaxis.set_tick_params(bottom=False, top=False, labelbottom=False, labeltop=False)
                        ax.ticklabel_format(useOffset=False, style='plain')
                        ax.grid(axis="y", b=False)
                        ax.set_ylabel(feature_type)

                        # Compute the non overlaping level where to plot the arrow
                        level = Level(offset=annotation_offset)
                        for n, feature in feature_df.iterrows():
                            fl = level(feature.ID, feature.start, feature.end, feature.strand)
                            if fl:
                                ax.add_patch( Arrow( posA=[fl.start, fl.level], posB=[fl.end, fl.level], linewidth=3,
                                    color=annotation_color, arrowstyle=fl.arrowstyle))
                                if annotation_label:
                                    text_end = fl.end if fl.end < end-annotation_offset else end-annotation_offset
                                    text_start = fl.start if fl.start > start+annotation_offset else start+annotation_offset
                                    ax.text (x=text_start+ (text_end-text_start)/2, y=fl.level, s=fl.ID, ha="center")

                        ax.set_ylim(level.min_level-0.5, level.max_level+0.5)

                        # First element exception
                        if first:
                            ax.set_title (track_name)
                            first = False
                    # Last elemet exception
                    ax.xaxis.set_tick_params(bottom=True, labelbottom=True)
