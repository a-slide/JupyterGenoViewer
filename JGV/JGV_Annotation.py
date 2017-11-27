# -*- coding: utf-8 -*-

# Strandard library imports
from collections import OrderedDict
from os import access, R_OK

# Third party import
import pandas as pd
from pycl.pycl import extensions_list, has_extension, file_basename, jprint, is_readable_file

#~~~~~~~CLASS~~~~~~~#
class Annotation(object):
    """
    Parse data from a file containing genomic annotation in GFF3, GTF or BED format.
    Can return the list of annotations for a given interval
    """

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, fp, name=None, min_len=None, max_len=None, refid_list=None, type_list=None, verbose=False, **kwargs):
        """
         * fp
            A path to a standard genomic file containing features annotations among the following format
            gff3: http://www.ensembl.org/info/website/upload/gff3.html
            gtf: http://www.ensembl.org/info/website/upload/gff.html
            bed:  http://www.ensembl.org/info/website/upload/bed.html
            Alternatively, one can use a python pickle file (.pkl) generated during a previous run.
            The file can eventually be compressed in ‘gzip’ format
        * min_len
            Minimal size (start to end) of a feature to be selected [default None]
        * max_len
            Maximal size (start to end) of a feature to be selected [default None]
        * refid_list
            List of reference id to select. Example: ["chr1", "chr2", "chr3"] [default None]
        * type_list
            List of feature type to select. Example: ["exon", "gene"] [default None]
        """
        if verbose: jprint ("Parse Annotation file")
        # Verify that the file is readable
        is_readable_file(fp)

        #Save self variable
        self.fp = fp
        self.name = name if name else file_basename(fp)

        # Find if gziped
        if has_extension (fp, pos=-1, ext=["gz","tgz"]):
            if verbose: jprint("\tFile is gziped")
            compression="gzip"
            ext_pos=-2
        else:
            if verbose: jprint("\tFile is not compressed")
            compression=None
            ext_pos=-1

        # Find extension type
        if has_extension (fp, pos=ext_pos, ext="gtf"):
            self.feature_df = self._gtf_parser(fp=fp, compression=compression, verbose=verbose)
        elif has_extension (fp, pos=ext_pos, ext="gff3"):
            self.feature_df = self._gff3_parser(fp=fp, compression=compression, verbose=verbose)
        elif has_extension (fp, pos=ext_pos, ext="bed"):
            self.feature_df = self._bed_parser(fp=fp, compression=compression, verbose=verbose)

        # Else try to import as a pickled file
        else:
            try:
                self.feature_df = self._pickle_parser(fp=fp, verbose=verbose)
            # If invalid file format
            except Exception as E:
                raise ValueError("Cannot open file or the file is not in a valid format")

        # Optional filterig steps
        if min_len or max_len:
            self.select_len (min_len=min_len, max_len=max_len, verbose=verbose)
        if refid_list:
            self.select_references (refid_list=refid_list, verbose=verbose)
        if type_list:
            self.select_types (type_list=type_list, verbose=verbose)

        # Sort the dataframe and reset index
        if verbose: jprint("Sorting and final cleanup")
        self.feature_df.sort_values(by=["refid","start","end"], inplace=True)
        self.feature_df.reset_index(drop=True, inplace=True)

        if verbose: jprint("\tNumber of features imported: {}".format(self.feature_count))

    def __repr__ (self):
        return ("{}: {} - Feature count {}".format(self.__class__.__name__, self.name, self.feature_count))

    #~~~~~~~PROPERTY METHODS~~~~~~~#
    @property
    def feature_count(self):
        """Number of features collected"""
        return len(self.feature_df)

    @property
    def refid_count(self):
        """Number of unique reference sequence ids found"""
        return self.feature_df["refid"].nunique()

    @property
    def type_count(self):
        """Number of unique feature type found"""
        return self.feature_df["type"].nunique()

    @property
    def refid_list(self):
        """List of unique reference sequence ids found"""
        return set(self.feature_df["refid"])

    @property
    def type_list(self):
        """List of unique feature type found"""
        return set(self.feature_df["type"])

    @property
    def refid_count_uniq(self):
        """List of unique reference sequence ids with count of associated features"""
        return pd.DataFrame(self.feature_df.groupby("refid").size().sort_values(ascending=False), columns=["count"])

    @property
    def type_count_uniq(self):
        """List of unique feature types with count of associated features"""
        return pd.DataFrame(self.feature_df.groupby("type").size().sort_values(ascending=False), columns=["count"])

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def interval_features (self, refid, start, end, feature_types=None, max_features_per_type=None, verbose=False, **kwargs):
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
            if verbose: jprint ("The reference {} is not in the list of references with alignment".format(refid))
            return pd.DataFrame(columns=["refid","start","end","strand","ID","type"])

        # Select the refid and coordinates
        df = self.feature_df[(self.feature_df["refid"] == refid)&(self.feature_df["end"] > start)&(self.feature_df["start"] < end)]
        if df.empty:
            if verbose: jprint ("No feature found in the requested interval")
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
                    if verbose: jprint ("No feature of type {} found in the requested interval".format(type_name))
                elif max_features_per_type and len(sdf)>max_features_per_type:
                    select_list.append(sdf.sample(max_features_per_type))
                else:
                    select_list.append(sdf)
        # Merge the selected features in a single df
        if select_list:
            df = pd.concat(select_list)
        else:
            if verbose: jprint ("No feature found in the requested interval")
            return pd.DataFrame(columns=["refid","start","end","strand","ID","type"])

        # Return a sorted copy of the  df
        df.sort_values(by=["refid","start","end"], inplace=True)
        return df.copy().reset_index(drop=True)

    def select_len (self, min_len=None, max_len=None, verbose=False, **kwargs):
        """ Select features longer or shorter that given values
        """
        if verbose:
            jprint ("Selecting features based on length")
            jprint ("\tFeatures before filtering: {}".format(self.feature_count))
        # Filter min len
        if min_len:
            self.feature_df = self.feature_df[((self.feature_df["end"]-self.feature_df["start"]) >= min_len)]
        if verbose: jprint ("\tFeatures after minimal length filtering: {}".format(self.feature_count))
        # Filter max len
        if max_len:
            self.feature_df = self.feature_df[((self.feature_df["end"]-self.feature_df["start"]) <= max_len)]
        if verbose: jprint ("\tFeatures after maximal length filtering: {}".format(self.feature_count))

    def select_references (self, refid_list, verbose=False, **kwargs):
        """ Select features which reference sequence id is in the given list or a single entry. Example: ["chr1", "chr2", "chr3"]
        """
        # Cast in list
        if type(refid_list) == str:
            refid_list = [refid_list]

        if verbose:
            jprint ("Selecting features based on reference id")
            jprint ("\tFeatures before filtering: {}".format(self.feature_count))
        self.feature_df = self.feature_df[(self.feature_df["refid"].isin(refid_list))]
        if verbose: jprint ("\tFeatures after filtering: {}".format(self.feature_count))

    def select_types (self, type_list, verbose=False, **kwargs):
        """ Select features which type is in the given list or a single entry. Example: ["exon", "gene"]
        """
        # Cast in list
        if type(type_list) == str:
            type_list = [type_list]

        if verbose:
            jprint ("Selecting features based on type")
            jprint ("\tFeatures before filtering: {}".format(self.feature_count))
        self.feature_df = self.feature_df[(self.feature_df["type"].isin(type_list))]
        if verbose: jprint ("\tFeatures after filtering: {}".format(self.feature_count))

    def to_pickle (self, fp=None, verbose=False, **kwargs):
        """
        Store the parsed file in a pickle file for further use.
        * fp
            Path to save the pickle file. By default original annotation file (- .gz/.tgz) + .pkl
        """
        if not fp:
            if self.fp.endswith(".gz") or self.fp.endswith(".tgz"):
                fp = self.fp.rpartition(".")[0]+".pkl"
            else:
                fp = self.fp+".pkl"

        if verbose: jprint ("Pickle dataframe in file {}".format(fp))
        self.feature_df.to_pickle(fp)
        return fp

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _bed_parser(self, fp, compression=None, verbose=False, **kwargs):
        """
        Parse a bed formated file
        """
        if verbose: jprint("\tUse BED parser to parse annotations")
        # try to import the file as a bed6 in a dataframe
        try:
            col_names = ["refid","start","end","ID","score","strand"]
            df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#", compression=compression)
            if verbose: jprint("\tSuccessfully imported as a bed6 file")

        # else try to import as a bed12
        except IndexError as E:
            col_names = ["refid","start","end","ID","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]
            df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#", compression=compression)
            if verbose: jprint("\tSuccessfully imported as a bed12 file")

        # Type is not available from bed files
        df['type'] = "."

        # Clean df
        df = self._clean_df(df, verbose=verbose)
        return df

    def _gff3_parser(self, fp, compression=None, verbose=False, **kwargs):
        """
        Parse a gff3 formated file
        """
        if verbose: jprint("\tUse GFF3 parser to parse annotations")
        # Import the file in a dataframe
        col_names = ["refid","source","type","start","end","score","strand","frame","attribute"]
        df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#", compression=compression)

        # Extract the ID field = first field of attribute
        df['ID'] = df["attribute"].str.split(';').str[0].str[3:]
        if verbose: jprint("\tSuccessfully imported as a gff3 file")

        # Clean df
        df = self._clean_df(df, verbose=verbose)
        return df

    def _gtf_parser(self, fp, compression=None, verbose=False, **kwargs):
        """
        Parse a gtf formated file
        """
        if verbose: jprint("\tUse GTF parser to parse annotations")
        # Import the file in a dataframe
        col_names = ["refid","source","type","start","end","score","strand","frame","attribute"]
        df = pd.read_csv(fp, sep="\t", names=col_names, index_col=False, comment="#", compression=compression)

        # Extract the ID field = first field of attribute=
        df['ID'] = df["attribute"].str.split('\"').str[1]
        if verbose: jprint("\tSuccessfully imported as a gtf file")

        # Clean df
        df = self._clean_df(df, verbose=verbose)
        return df

    def _clean_df (self, df, verbose=False, **kwargs):
        """
        Clean dataframe after parsing
        """
        # Select fields
        df = df[["refid","start","end","ID","score","strand","type"]].copy()

        # Drop column with NA values
        if verbose: jprint("\tRemove null values")
        l = len(df)
        df.dropna(inplace=True)
        if verbose: jprint("\tRemoved {} invalid lines".format(l-len(df)))

        # Cast the start and end field in integer
        if verbose: jprint("\tCast coordinates to integer and id to str")
        df[['start', 'end']] = df[['start', 'end']].astype(int)
        df[['ID']] = df[['ID']].astype(str)

        # Verify than the dataframe is not empty
        if df.empty:
            raise ValueError("No valid features imported. Is the file valid?")
        return df

    def _pickle_parser(self, fp, verbose=False, **kwargs):
        """
        Parse a pickle database
        """
        # Import the file in a dataframe
        if verbose: jprint ("\tTry to load as a pickle file")
        df = pd.read_pickle(fp)
        return df
