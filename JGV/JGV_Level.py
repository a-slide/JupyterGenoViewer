# -*- coding: utf-8 -*-

# Strandard library imports
from collections import OrderedDict, namedtuple, Counter

#~~~~~~~CLASS~~~~~~~#
class Level (object):
    """
    Compute the level of a given feature on the Annotation track to avoid annotation overlaping
    """

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#
    def __init__ (self,
                  max_depth=100,
                  offset=10,
                  filter_pos=False,
                  filter_neg=False,
                  filter_unstrand=False):
        """
        Define general options for Level class
        * max_depth
            Maximal total number of positive or negative levels. Safeguard value in case of highly
            overlapped annotations [ DEFAULT: 100 ]
        * offset
            Minimal distance between 2 contigous annotation features on the same level [ DEFAULT: 10 ]
        * filter_pos
            Filter-out annotation features on the positive strand [ DEFAULT: False ]
        * filter_neg
            Filter-out annotation features on the negative strand [ DEFAULT: False ]
        * filter_unstrand
            Filter-out annotation features with no strand specified [ DEFAULT: False ]
        """
        # Save general parameter
        self.max_depth = max_depth
        self.offset = offset
        self.pos_arrowstyle = "-|>,head_width=1,head_length=2"
        self.neg_arrowstyle = "<|-,head_width=1,head_length=2"
        self.unstrand_arrowstyle = "-"
        self.filter_pos = filter_pos
        self.filter_neg = filter_neg
        self.filter_unstrand = filter_unstrand

        # Create self containers
        self.level_dict={}
        self.count = Counter()
        self.enhanced_feature = namedtuple('enhanced_feature', ['ID','start', 'end', "arrowstyle", "level"])

    def __str__(self):
        """readable description of the object"""
        msg = "{} instance\n".format(self.__class__.__name__)

        # list all values in object dict in alphabetical order
        for k,v in OrderedDict(sorted(self.__dict__.items(), key=lambda t: t[0])).items():
            msg+="\t{}\t{}\n".format(k, v)
        return (msg)

    def __repr__ (self):
        return ("{}".format(self.__class__.__name__))

    @property
    def min_level(self):
        """Return the minimal level index"""
        return min(self.level_dict.keys())

    @property
    def max_level(self):
        """Return the minimal level index"""
        return max(self.level_dict.keys())

    @property
    def n_level(self):
        """Return the total number of levels index"""
        return len(self.level_dict)

    #~~~~~~~PUBLIC METHODS~~~~~~~#
    def __call__ (self, ID, start, end, strand):
        """
        Compute the level of an annnotation feature based on the instance options and the other features previously
        analysed to avoid overlapping. Iterative call of the function has to be done with annotation features sorted
        by start coordinates.
        * ID
            Name of the feature to fit in a level
        * start
            Start coordinate of the feature to fit in a level, on the positive strand
        * end
            End coordinate of the feature to fit in a level, on the positive strand
        * strand
            Strand of the feature. Can be + - or . if unknown
        """
        self.count["all_features"] +=1

        # For features on the positive strand
        if strand == "+" and not self.filter_pos:
            level=1
            self.count["positive_features"] +=1

            while level <= self.max_depth:
                # If level is empty or if the level is free at this position
                if level not in self.level_dict or (self.level_dict[level]+self.offset) < start:
                    self.level_dict[level] = end
                    return self.enhanced_feature (ID, start, end, self.pos_arrowstyle, level)
                level+=1

        # For features on the negative strand
        elif strand == "-" and not self.filter_neg:
            level=-1
            self.count["negative_features"] +=1

            while level >= -self.max_depth:
                # If level is empty or if the level is free at this position
                if level not in self.level_dict or (self.level_dict[level]+self.offset) < start:
                    self.level_dict[level] = end
                    return self.enhanced_feature (ID, start, end, self.neg_arrowstyle, level)
                level-=1

        elif strand == "." and not self.filter_unstrand:
            self.count["unstranded_features"] +=1
            self.level_dict[0] = end
            return self.enhanced_feature (ID, start, end, self.unstrand_arrowstyle, 0)

        return None
