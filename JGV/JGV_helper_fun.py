# -*- coding: utf-8 -*-

# Standard library imports
import os

# Third party imports
from matplotlib.cm import get_cmap
from IPython.core.display import display, HTML, Markdown

#~~~~~~~FUNCTIONS~~~~~~~#
def extensions (fp):
    """
    Return the extension of a file in lower-case.
    If the file is gziped the method will output the base extension + the archive extension
    """
    split_name = fp.split("/")[-1].split(".")
    # No extension ?
    if len (split_name) == 1:
        return []
    # Manage compressed files
    elif len (split_name) > 2 and split_name[-1].lower() == "gz":
        return [split_name[-2].lower(), split_name[-1].lower()]
    # Normal situation = return the last element of the list
    else:
        return [split_name[-1].lower()]

def file_basename (fp):
    """
    Return the base name of a file without extension nor path.
    If the file is gziped the method will output the name without extension + the archive extension
    """
    split_name = fp.split("/")[-1].split(".")
    # No extension ?
    if len (split_name) == 1:
        return split_name[0]
    # Manage compressed files
    elif len (split_name) > 2 and split_name[-1].lower() == "gz":
        return ".".join(split_name[0:-2])
    # Normal situation = return the last element of the list
    else:
        return ".".join(split_name[0:-1])

def dir_path (fp):
    """
    Return the directory path of a file
    """
    return fp.rpartition("/")[0]

def color_palette(n, colormap="brg"):
    """
    Return a list of n length with gradient colors from a given matplot lib colormap palette
    * n         Number of color scalar in the list
    * colormap  colormap color palette from matplotlib package.
                See http://matplotlib.org/examples/color/colormaps_reference.html
                example : inferno magma hot blues cool spring winter brg ocean hsv jet ... [DEFAULT: brg]
    """
    # Init variables
    cmap = get_cmap(colormap)
    index = 1 # skip the first value
    n_col_cmap = cmap.N-2 # remove border colors
    step = int(n_col_cmap/(n-1)) if n > 1 else n_col_cmap/2

    # Create the list of colors
    cl = []
    for i in range (n):
        cl.append(cmap(index))
        index+=step
    return cl

def jprint(*args, **kwargs):
    """
    Format a string in HTML and print the output. Equivalent of print, but highly customizable. Many options can be passed to the function.
    * args
        One or several objects that can be cast in str
    ** kwargs
        Formatting options to tweak the html rendering
        Boolean options : bold, italic, highlight, underlined, striked, subscripted, superscripted
        String options: font, color, size, align, background_color
    """

    # Join the different elements together and cast in string
    s =  " ".join([str(i) for i in args])

    # Replace new lines and tab by their html equivalent
    s = s.replace("\n", "<br>").replace("\t", "&emsp;")

    # For boolean options
    if "bold" in kwargs and kwargs["bold"]: s = "<b>{}</b>".format(s)
    if "italic" in kwargs and kwargs["italic"]: s = "<i>{}</i>".format(s)
    if "highlight" in kwargs and kwargs["highlight"]: s = "<mark>{}</mark>".format(s)
    if "underlined" in kwargs and kwargs["underlined"]: s = "<ins>{}</ins>".format(s)
    if "striked" in kwargs and kwargs["striked"]: s = "<del>{}</del>".format(s)
    if "subscripted" in kwargs and kwargs["subscripted"]: s = "<sub>{}</sub>".format(s)
    if "superscripted" in kwargs and kwargs["superscripted"]: s = "<sup>{}</sup>".format(s)

    # for style options
    style=""
    if "font" in kwargs and kwargs["font"]: style+= "font-family:{};".format(kwargs["font"])
    if "color" in kwargs and kwargs["color"]: style+= "color:{};".format(kwargs["color"])
    if "size" in kwargs and kwargs["size"]: style+= "font-size:{}%;".format(kwargs["size"])
    if "align" in kwargs and kwargs["align"]: style+= "text-align:{};".format(kwargs["align"])
    if "background_color" in kwargs and kwargs["background_color"]: style+= "background-color:{};".format(kwargs["background_color"])

    # Format final string
    if style: s = "<p style=\"{}\">{}</p>".format(style,s)
    else: s = "<p>{}</p>".format(s)

    display(HTML(s))

def jhelp(function, full=False):
    """
    Print a nice looking help string based on the name of a declared function. By default print the function definition and description 
    * full
        If True, the help string will included a description of all arguments
    """
    try:
        # For some reason signature is not aways importable. In these cases the build-in help in invoqued
        from inspect import signature, isfunction, ismethod
        if isfunction(function) or ismethod(function):
            name = function.__name__.strip()
            sig = str(signature(function)).strip()
            display(HTML ("<b>{}</b> {}".format(name, sig)))
            
            if function.__doc__:
                for line in function.__doc__.split("\n"):
                    line = line.strip()
                    if not full and line.startswith("*"):
                        break
                    display(Markdown(line.strip()))
        else:
            jprint("{} is not a function".format(function))

    except Exception:
        help(function)

def get_sample_file (package, path):
    """
    Verify the existence and return a file from the package data
    * package
        Name of the package
    * path
        Relative path to the file in the package. Usually package_name/data/file_name 
    """
    try:
        # Try to extract package with pkg_resources lib
        from pkg_resources import Requirement, resource_filename
        fp = resource_filename(Requirement.parse(package), path)
        if not os.access(fp, os.R_OK):
            raise IOError("Can not read {}".format(fp))
        else:
            return fp
        
        # Try local package instead
        fp = path
        if not os.access(fp, os.R_OK):
            raise IOError("Can not read {}".format(fp))
        else:
            return fp
                
    except Exception as E:
        jprint(E)
        jprint ("Please retrieve it from the github repository")
        return
