# -*- coding: utf-8 -*-

"""
  JGV_helper_fun.py
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

# Third party import
from matplotlib.cm import get_cmap

def extensions (fp):
    """ 
    Return the extension of a file in lower-case.
    If archived file "gz" the method will output the base extension + the archive extension
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
    If archived file "gz" the method will output the name without extension + the archive extension
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
    * colormap  colormap color palette from matplotlib package see http://matplotlib.org/examples/color/colormaps_reference.html
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
