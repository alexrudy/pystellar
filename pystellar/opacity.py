# 
#  opacity.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-03.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

from StringIO import StringIO

import numpy as np
import scipy as sp

from pkg_resources import resource_filename

# Alex's modules
from AstroObject.config import DottedConfiguration


def read_table_opal(fkey):
    """Read a given opacity file.
    
    """
    
    # Grab a file configuration system.
    c = DottedConfiguration({})
    c.load("data/OPAL.yml",silent=False)
    c.dn = DottedConfiguration
        
    with open(c[fkey]["filename"],'r') as f:
        flines = f.readlines()
        
        # Read Header Lines
        header = "".join(flines[c[fkey].get("header.begin",0):c[fkey]["header.end"]])
        
        # Read summary lines
        summary_txt = "".join(flines[c[fkey].get("summary.begin",0):c[fkey]["summary.end"]])
        sconv = lambda s : np.float(s.split("=")[1])
        table_summaries = np.genfromtxt(
            StringIO(summary_txt),
            skip_header=c[fkey].get("summary.head",0),
            skip_footer=c[fkey].get("summary.foot",0),
            comments="*",
            converters={
                0 : lambda v : 0,
                1 : lambda v : 0,
                3 : lambda v : 0,
                4 : lambda v : 0,
                5 : sconv,
                6 : sconv,
                7 : sconv,
                8 : sconv,
                9 : sconv,
            }
        )
        
        summary = np.array(table_summaries)
        # summary layout:
        # - Table Number
        # - X
        # - Y
        # - Z
        # - dXc
        # - dXo
        
        # Count Tables
        table_lines = len(flines) - c[fkey]["tables.begin"] - c[fkey]["tables.head"]
        table_num   = float(table_lines) / float(c[fkey]["table.end"] - c[fkey]["table.begin"])
        assert table_num == c[fkey]["tables.num"], "Number of Tables must align!"
        table_num = int(table_num)
        
        # Read tables
        tables = []
        rows = []
        heads = []
        begin = c[fkey]["tables.begin"] + c[fkey]["tables.head"]
        for i in xrange(table_num):
            begin += c[fkey]["table.begin"]
            end   = begin + c[fkey]["table.end"]
            
            # Normalize the table's text
            # Specifically, target short circuited lines.
            table_list = flines[begin:end] # List of lines in a given table
            table_txt = "".join(table_list[:c[fkey].get("table.head",0)]) # Text output of table normalization
            
            # Turn the header row into a numpy array
            head = np.genfromtxt(
                StringIO(table_list[c[fkey].get("table.header")]),
                )[1:]
            heads.append(head)
            
            # Normalize the data table line lengths
            for line in table_list[c[fkey].get("table.head",0):-c[fkey].get("table.foot",0)]:
                char_line = c[fkey]["table.cols"] * c[fkey]["table.chars"] # Character Length of a Line
                line = line.rstrip("\n")
                if len(line) < char_line:
                    line += "  ?????"*((char_line-len(line))/c[fkey]["table.chars"])
                table_txt += line + "\n"
            
            table_txt += "".join(table_list[-c[fkey].get("table.foot",0):]) #Tack the footer rows back on
            
            # Generate the data table from text
            table = np.genfromtxt(
                StringIO(table_txt),
                comments="*",
                filling_values=np.nan,
                missing_values=("?????","9.999"), #We inserted these as unknown values!
                skip_header=c[fkey].get("table.head",0),
                skip_footer=c[fkey].get("table.foot",0),
            )
            begin = end
            tables.append(table[:,1:])
            rows.append(table[:,0])
        
        tables = np.asarray(tables)
        rows   = np.asarray(rows)
        heads  = np.asarray(heads)
        
        nrows = c[fkey]["table.end"] - c[fkey]["table.begin"] - c[fkey]["table.head"] - c[fkey]["table.foot"]
        assert tables.shape == (c[fkey]["tables.num"],nrows,c[fkey]["table.cols"]-1), "Reading Table Shape Problem: %s" % repr(tables.shape)
        return tables, summary, rows, heads


def get_interpolator(tnum, tables, summary, rows, heads):
    """docstring for get_interpolator"""
    from scipy.interpolate import griddata
    
    logRs, logTs = np.meshgrid(rows[tnum],heads[tnum])
    rossk = tables[tnum]
    
    points = np.vstack((logRs.flat,logTs.flat)).T
    values = np.asarray(rossk.flat)
    return lambda p : griddata(points,values,p,method='nearest')
    
    