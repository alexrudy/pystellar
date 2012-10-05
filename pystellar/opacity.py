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
        tconv = lambda s : np.int(s.split("#")[1])
        table_summaries = np.genfromtxt(
            StringIO(summary_txt),
            skip_header=c[fkey].get("summary.head",0),
            skip_footer=c[fkey].get("summary.foot",0),
            comments="*",
#            usecols = [1,4,5,6,7,8],
            converters={
                1 : tconv,
                4 : sconv,
                5 : sconv,
                6 : sconv,
                7 : sconv,
                8 : sconv,
            }
        )
        summary = []
        for line in table_summaries:
            summary.append(np.array([line[i] for i in [1,4,5,6,7,8]]))
        summary = np.asarray(summary)
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

def prepare_interpolator(tables, summary, rows, heads):
    """docstring for get_interpolator"""
    from scipy.interpolate import griddata
    
    points = np.zeros((1,7))
    values = np.zeros((1))
    for tn in xrange(tables.shape[0]):
        X = summary[tn][1]
        Y = summary[tn][2]
        Z = summary[tn][3]
        dXc = summary[tn][4]
        dXo = summary[tn][5]
        logRs, logTs = np.meshgrid(rows[tn],heads[tn])
        rossk = tables[tn]
        Xs = np.ones(logRs.shape) * X
        Ys = np.ones(logRs.shape) * Y
        Zs = np.ones(logRs.shape) * Z
        dXcs = np.ones(logRs.shape) * dXc
        dXos = np.ones(logRs.shape) * dXo
        new_points = np.vstack((logRs.flat,logTs.flat,Xs.flat,Ys.flat,Zs.flat,dXcs.flat,dXos.flat)).T
        points = np.concatenate((points,new_points))
        values = np.concatenate((values,rossk.flat))
        
    return points, values

def reload_for_interpolator(fkey):
    """docstring for reload_for_interpolator"""
    c = DottedConfiguration({})
    c.load("data/OPAL.yml",silent=False)
    c.dn = DottedConfiguration
    
    points_filename = "%s.points.npy" % c[fkey]["filename"]
    values_filename = "%s.values.npy" % c[fkey]["filename"]
    
    points = np.load(points_filename)
    values = np.load(values_filename)
    
    return points, values

def save_for_interpolator(points,values,fkey):
    """docstring for save_for_interpolator"""
    
    c = DottedConfiguration({})
    c.load("data/OPAL.yml",silent=False)
    c.dn = DottedConfiguration
    
    points_filename = "%s.points.npy" % c[fkey]["filename"]
    values_filename = "%s.values.npy" % c[fkey]["filename"]
    
    np.save(points_filename,points)
    np.save(values_filename,values)
    
def prepare_point(logrho,logT,X,Y,Z,dXc=0.0,dXo=0.0):
    """Prepare a point for interpolation units:
    
    :param logrho: log(density) in cgs
    :param logT: log(temperature) in cgs
    :param X: fraction of Hydrogen
    :param Y: fraction of Helium
    :param Z: fraction of metals
    :param dXc: Excess fraction of carbon
    :param dXo: Excess fraction of oxygen
    
    """
    assert (X + Y + Z)-1.0 < 1e-10, "Abundances must add up to 1! X: %g, Y: %g, Z: %g Sum: %g" % (X,Y,Z,X+Y+Z)
    T6 = 1e-6 * (10**logT)
    rho = 10.0**logrho
    R = rho / (T6 ** 3)
    logR = np.log10(R)
    return np.array([logR,logT,X,Y,Z,dXc,dXo])
    
    
class OpacityTable(object):
    """This object can interpolate across opacity table information"""
    def __init__(self, fkey, method='nearest'):
        super(OpacityTable, self).__init__()
        self.fkey = fkey
        self.method = method
        try:
            self.load()
        except IOError:
            self.read()
            self.save()
        
    def read(self):
        """Read the opacity tables from an OPAL file"""
        self.points, self.values = prepare_interpolator(*read_table_opal(self.fkey))
        
    def load(self):
        """Load the interpolator values from a file."""
        self.points, self.values = reload_for_interpolator(self.fkey)
        
    def save(self):
        """Save this interpolator to a file"""
        save_for_interpolator(self.points,self.values,self.fkey)
        
    def lookup(self,points):
        """Lookup a group of points."""
        points = np.atleast_2d(points)
        assert points.shape[1] >= 5
        map(lambda p: prepare_point(*p),points)
        return self._lookup(points)
        
    def _lookup(self,point):
        """docstring for lookup"""
        from scipy.interpolate import griddata
        print self.points.shape, self.values.shape
        print point.shape
        return griddata(self.points,self.values,point)
        
    def lookup_single(self,logrho,logT,X,Y,Z,dXc=0.0,dXo=0.0):
        """lookup a single point"""
        point = prepare_point(logrho,logT,X,Y,Z,dXc=0.0,dXo=0.0)
        print "logR:%g,logT:%g,X:%g,Y:%g,Z:%g,dXc:%g,dXo:%g" % tuple(point)
        return self._lookup(point)
        