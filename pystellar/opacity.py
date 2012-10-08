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
from warnings import warn

import logging
# Alex's modules
from AstroObject.config import DottedConfiguration

log = logging.getLogger(__name__)

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
                filling_values=0.0,
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
            tables.append(table[1:,1:])
            rows.append(table[1:,0])
        
        tables = np.asarray(tables)
        rows   = np.asarray(rows)
        heads  = np.asarray(heads)
        
        nrows = c[fkey]["table.end"] - c[fkey]["table.begin"] - c[fkey]["table.head"] - c[fkey]["table.foot"] - 1
        assert tables.shape == (c[fkey]["tables.num"],nrows,c[fkey]["table.cols"]-1), "Reading Table Shape Problem: %s" % repr(tables.shape)
        return tables, summary, rows, heads

def prepare_interpolator(tables, summary, rows, heads):
    """docstring for get_interpolator"""
    from scipy.interpolate import griddata
    
    points = np.zeros((1,4))
    values = np.zeros((1))
    for tn in xrange(tables.shape[0]):
        X = summary[tn][1]
        Y = summary[tn][2]
        Z = summary[tn][3]
        dXc = summary[tn][4]
        dXo = summary[tn][5]
        assert X + Y + Z - 1.0 < 1e-10, "X: %g, Y: %g, Z: %g S: %g" % (X,Y,Z,X+Y+Z)
        logRs, logTs = np.meshgrid(heads[tn],rows[tn])
        rossk = tables[tn]
        Xs = np.ones(logRs.shape).astype(np.float) * X
        Ys = np.ones(logRs.shape).astype(np.float) * Y
        Zs = np.ones(logRs.shape).astype(np.float) * Z
        dXcs = np.ones(logRs.shape).astype(np.float) * dXc
        dXos = np.ones(logRs.shape).astype(np.float) * dXo
        new_points = np.vstack((logRs.flat,logTs.flat,Xs.flat,Ys.flat)).T
        points = np.concatenate((points,new_points))
        values = np.concatenate((values,rossk.flat))
        isnotnan = np.logical_not(np.isnan(values))
        values = values[isnotnan]
        points = points[isnotnan]
        
    
    return points[1:], values[1:]

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
    
def prepare_point(logrho,logT,X=None,Y=None,Z=None,dXc=0.0,dXo=0.0,nargs=7):
    """Prepare a point for interpolation units:
    
    :param logrho: log(density) in cgs
    :param logT: log(temperature) in cgs
    :param X: fraction of Hydrogen
    :param Y: fraction of Helium
    :param Z: fraction of metals
    :param dXc: Excess fraction of carbon
    :param dXo: Excess fraction of oxygen
    
    """
    if X is not None:
        assert (X + Y + Z)-1.0 < 1e-10, "Abundances must add up to 1! X: %g, Y: %g, Z: %g Sum: %g" % (X,Y,Z,X+Y+Z)
    T6 = 1e-6 * (10**logT)
    rho = 10.0**logrho
    R = rho / (T6 ** 3)
    logR = np.log10(R)
    return np.array([logR,logT,X,Y,Z,dXc,dXo])[:nargs]
    
    
class OpacityTable(object):
    """This object can interpolate across opacity table information"""
    def __init__(self, fkey, method='linear',load=True):
        super(OpacityTable, self).__init__()
        self.fkey = fkey
        self.method = method
        self.npts = 5
        log.info("Loading Tables...")
        if load:
            try:
                self.load()
            except IOError:
                self.read()
                self.save()
        else:
            self.read()
        log.info("Tables Loaded: ",self._points.shape)
        
        
    def read(self):
        """Read the opacity tables from an OPAL file"""
        self._points, self._values = prepare_interpolator(*read_table_opal(self.fkey))
        
    def load(self):
        """Load the interpolator values from a file."""
        self._points, self._values = reload_for_interpolator(self.fkey)
        
    def save(self):
        """Save this interpolator to a file"""
        save_for_interpolator(self._points, self._values,self.fkey)
        
    def composition(self,X,Y):
        """Set the composition limits for this sytem"""
        self._X = X
        self._Y = Y
        self._Z = 1.0 - X - Y
        tol = 1e-10
        keep_X = np.abs(self._points[:,2] - X) < tol
        keep_Y = np.abs(self._points[:,3] - Y) < tol
        use = np.logical_and(keep_X,keep_Y)
        if np.sum(use) < 1:
            raise ValueError("Invalid Composition Settings: X %g Y %g" % (X,Y))
        self.values = self._values[use]
        self.points = self._points[use,:2]
        self.npts = 2
        
    def check_range(self,points):
        """Check the range of this point"""
        for point in points:
            for ind,ele in enumerate(point):
                if ele > np.max(self.points[:,ind]):
                    raise ValueError("BOUNDS: p[%g]=%g > %g" % (ind,ele,np.max(self.points[:,ind])))
                elif ele < np.min(self.points[:,ind]):
                    raise ValueError("BOUNDS: p[%g]=%g < %g" % (ind,ele,np.min(self.points[:,ind])))
        
    def setup(self):
        """Set up the ND interpolation option."""
        if not hasattr(self,'points') and not hasattr(self,'values'):
            self.points = self._points
            self.values = self._values
        from scipy.interpolate import LinearNDInterpolator
        log.info("Initializing Interpolator...")
        log.info("Interpolating in %d dimensions on %d points" % (self.points.shape[1],self.points.shape[0]))
        log.debug("Input Shapes: %r, %r" % (repr(self.points.shape), repr(self.values.shape)))
        log.debug("Input NaNs: %g, %g" % (np.sum(np.isnan(self.points)), np.sum(np.isnan(self.values))))
        self._interpolator = LinearNDInterpolator(self.points,self.values)
        log.info("Initialized Interpolator")
        
    def _lookup(self,point):
        """docstring for lookup"""
        from scipy.interpolate import griddata
        point = np.atleast_2d(point)
        self.check_range(point)
        return self._interpolator(point)
        
    def lookup(self,points):
        """Lookup a group of points."""
        points = np.atleast_2d(points)
        assert points.shape[1] >= self.points.shape[1]
        map(lambda p: prepare_point(*p,nargs=self.points.shape[1]),points)
        return self._lookup(points)
        
    def lookup_single(self,logrho,logT,X=None,Y=None,Z=None,dXc=0.0,dXo=0.0):
        """lookup a single point"""
        if X is not None or Y is not None or Z is not None:
            if hasattr(self,'_X'):
                warn("Ignoring lookup_single(X=%g) as composition was fixed." % X)
        point = prepare_point(logrho,logT,X,Y,Z,dXc=0.0,dXo=0.0,nargs=self.points.shape[1])
        log.info("Single Point Requested: ", point)
        return self._lookup(point)
        