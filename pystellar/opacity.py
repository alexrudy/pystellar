# -*- coding: utf-8 -*-
# 
#  opacity.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-03.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 
from __future__ import division

from StringIO import StringIO

import numpy as np
import scipy as sp

from pkg_resources import resource_filename
from warnings import warn

from multiprocessing import Queue, Process

import logging
# Alex's modules
from AstroObject.config import DottedConfiguration

log = logging.getLogger(__name__)

def read_table_opal(fkey,cfg):
    """Read a given OPAL format opacity file.
    
    """
    # We use a configuration in a nested dictionary form to store all of the constants
    c = cfg[fkey]
    with open(c["filename"],'r') as f:
        flines = f.readlines()
        
    # SANITY CHECKS:
    # Count the number of OPAL tables to be read:
    table_lines = len(flines) - c["tables.begin"] - c["tables.head"]
    table_num   = table_lines / (c["table.end"] - c["table.begin"])
    table_line  = table_lines / table_num
    assert table_num == c["tables.num"], "Number of Tables must add up to configured value! Counted: %f (%f), CFG: %f" % (table_num,table_line,c["tables.num"])
    table_num = int(table_num)
    
    # Read Header Lines
    header = "".join(flines[c.get("header.begin",0):c["header.end"]])
    
    # Read summary lines
    # Summary Table lines show the compositional fractions X,Y and Z for each table.
    summary_text = "".join(flines[c.get("summary.begin",0):c["summary.end"]])
    
    # Conversion functions read ASCII data and automatically convert it into numbers.
    sconv = lambda s : np.float(s.split("=")[-1]) # Converts text which is of the form X=0.0 to a float
    tconv = lambda s : np.float(s.split("#")[-1])   # Converts text which is of the form TABLE #  1 to a float.
    summary = np.genfromtxt(StringIO(summary_text),
        skip_header=c.get("summary.head",0),
        skip_footer=c.get("summary.foot",0),
        # dtype=np.float,
        # dtype=str,
        comments="*",
        delimiter=c["summary.cols"],
        autostrip=True,
        converters={
            0 : tconv,
            6 : sconv, 3 : sconv, 4 : sconv, 5 : sconv, 2 : sconv,
            }
    )[:,[0,2,3,4,5,6]]
    # summary layout:
    # - Table Number - X - Y - Z - dXc - dXo
    # We have cropped out the 'source' column which 'cites' each OPAL table.
    
    # Read tables
    nrows = c["table.end"] - c["table.begin"] - c["table.head"] - c["table.foot"]
    tables = np.empty((c["tables.num"],nrows,c["table.cols"]-1))
    rows = np.empty((c["tables.num"],nrows))
    heads = np.empty((c["tables.num"],c["table.cols"]-1))
    begin = c["tables.begin"] + c["tables.head"]
    char_line = c["table.cols"] * c["table.chars"] # Character Length of a Line
    for i in xrange(table_num):
        begin += c["table.begin"]
        end   = begin + c["table.end"]
        # Normalize the table's text
        # Specifically, target short circuited lines.
        table_list = flines[begin:end] # List of lines in a given table
        table_txt = "".join(table_list[:c.get("table.head",0)]) # Text output of table normalization
        
        # Turn the header row into a numpy array
        head = np.genfromtxt(
            StringIO(table_list[c.get("table.header")]),
            filling_values=0.0,
            )[1:]
        
        if -c.get("table.foot",0) == 0:
            table_data = table_list[c.get("table.head",0):]
        else:
            table_data = table_list[c.get("table.head",0):-c.get("table.foot",0)]
            
        # Normalize the data table line lengths
        for line in table_data:
            line = line.rstrip("\n")
            if len(line) < char_line:
                line += "  9.999"*((char_line-len(line))//c["table.chars"])
            table_txt += line + "\n"
        table_txt += "".join(table_list[-c.get("table.foot",0):0]) #Tack the footer rows back on
        # Generate the data table from text
        table = np.genfromtxt(
            StringIO(table_txt),
            comments="*",
            filling_values=np.nan,
            missing_values="9.999", #We inserted these as unknown values!
            skip_header=c.get("table.head",0),
            skip_footer=c.get("table.foot",0)-1,
        )
        
        table[table == 9.999] = np.nan
        
        begin = end
        
        tables[i,:,:] = table[:,1:]
        rows[i,:] = table[:,0]
        heads[i,:] = head
        
    tables = np.asarray(tables)
    rows   = np.asarray(rows)
    heads  = np.asarray(heads)
    # Test that we got the right number of table cellls out in the end
    assert tables.shape == (c["tables.num"],nrows,c["table.cols"]-1), "Reading Table Shape Problem: %s" % repr(tables.shape)
    
    # Now we want to convert the tables into co-ordinate format.
    # For tables, this will be [logR,logT,kappa]
    # For composition, this will be [X,Y,dXc,dXo,Table Number]
    
    # Co-ordinate arrays for tables
    tbls = np.empty((table_num,rows.shape[1]*heads.shape[1],3))
    
    # Co-ordinate arrays for composition
    comp = np.empty((table_num,5))
    
    # Iterate through each table to make the table arrays
    for tn in xrange(tables.shape[0]):
        X = summary[tn][1]
        Y = summary[tn][2]
        Z = summary[tn][3]
        dXc = summary[tn][4]
        dXo = summary[tn][5]
        assert X + Y + Z - 1.0 < 1e-10, "X: %g, Y: %g, Z: %g S: %g" % (X,Y,Z,X+Y+Z)
        logRs, logTs = np.meshgrid(heads[tn],rows[tn])
        rossk = tables[tn]
        
        # Add this table's composition co-ordinates to global array
        comp[tn] = np.array([X,Y,dXc,dXo,tn])
        
        # Make our co-ordinate arrays for this table.
        tbls[tn] = np.vstack((logRs.flat,logTs.flat,rossk.flat)).T
        
    
    return comp, tbls


class OpacityTable(object):
    """This object can interpolate across opacity table information.
    
    We use a nearest point interpolator to find the composition table that best matches the requested composition.
    A linear interpolation is then used to find requested opacity value.
    
      
    """
    def __init__(self, fkey,load=True, filename="OPAL.yml", X=None, Y=None):
        super(OpacityTable, self).__init__()
        
        # Initialize our attribute values
        self._X = None
        self._Y = None
        self._dXc = None
        self._dXo = None
        
        # Set up the configuration
        self.fkey = fkey
        self.cfg = DottedConfiguration({})
        self.cfg.load(filename,silent=False)
        self.cfg.dn = DottedConfiguration
        
        # Load the tables (from cached .npy files if appropriate)
        log.info("Loading Tables...")
        if load:
            try:
                self.load()
            except IOError:
                self.read()
                self.save()
        else:
            self.read()
        log.info("Tables Loaded: %s",repr(self._tbls.shape))
        
        # Make ourselves a nearest-neighbor composition interpolator
        self._composition_interpolator()
        log.info("Opacity Interpolator Initialzied")
        
        # If we 
        if X is not None and Y is not None:
            self.composition(X,Y)
        
    def _composition_interpolator(self):
        """Creates the compositional (X,Y,Z) interpolator"""
        from scipy.interpolate import NearestNDInterpolator
        points = self._comp[:,:4]
        values = self._comp[:,4]
        log.info("Interpolating in %d dimensions on %d points" % (points.shape[1],points.shape[0]))
        log.debug("Input Shapes: [x,y]=%r, [v]=%r" % (repr(points.shape), repr(values.shape)))
        nans = (np.sum(np.isnan(points)), np.sum(np.isnan(values)))
        if nans[0] > 0 or nans[1] > 0:
            log.debug("Input NaNs: [x,y]=%g, [v]=%g" % nans)
        self._table_number = NearestNDInterpolator(points,values)
        
    def _temperature_density_interpolator(self):
        """Creates the temperature-density interpolator"""
        from scipy.interpolate import LinearNDInterpolator
        log.info("Initializing Interpolator...")
        points = self._tbls[self.n,:,:2]
        values = self._tbls[self.n,:,2]
        points = points[np.isfinite(values)]
        values = values[np.isfinite(values)]
        log.info("Interpolating in %d dimensions on %d points" % (points.shape[1],points.shape[0]))
        log.debug(u"Input Shapes: [logR,logT]=%r, [κ]=%r" % (repr(points.shape), repr(values.shape)))
        
        nans = (np.sum(np.isnan(points)), np.sum(np.isnan(values)))
        if nans[0] > 0 or nans[1] > 0:
            log.debug(u"Input NaNs: [logR,logT]=%g, [κ]=%g" % nans)
        self._interpolator = LinearNDInterpolator(points,values)
        log.info("Created Interpolator")
        
    def read(self):
        """Read the opacity tables from an OPAL file"""
        self._comp, self._tbls = read_table_opal(self.fkey,cfg=self.cfg)
        
    def load(self):
        """Load the interpolator values from a file."""
        c = self.cfg[self.fkey]
        self._comp = np.load("%s.comp.npy" % c["filename"])
        self._tbls = np.load("%s.tbls.npy" % c["filename"])
        
    def save(self):
        """Save this interpolator to a file"""
        c = self.cfg[self.fkey]
        np.save("%s.comp.npy" % c["filename"],self._comp)
        np.save("%s.tbls.npy" % c["filename"],self._tbls)
        
    @property
    def X(self):
        """Fraction of H atoms"""
        return self._X
    
    @property
    def Y(self):
        """Fraction of He atoms"""
        return self._Y
        
    @property
    def Z(self):
        """Fraction of 'metal' atoms"""
        return (1.0 - self.X - self.Y - self.dXo - self.dXc)
        
    @property
    def dXo(self):
        """Fraction of oxygen greater than Z"""
        return self._dXo
        
    @property
    def dXc(self):
        """Fraction of carbon greater than Z"""
        return self._dXc
        
    @property
    def n(self):
        """Table number"""
        return self._n
        
    def properties(self):
        """Return the properties of this object as a tuple"""
        return (self.n, self.X, self.Y, self.Z, self.dXc, self.dXo)
        
    def composition(self,X,Y,dXc=0.0,dXo=0.0):
        """Set the composition for this sytem. Composition must b"""
        assert X + Y + dXc + dXo <= 1.0, u"Composition must be less than or equal to 1.0! ∑X=%0.1f" % (X + Y + dXc + dXo)
        if X == self.X and Y == self.Y and dXc == self.dXc and dXo == self.dXo:
            log.debug("Values are unchanged")
            return
        point = np.atleast_2d(np.array([X,Y,dXc,dXo]))
        self._n = self._table_number(point)[0]
        self._X = self._comp[self.n,0]
        self._Y = self._comp[self.n,1]
        self._dXc = self._comp[self.n,2]
        self._dXo = self._comp[self.n,3]
        self._temperature_density_interpolator()
        
    def invert_points(self,logR,logT):
        """Return a log(rho) and log(T) for us."""
        logR = np.asarray(logR)
        logT = np.asarray(logT)
        R = np.power(10,logR)
        T = np.power(10,logT)
        T6 = 1e-6 * np.power(10,logT)
        rho = R * np.power(T6,3)
        logrho = np.log10(rho)
        return np.vstack((logrho,logT)).T
        
    def _make_points(self,logrho,logT):
        """Convert the units for a point into proper table units"""
        logrho = np.asarray(logrho)
        logT = np.asarray(logT)
        T6 = 1e-6 * (np.power(10,logT))
        rho = np.power(10,logrho)
        R = rho / (np.power(T6,3))
        logR = np.log10(R)
        return np.vstack((logR,logT)).T
        
    def check_range(self,points):
        """Check the range of this point"""
        maxes = np.array([(points[:,i] <= np.max(self._tbls[self.n,:,i])).any() for i in xrange(2)]).all()
        mines = np.array([(points[:,i] >= np.min(self._tbls[self.n,:,i])).any() for i in xrange(2)]).all()
        if (mines and maxes):
            log.debug("Passed Tests: %r, %r" % (mines,maxes))
            return
        
        for point in points:
            for ind,ele in enumerate(point):
                if ele > np.max(self._tbls[self.n,:,ind]):
                    raise ValueError("BOUNDS: p[%g]=%g > %g" % (ind,ele,np.max(self._tbls[self.n,:,ind])))
                elif ele < np.min(self._tbls[self.n,:,ind]):
                    raise ValueError("BOUNDS: p[%g]=%g < %g" % (ind,ele,np.min(self._tbls[self.n,:,ind])))
                    
        raise RuntimeError("BOUNDS: Error Index Unknown!!!, %r" % points)
        
    def lookup(self,points=None,logrho=None,logT=None):
        """A lookup function for our interpolator."""
        if points is None:
            assert logrho is not None, u"Must provide a log(ρ)=?"
            assert logT is not None, u"Must provide a log(T)=?"
            logrho = np.asarray(logrho)
            logT = np.asarray(logT)
            points = self._make_points(logrho=logrho,logT=logT)
        else:
            points = self._make_points(logrho=points[:,0],logT=points[:,1])
        
        self.check_range(points)
        kappa = self._interpolator(points)
        if np.isnan(kappa).any():
            raise ValueError("BOUNDS: Interpolator returned NaN")
        return kappa
    