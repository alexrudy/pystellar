# -*- coding: utf-8 -*-
# 
#  opacity.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-10-03.
#  Copyright 2012 Alexander Rudy. All rights reserved.
#
"""
:mod:`opacity` - Rosseland Mean Opacity Interpolation
=====================================================

.. autoclass::
    OpacityTable
    :members:
    
Helper Functions
****************

.. autofunction::
    read_table_opal


"""

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

# Internal Modules
from .errors import CodedError

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

class OpacityTableError(CodedError):
    """Error raised due to a poor configuration of the opacity table."""
    
    def __init__(self,*args,**kwargs):
        """Initialization pass-through"""
        super(OpacityTableError, self).__init__(*args, **kwargs)
        log.error(str(self))
    
    codes = {
        2**0 : "Unkown Error!",
        2**1 : "Interpolator was not initialized with a composition",
        2**2 : "Interpolation value for LogR was too large",
        (2**2)+1 : "Interpolation value for LogR was too small",
        2**3 : "Interpolation value for LogT was too large",
        (2**3)+1 : "Interpolation value for LogT was too small",
        2**4 : "Interpolation value failed for unkown reason",
    }
    
    def __str__(self):
        """Set the string representation of this object."""
        if self.msg is None and self.code in self.codes:
            self.msg = self.codes[self.code]
        return super(OpacityTableError, self).__str__()
        
    def __repr__(self):
        """Representation of this error"""
        return "<" + str(self) + ">"
        

class OpacityTable(object):
    """This object can interpolate across opacity table information.
    
    We use a nearest point interpolator to find the composition table that best matches the requested composition.
    A linear interpolation is then used to find requested opacity value.
    
      
    """
    def __init__(self, fkey,load=True, filename="OPAL.yml", X=None, Y=None, snap=False, error=True, wmax=100):
        super(OpacityTable, self).__init__()
        
        # Initialize our attribute values
        
        # Compositional Values
        self._X = None
        self._Y = None
        self._dXc = None
        self._dXo = None
        
        self._interpolator = None # Object for the interpolator
        self._snap = snap   # Snap out-of-bounds objects to the grid.
        self._error = error # Use python errors, or return np.nan
        self._warnings = {
            "NaNs" : 0
        }
        self._warnings_max = wmax
        
        # Logging Values:
        self.log = logging.getLogger(__name__)
        
        # Set up the configuration
        self.fkey = fkey
        self.cfg = DottedConfiguration({})
        self.cfg.load(filename,silent=False)
        self.cfg.dn = DottedConfiguration
        
        # Load the tables (from cached .npy files if appropriate)
        self.log.debug("Loading Tables...")
        if load:
            try:
                self.load()
            except IOError:
                self.read()
                self.save()
        else:
            self.read()
        self.log.debug("Tables Loaded: %s",repr(self._tbls.shape))
        
        # Make ourselves a nearest-neighbor composition interpolator
        self._composition_interpolator()
        
        # If we have a composition, we should use it.
        if X is not None and Y is not None:
            self.composition(X,Y)
        
        self.log.info("Opacity Interpolator Initialzied")
        
    def _composition_interpolator(self):
        """Creates the compositional (X,Y,Z) interpolator"""
        from scipy.interpolate import NearestNDInterpolator
        points = self._comp[:,:4]
        values = self._comp[:,4]
        self.log.debug("Composition Interpolating in %d dimensions on %d points" % (points.shape[1],points.shape[0]))
        self.log.debug("Input Shapes: [x,y]=%r, [v]=%r" % (repr(points.shape), repr(values.shape)))
        nans = (np.sum(np.isnan(points)), np.sum(np.isnan(values)))
        if nans[0] > 0 or nans[1] > 0:
            log.debug("Input NaNs: [x,y]=%g, [v]=%g" % nans)
        self._table_number = NearestNDInterpolator(points,values)
        
    def _temperature_density_interpolator(self):
        """Creates the temperature-density interpolator"""
        from scipy.interpolate import LinearNDInterpolator
        self.log.debug("Initializing Main Interpolator...")
        points = self._tbls[self.n,:,:2]
        values = self._tbls[self.n,:,2]
        points = points[np.isfinite(values)]
        values = values[np.isfinite(values)]
        self.log.debug("Interpolating Opacity in %d dimensions on %d points" % (points.shape[1],points.shape[0]))
        self.log.debug(u"Input Shapes: [logR,logT]=%r, [κ]=%r" % (repr(points.shape), repr(values.shape)))
        
        nans = (np.sum(np.isnan(points)), np.sum(np.isnan(values)))
        if nans[0] > 0 or nans[1] > 0:
            log.debug(u"Input NaNs: [logR,logT]=%g, [κ]=%g" % nans)
        self._interpolator = LinearNDInterpolator(points,values)
        self.log.debug("Created Opacity Interpolator")
        
    def read(self):
        """Read the opacity tables from an OPAL file."""
        self._comp, self._tbls = read_table_opal(self.fkey,cfg=self.cfg)
        self.log.debug("Read Opacity Tables from %(filename)s" % self.cfg[self.fkey])
        
    def load(self):
        """Load the interpolator values from a pair of files. Will load from two numpy files the composition table and the opacity tables."""
        c = self.cfg[self.fkey]
        self._comp = np.load("%s.comp.npy" % c["filename"])
        self._tbls = np.load("%s.tbls.npy" % c["filename"])
        self.log.debug("Loaded Opacity Tables from Numpy Files: %(filename)s.*.npy" % self.cfg[self.fkey])
        
    def save(self):
        """Save this interpolator to a numpy file. The composition and tables will be saved to separate tables."""
        c = self.cfg[self.fkey]
        np.save("%s.comp.npy" % c["filename"],self._comp)
        np.save("%s.tbls.npy" % c["filename"],self._tbls)
        self.log.debug("Saved Opacity Tables to Numpy Files: %(filename)s.*.npy" % self.cfg[self.fkey])
        
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
        """Table number currently in use."""
        return self._n
        
    def properties(self):
        """Return the properties of this object as a tuple:
        
        (n, X, Y, Z, dXc, dXo)
        
        """
        return (self.n, self.X, self.Y, self.Z, self.dXc, self.dXo)
        
    def composition(self,X,Y,dXc=0.0,dXo=0.0):
        ur"""Set the composition for this sytem. Composition must total to 1.0.
        
        Assumes that given X, Y, dXc, dXo, that
        
        .. math::
            X + Y + dXc + dXo = 1 - Z
            
        
        :param X: Fractional content of hydrogen.
        :param Y: Fractional content of helium.
        :param dXc: Fractional content of carbon.
        :param dXo: Fractional content of oxygen.
        
        This function will interpolate to the nearest composition value that it can find. As such, beware that the requested composition will be moved to the nearest acceptable composition in the current opacity tables.
        """
        assert X + Y + dXc + dXo <= 1.0, u"Composition must be less than or equal to 1.0! ∑X=%0.1f" % (X + Y + dXc + dXo)
        if X == self.X and Y == self.Y and dXc == self.dXc and dXo == self.dXo:
            self.log.debug("Values are unchanged")
            return
        point = np.atleast_2d(np.array([X,Y,dXc,dXo]))
        self._n = self._table_number(point)[0]
        self._X = self._comp[self.n,0]
        self._Y = self._comp[self.n,1]
        self._dXc = self._comp[self.n,2]
        self._dXo = self._comp[self.n,3]
        if X != self.X or Y != self.Y or dXc != self.dXc or dXo != self.dXo:
            self.log.warning("Interoplation reached a nearby value, not exact requested composition: X=%.4f, Y=%.4f, dXc=%.4f, dXo=%.4f" % (self.X, self.Y, self.dXc, self.dXo))
        
        self._temperature_density_interpolator()
        
    @classmethod
    def invert_points(cls,logR,logT):
        ur"""Return a log(rho) and log(T) for us.
        
        :param logR: An array (or single value) of :math:`\log(R)` (Table Units).
        :param logT: An array (or single value) of :math:`\log(T)` (Table Units).
        :returns: An array of [:math:`\log(\rho)`, :math:`\log(T)` ].
        
        """
        logR = np.asarray(logR)
        logT = np.asarray(logT)
        R = np.power(10,logR)
        T = np.power(10,logT)
        T6 = 1e-6 * np.power(10,logT)
        rho = R * np.power(T6,3)
        logrho = np.log10(rho)
        return np.vstack((logrho,logT)).T
        
    @classmethod
    def make_points(cls,logrho,logT):
        ur"""Convert the units for a point into proper table units.
        
        :param logrho: An array (or single value) of :math:`\log(\rho)`.
        :param logT: An array (or single value) of :math:`\log(T)`.
        :returns: An array of [:math:`\log(R)`, :math:`\log(T)`].
        
        """
        logrho = np.asarray(logrho)
        logT = np.asarray(logT)
        T6 = 1e-6 * (np.power(10,logT))
        rho = np.power(10,logrho)
        R = rho / (np.power(T6,3))
        logR = np.log10(R)
        return np.vstack((logR,logT)).T
        
    def validate(self,points,RMode=False,Describe=False):
        ur"""Return a boolean if the points guessed are in the table at all.
        
        :param np.array points: An array of :math:`\log(\rho)` and :math:`\log(T)` to be checked.
        :param bool RMode: Whether to assume [logR,logT] (else, assume [logrho,logT])
        :returns: True if the points fit within the table.
        """
        if not RMode:
            points = self.make_points(*points)
        try:
            self.__valid__(points)
        except OpacityTableError as e:
            if Describe:
                return False, e.code, e.codes[e.code]
            else:
                return False
        else:
            if Describe:
                return True, 0, "No Error"
            else:
                return True
        
    def bounds(self,logR=None,logT=None):
        """Return the bounds of the selected opacity table."""
        top = np.array([np.max(self._tbls[self.n,:,i]) for i in xrange(2)])
        bot = np.array([np.min(self._tbls[self.n,:,i]) for i in xrange(2)])
        return np.vstack((top,bot)).T
        
    def snap(self,points):
        """Take a pair of points and place them back on the valid area."""
        maxes = np.array([(points[:,i] <= np.max(self._tbls[self.n,:,i])).any() for i in xrange(2)]).all()
        mines = np.array([(points[:,i] >= np.min(self._tbls[self.n,:,i])).any() for i in xrange(2)]).all()
        if (mines and maxes):
            return points
        
        for point in points:
            for ind,ele in enumerate(point):
                vmax, vmin = np.max(self._tbls[self.n,:,ind]),np.min(self._tbls[self.n,:,ind])
                if ele > vmax:
                    point[ind] = vmax
                elif ele < vmin:
                    point[ind] = vmin
                    
        return points
        
    def __valid__(self,points):
        ur"""Check the range of this point compared to the opacity table range.
        
        :param np.array points: An array of :math:`\log(\rho)` and :math:`\log(T)` to be checked.
        :raises: :exc:`ValueError` when points are out of bounds.
        
        """
        maxes = np.array([(points[:,i] <= np.max(self._tbls[self.n,:,i])).any() for i in xrange(2)]).all()
        mines = np.array([(points[:,i] >= np.min(self._tbls[self.n,:,i])).any() for i in xrange(2)]).all()
        if (mines and maxes):
            return True
        
        cols = {0:"logR",1:"logT"}
        
        for point in points:
            for ind,ele in enumerate(point):
                vmax, vmin = np.max(self._tbls[self.n,:,ind]),np.min(self._tbls[self.n,:,ind])
                if ele > vmax:
                    raise OpacityTableError(msg="BOUNDS: %s=%g > %g" % (cols[ind],ele,vmax),code=2**(ind+2),val=vmax)
                elif ele < vmin:
                    raise OpacityTableError(msg="BOUNDS: %s=%g < %g" % (cols[ind],ele,vmin),code=2**(ind+2)+1,val=vmin)
                    
        raise OpacityTableError(msg="BOUNDS: Error Index Unknown!!!, %r" % points,code=2**4)
        
    def lookup(self,points=None,logrho=None,logT=None):
        ur"""
        A lookup function for our interpolator. Does the pure lookup.
        
        :param points: An array of :math:`\log(\rho)` and :math:`\log(T)`
        :param logrho: An array (or single value) of :math:`\log(\rho)` used only when `points` is not provided.
        :param logT: An array (or single value) of :math:`\log(T)` used only when `points` is not provided.
        :returns: κ, an array of the rosseland mean opacity.
        
        """
        if self._interpolator is None:
            raise OpacityTableError(code=2**1)
        
        if points is None:
            assert logrho is not None, u"Must provide a log(ρ)=?"
            assert logT is not None, u"Must provide a log(T)=?"
            logrho = np.asarray(logrho)
            logT = np.asarray(logT)
            points = self.make_points(logrho=logrho,logT=logT)
        else:
            points = self.make_points(logrho=points[:,0],logT=points[:,1])
        
        if self._snap:
            points = self.snap(points)
        elif self._error:
            self.__valid__(points)
        kappa = self._interpolator(points)
        if np.isnan(kappa).any() and self._error:
            raise ValueError("BOUNDS: Interpolator returned NaN")
        return kappa
        
    def kappa(self,logrho=None,logT=None,rho=None,T=None):
        ur"""Return a rosseland mean opacity at a temperature and density.
        
        :param logrho: Logarithmic Density, base 10, :math:`\log(\rho)`
        :param rho: Density. Accepts `logrho` or `rho`.
        :param logT: Logarithmic Temperature, :math:`\log(T)`
        :param T: Temperature. Accepts `logT` or `T`.
        :returns: κ, the rosseland mean opacity.
        
        """
        assert (logrho is not None) ^ (rho is not None), u"Must provide one and only one value for ρ."
        assert (logT is not None) ^ (T is not None), u"Must provide one and only one value for T"
        
        logT = logT if logT is not None else np.log10(T)
        logrho = logrho if logrho is not None else np.log10(rho)
        kappa = np.power(10,self.lookup(logT=logT,logrho=logrho))
        knans = np.sum(np.isnan(kappa))
        if knans > 0:
            self._warnings["NaNs"] += 1
        if knans > 0 and self._warnings["NaNs"] < self._warnings_max:
            if len(kappa) == 1:
                self.log.warning("Opacity Table Returned NaN: Kappa: %g logT: %g, logrho: %g" % (kappa,logT,logrho))
            else:
                inans = np.sum(np.isnan(logT)) + np.sum(np.isnan(logrho))
                inputs = logT.size + logrho.size
                self.log.warning("Opacity Table Returned NaNs: Kappas: %d/%d, Inputs: %d/%d" % (knans,kappa.size,inans,inputs))
            if T is not None and rho is not None:
                self.log.debug("T: %r Rho: %r" % (T,rho))
        elif knans > 0 and self._warnings["NaNs"] == self._warnings_max:
            self.log.warning("Caught %d NaN Warnings. Future warnings will be suppressed." % self._warnings["NaNs"])
            
        return kappa
    