# -*- coding: utf-8 -*-
# 
#  results.py
#  pystellar
#  
#  Created by Alexander Rudy on 2012-12-02.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 
from __future__ import division
import numpy as np
import scipy as sp
"""
Functions for outputting pretty results at the end (or from saved data files.)
"""

def compare(end_guesses,comp_guesses):
    """Comparison to a ZAMS model."""
    names = ["R_s","L_s","P_c","T_c"]
    units = [r"\unit{cm}",r"\unit{ergs s}^{-1}",r"\unit{dyne cm}^{-2}",r"\unit{K}"]
    header = r"& This Model & ZAMS Model \\"
    template = r"""
    \begin{tabular}{l|l|l}
    %(header)s \hline \hline
    %(table)s
    \end{tabular}
    """
    rowTemplate = r"$%(name)s$ & $%(start)s$ & $%(end)s$ \\"
    rows = []
    for name,unit,start,end in zip(names,units,end_guesses,comp_guesses):
        fdict = {
            'name': name,
            'start': make_sci_tex(start,unit),
            'end': make_sci_tex(end,unit),
        }
        rows.append(rowTemplate % fdict)
    fdict = {
        'header' : header,
        'table' : "\n".join(rows),
    }
    return template % fdict
    

def guess_table(start_guesses,end_guesses):
    """Print out a guess table in TeX format"""
    names = ["R_s","L_s","P_c","T_c"]
    units = [r"\unit{cm}",r"\unit{ergs s}^{-1}",r"\unit{dyne cm}^{-2}",r"\unit{K}"]
    header = r"& Initial & Converged \\"
    template = r"""
    \begin{tabular}{l|l|l}
    %(header)s \hline \hline
    %(table)s
    \end{tabular}
    """
    rowTemplate = r"$%(name)s$ & $%(start)s$ & $%(end)s$ \\"
    rows = []
    for name,unit,start,end in zip(names,units,start_guesses,end_guesses):
        fdict = {
            'name': name,
            'start': make_sci_tex(start,unit),
            'end': make_sci_tex(end,unit),
        }
        rows.append(rowTemplate % fdict)
    fdict = {
        'header' : header,
        'table' : "\n".join(rows),
    }
    return template % fdict
    
def fp_table(start_fp,end_fp,start_c,end_c):
    """Fitting points table."""
    names = ["\Delta r(m_{fp})","\Delta l(m_{fp})","\Delta P(m_{fp})","\Delta T(m_{fp})"]
    units = [r"\unit{cm}",r"\unit{ergs s}^{-1}",r"\unit{dyne cm}^{-2}",r"\unit{K}"]
    header = r"& Initial & Initial $c$ & Converged & Final $c$ \\"
    template = r"""
    \begin{tabular}{l|l|l|l|l}
    %(header)s \hline \hline
    %(table)s
    \end{tabular}
    """
    rowTemplate = r"$%(name)s$ & $%(start)s$ & $%(sconv)s$ & $%(end)s$ & $%(econv)s$\\"
    rows = []
    for name,unit,start,end,sconv,econv in zip(names,units,start_fp,end_fp,start_c,end_c):
        fdict = {
            'name': name,
            'start': make_sci_tex(start,unit),
            'end': make_sci_tex(end,unit),
            'sconv': make_sci_tex(sconv),
            'econv': make_sci_tex(econv)
        }
        rows.append(rowTemplate % fdict)
    fdict = {
        'header' : header,
        'table' : "\n".join(rows),
    }
    return template % fdict
    
    
def macros(config):
    """Simulation paramter macros"""
    from astropysics.constants import Rs, Lsun, Ms
    template = r"\newcommand{\%(macro)s}{%(value)s}"
    macros = {
        'fittingpoint': (config["Star.Integration.fp"] * config["Star.Properties.M"],"M_\odot"),
        'starmass': (config["Star.Properties.M"],"M_\odot"),
        'tolparam': (config["System.NewtonRapson.tol"],""),
        'initaldm': (config["Star.Initial.dm"] * config["Star.Properties.M"],"M_\odot"),
        'fdjacdx' : (config["System.NewtonRapson.Jac.dx"],"")
    }
    templates = [ template % dict(macro=macro,value=make_sci_tex(*args)) for macro,args in macros.iteritems() ]
    return "\n".join(templates)
        
    
    
    
def make_sci_tex(value,unit=None):
    """Return the math-mode TeX for a number in scientific notaiton."""
    vsign = np.sign(value)
    esign = np.sign(np.log10(np.abs(value)))
    exp = np.floor(np.log10(np.abs(value)))
    
    val = value / (np.power(10,exp))
    
    sciNotation = r"%(val).3f \times 10^{%(exp)d}"
    
    if unit is not None:
        sciNotation += unit
    else:
        unit = ""
    if np.abs(exp) < 3:
        return "%.3f %s" % (value,unit)
    else:
        return sciNotation % dict(exp=exp,val=val)