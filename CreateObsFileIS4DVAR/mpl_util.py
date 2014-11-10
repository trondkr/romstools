import matplotlib
from numpy import ma
import pylab as pl



"""A set of utility functions for matplotlib"""

# ------------------
# Plot land mask
# ------------------

def landmask(M, color='0.8'):

   # Make a constant colormap, default = grey
   constmap = pl.matplotlib.colors.ListedColormap([color])

   jmax, imax = M.shape
   # X and Y give the grid cell boundaries,
   # one more than number of grid cells + 1
   # half integers (grid cell centers are integers)
   X = -0.5 + pl.arange(imax+1)
   Y = -0.5 + pl.arange(jmax+1)

   # Draw the mask by pcolor
   M = ma.masked_where(M > 0, M)
   pl.pcolor(X, Y, M, shading='flat', cmap=constmap)

# -------------
# Colormap
# -------------

# Colormap, smlgn. med Rob Hetland

def LevelColormap(levels, cmap=None):
    """Make a colormap based on an increasing sequence of levels"""
    
    # Start with an existing colormap
    if cmap == None:
        cmap = pl.get_cmap()

    # Spread the colours maximally
    nlev = len(levels)
    S = pl.arange(nlev, dtype='float')/(nlev-1)
    A = cmap(S)

    # Normalize the levels to interval [0,1]
    levels = pl.array(levels, dtype='float')
    L = (levels-levels[0])/(levels[-1]-levels[0])

    # Make the colour dictionary
    R = [(L[i], A[i,0], A[i,0]) for i in xrange(nlev)]
    G = [(L[i], A[i,1], A[i,1]) for i in xrange(nlev)]
    B = [(L[i], A[i,2], A[i,2]) for i in xrange(nlev)]
    cdict = dict(red=tuple(R),green=tuple(G),blue=tuple(B))

    # Use 
    return matplotlib.colors.LinearSegmentedColormap(
        '%s_levels' % cmap.name, cdict, 256)
