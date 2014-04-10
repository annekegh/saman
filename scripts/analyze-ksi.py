#!/bin/python

"""Analyze ksi-coordinate
AG, 19 March, 2014"""

import numpy as np
# for matplotlib figure generation...
import matplotlib
matplotlib.use('Agg')


###################################################################################

from saman.ringfunc import read_coords,read_rings,read_unitcell
from saman.ringfunc import read_ksi_truncated, passing_average
from saman.ringfunc import Passing, RingGroup
import sys,glob

system = sys.argv[1]
temp = sys.argv[2]
datadir = sys.argv[3]
#list_datadir = glob.glob("../data-channels/%s_%sK/%s_*ethylene_%sK/" %(system,temp,system,temp))
list_datadir = glob.glob("%s/"%datadir)
print "list_datadir",list_datadir
assert len(list_datadir) == 1
datadir = list_datadir[0]

# ring definitions
#somedir = "../channels/ringdefinitions/"
somedir = sys.argv[4]
fn_xyz = "%s/config.xyz" %(somedir,)
fn_rings = "%s/ringatoms" %(somedir,)
#fn_rings_atomtypes = "../%s/atomtypes" %(somedir,)
fn_unitcell = "%s/cellvectors" %(somedir,)


## set datadir
#system = "sapo34"
#temp = "450"
#temp = "300"
#import glob
#list_datadir = glob.glob("../data-channels/%s_%sK/%s_*ethylene_%sK/" %(system,temp,system,temp))
#list_datadir = glob.glob("/home/an/%s_*ethylene_%sK/" %(system,temp))
#assert len(list_datadir) == 1
#datadir = list_datadir[0]
##datadir = "../data-channels/%s_450K/%s_26ethylene_450K/" %(system,system)
#
## ring definitions
#fn_xyz = "../zeoTsites/results_sam/%s/config.xyz" %system
#fn_rings = "../zeoTsites/results_sam/%s/ringatoms" %system
#fn_rings_atomtypes = "../zeoTsites/results_sam/%s/atomtypes" %system

#################################################################################

natom,atomtypes,pos = read_coords(fn_xyz)
list_rings_indices = read_rings(fn_rings)
unitcell0 = read_unitcell(fn_unitcell)

#### for testing!!!!!!!!!!!!1
#list_rings_indices = list_rings_indices[:2]

rg = RingGroup(list_rings_indices,pos,atomtypes,unitcell=unitcell0)
rg.set_ellips()
rg.set_passing()

#################################################################################
print "*"*20
print "Analyze composition"
print "*"*20
#ingredients = None
ingredients = ["Al","O1","O2","P","Si",]
rg.analyze_composition(ingredients=ingredients)
rg.print_composition()

print "*"*20


#################################################################################
print "*"*20
print "Analyze Ring Neighbours"
print "*"*20

fn_all_rings = "%s/all_ringatoms" %somedir
list_all_rings_indices = read_rings(fn_all_rings)
all_rg = RingGroup(list_all_rings_indices,pos,atomtypes)
all_rg.set_ellips()

print "nring",rg.nring
for i,ring in enumerate(rg.list_rings):
    print "ring ",i
    ring.set_neighbors(all_rg)


print "*"*20
for i,ring in enumerate(rg.list_rings):
    print "ring %i %8.3f %8.3f"%(i,min(rg.list_rings[i].ax1,rg.list_rings[i].ax2),max(rg.list_rings[i].ax1,rg.list_rings[i].ax2))

print "*"*20


#################################################################################
from saman.ringfunc import create_summarytransitions,write_summarytransitions,write_averagetransitions
from saman.ringfunc import plot_Fprofiles,plot_Fprofiles_perringtype,plot_Fprofiles_ringtypeidentical
from saman.ringfunc import write_Fprofiles

create_summarytransitions(datadir,rg,runs=np.arange(1,21))
logfile = "%s.%s.transitions.dat" % (system,temp)
write_summarytransitions(logfile,rg)

logfile = "%s.%s.transitions.average.dat" % (system,temp)
write_averagetransitions(logfile,rg)
print "-"*20
write_averagetransitions(sys.stdout,rg)
print "-"*20

#################################################################################

#plot_Fprofiles("histogram.ksi_ma.%s_%s"%(system,temp),rg,)
#plot_Fprofiles_perringtype("histogram.perringtype1.ksi_ma.%s_%s"%(system,temp),rg,)
#plot_Fprofiles_ringtypeidentical("identical.histogram.ksi_ma.%s_%s"%(system,temp),rg,)
write_Fprofiles("histogram.ksi_ma.%s_%s"%(system,temp),rg,)

    
#################################################################################
      #  import numpy.ma as ma
      #  dist_ma = ma.masked_greater(dist,5.)
      #  ksi_ma = ma.masked_array(ksi,mask=dist_ma.mask)
      #  
      #  print "masked dist:",sum(dist_ma.mask)  # this is the number of "True" values = the number of masked values
      #  print "masked ksi:",sum(ksi_ma.mask)  # this is the number of "True" values = the number of masked values
      #  
      #  # detect ring passage
      #  dksi_ma = ma.masked_array(signchange,mask=dist_ma.mask[:-1])
      #  hist,edges = np.histogram(ksi_ma[~ksi_ma.mask].ravel(),bins=50)

