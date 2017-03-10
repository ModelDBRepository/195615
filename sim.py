# written by sam neymotin, modified by ernie forzano
from neuron import h
h.load_file("stdrun.hoc")
from pylab import *
import time 
from time import time, clock
import datetime # to format time of run
import sys
import pickle
import numpy
h.install_vecst() # for samp and other NQS/vecst functions
from conf import *
import os

# determine config file name
def setfcfg ():
  fcfg = "PTcell.BS0284.cfg" # default config file name
  for i in xrange(len(sys.argv)):
    if sys.argv[i].endswith(".cfg") and os.path.exists(sys.argv[i]):
      fcfg = sys.argv[i]
  print "config file is " , fcfg
  return fcfg

fcfg=setfcfg() # config file name
dconf = readconf(fcfg)
dprm = dconf['params']
dfixed = dconf['fixed']
sampr = dconf['sampr'] # sampling rate

exec('from ' + dconf['cellimport'] + ' import ' + dconf['cellfunc']) # import cell creation function
if fcfg.startswith('PTcell'):
  exec('cell = ' + dconf['cellfunc'] + '(' + str(dconf['cellfuncargs']) + ')') # create the cell - can use different morphologies)
else:
  exec('cell = ' + dconf['cellfunc'] + '()') # create the cell

exec('import ' + dconf['cellimport']) # import the file's variables too, so can access them

for p in dfixed.values(): exec(p.assignstr(p.origval)) # fixed values
try:
  ic = h.IClamp(0.5,sec=cell.soma[0])
except:
  ic = h.IClamp(0.5,sec=cell.soma)


voltLoc = None # voltage location used for optimization
voltLocInterp = None

# record variables (includes voltage, spike times, etc.)
def setrec (cell):
  global voltLoc,voltLocInterp
  rd = {}
  rd['vt'] = h.Vector(); rd['vt'].record(h._ref_t) # Record simulation time
  # Record cell voltages
  for i,loc in enumerate(dconf['recordV']):
    sloc = loc + '.v' # easier to read string representation
    ploc = loc + '._ref_v' # pointer location
    if i == 0:
      voltLoc = sloc # v
      voltLocInterp = voltLoc + '_interp'
    rd[sloc] = h.Vector()
    cmd = 'rd[sloc].record(' + ploc + ')'#
    exec(cmd)
  # record spikes at specified Section(location)
  vspike = h.Vector()
  sloc = dconf['recordSpike']
  cmd = 'spikerec = h.NetCon(' + sloc + '._ref_v, None, sec=' + sloc.split('(')[0] +')'
  exec(cmd)
  spikerec.record(vspike)
  spikerec.threshold = 0 # 0 mV threshold
  rd['vspike'] = vspike
  rd['spikerec'] = spikerec
  return rd

rd = setrec(cell)
vt = rd['vt']

# look at output
def plotout (rd,lk = [voltLoc],clr = ['r','g'],drdot=False):
  tt = rd['vt'].as_numpy()
  for gn,k in enumerate(lk):
    if k.count('interp') > 0: plot(rd['vt_interp'],rd[k],clr[gn],linewidth=1);
    else: plot(tt,rd[k],clr[gn],linewidth=1);
    if drdot: plot(tt,rd[k],clr[gn]+'o',markersize=5)
    xlabel('Time (ms)',fontsize=16); ylabel('Vm',fontsize=16); xlim((0,mytstop-mybase))

if dconf['usecvode']: h.cvode.active(1) #much faster

mytstop = dconf['tstop'] 
mybase = dconf['baset'] 
stimdel = dconf['stimdel']
stimdur = dconf['stimdur']
I = numpy.load(dconf['lstimamp']) # somatic current injection levels

# interpolate voltage recorded in simulation to a fixed grid (dt millisecond spacing)
# output voltage,time is stored in rd (dictionary) with keys voltLocInterp, vt_interp
def interpvolt (rd,dt):
  tsrc,vsrc = rd['vt'], rd[voltLoc]
  tdest = h.Vector(); tdest.indgen(mybase,mytstop,dt)
  vdest = h.Vector(); vdest.interpolate(tdest,tsrc,vsrc)
  tdest.sub(mybase)
  rd['vt_interp'] = tdest
  rd[voltLocInterp] = vdest
 
# 
def myrun (tstop=mytstop,inj=0.75,draw=True,prtime=False):
  if prtime: clockStart = time()
  h.tstop = tstop 
  ic.amp = inj
  ic.delay = mybase + stimdel
  ic.dur = stimdur # BS0284 experiment has only 1 s of stim
  h.run()
  interpvolt(rd,1e3/sampr) # interpolate recorded voltage to fixed temporal grid
  if draw: plotout(rd)
  if prtime:
    print rd['vspike'].size()*1e3/ic.dur , 'Hz, during stim.'
    clockEnd = time()
    print '\nsim runtime:',str(round(clockEnd-clockStart,2)),'secs'   

# run all current injections and return output voltage 
def gathertraces ():
  dvec = {}
  for i,inj in enumerate(I):
    print i,inj
    myrun(tstop=mytstop,inj=inj,draw=False,prtime=True)
    if i == 0: dvec['vt'] = rd['vt_interp'].to_python()
    dvec[inj] = rd[voltLocInterp].to_python()
  return dvec

evolts = numpy.load(dconf['evolts']) # experimental voltage traces
  
if __name__ == "__main__":
  if len(sys.argv) < len(dprm.keys()):
    print 'usage: sim.py [cfg file] [params] [fout]'
  else:
    i = 1
    if sys.argv[i].count('python') > 0: i += 1 # skip -python or python
    if sys.argv[i].count('sim.py') > 0: i += 1 # skip sim.py
    if sys.argv[i].count('.cfg') > 0: i += 1 # skip cfg file
    for j,p in enumerate(dprm.keys()): # parameter values to evaluate - optimized by evolution
      print i,j,p,sys.argv[i]
      exec(dprm[p].assignstr(float(sys.argv[i])))
      i += 1
    if 'postassign' in dconf: exec(dconf['postassign'])
    while i < len(sys.argv): # check optional args
      if sys.argv[i].count('.cfg') > 0:
        pass # config file - skip it since read above
      elif sys.argv[i].count('python') > 0 or sys.argv[i].count('sim.py') > 0:
        pass # pass, it's call from commandline
      i += 1
  print 'running sim ... '
  dvec = gathertraces()
  pickle.dump(dvec,open(os.path.join('data',dconf['cellimport'].split('.cfg')[0]+'.pkl'),'w'))
