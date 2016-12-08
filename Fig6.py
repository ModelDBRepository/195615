# code by sam neymotin & ernie forzano
from neuron import h
h.load_file("stdrun.hoc")
from pylab import *
import sys
import pickle
import numpy
h.install_vecst() # for samp and other NQS/vecst functions
from conf import *
import os

# check if data exists - if not, run simulations, saving output data
dataDict = {'SPI6': 'SPI6','morph':'PTcell'}
for k, v in dataDict.items() :
  if not os.path.exists('data/'+ k +'.pkl'):
    os.system('python sim.py ' + v +'.cfg')

ion()

rcParams['lines.markersize'] = 15
rcParams['lines.linewidth'] = 4
rcParams['font.size'] = 45
tl = tight_layout

#
def drtxt (ax,lett,tx=-0.075,ty=1.03,fsz=45): text(tx,ty,lett,fontweight='bold',transform=ax.transAxes,fontsize=fsz)

def naxbin (ax,nb): ax.locator_params(nbins=nb);

lmodel = ['Detailed', 'Simple']
lmodelpath = ['data/morph.pkl','data/SPI6.pkl']
lclr = ['r', 'b']

dmod = {}
for m,p in zip(lmodel,lmodelpath): dmod[m] = pickle.load(open(p))

# determine config file name
def setfcfg ():
  fcfg = "PTcell.cfg" # default config file name
  for i in xrange(len(sys.argv)):
    if sys.argv[i].endswith(".cfg") and os.path.exists(sys.argv[i]):
      fcfg = sys.argv[i]
  #print "config file is " , fcfg
  return fcfg

fcfg=setfcfg() # config file name
dconf = readconf(fcfg)
dprm = dconf['params']
dfixed = dconf['fixed']
sampr = dconf['sampr'] # sampling rate
I = numpy.load(dconf['lstimamp'])
evolts = numpy.load(dconf['evolts']) # experimental voltage traces
tt = numpy.array(dmod[lmodel[0]]['vt'])
tte = linspace(0, 1e3*evolts.shape[0]/sampr, evolts.shape[0])

#
def indexof (a,f):
  for i,val in enumerate(a):
    if abs(val-f) < 0.01: return i
  return -1

mytstop = dconf['tstop'] 
mybase = dconf['baset'] 
stimdel = dconf['stimdel']
stimdur = dconf['stimdur']
evolts = numpy.load(dconf['evolts']) # experimental voltage traces

ISubth = I[0:6]
ISup = I[6:] # subthresh right before threshold & superthreshold traces 

#
def drawsuptraces ():
  tx,ty=-.05,1.02; offy = amin(tt[0]) - 30
  cdx = 1
  ncol = len(lmodel) + 1
  ax=subplot(1,ncol,cdx); ax.set_xticks([]); ax.set_yticks([]);
  plot([1420,1520],[590,590],'k',linewidth=12)
  plot([1520,1520],[580,590],'k',linewidth=12)
  ypos = offy 
  ltxt = ['','a','b','c','d']
  for j,i in enumerate(ISup):
    idx = indexof(I,i)
    plot(tte,evolts[:,idx] + ypos,'k')
    if j > 0: ypos += 95
    else: ypos += 15  
  text(tx,ty,ltxt[cdx],fontweight='bold',transform=ax.transAxes); 
  title('Experiment',fontweight='bold')
  cdx += 1
  for mdx,m in enumerate(lmodel):
    ax=subplot(1,ncol,cdx); title(m,fontweight='bold')
    text(tx,ty,ltxt[cdx],fontweight='bold',transform=ax.transAxes);
    ypos = offy 
    for j,i in enumerate(ISup):
      plot(tt, dmod[m][i] + ypos,lclr[mdx])
      if j > 0: ypos += 95
      else: ypos += 15
    cdx += 1
    ax.set_xticks([]); ax.set_yticks([]);
  for i in xrange(1,ncol+1,1):
    ax=subplot(1,ncol,i)
    xlim((400,1600)); 
    ylim((-125,600)); 

drawsuptraces()
