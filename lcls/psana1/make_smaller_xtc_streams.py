from __future__ import absolute_import, division, print_function
import glob
import os
import shutil
import stat

import numpy as np
import os
from xfel.cxi.cspad_ana import cspad_tbx

message = """ Script to write out a smaller xtc file based on provided list of timestamps
              Initial script courtesy of Christopher O'Grady [SLAC]"""

class Dgram:
  def __init__(self,f):
    headerwords = 10
    self._header = np.fromfile(f,dtype=np.uint32,count=headerwords)
    self._xtcsize = 20
    self._payload = np.fromfile(f,dtype=np.uint8,count=self.extent()-self._xtcsize)
    #print ('payload',self.extent(),len(self._payload))
  def clocklow(self): return self._header[0]
  def clockhigh(self): return self._header[1]
  def tslow(self): return self._header[2]&0xffffff
  def transitionId(self): return (self._header[2]>>24)&0x1f
  def tshigh(self): return self._header[3]
  def env(self): return self._header[4]
  def dmg(self): return self._header[5]
  def srclog(self): return self._header[6]
  def srcphy(self): return self._header[7]
  def contains(self): return self._header[8]
  def extent(self): return self._header[9]
  def next(self): return self.extent()+self._xtcsize
  def data(self): return self._header
  def write(self,outfile):
      self._header.tofile(outfile)
      self._payload.tofile(outfile)


if __name__=='__main__':
  '''
  Note these are cctbx timestamps. Look at corresponding high/low timestamps below
  that are read by psana
  timestamps_indexe=[
  '20180501143548650',
  '20180501143551715',
  '20180501143555114',
  '20180501143602713',
  '20180501143625171',
  '20180501143628702',
  '20180501143632300',
  '20180501143643662',
  '20180501143701853']
  '''
  # Timestamps indexed by FFT1D using FEE wavelength for run 222
  timestamps_indexed = [
    (1525185348, 650044434),
    (1525185351, 715818581),
    (1525185355, 114751564),
    (1525185362, 713026761),
    (1525185385, 171252956),
    (1525185388, 702989123),
    (1525185392, 300990237),
    (1525185403, 662398118),
    (1525185421, 853941994)
  ]

  outfile = open('junk.xtc','w')
  for ii in range(1):
    infile = open('/reg/d/psdm/mfx/mfxls4916/xtc/e1182-r0222-s8%d-c00.xtc'%ii,'r')
    counter = 0
    while True:
      try:
        dg = Dgram(infile)
        counter +=1
        if dg.transitionId() !=12:
          print ('DG stuff:: not transitionId 12  = ',dg.transitionId(),dg.clocklow(),dg.clockhigh())
          dg.write(outfile)
        if (dg.clockhigh(),dg.clocklow()) not in timestamps_indexed: continue
        print ('DG stuff = ',dg.transitionId(),dg.clocklow(),dg.clockhigh())
        dg.write(outfile)
      except Exception as e:
        print ('Stream %d has %d frames = '%(ii, counter))
        break
    infile.close()
  outfile.close()
