from __future__ import absolute_import, division, print_function
import glob
import os
import shutil
import stat

import numpy as np
import os
from xfel.cxi.cspad_ana import cspad_tbx

message = """ Script to write out a combined xtc file based on provided list of runs 
              Initial script courtesy of Christopher O'Grady [SLAC]

              Make sure you always save the stuff with datagram id not equal to 12
              as that has the configuration info
              Once xtc file is created use the following for smd/idx files
              a) xtcindex -f junk.xtc -o junk.xtc.idx
              b) smldata -f junk.xtc -o junk.smd.xtc"""

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
  def fiducials(self): return self.tshigh()&0x1ffff
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


def should_write_out_fee_data(seconds, fiducials, event_info):
  ''' A function that evaluates if the seconds, fiducials are close enough 
      to one of the entries in event_info
      Close enough criterion:
      a. difference in seconds < 5
      b. fiducials should be equal
      Function should return a bool '''
   
  for iid, evt in enumerate(event_info):
    if evt[1] !=fiducials:
      continue
    if abs(evt[0]-seconds) >= 5:
      continue
    return True
  return False

if __name__=='__main__':

  write_xrd_data=True
  write_fee_data=True
  runs = range(222, 222+1)
  all_event_info = {}
  if write_xrd_data:
    print ('Writing out XRD data')
    outfile = open('junk.xtc','w')
    for run in runs:
      local_event_info = []
      if run == 247:
        continue
      else:
        for ii in range(4):
          infile = open('/reg/d/psdm/mfx/mfxls4916/xtc/e1182-r%04d-s0%d-c00.xtc'%(run,ii),'r')
          counter = 0
          while True:
            try:
              dg = Dgram(infile)
              counter +=1
              if dg.transitionId() !=12:
                print ('DG stuff:: not transitionId 12  = ',dg.transitionId(),dg.clocklow(),dg.clockhigh())
                dg.write(outfile)
              else:
                #seconds=dg.clockhigh() 
                #fiducials=dg.fiducials() 
                local_event_info.append((dg.clockhigh(), dg.fiducials())) 
              #print ('DG stuff = ',dg.transitionId(),dg.clocklow(),dg.clockhigh(), dg.fiducials())
              dg.write(outfile)
            except Exception as e:
              print ('Stream %d has %d frames = '%(ii, counter))
              break
          infile.close()
      all_event_info[run]=local_event_info
    outfile.close()
    print ('Done writing out junk xtc file')


  # Now do the same with the s80 streams for the fee spec data
  if write_fee_data:
    print ('Writing out FEE data')
    print ('First read in junk.xtc to read in the seconds,fiducials to be printed out')
    infile=open('junk.xtc')
    event_info = []
    while False:
      try:
        dg=Dgram(infile)
        if dg.transitionId() !=12: continue 
        seconds=dg.clockhigh()
        fiducials=dg.fiducials()
        event_info.append((seconds, fiducials))
      except Exception as e:
        infile.close()
        break

    print ('Now read in the s80 stream and see which events should be written out')
    outfile = open('fee.xtc','w')
    for run in runs:
      event_info = all_event_info[run]
      for ii in range(1):
        infile = open('/reg/d/psdm/mfx/mfxls4916/xtc/e1182-r%04d-s8%d-c00.xtc'%(run,ii),'r')
        counter = 0
        while True:
          try:
            dg = Dgram(infile)
            counter +=1
            if dg.transitionId() !=12:
              print ('DG stuff:: not transitionId 12  = ',dg.transitionId(),dg.clocklow(),dg.clockhigh())
              dg.write(outfile)
              continue
            seconds=dg.clockhigh()
            fiducials=dg.fiducials()
            if not should_write_out_fee_data(seconds, fiducials, event_info): continue
            print ('DG stuff = ',dg.transitionId(),dg.clocklow(),dg.clockhigh(), dg.fiducials())
            dg.write(outfile)
          except Exception as e:
            print ('Stream %d has %d frames = '%(ii, counter))
            break
        infile.close()
    outfile.close()
