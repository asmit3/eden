# source psana 
from psana import *
import matplotlib.pyplot as plt
import numpy as np
from xfel.cxi.cspad_ana import cspad_tbx

ds = DataSource('exp=mfxls4916:run=260:smd')
#det=Detector('FEE_Spec')
det_jun = Detector('MfxEndstation.0:Jungfrau.0')
det_ray = Detector('MfxEndstation.0:Rayonix.0')

all_waves = []

ts_ray=[]
ts_jun=[]

counter_ray = 0
counter_jun = 0


ts_to_save = []

if False:
  with open('timestamps_to_dump.dat','r') as fin:
    for line in fin:
      if line !='\n':
        ts_to_save.append(line.split()[0].strip())


for nevent,evt in enumerate(ds.events()):
    img_ray = det_ray.image(evt)
    img_jun = det_jun.image(evt)
    ts=cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt))
    xfel_ts = ts[0:4] + ts[5:7] + ts[8:10] + ts[11:13] + ts[14:16] + ts[17:19] + ts[20:23]
    #if xfel_ts in ts_to_save:
    print ('EVENT_NUMBER=',nevent)
    if True:
      
      if img_ray is not None:
        counter_ray +=1
        ts_ray.append(ts)
        #if np.max(img_ray) >0:
        #  from IPython import embed; embed(); exit()
      if img_jun is not None:
        from IPython import embed; embed(); exit()
        counter_jun +=1
        ts_jun.append(ts)
      if True:
        wavelength = cspad_tbx.evt_wavelength(evt)
        if wavelength is None:
          continue
        else:
        #print "Energy:", 12398.4/wavelength
        #from IPython import embed; embed(); exit()
        #print "Wavelength:", nevent, wavelength
          all_waves.append(12398.4187/wavelength)
#    from IPython import embed; embed(); exit()
    #plt.imshow(img,vmin=-2,vmax=2)
    #plt.show()

from IPython import embed; embed(); exit()

print 'Done.'
