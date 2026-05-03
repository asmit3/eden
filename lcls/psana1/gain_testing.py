import psana
import numpy as np
# Calculates gain and then prints the sum of the gain for Jungfrau detector.


def get_gain(ds,det):
  for nevent,evt in enumerate(ds.events()):
    if nevent>=1: break
    env = ds.env()
    src = 'MfxEndstation.0:Jungfrau.0'
    det = psana.Detector(src,env)
    gain = det.gain(evt)
    return gain

ds1 = psana.DataSource('exp=xpptut15:run=430:smd')
det1 = psana.Detector('MfxEndstation.0:Jungfrau.0')
print 'Sum of all gain from XPPTUT15, run=430 is ',np.sum(get_gain(ds1,det1))
ds2 = psana.DataSource('exp=mfxls4916:run=72:smd')
det2 = psana.Detector('MfxEndstation.0:Jungfrau.0')
print 'Sum of all gain from MFXLS4916, run=72 is ',np.sum(get_gain(ds2,det2))
