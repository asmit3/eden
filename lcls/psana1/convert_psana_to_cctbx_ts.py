# source psana 
# use this command: source /reg/g/psdm/etc/ana_env.sh 
#import matplotlib.pyplot as plt
#import numpy as np

import psana
import time

#ds = DataSource('exp=xpptut15:run=430:smd')
ds = psana.DataSource('exp=xpplv6818:run=93:smd')

for nevent,evt in enumerate(ds.events()):
    if nevent>=5: break
    evt_id=evt.get(psana.EventId)
    sec, nsec = evt_id.time()
    nsec=nsec/1e6
    cctbx_ts = time.strftime("%Y-%m-%dT%H:%MZ%S", time.gmtime(sec)) + \
      (".%03d" % nsec) 
    print (cctbx_ts)
    #from IPython import embed; embed(); exit()

#evt_id.time()[0] << 32 | evt_id.time()[1]) evt_id = evt.get(psana.EventId)
#        evt = self._get_event(index)
#        time = evt.get(psana.EventId).time()
#        # fid = evt.get(psana.EventId).fiducials()
#
#        sec = time[0]
#        nsec = time[1]
#
#        return cspad_tbx.evt_timestamp((sec, nsec / 1e6))
#  if t is None:
#    t = evt_time(evt=None)
#    if t is None:
#      return None
#  return time.strftime("%Y-%m-%dT%H:%MZ%S", time.gmtime(t[0])) + \
#      (".%03d" % t[1])
