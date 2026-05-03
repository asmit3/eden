# source psana 
# some code taken from ~philiph/psana/jungfrau/V2/mfx11116/makeTuple2018.py
from psana import *
import matplotlib.pyplot as plt
import numpy as np
import copy

g0cut = 1<<14
g1cut = 2<<14
g2cut = 3<<14
ds = DataSource('exp=mfxls1016:run=350:smd')
det = Detector('MfxEndstation.0:Rayonix.0')
src=Source('DetInfo(MfxEndstation.0:Rayonix.0)')
for nevent,evt in enumerate(ds.events()):
    if nevent>=3771: break
    if nevent != 3770: continue
    print 'nevent = ','\n'.join([str(k) for k in evt.keys()])
    geom = det.pyda.geoaccess(evt)
    print 'Fetching event number',nevent
#    img = det.image(evt)
    rawFrame = det.raw(evt)
    #from IPython import embed; embed(); exit()
#    calibFrame=det.calib(evt).astype(np.float64) # avoid negative values
    calibFrame=evt.get(Camera.FrameV1,src).data16().astype(np.float64)
    plt.figure()
    #from IPython import embed; embed();exit()
    plt.hist(calibFrame,bins=50,alpha=0.2,facecolor='red')
    #plt.xlim(0,5000)
    #plt.hist(calibFrame[fG0])
    #plt.title('Gain State 0')
    #plt.figure()
    ##plt.hist(calibFrame[fG1])
    #plt.title('Gain State 1')
    #plt.figure()
    #plt.hist(calibFrame[fG2])
    #plt.title('Gain State 2')
    
#    plt.imshow(img,vmin=-2,vmax=2)
    plt.show()

#get_num_data(ds,1)

print 'Done.'
