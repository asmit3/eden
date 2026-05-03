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
det = Detector('MfxEndstation.0:Jungfrau.0')
for nevent,evt in enumerate(ds.events()):
    if nevent>=3771: break
    if nevent !=3770: continue
    print 'nevent = ','\n'.join([str(k) for k in evt.keys()])
    geom = det.pyda.geoaccess(evt)
    print 'Fetching event number',nevent
#    img = det.image(evt)
    rawFrame = det.raw(evt)
    offset=copy.deepcopy(det.offset(evt))
    offset[1,:]*=0.01
    offset[2,:]*=-1.0
    #from IPython import embed; embed(); exit()
    calibFrame=det.calib(evt) # avoid negative values
    calibFrame_offset_fudge=det.calib(evt,offs=offset) # avoid negative values
    fG0 = rawFrame<g0cut
    fG1 = (rawFrame>=g0cut) & (rawFrame<g1cut)
    fG2 = rawFrame>=g2cut
    fGval = fG0*1 + fG1*2 + fG2*3
    #from IPython import embed; embed(); exit()
    plt.figure()
    plt.hist(calibFrame.reshape(1024*1024),2000, alpha=0.2,facecolor='red')
    plt.hist(calibFrame_offset_fudge.reshape(1024*1024),2000, alpha=0.5,facecolor='green')
    #plt.xlim(-10,2000)
    plt.legend(['Default constants from PSANA', 'Fudged offset: offset[1,:]*=0.01,offset[2,:]*=-1.0'])
    plt.xlabel('Calibrated Intensity')
    plt.ylabel('Occurance')
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
