# source psana 
from psana import *
import matplotlib.pyplot as plt
import numpy as np
from xfel.cxi.cspad_ana import cspad_tbx

ds = DataSource('exp=mfxlx4519:run=196:smd')
#det=Detector('FEE_Spec')
det_ray = Detector('MfxEndstation.0:Rayonix.0')

for nevent,evt in enumerate(ds.events()):
    img = det_ray.image(evt)
    plt.imshow(img,vmin=-100,vmax=100)
    plt.show()

from IPython import embed; embed(); exit()

