import numpy as np
import matplotlib.pyplot as plt

from cctbx import geometry_restraints

def get_residual(b):
  return b.residual()
  # assert approx_equal(b.sites, [(1,2,3),(2,4,6)])
  # assert approx_equal(b.distance_ideal, 3.5)
  # assert approx_equal(b.weight, 1)
  # assert approx_equal(b.slack, 0)
  # assert approx_equal(b.limit, -1)
  # assert not b.top_out
  # assert approx_equal(b.origin_id, 0)
  # assert approx_equal(b.distance_model**2, 14)
  # assert approx_equal(b.delta, -0.241657386774)
  # assert approx_equal(b.residual(), 0.0583982925824)
  # assert approx_equal(b.gradients(),
  #   ((-0.12917130661302928, -0.25834261322605856, -0.38751391983908784),
  #    ( 0.12917130661302928,  0.25834261322605856,  0.38751391983908784)))

def cctbx_slack():
  mean = 1.95
  flat_top = 0.35
  sigma=0.05
  sigma0p02=0.02
  k0 = 1/(sigma*sigma)
  k0p02 = 1/(sigma0p02*sigma0p02)
  offset = 0.0
  distances = np.arange(0, 4, 0.01)
  penalty_3b = []
  penalty_3b0p02 = []
  penalty_3c = []

  for x in distances:
    b = geometry_restraints.bond(
      sites=[(x,0.,0.),(0.,0.,0.)],
      distance_ideal=mean,
      weight=k0)
    b0p02 = geometry_restraints.bond(
      sites=[(x,0.,0.),(0.,0.,0.)],
      distance_ideal=mean,
      weight=k0p02)
    s = geometry_restraints.bond(
      sites=[(x,0.,0.),(0.,0.,0.)],
      distance_ideal=mean,
      slack=flat_top,
      weight=k0)
    print('%5.2f %4.1f %4.1f' % (x, b.residual(), s.residual()))
    penalty_3b.append(b.residual())
    penalty_3b0p02.append(b0p02.residual())
    penalty_3c.append(s.residual())


  if False:
    penalty_3b = np.log(penalty_3b)
    penalty_3c = np.log(penalty_3c)


  plt.scatter(distances, penalty_3b, color='skyblue',label='Harmonic Restraints $\mu=1.95, \sigma=0.05$, this work',s=50)
  plt.scatter(distances, penalty_3b0p02, color='grey',label='Harmonic Restraints $\mu=1.95, \sigma=0.02$,', alpha=0.5,s=50)
  plt.scatter(distances, penalty_3c, color='wheat', label='Flat-bottom Harmonic Restraints (FBHR), $\mu=1.95, flat=0.35, \sigma=0.05$',s=50)
  #plt.yscale('log')
  plt.xlim([1.4, 2.5])
  plt.ylim([-5, 150])
  plt.xticks(fontsize=20)
  plt.yticks(fontsize=20)
  # Z=2 axvlines 
  #plt.axvline(1.99)
  #plt.axvline(1.91)
  #plt.axvline(1.85, c='blue')
  #plt.axvline(2.05, c='blue')
  #plt.legend(fontsize=16) 
  #plt.tight_layout()
  #plt.show()
  #plt.ylim([0, 300])
  #plt.savefig('restraints.pdf', format='pdf')
  plt.legend(loc='best', bbox_to_anchor=(1.04, 0.97), fontsize=20)
  plt.xlabel('Distance($\AA$)', fontsize=20)
  plt.ylabel('Residual', fontsize=20)
  #plt.savefig('SI_fig13_restraints_detail.pdf', format='pdf')
  plt.show()


if __name__ == '__main__':
  cctbx_slack()
