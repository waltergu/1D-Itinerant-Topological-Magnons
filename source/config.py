from HamiltonianPy import *

__all__=['name','nneighbour','parametermap','idfmap','tbaidfa','tbaidfb','t','sd','dt','Us','Ud','S1','S2x']

# The configs of the model
name='1DIF'
nneighbour=2

# parametermap
parametermap=None

# idfmap
idfmap=lambda pid: Fermi(atom=(pid.site+1)%2,norbital=1,nspin=2,nnambu=1)
tbaidfa=lambda pid: Fermi(atom=(pid.site+0)%2,norbital=1,nspin=1,nnambu=1)
tbaidfb=lambda pid: Fermi(atom=(pid.site+1)%2,norbital=1,nspin=1,nnambu=1)

# terms
t=lambda **parameters: Hopping('t',parameters['t'],neighbour=2,atoms=[0,0])
sd=lambda **parameters: Hopping('sd',parameters['sd'],neighbour=1,modulate=True)
dt=lambda **parameters: Onsite('dt',parameters['dt'],atoms=[1,1],modulate=True)
Us=lambda **parameters: Hubbard('Us',parameters['Us'],atom=0,modulate=True)
Ud=lambda **parameters: Hubbard('Ud',parameters['Ud'],atom=1,modulate=True)

# cluster
S1,S2x=Square('S1'),Square('S2x')
