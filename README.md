# MAGSPub

## Description
Public version of MAGSPEED project. MAGSPEED is a magnetic anomaly (MA)
navigation system that derives a velocity estimate from anomaly field 
gradients. MAGSPEED development was funded by the NAVSEA contract 
N68335-22-C-0599.

## Dependencies

* [NumPy](https://www.numpy.org)
* [Matplotlib](https://www.matplotlib.org) - if you want to run notebook(s) 
  or visualization.

## Installation
To do ...

## Quickstart
See the [HOWTO](https://github.com/navigationis/MAGSPub/blob/main/MAGSim_HOWTO.ipynb)
Jupyter notebook for a few more details. Very basic use is shown below:
    
    >>> import magsim
    >>> field_model = magsim.TabulatedFieldModel('data/field-3m.npy')
    >>> sim = magsim.MAGSim(field_model=field_model)


