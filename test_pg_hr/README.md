
# mitgcm_configs
## High resolution idealized North atlantic

This configuration is an idealized mid-latitude ocean model. Initial conditions are derived from Samelson and Vallis (1997, jmr); with a slightly different scaling (w = 6e-4 m/s instead of w = 1e-4 m/s).

## Create initial condition

The library spoisson is required to create the initial condition (needed to get the sea surface elevation from the velocity field).

in any directory:
```
git clone https://github.com/bderembl/spoisson
cd spoisson
python setup.py install --user
```
