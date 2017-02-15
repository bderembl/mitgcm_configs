
# mitgcm_configs
## High resolution idealized North atlantic

This configuration is an idealized mid-latitude ocean model. Initial conditions are derived from Samelson and Vallis (1997, jmr); with a slightly different scaling (w = 6e-4 m/s instead of w = 1e-4 m/s).

## Setup

```
export PATH=${PATH}:<mitgcmpath>/tools
mkdir build
mkdir run
```

Build mitgcm
```
cd build
genmake2 -mods=../code
make depend
make -j4
```

Create init files 
```
cd ../input
python mygendata.py
```

Run the model
```
cd ../run
ln -s ../input/* .
cp ../build/mitgcmuv .
./mitgcmuv
```
