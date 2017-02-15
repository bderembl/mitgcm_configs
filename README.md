# mitgcm_configs

Set of MITgcm runs

## install MITgcm from scratch

* Visit http://mitgcm.org/download/ and copy the name of the latest available version: MITgcm_XXXX.tar.gz

* Dowload MITgcm in a local directory (change XXXX):
```
wget http://mitgcm.org/download/MITgcm_XXXX.tar.gz
```

* untar the source
```
tar xvf MITgcm_XXXX.tar.gz
mv MITgcm_XXXX MITgcm
cd MITgcm
```

* get mitgcm_configs
```
git clone https://github.com/bderembl/mitgcm_configs
cd mitgcm_configs
```

## Running a configuration

* In each configuration (e.g. for 'corner' here)

```
cd corner	
mkdir build
mkdir run
```

* Build mitgcm

At this point you need a fortran compiler (and openmpi for big runs)
Iin most HPC configurations, you will need to load the appropriate module. Check
```
module avail
```
and load the compiler you want (something like):
```
module load gnu
module load openmpi
```



```
cd build
../../../tools/genmake2 -mods=../code (-mpi)
make depend
make -j4
```
(Use the '-mpi' option for parallel runs)


* Create init files 
```
cd ../input
python mygendata.py
```

* Run the model
```
cd ../run
ln -s ../input/* .
cp ../build/mitgcmuv .
./mitgcmuv
```
