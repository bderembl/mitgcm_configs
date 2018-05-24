# mitgcm_configs

Set of MITgcm runs

## install MITgcm from scratch



* Dowload MITgcm in a local directory 
```
git clone https://github.com/MITgcm/MITgcm.git
cd MITgcm
echo "export MITGCMDIR=$PWD" >> ~/.bashrc
echo "export PATH=\$PATH:\$MITGCMDIR/tools" >> ~/.bashrc
```

* In a separate repository get mitgcm_configs
```
git clone https://github.com/bderembl/mitgcm_configs
cd mitgcm_configs
```

## Run a configuration

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
genmake2 -mods=../code (-mpi)
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
