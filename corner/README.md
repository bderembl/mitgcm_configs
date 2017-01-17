# mitgcm_configs
## corner experiment

'''
export PATH=${PATH}:<mitgcmpath>/tools
mkdir build
mkdir run
'''

Build mitgcm
'''
cd build
genmake2 -mods=../code
make depend
make -j4
'''

Create init files 
'''
cd ../input
python mygendata.py
'''

Run the model
'''
cd ../run
ln -s ../input/* .
cp ../build/mitgcmuv .
./mitgcmuv
'''
