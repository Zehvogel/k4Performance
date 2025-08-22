# k4Performance


This repository hosts Gaudi functionals useful for performance studies


## Dependencies

* ROOT

* PODIO

* Gaudi

* k4FWCore

## Installation

Run, from the `k4Performance` directory:

``` bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
k4_local_repo
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
```

## Execute Examples

Make sure that `../install/lib` and `../install/python` are in `LD_LIBRARY_PATH`
and `PYTHONPATH` respectively (`k4_local_repo` should take care of this).
If they are not, they can be added by running:
``` bash
export LD_LIBRARY_PATH=$PWD/../install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/../install/python:$PYTHONPATH
```
and then run the examples like this:

``` bash
k4run ../RecoMCTruthLinkers/options/runCaloClusterMCParticleLinker.py
```


## References:
These could perhaps be usefule for newcomers.
1. [lhcb-98-064 COMP](https://cds.cern.ch/record/691746/files/lhcb-98-064.pdf)
2. [Hello World in the Gaudi Framework](https://lhcb.github.io/DevelopKit/02a-gaudi-helloworld)
