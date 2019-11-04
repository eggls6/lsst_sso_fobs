#!/bin/bash

cd src
make clean
wait
make
wait
cp lsst_mba ../test/.
wait
cd ../test
wait
./lsst_mba lsst.dat amor.com amor.obs
wait
echo 'all done'
exit 0

