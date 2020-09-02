#!/bin/bash

# find 2H, 1T  structure inside the fracture data
echo "========find 2H, 1T  structure inside the fracture data========="
cd find_structure
make clean
make
./2d_phase ../fracture.xyz
mv output.xyz ../structure.xyz
make clean
cd ..

# find defect ring structure insde the 1T phase of the fracture data
echo "=========find defect ring structure insde the 1T phase of the fracture data========="
cd find_ring
python3.6 parallelmain.py ../structure.xyz
mv MoWSe2ring.xyz ../structure_ring.xyz
cd ..
make clean
rm structure.xyz

#find interfact between 2H and 1T regions insde the fracture data
echo "==========find interfact between 2H and 1T regions insde the fracture data==========="
cd find_interface
make clean
make 
./2d_phase ../structure_ring.xyz
mv output.xyz ../interface.xyz
make clean
cd .. 
rm structure_ring.xyz

#find alpha and beta structures at the interface region
echo "=======find alpha and beta structures at the interface region================="
cd find_interface_structure
make clean
make
./2d_phase ../interface.xyz
mv output.xyz ../train.xyz
make clean
cd ..
rm interface.xyz

#construct a n*64*64*3 tensor using the processed fracture frame
echo "=======construct  N*64*64*3 tensor using the processed fracture frame in numpy format========"
cd construct_tensor
limit=2000
tag1=1
tag2=2
tag3=3
tag4=4
tag5=5
make clean
make
./c_feature ../train.xyz $limit $tag1
mv dcnn_feature.txt 2H_1T.txt
./c_feature ../train.xyz $limit $tag2
mv dcnn_feature.txt defect.txt
./c_feature ../train.xyz $limit $tag3
mv dcnn_feature.txt interface_defect.txt
./c_feature ../train.xyz $limit $tag4
mv dcnn_feature.txt alpha.txt
./c_feature ../train.xyz $limit $tag5
mv dcnn_feature.txt beta.txt
echo "=======Information of the generated training dataset========"
python3.6 create_tensor.py
rm defect.txt 2H_1T.txt
mv train_XX.p train_YY.npy train_pos.npy ../data/.
make clean
cd ..

