from readinput import *
import matplotlib.pyplot as plt
import pickle

#dir_name = ''
#file_2H_1T = '../10/H_T000182.xyz'
file_2H_1T = '2H_1T.txt'
interface_d = ['interface_defect.txt']
interface_a = ['alpha.txt']
interface_b = ['beta.txt']
surface = ['defect.txt']

#2H and 1T data set
Natoms,ngrids,position,property,property2,tensor_tot= readfile(file_2H_1T)
print(tensor_tot.shape)
print("Total number of atoms: ",Natoms)

#read surface defects
for fname in surface:
    N_temp, _, position_temp, property_temp, property2_temp, tensor_tot_temp = readfile(fname)
    Natoms += N_temp
    position = np.concatenate((position, position_temp), axis=0)
    property = np.concatenate((property,property_temp),axis=0)
    property2 = np.concatenate((property2, property2_temp), axis=0)
    tensor_tot = np.concatenate((tensor_tot, tensor_tot_temp), axis=0)
print("Total number of surface: ",N_temp)

#read interface defects
for fame in interface_d:
    N_temp, _, position_temp, property_temp, property2_temp, tensor_tot_temp = readfile(fname)
    Natoms += N_temp
    position = np.concatenate((position, position_temp), axis=0)
    property = np.concatenate((property,property_temp),axis=0)
    property2 = np.concatenate((property2, property2_temp), axis=0)
    tensor_tot = np.concatenate((tensor_tot, tensor_tot_temp), axis=0)
print("Total number of interface defect: ",N_temp)

Nalpha=0
#read interface alpha
for fname in interface_a:
    N_temp, _, position_temp, property_temp, property2_temp, tensor_tot_temp = readfile(fname)
    Natoms += N_temp
    Nalpha +=N_temp
    position = np.concatenate((position, position_temp), axis=0)
    property = np.concatenate((property,property_temp),axis=0)
    property2 = np.concatenate((property2, property2_temp), axis=0)
    tensor_tot = np.concatenate((tensor_tot, tensor_tot_temp), axis=0)

print("Total number of alpha: ",Nalpha)

Nbeta=0
#read interface beta
for fname in interface_b:
    N_temp, _, position_temp, property_temp, property2_temp, tensor_tot_temp = readfile(fname)
    Natoms += N_temp
    Nbeta += N_temp
    position = np.concatenate((position, position_temp), axis=0)
    property = np.concatenate((property,property_temp),axis=0)
    property2 = np.concatenate((property2, property2_temp), axis=0)
    tensor_tot = np.concatenate((tensor_tot, tensor_tot_temp), axis=0)

print("Total number of beta: ",Nbeta)
print("Total number of atoms: ",Natoms)
print("position",position.shape)
print("property",property.shape)
print("property2",property2.shape)
print("tensor_tot",tensor_tot.shape)

features={'data':tensor_tot}
with open('train_XX.p', 'wb') as infile:
    pickle.dump(features, infile, -1)

total_atoms = np.concatenate((position,property2),axis=1)
np.save("train_pos",total_atoms)
np.save("train_YY",property)
for keys in features:
    print(keys,features[keys].shape)

#plt.imshow(tensor_tot[5000][1][:][:],cmap=plt.cm.cool)
#plt.show()

