import numpy as np
import pickle 

class dataset():
    def __init__(self,feature_vec,labels,coordinate):
        self.feature_vec = feature_vec
        self.labels = labels
        self.coordinate = coordinate
        self.elements = feature_vec.shape[0]
        self.Ngrids = feature_vec.shape[1]
        self.Nchannels = feature_vec.shape[3]
    def print_info(self):
        print("Number of elements: ",self.elements)
        print("Ngrids: ",self.Ngrids)
        print("Channels: ",self.Nchannels)


#read training data
def readtraining_data(file1,file2,file3):
    print("files to read: ",file1,file2,file3)
    train_label=np.load(file2)
    train_corr=np.load(file3)
    for ii in range(0,len(train_label)):
        if(train_label[ii] == 4):
            train_label[ii] = 3
        elif(train_label[ii] == 5):
            train_label[ii] = 4
    #    elif(train_label[ii] == 2):  
    #        train_label[ii] = 3
    #    elif(train_label[ii] == 5): 
    #        train_label[ii] = 3
    with open(file1, 'rb') as infile:
        train_dXX = pickle.load(infile,encoding='bytes')
    for keys in train_dXX:
        print(keys,train_dXX[keys].shape)
    train_XX=train_dXX['data'].transpose((0,2,3,1)).astype('float')
    del train_dXX
    print(type(train_XX),type(train_label),type(train_corr)) 
    print("rr",train_XX.shape)
    print("rr",train_corr.shape)
    print("rr",train_label.shape)
    train_data = dataset(train_XX,train_label,train_corr)
    return train_data

def mergedata(data1,data2):
    train_XX=np.concatenate((data1.feature_vec[:][:],data2.feature_vec[:][:]),axis=0)
    train_label=np.concatenate((data1.labels[:][:],data2.labels[:][:]),axis=0)
    train_corr=np.concatenate((data1.coordinate[:][:],data2.coordinate[:][:]),axis=0)
    train_data = dataset(train_XX,train_label,train_corr)
    return train_data
    
def count_stats(total_class,labels,N_atoms):
    stats=[0.0]*total_class
    for val in labels:
        stats[int(val)]+=1.0
    stats=np.asarray(stats)
    print("stats: ",stats)
    return stats,np.min(stats)

def select_data(input_data,min_frac,stats):
    stat_frac=min_frac/stats
    print(stats)
    print(stat_frac)
    iselect_val=list()
    for ii in range(0,input_data.elements):
        val = np.random.random()
        itype=int(input_data.labels[ii])
        if val <= stat_frac[itype]:
            iselect_val.append(True)
        else:
            iselect_val.append(False)
    iselect_val=np.asarray(iselect_val)
    print(iselect_val.shape,input_data.labels.shape)
    train_XX=input_data.feature_vec[iselect_val][:]
    train_label=input_data.labels[iselect_val]
    train_corr=input_data.coordinate[iselect_val][:]
    red_data = dataset(train_XX,train_label,train_corr)
    print("reduced array",input_data.elements,red_data.elements)
    return red_data

def write_xyz(coordinate,label,prediction,Natoms):
    output = open("viz.xyz",'w')
    output.write(str(Natoms)+"\n")
    output.write(str(Natoms)+"\n")
    for ii in range(0,Natoms):
        itype=int(coordinate[ii][3])
        output.write("%6d %12.6f %12.6f %12.6f %6d %6d \n" % 
        (itype,coordinate[ii][0],coordinate[ii][1],coordinate[ii][2],prediction[ii],label[ii]))

