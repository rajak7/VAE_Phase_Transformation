import numpy as np

def readfile(filename):
    with open(filename,'r') as inputfile:
        Natoms=int(inputfile.readline().strip())
        ngrids = int(inputfile.readline().strip())
        nfeatures = inputfile.readline()
        position = np.empty((Natoms, 3), dtype='float32')
        property1 = np.empty(Natoms, dtype='float32')
        property2 = np.empty((Natoms,1), dtype='float32')
        tensor_tot = np.empty((Natoms,3,ngrids,ngrids),dtype='float32')
        ii=0; t_mo=-1; t_set=-1; t_seb=-1;
        count =0
        print("ngrids: ",ngrids,"Natoms: ",Natoms)
        tt=0
        for val in inputfile:
            tt+=1
            val = val.strip().split()
            #print("count: ",count)
            if count == 0:
                position[ii,:] = np.array(val[1:4]).astype('float32')
                property1[ii] = np.array(val[5]).astype('float32')
                property2[ii] = np.array(val[0]).astype('float32')
                count += 1
            elif count >0 and count <= ngrids:
                count+=1
                t_set+=1
                tensor_tot[ii][0][t_set][:]= np.asarray(val,dtype=float)
            elif count > ngrids and count <=2*ngrids:
                count+=1
                t_mo +=1
                tensor_tot[ii][1][t_mo][:] = np.asarray(val, dtype=float)
            elif count > 2*ngrids and count <=3*ngrids:
                count+=1
                t_seb+=1
                tensor_tot[ii][2][t_seb][:] = np.asarray(val, dtype=float)
                if count == 3*ngrids+1:
                    #print("done: ", ii)
                    count = 0
                    ii+=1; t_mo=-1; t_set=-1; t_seb=-1;
                    #break
        print("total lines",tt)
    return Natoms,ngrids,position,property1,property2,tensor_tot

