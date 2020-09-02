from readwrite import *
from queue import Queue
import multiprocessing
import time
import sys
# Ring Analysis Code for Ceramic Sysytems

# Decleration of Variables
filename=sys.argv[1]        # Input file name
print("filename: ",filename)
Natom=0                        # Total number of atoms in the system
Rcut_al_O=3.3                 # Al-O Bond cutoff distance
Rcutsq=Rcut_al_O*Rcut_al_O     # Square of the cutoff distance
boxsize=list()                 # Simulation Box Size
itype=list()                   # Atom Type
igd=list()                     # Global Id
pos=list()                     # coordinate of the atom
deform=list()                  # local deformation parameter
lshd=list()                    # Link List of the atoms
llst=list()                    # Link list of the atoms
nlist=list()                   # Link list of the atoms
cellsize=[Rcut_al_O,Rcut_al_O,Rcut_al_O]         # Link List size
Neighbor=list()                # Neighborlist for each atom
Ringstat=list()                # Ring stastistics for each atom
MaxRingsize=4                 # Maximum size of the ring
Ringlistsize=MaxRingsize+3     # size of the list which keeps track of ring statistics
Neighborinfo=dict()            # Stores all the ids of the node which makes ring with a given atom
Alcord=6                       # Coordination number for Al
Wcord=6                         # Coordination number for W
Ocord=3                        # Coordination number for O
Alring_ref = [0, 0, 0, 0, 6, 0, 24, 0, 0, 0, 48, 0, 24]   # refrence value of rings of different size in Al
Oring_ref = [0, 0, 0, 0, 3, 0, 6, 0, 6, 0, 6, 0, 6]  # refrence value of rings of different size in O
nprocessor=8                  #total number of processors used

def ringanalysis(startid):
    exploredA=set()        #Visited neighbours from side A
    exploredB=set()        #Visited neighbours from side B
    iringtag=False
    endid=startid+1
    for ii in range(startid,endid):
#        print("atomid: ",ii)
        totneigh=Neighbor[ii][0]+1
        Neighborinfo[ii]=[]
        for Vii in range(1,totneigh-1):
            for Vjj in range(Vii+1,totneigh):
#                print("Vii,Vjj",Vii,Vjj,Neighbor[ii][Vii],Neighbor[ii][Vjj])
                exploredA = set()                   # Initialize the explored set for neighbor Vii & Vjj pair
                exploredB = set()                   # Initialize the explored set for neighbor Vii & Vjj pair
                exploredA.add(ii)                   # Add atom ii into the explored setA
                exploredB.add(ii)                   # Add atom ii into the explored setB
                exploredA.add(Neighbor[ii][Vii])    # Add Neighbor Vii into explored setA
                exploredB.add(Neighbor[ii][Vjj])    # Add Neighbor Vjj into explored setB
                sideA=Queue()                       # Initialize the Queue from side A
                sideB=Queue()                       # Initialize the Queue from side B
                sideA.put(Neighbor[ii][Vii])
                sideB.put(Neighbor[ii][Vjj])
                Rcurrent=0
                while Rcurrent < MaxRingsize:
#                    iringtag = False
                    Rcurrent+=2
                    # BFS from side A
                    sideAnew = set()
                    while sideA.empty()==False:
                        tempa1=sideA.get()                     # Take a neighbor node of ii out from sideA
                        NmychildA = Neighbor[tempa1][0] + 1    # number of neighbor node this tempa1 node has
                        for kk in range(1,NmychildA):
                            mychild = Neighbor[tempa1][kk]
                            if mychild not in exploredA:  sideAnew.add(mychild)
                    # BFS from side B
                    sideBnew = set()
                    while sideB.empty()==False:
                        tempa2 = sideB.get()                        # Take a neighbor node of ii out from sideB
                        if tempa2 in exploredA: continue            # Making sure same node is not present on both side
                        NmychildB = Neighbor[tempa2][0] + 1         # number of neighbor node this tempa2 node has
                        for kk in range(1,NmychildB):
                            mychild=Neighbor[tempa2][kk]
                            if mychild not in exploredB:
                                sideBnew.add(mychild)
                                if mychild in sideAnew:
#                                    iringtag = False
#                                    print("ring", ii, Rcurrent, Vii, Vjj)
                                    Ringstat[ii][Rcurrent + 2] += 1
                    for val in sideAnew:                    # Val is thr key of the dictionary sideAnew
                        sideA.put(val)                      # Adding the node: val into the queue sideA
                        exploredA.add(val)
                    for val in sideBnew:                   # Val is thr key of the dictionary sideBnew
                        sideB.put(val)                     # Adding the node: val into the queue sideB
                        exploredB.add(val)
    return Ringstat[startid:endid]


def writeringstat(Natom,Neighbor,pos,itype,Ringstat):
    outputfile = open('MoWSe2ring.xyz', 'w')
    outputfile.write(str(Natom) + "\n")
    outputfile.write(str(boxsize[0]) + " "+str(boxsize[1])+" 100.00 \n")
    for ii in range(0, Natom):
        if itype[ii] == 1 or itype[ii] ==3:
            ring4val = Ringstat[ii][4] - Alring_ref[4]
 #           ring6val = Ringstat[ii][6] - Alring_ref[6]
            if ring4val == 0 : class1 = deform[ii]
            else: class1 = 2
        if itype[ii] == 2:
            ring4val = Ringstat[ii][4] - Oring_ref[4]
#            ring6val = Ringstat[ii][6] - Oring_ref[6]
            if ring4val == 0 : class1 = deform[ii]
            else: class1 = 2
        if itype[ii] == 1: outputfile.write("1 %12.6f %12.6f %12.6f %3d %3d \n" % (pos[ii][0], pos[ii][1], pos[ii][2],class1,class1))
        if itype[ii] == 2: outputfile.write("2 %12.6f %12.6f %12.6f %3d %3d \n" % (pos[ii][0], pos[ii][1], pos[ii][2],class1,class1))
        if itype[ii] == 3: outputfile.write("1 %12.6f %12.6f %12.6f %3d %3d \n" % ( pos[ii][0], pos[ii][1], pos[ii][2],class1,class1))



#Program Starts from here

Natom=readcoordinate(filename,boxsize,itype,igd,pos,deform)     # Read-Coordinate
#writexyz(Natom,itype,igd,pos)

lshd,llst,nlist=makelinklist(Natom,pos,boxsize,cellsize)                          # making link list
Neighbor=makeneighbourlist(lshd,llst,pos,boxsize,Natom,itype,cellsize,nlist,Rcutsq)     # Making Neighbor list for each atom

print("Bozsize:  ", boxsize)
print("Cellsize: ", cellsize)
print("Bin Size: ", nlist)

#-------Calculation for ring statistics start from here
startid=0
ichunk=int(Natom/nprocessor)
print("Chunk size: ",ichunk,"total processers: ",nprocessor)
Ringstat=[[0 for j in range(0,Ringlistsize)] for i in range(0,Natom)]

if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=nprocessor)
    start_time = time.time()
    result=pool.imap(ringanalysis,range(0,Natom),ichunk)
    pool.close()
    pool.join()
    end_time = time.time()

print("total time taken:",end_time-start_time)
count=-1
for val in result:
    count+=1
    for val2 in val:
        Ringstat[count] = val2

writeringstat(Natom,Neighbor,pos,itype,Ringstat)
#--------------------****************-----------------------------



