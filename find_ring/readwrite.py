import re


## Read input coordinate
def readcoordinate(filename,boxsize,itype,igd,pos,deform):
    inputfile = open(filename, 'r')
    count=0
    for line in inputfile:
        count+=1
        line=line.replace('\n',"")
        if count==1: Natom=int(line)
        if count==2:
            temp=re.findall(r'\S+',line)
            for val in temp:
                boxsize.append(float(val))
        if count > 2:
            temp = re.findall(r'\S+', line)
            itype.append(int(temp[0]))
            igd.append(count-2)
            pos.append([float(temp[1]),float(temp[2]),float(temp[3])])
            deform.append(int(temp[5]))
    return Natom


def makelinklist(Natom,pos,boxsize,cellsize):

    ng=[0,0,0]
    ig=[0,0,0]

    for i in range(0,3):
        ng[i]=int(boxsize[i]/cellsize[i])
        cellsize[i]=boxsize[i]/float(ng[i])
    print("Bozsize:  ",boxsize)
    print("Cellsize: ",cellsize)
    print("Bin Size: ",ng)

    lshd=[[[-1 for k in range(ng[2])] for j in range(ng[1])] for i in range(ng[0])]
    llst=[None for i in range(0,Natom)]
    for ii in range(0,Natom):
        for jj in range(0,3):
            ig[jj]=int(pos[ii][jj]/float(cellsize[jj]))
            if (ig[jj] < 0 or ig[jj] >= ng[jj]): print("ig out of bound for ",jj,ig[jj],ng[jj])
        llst[ii]=lshd[ig[0]][ig[1]][ig[2]]
        lshd[ig[0]][ig[1]][ig[2]]=ii
    return lshd,llst,ng

def makeneighbourlist(lshd,llst,pos,boxsize,Natom,itype,cellsize,nlist,Rcutsq):

    ig = [0,0,0]
    ic1 = [0,0,0]
    ndir=[-1,0,1]              # neighbor directions
    dr=[0.0,0.0,0.0]          # distance between atoms
    halfbox=[0.5*val for val in boxsize]
    Neighbor=[[0] for j in range(0,Natom)]
    for ii in range(0,Natom):
        iatom=ii
        for jj in range(0, 3):
            ig[jj] = int(pos[iatom][jj] / float(cellsize[jj]))
        for ix in ndir:
            for iy in ndir:
                for iz in ndir:
                     ic1[0] = (ig[0] + ix + nlist[0]) % nlist[0]
                     ic1[1] = (ig[1] + iy + nlist[1]) % nlist[1]
                     ic1[2] = (ig[2] + iz + nlist[2]) % nlist[2]
                     jatom=lshd[ic1[0]][ic1[1]][ic1[2]]
                     while jatom >= 0:
                         rsq=0.0
                         if not jatom==iatom:
                             for kk in range(0,3):
                                 dr[kk]=pos[iatom][kk]-pos[jatom][kk]
                                 if (dr[kk] > halfbox[kk]) : dr[kk] -= boxsize[kk]
                                 if (dr[kk] < -halfbox[kk]): dr[kk] += boxsize[kk]
                                 rsq+=dr[kk]*dr[kk]
                             if rsq <= Rcutsq and (itype[iatom] != itype[jatom]):
                                 if(itype[iatom] == 1 and itype[jatom] ==3) or (itype[iatom] == 3 and itype[jatom] ==1):
                                     pass
                                 else:
                                    Neighbor[iatom][0]+=1
                                    Neighbor[iatom].append(jatom)
                         jatom=llst[jatom]
    return Neighbor


# Write XYZ File

def writexyz(Natom,itype,igd,pos):
    outputfile = open('ring.xyz', 'w')
    outputfile.write(str(Natom)+"\n")
    outputfile.write(str(Natom) + "\n")
    for ii in range(0,Natom):
        if itype[ii] == 1: outputfile.write("Al %12.6f %12.6f %12.6f \n" % ( pos[ii][0] , pos[ii][1] , pos[ii][2]))
        if itype[ii] == 2: outputfile.write("O  %12.6f %12.6f %12.6f \n" % ( pos[ii][0] , pos[ii][1] , pos[ii][2]))


# Write Single atom information into the file
def writesingleatom(atomid,Neighborinfo,pos,itype,ringlength):
    outputfile=open('singleatom.xyz','w')
    itag = 0
    if itype[atomid]==1:
        outputfile.write("Al %12.6f %12.6f %12.6f %3d\n" % ( pos[atomid][0] , pos[atomid][1] , pos[atomid][2],itag))
    else:
        outputfile.write("O  %12.6f %12.6f %12.6f %3d\n" % (pos[atomid][0], pos[atomid][1], pos[atomid][2],itag))
    for val in Neighborinfo[atomid]:
        if(len(val)==(ringlength-1)):
            itag+=1
            print("list: ",val,"length",len(set(val)))
            for ii in val:
                if itype[ii]==1:
                    outputfile.write("Mo %12.6f %12.6f %12.6f %3d\n" % ( pos[ii][0],pos[ii][1],pos[ii][2],itag))
                elif type[ii]==2:
                    outputfile.write("Se  %12.6f %12.6f %12.6f %3d\n" % (pos[ii][0],pos[ii][1],pos[ii][2],itag))
                else:
                    outputfile.write("W  %12.6f %12.6f %12.6f %3d\n" % (pos[ii][0], pos[ii][1], pos[ii][2], itag))

