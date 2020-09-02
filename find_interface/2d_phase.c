#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atom_property.h"

int main(int argc, char* argv[]){
    char *filename;
    a_systeminfo mdatom_info;
    a_coodrinates *input_atoms;
    int *nlist;
    int H_phase=0,T_phase=0,jatom,thresh=3;

    filename=argv[1];
    read_input(filename,&mdatom_info,&input_atoms);
    makelinkedlist(input_atoms,mdatom_info.Natoms,mdatom_info.cellsize,mdatom_info.ng,mdatom_info.llst,mdatom_info.lshd);
    fprintf(stdout,"Total number of atoms %10d \n",mdatom_info.Natoms);
    fprintf(stdout,"Box size %12.6f %12.6f %12.6f \n",mdatom_info.boxmd[0],mdatom_info.boxmd[1],mdatom_info.boxmd[2]);

    //find interface atoms 

    for(int i=0;i<mdatom_info.Natoms;i++){
	     if (input_atoms[i].property2 ==2){
		     makeneighboutlist_peratom(input_atoms,i,mdatom_info.Natoms,mdatom_info.boxmd,mdatom_info.halfboxmd,
				     mdatom_info.cellsize,mdatom_info.ng,mdatom_info.rcutoffsq,
				     mdatom_info.llst,mdatom_info.lshd,&nlist);
		     H_phase=0; T_phase =0;
		     for (int j=1; j <= nlist[0]; j++){
			  jatom = nlist[j];
			  if ( input_atoms[jatom].property2 == 0) {H_phase +=1;}
			  else if ( input_atoms[jatom].property2 == 1) {T_phase +=1;}
		     }
		    if (H_phase > thresh && T_phase > thresh) {input_atoms[i].property2 = 3;} 
		    free(nlist);
	     }
    }
    write_coordinate(input_atoms,"output.xyz",mdatom_info.Natoms,mdatom_info.boxmd);
    return(0);
}
