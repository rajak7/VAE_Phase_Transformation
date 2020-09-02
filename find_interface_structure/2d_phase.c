#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atom_property.h"

int cal_distance2d(int iatom, int jatom,a_coodrinates *atoms,float *boxmd,float *halfboxmd,float *retval){
    float rsq=0.0;
    float dr[2]={0.0,0.0};

    for(int kk=0;kk<2;kk++){
	dr[kk] = atoms[iatom].loc[kk] - atoms[jatom].loc[kk];
	if (dr[kk] > halfboxmd[kk]) {dr[kk] -= boxmd[kk];}
	if (dr[kk] < -halfboxmd[kk]) {dr[kk] += boxmd[kk];}
	rsq+=dr[kk]*dr[kk];
	retval[kk]=dr[kk];
    }
    retval[2]=sqrt(rsq);
    return(0);
}

int cmpfunc (const void * a, const void * b) {
     return ( *(float*)a - *(float*)b );
}

int cal_angle(a_coodrinates *atoms,int iatom,int *nlist,float *boxmd,float *halfboxmd){

    int jatom,katom,count=0;
    int N_top=0,N_bottom=0;
    int toplist[3],botlist[3];
    float dis2d_ij[3],dis2d_ik[3];
    float angleval[9],costheta_ijk,top_val;
    float threshold = 20.0;

    if (nlist[0] == 6){
	 //calculate top and bottom list
	 for(int ii=1;ii<=nlist[0];ii++){
	    jatom = nlist[ii];
	    if(atoms[jatom].loc[2] > atoms[iatom].loc[2]){
	        N_top +=1;
		if(N_top >3){ return 2; }
		toplist[N_top-1] = jatom;
	    }else{
	        N_bottom+=1;
		if(N_bottom >3){ return 2; }
		botlist[N_bottom-1] = jatom;
	    }
	 }
	 // calculate angle matrix
	 for(int ii=0;ii<3;ii++){
	    jatom = toplist[ii];
	    cal_distance2d(iatom,jatom,atoms,boxmd,halfboxmd,dis2d_ij);
	    for(int jj=0;jj<3;jj++){
		katom = botlist[jj];
		cal_distance2d(iatom,katom,atoms,boxmd,halfboxmd,dis2d_ik);
		top_val = dis2d_ij[0]*dis2d_ik[0] + dis2d_ij[1]*dis2d_ik[1];
		costheta_ijk = top_val/(dis2d_ij[2]*dis2d_ik[2]);
		if(costheta_ijk <= -1.00) {costheta_ijk += 0.001;}
		if (costheta_ijk >= 1.00) {costheta_ijk -= 0.001;}
		angleval[count++] = 180.0*acos(costheta_ijk)/3.14;
	    }
	 }
	 qsort(angleval,sizeof(angleval)/sizeof(angleval[0]),sizeof(angleval[0]),cmpfunc);
	 if (angleval[0]+angleval[1]+angleval[2] < threshold) {
		 return 0;  //2H
	 }else {
		 return 1; //1T
	 }
    }else{
         return 2;  //defect
    }
}
int main(int argc, char* argv[]){
    char *filename;
    a_systeminfo mdatom_info;
    a_coodrinates *input_atoms;
    int *nlist,*nlist2;
    int jatom,itag=0;
    int nbeta=0,nalpha=0,ndefect=0;

    filename=argv[1];
    read_input(filename,&mdatom_info,&input_atoms);
    makelinkedlist(input_atoms,mdatom_info.Natoms,mdatom_info.cellsize,mdatom_info.ng,mdatom_info.llst,mdatom_info.lshd);
    fprintf(stdout,"Total number of atoms %10d \n",mdatom_info.Natoms);
    fprintf(stdout,"Box size %12.6f %12.6f %12.6f \n",mdatom_info.boxmd[0],mdatom_info.boxmd[1],mdatom_info.boxmd[2]);

    //alpha boundary
    for(int i=0;i<mdatom_info.Natoms;i++){
	     if (input_atoms[i].atype ==1  && input_atoms[i].property2 == 3){
		     makeneighboutlist_peratom(input_atoms,i,mdatom_info.Natoms,mdatom_info.boxmd,mdatom_info.halfboxmd,
				     mdatom_info.cellsize,mdatom_info.ng,mdatom_info.rcutoffsq,
				     mdatom_info.llst,mdatom_info.lshd,&nlist);
		     if (nlist[0] == 5 || nlist[0] == 7) {
			 itag=0;
			 for (int j=1; j <= nlist[0]; j++){
			      jatom = nlist[j];
			      makeneighboutlist_peratom(input_atoms,jatom,mdatom_info.Natoms,mdatom_info.boxmd,mdatom_info.halfboxmd,
					       mdatom_info.cellsize,mdatom_info.ng,mdatom_info.rcutoffsq,
					       mdatom_info.llst,mdatom_info.lshd,&nlist2);
			      if (nlist2[0] == 3){
				 itag += 1;
			      }
			      free(nlist2);
			 }
			 if (itag >=4){  // atleast 4 se are 3 coodrinated
			     input_atoms[i].property2 = 4;
			     for (int j=1; j <= nlist[0]; j++){
				 jatom = nlist[j];
				 input_atoms[jatom].property2 = 4;
			     }
			 }
		     }
		     free(nlist);
	     }
    }
    //beta boundary
    for(int i=0;i<mdatom_info.Natoms;i++){
        if (input_atoms[i].atype ==2  && input_atoms[i].property2 == 3){
		makeneighboutlist_peratom(input_atoms,i,mdatom_info.Natoms,mdatom_info.boxmd,mdatom_info.halfboxmd,
				mdatom_info.cellsize,mdatom_info.ng,mdatom_info.rcutoffsq,
				mdatom_info.llst,mdatom_info.lshd,&nlist);
		if (nlist[0] == 4) {
		    itag=0;
		    for (int j=1; j <= nlist[0]; j++){
			 jatom = nlist[j];
			 makeneighboutlist_peratom(input_atoms,jatom,mdatom_info.Natoms,mdatom_info.boxmd,mdatom_info.halfboxmd,
					 mdatom_info.cellsize,mdatom_info.ng,mdatom_info.rcutoffsq,
					 mdatom_info.llst,mdatom_info.lshd,&nlist2);
			 if (nlist2[0] == 6 || nlist2[0] == 7 ){
		             itag += 1;
			 }
			 free(nlist2);
		    }
		    if (itag ==4){
			input_atoms[i].property2 = 5;
			for (int j=1; j <= nlist[0]; j++){
			     jatom = nlist[j];
			     input_atoms[jatom].property2 = 5;
			}
		    }
		}
		free(nlist);
	}
    }
    for(int i=0;i<mdatom_info.Natoms;i++){
	 if (input_atoms[i].atype ==1  && input_atoms[i].property2 == 3){
            ndefect +=1;
	 }else if (input_atoms[i].atype ==1  && input_atoms[i].property2 == 4){
	    nalpha+=1;
	 }else if (input_atoms[i].atype ==1  && input_atoms[i].property2 == 5){
	    nbeta+=1;
	 }
    }
    fprintf(stdout,"alpha %6d beta %6d defects %6d total %6d \n",nalpha,nbeta,ndefect,nalpha+nbeta+ndefect);
    write_coordinate(input_atoms,"output.xyz",mdatom_info.Natoms,mdatom_info.boxmd);
    return(0);
}
