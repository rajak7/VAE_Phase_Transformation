#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atom_property.h"
#include "createfeature.h"


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


int cal_tensor(a_coodrinates *atoms,int iatom,int *nlist,float *boxmd,float *halfboxmd,float *gridsize){

    int itype,jtype,jatom;
    float dis2d_ij[3];
    float **tensor_mo,**tensor_set,**tensor_seb;
    float temp_xx,temp_yy,temp;
    itype=atoms[iatom].atype;

    //fprintf(stdout,"grid size %12.6f \n",*gridsize);
    //fprintf(stdout,"atom info %3d \n",itype);

    tensor_mo=(float **)malloc(ngrids*sizeof(float *));
    tensor_set=(float **)malloc(ngrids*sizeof(float *));
    tensor_seb=(float **)malloc(ngrids*sizeof(float *));
    for(int ii=0;ii<ngrids;ii++){
        tensor_mo[ii]=(float *)malloc(ngrids*sizeof(float));
        tensor_set[ii]=(float *)malloc(ngrids*sizeof(float));
        tensor_seb[ii]=(float *)malloc(ngrids*sizeof(float));
        temp_yy=ii*(*gridsize)-t_length;
        for(int jj=0;jj<ngrids;jj++){
            temp_xx=jj*(*gridsize)-t_length;
            temp = temp_xx*temp_xx + temp_yy*temp_yy;
            tensor_mo[ii][jj] = exp(-temp/gamma_2);
            tensor_set[ii][jj] = 0.0;
            tensor_seb[ii][jj] = 0.0;
//          tensor_mo[ii][jj] = 0.0;
        }
    }
    // add central atom information to the tensor_mo

    for(int ii=1;ii<=nlist[0];ii++){
        jatom=nlist[ii];
        jtype=atoms[jatom].atype;
        //Mo feature
        if (itype == jtype && itype ==1){
           cal_distance2d(iatom,jatom,atoms,boxmd,halfboxmd,dis2d_ij);
           if(fabsf(dis2d_ij[0]) < t_length && fabsf(dis2d_ij[1]) < t_length){
              updatetensor(tensor_mo,dis2d_ij,gridsize);
           }
        }
        // Se features
        if (itype != jtype && itype ==1){
           cal_distance2d(iatom,jatom,atoms,boxmd,halfboxmd,dis2d_ij);
           if(fabsf(dis2d_ij[0]) < t_length && fabsf(dis2d_ij[1]) < t_length){
                 if (atoms[jatom].loc[2] > atoms[iatom].loc[2]){
                    updatetensor(tensor_set,dis2d_ij,gridsize);
                 } else if (atoms[jatom].loc[2] < atoms[iatom].loc[2]){
                   updatetensor(tensor_seb,dis2d_ij,gridsize);
                 } else{
                   fprintf(stderr,"Invalid input type");
                   exit(1);
                 }
           }
        }
    }
    atoms[iatom].coordination = nlist[0];
    atoms[iatom].tensor_mo = tensor_mo;
    atoms[iatom].tensor_set = tensor_set;
    atoms[iatom].tensor_seb = tensor_seb;
    return(0);
}

int write_features(a_coodrinates *input_atoms,char *filename,int Natoms_tot,float *box,int Natoms){
    FILE *fp;
    int count=0;

    fp=fopen(filename,"w");
    fprintf(fp,"%d \n",Natoms);
    fprintf(fp,"%d \n",ngrids);
    fprintf(fp,"%d \n",ngrids*ngrids*3);

    for (int i=0;i<Natoms_tot;i++){
	if (input_atoms[i].atype ==1 && input_atoms[i].i_include == 1){
	   fprintf(fp,"%d \t %12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \t",input_atoms[i].atype,input_atoms[i].loc[0],
	   input_atoms[i].loc[1],input_atoms[i].loc[2],input_atoms[i].property,input_atoms[i].property2);
	   fprintf(fp,"\n");
	   //write se top tensor
	   for (int jj=0;jj<ngrids;jj++){
	       for (int kk=0;kk<ngrids;kk++){
		   fprintf(fp,"%12.6f ",input_atoms[i].tensor_set[jj][kk]);
	       }
	       fprintf(fp,"\n");
	   }
	   //write mo tensor
	   for (int jj=0;jj<ngrids;jj++){
	        for (int kk=0;kk<ngrids;kk++){
	            fprintf(fp,"%12.6f ",input_atoms[i].tensor_mo[jj][kk]);
		}
		fprintf(fp,"\n");
	   }
	   // write Se bot tensor
	   for (int jj=0;jj<ngrids;jj++){
	       for (int kk=0;kk<ngrids;kk++){
		    fprintf(fp,"%12.6f ",input_atoms[i].tensor_seb[jj][kk]);
	       }
	       fprintf(fp,"\n");
	   }
	   free(input_atoms[i].tensor_set);
	   free(input_atoms[i].tensor_mo);
	   free(input_atoms[i].tensor_seb);
	}
    }
    fclose(fp);
    return(0);
}

int viz_tensor(float **tensor_mo,float **tensor_set,float **tensor_seb){
    FILE  *fp1=fopen("image_mo.txt","w");
    FILE  *fp2=fopen("image_set.txt","w");
    FILE  *fp3=fopen("image_seb.txt","w");

    for (int ii=0;ii<ngrids;ii++){
	for (int jj=0;jj<ngrids;jj++){
	    fprintf(fp1,"%12.6f ",tensor_mo[ii][jj]);
	}
	fprintf(fp1,"\n");
    }
    for (int ii=0;ii<ngrids;ii++){
	for (int jj=0;jj<ngrids;jj++){
	    fprintf(fp2,"%12.6f ",tensor_set[ii][jj]);
	}
	fprintf(fp2,"\n");
    }
    for (int ii=0;ii<ngrids;ii++){
	for (int jj=0;jj<ngrids;jj++){
	    fprintf(fp3,"%12.6f ",tensor_seb[ii][jj]);
	}
	fprintf(fp3,"\n");
    }
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    return(0);
}
int updatetensor(float **tensor,float *dis2d_ij,float *gridsize){
    float temp_xx,temp_yy,temp;
    for(int k1=0; k1 < ngrids;k1++){
        temp_yy=k1*(*gridsize)-t_length;
        for(int k2=0; k2 < ngrids; k2++){
           temp_xx=k2*(*gridsize)-t_length;
           temp=(temp_xx-dis2d_ij[0])*(temp_xx-dis2d_ij[0])+(temp_yy-dis2d_ij[1])*(temp_yy-dis2d_ij[1]);
           tensor[k1][k2] += exp(-temp/gamma_2);
        }
    }
    return(0);
}


int main(int argc, char* argv[]){
    char *filename;
    a_systeminfo mdatom_info;
    a_coodrinates *input_atoms;
    int *nlist,N_moatoms=0;
    float gridsize=(2.0*t_length)/(float)ngrids;
    int H_frac=0,T_frac=0,maxlimit=10000,itag;
    float prob_val,prob=0.5;
    int cur_element=3;

    filename=argv[1];
    maxlimit=atoi(argv[2]);
    cur_element=atoi(argv[3]);
    fprintf(stdout,"Input parameters: %s %d %d\n",filename,maxlimit,cur_element);
    read_input(filename,&mdatom_info,&input_atoms);
    makelinkedlist(input_atoms,mdatom_info.Natoms,mdatom_info.cellsize,mdatom_info.ng,mdatom_info.llst,mdatom_info.lshd);
    fprintf(stdout,"Total number of atoms %10d \n",mdatom_info.Natoms);
    fprintf(stdout,"Box size %12.6f %12.6f %12.6f \n",mdatom_info.boxmd[0],mdatom_info.boxmd[1],mdatom_info.boxmd[2]);
    //write_coordinate(input_atoms,"output.xyz",mdatom_info.Natoms,mdatom_info.boxmd);

    //create tensor for all Mo atoms 

    for(int i=0;i<mdatom_info.Natoms;i++){
	   input_atoms[i].i_include = 0;
	   itag = 0;
       	   if(input_atoms[i].atype == 1){
//		itag = 1;
		if ( cur_element <= 1){
		   if (H_frac <= maxlimit && input_atoms[i].property2 == 0){
		       prob_val =((double)rand())/((double)RAND_MAX);
		       if (prob_val < prob) {H_frac += 1; itag =1;}
		   } else if (T_frac <= maxlimit && input_atoms[i].property2 == 1){
		       prob_val =((double)rand())/((double)RAND_MAX);
		       if (prob_val < prob) {T_frac += 1; itag =1;}
		   }
		} 
		else if (input_atoms[i].property2 == cur_element && cur_element >= 2){
		     prob_val =((double)rand())/((double)RAND_MAX);
		     if (prob_val <= 0.9 && N_moatoms <= maxlimit) {itag = 1;}
		}
	   }
	   if (itag == 1 && input_atoms[i].property2 == cur_element) {
	       makeneighboutlist_peratom(input_atoms,i,mdatom_info.Natoms,mdatom_info.boxmd,mdatom_info.halfboxmd,
			       mdatom_info.cellsize,mdatom_info.ng,mdatom_info.rcutoffsq,
			       mdatom_info.llst,mdatom_info.lshd,&nlist);
	        N_moatoms+=1;
		input_atoms[i].i_include = 1;
		cal_tensor(input_atoms,i,nlist,mdatom_info.boxmd,mdatom_info.halfboxmd,&gridsize);
		free(nlist);
//		viz_tensor(input_atoms[i].tensor_mo,input_atoms[i].tensor_set,input_atoms[i].tensor_seb);
//		break;
	   }
    }
    fprintf(stdout,"Total number of 2H %d and 1T %d \n",H_frac,T_frac);
    fprintf(stdout,"Total number of atoms written %d \n",N_moatoms);
    write_features(input_atoms,"dcnn_feature.txt",mdatom_info.Natoms,mdatom_info.boxmd,N_moatoms);
    return(0);
}
