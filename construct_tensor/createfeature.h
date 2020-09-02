int ngrids = 64;
float t_length=10.0; //length of the tensor
float gamma_2 = 2.0*0.2*0.2;  // width of exponetial kernel

extern int cal_distance2d(int iatom, int jatom,a_coodrinates *atoms,float *boxmd,float *halfboxmd,float *retval);
extern int updatetensor(float **tensor,float *dis2d_ij,float *gridsize);
extern int cal_tensor(a_coodrinates *atoms,int iatom,int *nlist,float *boxmd,float *halfboxmd,float *gridsize);
extern int viz_tensor(float **tensor_mo,float **tensor_set,float **tensor_seb);
extern int write_features(a_coodrinates *input_atoms,char *filename,int Natoms_tot,float *box,int Natoms);
