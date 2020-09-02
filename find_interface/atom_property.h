//function defination for atom properties

#define BINSIZE_DEFALUT 12.0
#define RCUTOFF_DEFAULT 12.0
#define max_neighbour 500
#define NTYPE 3

typedef struct systeminfo{
    float boxmd[3],halfboxmd[3],cellsize[3];
    int Natoms;
    int ng[3];
    float binsize;              //binsize is minimum linklist cellsize
    float rcutoff,rcutoffsq;   //cutoff distance for atoms
    long int *llst;
    long int ***lshd;
    int **neigh_info;
} a_systeminfo;

typedef struct coodrinates{
    int atype;
    char aname[2];
    float loc[3];            //x,y,z coordinate
    int coordination;
    int property;
    int property2;
} a_coodrinates;

extern int read_input(char *filename,a_systeminfo *mdatom_info,a_coodrinates **retval);
extern int write_coordinate(a_coodrinates *input_atoms,char *filename,int Natoms,float *box);
extern int setsystemparameter(a_systeminfo *mdatom_info);
extern int makelinkedlist(a_coodrinates *input_atoms,int Natoms,float *cellsize,int *ng,long int *llst,long int ***lshd);
extern int makeneighboutlist(a_coodrinates *input_atoms,int Natoms,float *boxmd,float *halfboxmd,
           float *cellsize,int *ng,float rcutoffsq,long int *llst,long int ***lshd,int ***retval);
extern int makeneighboutlist_peratom(a_coodrinates *input_atoms,int iatom,int Natoms,float *boxmd,float *halfboxmd,
                      float *cellsize,int *ng,float rcutoffsq,long int *llst,long int ***lshd,int **retval);

