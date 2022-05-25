typedef struct {
        t_grp   grp;
        t_vecR  *ref;
        t_vecR  *last;
        real    msf;
        int     nTimes;
} t_msf;

typedef struct {
        real    msf;
        int     cnt;
} t_msfTotal;

int allocMSF(t_msf *msf,int nTimes);

int initMSF1(t_atom *atoms,t_grp grp,t_msf *msf);

int initMSF2(FILE *ref,t_atom *atoms,t_grp grp,t_msf *msf);

int initMSFtotal(t_msfTotal *msfTotal);

int reinitMSF1(FILE *ref,t_atom *atoms,t_grp grp,t_msf *msf);

int reinitMSF2(FILE *ref,t_atom *atoms,t_grp grp,t_msf *msf,t_msfTotal *msfTotal);

int doMSF1(t_atom *atoms,t_box box,t_msf *msf);

int doMSF2(t_atom *atoms,t_box box,t_msf *msf);

