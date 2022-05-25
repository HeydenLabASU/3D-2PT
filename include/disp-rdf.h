typedef struct {
        t_grp grp1;
        t_grp grp2;
        t_vecR *ref1;
        t_vecR *ref2;
        t_vecR *last1;
        t_vecR *last2;
        real vol;
} t_dispRDF;

typedef struct {
        real dispTot1;
        int n1;
        real dispTot2;
        int n2;
        real dens2;
        real *dispRDF;
        real binSize;
        int nBins;
        int cnt;
} t_dispRDFtotal;

int getBoxVol(t_box box,real *vol);

int initDispRDF(t_atom *atoms,t_box box,t_grp grp1,t_grp grp2,t_dispRDF *dispRDF);

int initDispRDFtotal(t_dispRDFtotal *dispRDFtotal,real binSize,int nBins);

int reinitDispRDF(t_atom *atoms,t_box box,t_grp grp1,t_grp grp2,t_dispRDF *dispRDF,t_dispRDFtotal *dispRDFtotal);

int doDispRDF(t_atom *atoms,t_box box,t_dispRDF *dispRDF);

