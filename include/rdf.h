typedef struct {
	int grp1Idx;
	int grp2Idx;
	real binSize;
	int nBins;
	real *bins;
	real *bins2;
	int cnt;
	real averDens; 
} t_rdf;

int initRDF(int grp1Idx,int grp2Idx,real binSize,int nBins,t_rdf *rdf);

int getBoxVol(t_box box,real *vol);

int doRDF(t_atom *atoms,int nAtoms,t_grp *grps,t_rdf *rdf,t_box box);

int finalizeRDF(t_rdf *rdf);

