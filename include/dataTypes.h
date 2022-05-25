#include "fftw3.h"

typedef float real;

typedef int integer;

typedef struct {
        real x,y,z;
} t_vecR;

typedef struct {
	t_vecR	c1;
	t_vecR	c2;
	t_vecR	c3;
} t_matR;

typedef struct {
        int      resNr;
        char     resName[7];
        char     atomName[7];
        int      atomNr;
        t_vecR   crd;
        t_vecR   vel;
	t_vecR   force;
        char     elem[5];
        real     mass;
        real     charge;
	real     pol;
	real     sigmaLJ;
	real     epsLJ;
	t_vecR   DM;
	int	 HBinfo[2];
} t_atom;

typedef struct {
	char     siteName[5];
	t_vecR   crd;
	real     charge;
        real     mass;
	real     pol;
	real     sigmaLJ;
	real     epsLJ;
	t_vecR   DM;
} t_site;

typedef struct {
	t_vecR	a;
	t_vecR	b;
	t_vecR	c;
	int	ortho;
	t_matR trans;
} t_box;

/* integers to be stored in the following structures
   are supposed to refer to the index (starting with 0)
   of the involved atoms
*/
typedef int t_bond[2];

/* example: H O H (water) */
typedef int t_angle[3];

/* example: H C C H  (ethan) */
typedef int t_dihed[4];

/* example: N H H H (ammonia) */
typedef int t_improp[4];

typedef struct {
	char      molName[7];
	char	  chain[5];
	int       nAtoms;
	int       *atoms;
	int       nSites;
	t_site    *sites;
	int       nBonds;
	t_bond    *bonds;
	real      *bondLen;
	int       nAngles;
	t_angle   *angles;
	real      *angleValues;
	int       nDihed;
	t_dihed   *dihed;
	real      *dihedValues;
	int       nImprop;
	t_improp  *improp;
	real      *impropValues;
	t_vecR    DM;
	t_vecR    dDM;
	t_vecR    COM;
	t_vecR    dCOM;
	t_vecR    axes[3];
	real      mass;
} t_mol;

typedef struct {
	int   acc;
	int   don;
	int   donH;
	real  distDonAcc;
	real  angle;
	int   intact;
	int   cont;
} t_hb;

typedef struct {
	int   *list;
	int   len;
	char  name[100];
} t_listI;

typedef struct {
	real  *list;
	int   len;
} t_listR;

typedef struct {
	char get[100][10];
	int  nGet;
	char getNot[100][10];
	int  nGetNot;
	char getBeginWC[100][10];
	int  getBeginWClen[100];
	int  nGetBeginWC;
	char getNotBeginWC[100][10];
	int  getNotBeginWClen[100];
	int  nGetNotBeginWC;
	char getEndWC[100][10];
        int  getEndWClen[100]; 
        int  nGetEndWC;
        char getNotEndWC[100][10];
        int  getNotEndWClen[100];
        int  nGetNotEndWC;
} t_nameRule;

typedef struct {
	char 		name[100];
	char		selChain[100][10];
	int		nSelChain;
	t_nameRule	chainRule;
	char		selResName[100][10];
	int		nSelResName;
	t_nameRule	resRule;
	char		selAtomName[100][10];
	int		nSelAtomName;
	t_nameRule	atomRule;
	int		selResNr[1000];
	int		nResNr;
	int		resNrRange[2];
	int		doResNrRange;
	int		selAtomNr[1000];
	int		nAtomNr;
	int		atomNrRange[2];
	int		doAtomNrRange;
} t_statGrpDef;

typedef struct {
	char	name[100];
	char	refGrpName[100];
	int	refGrpIdx;
	int	doRef;
	char	preGrpName[100];
	int	preGrpIdx;
	real	distRange[2];
	real    distSqRange[2];
	int	doDist;
	real	xRange[2];
	int	doX;
	real	yRange[2];
	int	doY;
	real	zRange[2];
	int	doZ;
	t_vecR  crd;
	int	doCrd;
	int	nSel;       /* nSel overrides DistRange */
	int	doN;
	char	refMode[100];  /* COM, closest_atom */
	char	selMode[100];  /* atoms, mol_by_COM, mol_by_atom */
} t_dynGrpDef;

typedef struct {
	char	name[100];
	char    grpNames[10][100];
	int	grpIdx[10];
	int	nGrpNames;
} t_combGrpDef;

typedef struct {
	t_statGrpDef	*statGrpDef;
	int		nStat;
	t_dynGrpDef	*dynGrpDef;
	int		nDyn;
	t_combGrpDef	*combGrpDef;
	int		nComb;
} t_grpDef;

typedef struct {
	char	name[100];
	int	*atoms;
	int	nAtoms;
	int	*mols;
	int	nMols;
	t_vecR	COM;
	t_vecR  COMvel;
	t_matR  T;         /* inertiaTensor */
	t_vecR  L;         /* angular momentum */
	real	spread;
} t_grp;

typedef struct {
	t_vecR  *coords;
	real    *values;
	int     nx,ny,nz;
} t_gridR;

typedef struct {
        t_vecR         *coords;
        fftwf_complex  *values;
        int            nx,ny,nz;
} t_gridCf;

