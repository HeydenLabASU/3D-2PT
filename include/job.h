int getLineFromJOB(FILE *in,char *buffer,int *lineCnt);

int getIntFromJOB(char *fn,FILE *in,int *dest,char *keyword,int *lineCnt);

int getIntFromJOB2(char *fn,FILE *in,int *dest,int *activator,char *keyword,int *lineCnt);

int getStringFromJOB(char *fn,FILE *in,char *dest,char *keyword,int *lineCnt);

int getMultStringsFromJOB(char *fn,FILE *in,char (*dest)[6],int *nRead,char *keyword,int *lineCnt);

int getMultStringsFromJOB2(char *fn,FILE *in,char (*dest)[100],int *nRead,char *keyword,int *lineCnt);

int getMultIntsFromJOB(char *fn,FILE *in,int *dest,int *nRead,char *keyword,int *lineCnt);

int getIntRangeFromJOB(char *fn,FILE *in,int *dest,int *activate,char *keyword,int *lineCnt);

int getRealRangeFromJOB(char *fn,FILE *in,real *dest,int *activate,char *keyword,int *lineCnt);

int getVecRFromJOB(char *fn,FILE *in,t_vecR *dest,int *activate,char *keyword,int *lineCnt);

int makeNameRule(char (*names)[6],int nNames,t_nameRule *nameRule);

int getStatGrpFromJOB(char *fn,FILE *in,t_statGrpDef *statGrpDef,int *lineCnt);

int getDynGrpFromJOB(char *fn,FILE *in,t_dynGrpDef *dynGrpDef,int *lineCnt);

int getCombGrpFromJOB(char *fn,FILE *in,t_combGrpDef *combGrpDef,int *lineCnt);

int getGrpDefFromJOB(char *fn,FILE *in,t_grpDef *grpDef,int *lineCnt);

int readJob(char *fn,t_grpDef *grpDef);


