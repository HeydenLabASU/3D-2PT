#include <stdio.h>
#include <string.h>
#include "../../include/fatal.h"
#include "../../include/dataTypes.h"
#include "../../include/alloc.h"
#include "../../include/io.h"

int getLineFromJOB(FILE *in,char *buffer,int *lineCnt)
{
        fgets(buffer,5000,in);
        lineCnt[0]++;
        while(strncmp(buffer,"#",1)==0)
        {
                fgets(buffer,5000,in);
                lineCnt[0]++;
        }

        return 0;
}


int getIntFromJOB(char *fn,FILE *in,int *dest,char *keyword,int *lineCnt)
{
        char buffer[5000];
        int readCnt;
        char object[100];
        char err[5000];

        getLineFromJOB(in,buffer,lineCnt);
        readCnt=sscanf(buffer,"%s %d",object,dest);
        if(strcmp(object,keyword)!=0 || readCnt!=2) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getIntFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
        return 0;
}

int getIntFromJOB2(char *fn,FILE *in,int *dest,int *activator,char *keyword,int *lineCnt)
{
        char buffer[5000];
	char tmp[15];
        int readCnt;
        char object[100];
        char err[5000];

        getLineFromJOB(in,buffer,lineCnt);
        readCnt=sscanf(buffer,"%s %s",object,tmp);
        if(strcmp(object,keyword)!=0 || readCnt!=2) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getIntFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
	if(strcmp(tmp,"*")==0) {
		activator[0]=0;
		dest[0]=0;
	} else {
		readCnt=sscanf(tmp,"%d",dest);
		activator[0]=1;
		if(readCnt!=1) {
			sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getIntFromJOB2\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
		}
	}
        return 0;
}

int getStringFromJOB(char *fn,FILE *in,char *dest,char *keyword,int *lineCnt)
{
        char buffer[5000];
        int readCnt;
        char object[100];
        char err[5000];

        getLineFromJOB(in,buffer,lineCnt);
        readCnt=sscanf(buffer,"%s %s",object,dest);
        if(strcmp(object,keyword)!=0 || readCnt!=2) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getStringFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
        return 0;
}

int getMultStringsFromJOB(char *fn,FILE *in,char (*dest)[10],int *nRead,char *keyword,int *lineCnt)
{
	char buffer[5000];
	char truncBuffer[5000];
	int readCnt;
	char object[100];
	char err[5000];
	int i;
	int wildCard=0;
	char *ptr;

	readCnt=0;
	getLineFromJOB(in,buffer,lineCnt);
	readCnt=sscanf(buffer,"%s %s",object,dest[0]);
	if(strcmp(object,keyword)!=0 || readCnt!=2) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getMultStringsFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
	if(strcmp(dest[0],"*")==0) wildCard=1;
	strcpy(truncBuffer,&buffer[strlen(object)]);
	ptr=strstr(truncBuffer,dest[0]);
	strcpy(buffer,ptr);
	strcpy(truncBuffer,&buffer[strlen(dest[0])]);
	strcpy(buffer,truncBuffer);
	i=1;
	while(sscanf(buffer,"%s",dest[i])==1 && i<100)
	{
		if(strcmp(dest[i],"*")==0) wildCard=1;
		ptr=strstr(buffer,dest[i]);
		strcpy(truncBuffer,ptr);
		strcpy(buffer,&truncBuffer[strlen(dest[i])]);
		i++;
	}

	if(wildCard==0) nRead[0]=i;
        else nRead[0]=0;

	if(i==100) printf("WARNING: NUMBER OF STRINGS EXCEEDED LIMIT (100) ON LINE %d OF JOB FILE %s\n",lineCnt[0],fn);

	return 0;
}

int getMultStringsFromJOB2(char *fn,FILE *in,char (*dest)[100],int *nRead,char *keyword,int *lineCnt)
{
        char buffer[5000];
        char truncBuffer[5000];
        int readCnt;
        char object[100];
        char err[5000];
        int i;
        int wildCard=0;
        char *ptr;

        readCnt=0;
        getLineFromJOB(in,buffer,lineCnt);
        readCnt=sscanf(buffer,"%s %s",object,dest[0]);
        if(strcmp(object,keyword)!=0 || readCnt!=2) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getMultStringsFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
        if(strcmp(dest[0],"*")==0) wildCard=1;
        strcpy(truncBuffer,&buffer[strlen(object)]);
        ptr=strstr(truncBuffer,dest[0]);
        strcpy(buffer,ptr);
        strcpy(truncBuffer,&buffer[strlen(dest[0])]);
        strcpy(buffer,truncBuffer);
        i=1;
        while(sscanf(buffer,"%s",dest[i])==1 && i<10)
        {
                if(strcmp(dest[i],"*")==0) wildCard=1;
                ptr=strstr(buffer,dest[i]);
                strcpy(truncBuffer,ptr);
                strcpy(buffer,&truncBuffer[strlen(dest[i])]);
                i++;
        }

        if(wildCard==0) nRead[0]=i;
        else nRead[0]=0;

        if(i==10) printf("WARNING: NUMBER OF STRINGS EXCEEDED LIMIT (10) ON LINE %d OF JOB FILE %s\n",lineCnt[0],fn);

        return 0;
}

int getMultIntsFromJOB(char *fn,FILE *in,int *dest,int *nRead,char *keyword,int *lineCnt)
{
        char buffer[5000];
	char tmp[15];
        char truncBuffer[5000];
        int readCnt;
        char object[100];
        char err[5000];
        int i;
	int wildCard=0;
        char *ptr;

        readCnt=0;
        getLineFromJOB(in,buffer,lineCnt);
        readCnt=sscanf(buffer,"%s %s",object,tmp);
        if(strcmp(object,keyword)!=0 || readCnt!=2) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getMultIntsFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
	if(strcmp(tmp,"*")==0) wildCard=1;
	else {
		readCnt=sscanf(tmp,"%d",&dest[0]);
		if(readCnt==0) {
			sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getMultIntsFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                	fatal(err);
        	}
	}
        strcpy(truncBuffer,&buffer[strlen(object)]);
        ptr=strstr(truncBuffer,tmp);
        strcpy(buffer,ptr);
        strcpy(truncBuffer,&buffer[strlen(tmp)]);
        strcpy(buffer,truncBuffer);
        i=1;
        while(sscanf(buffer,"%s",tmp)==1 && i<1000)
        {
		if(strcmp(tmp,"*")==0) wildCard=1;
		else {
			readCnt=sscanf(tmp,"%d",&dest[i]);
			if(readCnt==0) {
                		sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getMultIntsFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                		fatal(err);
        		}
		}
                ptr=strstr(buffer,tmp);
                strcpy(truncBuffer,ptr);
                strcpy(buffer,&truncBuffer[strlen(tmp)]);
                i++;
        }

	if(wildCard==0) nRead[0]=i;
        else nRead[0]=0;

        if(i==1000) printf("WARNING: NUMBER OF INTS EXCEEDED LIMIT (1000) ON LINE %d OF JOB FILE %s\n",lineCnt[0],fn);

        return 0;
}

int getIntRangeFromJOB(char *fn,FILE *in,int *dest,int *activate,char *keyword,int *lineCnt)
{
	char buffer[5000];
        char tmp[15],tmp2[15];
        int readCnt;
        char object[100];
        char err[5000];
	int wildCard1=0;
	int wildCard2=0;

        readCnt=0;
        getLineFromJOB(in,buffer,lineCnt);
        readCnt=sscanf(buffer,"%s %s %s",object,tmp,tmp2);
        if(strcmp(object,keyword)!=0) {
                sprintf(err,"1 WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getIntRangeFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
	if(strcmp(tmp,"*")!=0 && readCnt!=3) {
                sprintf(err,"2 WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getIntRangeFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }

	if(strcmp(tmp,"*")==0) wildCard1=1;
	if (readCnt==3)
	{
		if(strcmp(tmp2,"*")==0) wildCard2=1;
	} else wildCard2=1;
	if(wildCard1==1 && wildCard2==1) activate[0]=0;
	else {
		activate[0]=1;
		if(wildCard1==1) dest[0]=-1;
		else {
			readCnt=sscanf(tmp,"%d",&dest[0]);
			if(readCnt!=1) {
				sprintf(err,"3 WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getIntRangeFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
				fatal(err);
        		}
		}
		if(wildCard2==1) dest[1]=-1;
                else {
                        readCnt=sscanf(tmp2,"%d",&dest[1]);
                        if(readCnt!=1) {
                                sprintf(err,"4 WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getIntRangeFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                                fatal(err);
                        }
                }
	}

	return 0;
}

int getRealRangeFromJOB(char *fn,FILE *in,real *dest,int *activate,char *keyword,int *lineCnt)
{
        char buffer[5000];
        char tmp[15],tmp2[15];
        int readCnt;
        char object[100];
        char err[5000];
	float dummy;
 
        readCnt=0; 
        getLineFromJOB(in,buffer,lineCnt);
        readCnt=sscanf(buffer,"%s %s %s",object,tmp,tmp2);
        if(strcmp(object,keyword)!=0) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getRealRangeFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
        if(strcmp(tmp,"*")==0 && readCnt!=2) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getRealRangeFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);      
                fatal(err);
        }
 
        if(strcmp(tmp,"*")==0) activate[0]=0;
	else {
		activate[0]=1;
		readCnt=sscanf(tmp,"%f",&dummy);
		if(readCnt!=1) {
			sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getRealRangeFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
			fatal(err);
                }
		dest[0]=(real)dummy;
		readCnt=sscanf(tmp2,"%f",&dummy);
                if(readCnt!=1) {
                        sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getRealRangeFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
			fatal(err);
                }
                dest[1]=(real)dummy;
        }

	return 0;
}

int getVecRFromJOB(char *fn,FILE *in,t_vecR *dest,int *activate,char *keyword,int *lineCnt)
{
        char buffer[5000];
        char tmp[30],tmp2[30],tmp3[30];
        int readCnt;
        char object[100];
        char err[5000];
        float dummy;

        readCnt=0;
        getLineFromJOB(in,buffer,lineCnt);
        readCnt=sscanf(buffer,"%s %s %s %s",object,tmp,tmp2,tmp3);
        if(strcmp(object,keyword)!=0) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getVecRFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
        if(strcmp(tmp,"*")==0 && readCnt!=2) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getVecRFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }

        if(strcmp(tmp,"*")==0) activate[0]=0;
        else {
                activate[0]=1;
                readCnt=sscanf(tmp,"%f",&dummy);
                if(readCnt!=1) {
                        sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getVecRFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                        fatal(err);
                }
                dest[0].x=(real)dummy;
                readCnt=sscanf(tmp2,"%f",&dummy);
                if(readCnt!=1) {
                        sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getVecRFromJOB\n -> job.c\n",fn,lineCnt[0],buffer);
                        fatal(err);
                }
                dest[0].y=(real)dummy;
		readCnt=sscanf(tmp3,"%f",&dummy);
                if(readCnt!=1) {
                        sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getVecRFromJOB\n -> job.c\n",fn,lineCnt[
0],buffer);
                        fatal(err);
                }
                dest[0].z=(real)dummy;
        }

        return 0;
}

int makeNameRule(char (*names)[10],int nNames,t_nameRule *nameRule)
{
	int i;
	char tmp[10];

	nameRule[0].nGet=0;
	nameRule[0].nGetNot=0;
	nameRule[0].nGetBeginWC=0;
	nameRule[0].nGetNotBeginWC=0;
	nameRule[0].nGetEndWC=0;
	nameRule[0].nGetNotEndWC=0;

	for(i=0;i<nNames;i++)
	{
		if(strncmp(names[i],"!",1)!=0 && strstr(names[i],"*")==NULL)
		{
			strcpy(nameRule[0].get[nameRule[0].nGet],names[i]);
			nameRule[0].nGet++;

		} else if(strncmp(names[i],"!",1)==0 && strstr(names[i],"*")==NULL)
		{
			strcpy(nameRule[0].getNot[nameRule[0].nGetNot],&names[i][1]);
			nameRule[0].nGetNot++;

		} else if(strncmp(names[i],"*",1)==0)
		{
			strcpy(nameRule[0].getBeginWC[nameRule[0].nGetBeginWC],&names[i][1]);
			nameRule[0].getBeginWClen[nameRule[0].nGetBeginWC]=strlen(names[i])-1;
			nameRule[0].nGetBeginWC++;

		} else if(strncmp(names[i],"!",1)==0 && strncmp(&names[i][1],"*",1)==0)
		{
			strcpy(nameRule[0].getNotBeginWC[nameRule[0].nGetNotBeginWC],&names[i][2]);
			nameRule[0].getNotBeginWClen[nameRule[0].nGetNotBeginWC]=strlen(names[i])-2;
			nameRule[0].nGetNotBeginWC++;

		} else if(strncmp(names[i],"!",1)!=0 && strncmp(&names[i][strlen(names[i])-1],"*",1)==0)
		{
			strcpy(tmp,names[i]);
			tmp[strlen(tmp)-1]=(char)0;
			strcpy(nameRule[0].getEndWC[nameRule[0].nGetEndWC],tmp);
			nameRule[0].getEndWClen[nameRule[0].nGetEndWC]=strlen(tmp);
			nameRule[0].nGetEndWC++;
		} else if(strncmp(names[i],"!",1)==0 && strncmp(&names[i][strlen(names[i])-1],"*",1)==0)
		{
			strcpy(tmp,&names[i][1]);
			tmp[strlen(tmp)-1]=(char)0;
			strcpy(nameRule[0].getNotEndWC[nameRule[0].nGetNotEndWC],tmp);
			nameRule[0].getNotEndWClen[nameRule[0].nGetNotEndWC]=strlen(tmp);
			nameRule[0].nGetNotEndWC++;
		} else fatal("UNEXPECTED BEHAVIOUR WHEN EXECUTING makeNameRule\n -> makeNameRule\n -> job.c\n");

	}

	return 0;
}

int getStatGrpFromJOB(char *fn,FILE *in,t_statGrpDef *statGrpDef,int *lineCnt)
{
	getStringFromJOB(fn,in,statGrpDef[0].name,"name",lineCnt);

	getMultStringsFromJOB(fn,in,statGrpDef[0].selChain,&statGrpDef[0].nSelChain,"chain(s)",lineCnt);
	makeNameRule(statGrpDef[0].selChain,statGrpDef[0].nSelChain,&statGrpDef[0].chainRule);

	getMultStringsFromJOB(fn,in,statGrpDef[0].selResName,&statGrpDef[0].nSelResName,"residue(s)",lineCnt);
	makeNameRule(statGrpDef[0].selResName,statGrpDef[0].nSelResName,&statGrpDef[0].resRule);

	getMultStringsFromJOB(fn,in,statGrpDef[0].selAtomName,&statGrpDef[0].nSelAtomName,"atom(s)",lineCnt);
	makeNameRule(statGrpDef[0].selAtomName,statGrpDef[0].nSelAtomName,&statGrpDef[0].atomRule);

	getMultIntsFromJOB(fn,in,statGrpDef[0].selResNr,&statGrpDef[0].nResNr,"residue_number(s)",lineCnt);

	getIntRangeFromJOB(fn,in,statGrpDef[0].resNrRange,&statGrpDef[0].doResNrRange,"residue_range",lineCnt);

	getMultIntsFromJOB(fn,in,statGrpDef[0].selAtomNr,&statGrpDef[0].nAtomNr,"atom_number(s)",lineCnt);

	getIntRangeFromJOB(fn,in,statGrpDef[0].atomNrRange,&statGrpDef[0].doAtomNrRange,"atom_range",lineCnt);

	return 0;
}

int getDynGrpFromJOB(char *fn,FILE *in,t_dynGrpDef *dynGrpDef,int *lineCnt)
{
	char err[300];

	getStringFromJOB(fn,in,dynGrpDef[0].name,"name",lineCnt);

	getStringFromJOB(fn,in,dynGrpDef[0].refGrpName,"reference_group",lineCnt);
	if(strncmp(dynGrpDef[0].refGrpName,"*",1)==0) dynGrpDef[0].doRef=0;
	else dynGrpDef[0].doRef=1;

	getStringFromJOB(fn,in,dynGrpDef[0].preGrpName,"select_from",lineCnt);

	getRealRangeFromJOB(fn,in,dynGrpDef[0].distRange,&dynGrpDef[0].doDist,"distance_range",lineCnt);
	dynGrpDef[0].distSqRange[0]=dynGrpDef[0].distRange[0]*dynGrpDef[0].distRange[0];
	dynGrpDef[0].distSqRange[1]=dynGrpDef[0].distRange[1]*dynGrpDef[0].distRange[1];

	getRealRangeFromJOB(fn,in,dynGrpDef[0].xRange,&dynGrpDef[0].doX,"x_range",lineCnt);

	getRealRangeFromJOB(fn,in,dynGrpDef[0].yRange,&dynGrpDef[0].doY,"y_range",lineCnt);

	getRealRangeFromJOB(fn,in,dynGrpDef[0].zRange,&dynGrpDef[0].doZ,"z_range",lineCnt);

	getVecRFromJOB(fn,in,&dynGrpDef[0].crd,&dynGrpDef[0].doCrd,"dist_to_vec",lineCnt);

	getIntFromJOB2(fn,in,&dynGrpDef[0].nSel,&dynGrpDef[0].doN,"closest_n",lineCnt);

	getStringFromJOB(fn,in,dynGrpDef[0].refMode,"ref_mode",lineCnt);
	if(strcmp(dynGrpDef[0].refMode,"COM")!=0 && strcmp(dynGrpDef[0].refMode,"closest_atom")!=0)
	{
		sprintf(err,"UNKNOWN REFERENCE MODE %s SPECIFIED IN FILE %s LINE %d\nPOSSIBLE CHOICES ARE 'COM' OR 'closest_atom'\n -> getDynGrpFromJOB\n -> job.c\n",dynGrpDef[0].refMode,fn,lineCnt[0]);
		fatal(err);
	}

	getStringFromJOB(fn,in,dynGrpDef[0].selMode,"sel_mode",lineCnt);
	if(strcmp(dynGrpDef[0].selMode,"atoms")!=0 && strcmp(dynGrpDef[0].selMode,"mol_by_COM")!=0 && strcmp(dynGrpDef[0].selMode,"mol_by_atom")!=0)
        {
		sprintf(err,"UNKNOWN SELECTION MODE %s SPECIFIED IN FILE %s LINE %d\nPOSSIBLE CHOICES ARE 'atoms' OR 'mol_by_COM' OR 'mol_by_atom'\n -> getDynGrpFromJOB\n -> job.c\n",dynGrpDef[0].refMode,fn,lineCnt[0]);
                fatal(err);
        }

	return 0;
}

int getCombGrpFromJOB(char *fn,FILE *in,t_combGrpDef *combGrpDef,int *lineCnt)
{
	getStringFromJOB(fn,in,combGrpDef[0].name,"name",lineCnt);

	getMultStringsFromJOB2(fn,in,combGrpDef[0].grpNames,&combGrpDef[0].nGrpNames,"groups",lineCnt);

	return 0;
}

int getGrpDefFromJOB(char *fn,FILE *in,t_grpDef *grpDef,int *lineCnt)
{
	int i,j,k;
	char err[300];

	getIntFromJOB(fn,in,&grpDef[0].nStat,"[static_groups]",lineCnt);
	allocStatGrpDefs(&grpDef[0].statGrpDef,grpDef[0].nStat);
	for(i=0;i<grpDef[0].nStat;i++)
	{
		getStatGrpFromJOB(fn,in,&grpDef[0].statGrpDef[i],lineCnt);
	}

	getIntFromJOB(fn,in,&grpDef[0].nDyn,"[dynamic_groups]",lineCnt);
        allocDynGrpDefs(&grpDef[0].dynGrpDef,grpDef[0].nDyn);
        for(i=0;i<grpDef[0].nDyn;i++)
        {
                getDynGrpFromJOB(fn,in,&grpDef[0].dynGrpDef[i],lineCnt);
		if(grpDef[0].dynGrpDef[i].doRef==1)
		{
			grpDef[0].dynGrpDef[i].refGrpIdx=-1;
			for(j=0;j<grpDef[0].nStat;j++)
			{
				if(strcmp(grpDef[0].dynGrpDef[i].refGrpName,grpDef[0].statGrpDef[j].name)==0) grpDef[0].dynGrpDef[i].refGrpIdx=j;
			}
			for(j=0;j<i;j++)
			{
				if(strcmp(grpDef[0].dynGrpDef[i].refGrpName,grpDef[0].dynGrpDef[j].name)==0) grpDef[0].dynGrpDef[i].refGrpIdx=j+grpDef[0].nStat;
			}
			if(grpDef[0].dynGrpDef[i].refGrpIdx==-1)
			{
				sprintf(err,"REFERENCE GROUP %s FOR DYNAMIC GROUP %s NOT FOUND!\n -> getGrpDefFromJOB\n -> job.c\n",grpDef[0].dynGrpDef[i].refGrpName,grpDef[0].dynGrpDef[i].name);
				fatal(err);
			}
		}
		grpDef[0].dynGrpDef[i].preGrpIdx=-1;
		for(j=0;j<grpDef[0].nStat;j++)
		{
			if(strcmp(grpDef[0].dynGrpDef[i].preGrpName,grpDef[0].statGrpDef[j].name)==0) grpDef[0].dynGrpDef[i].preGrpIdx=j;
		}
		for(j=0;j<i;j++)
                {
                        if(strcmp(grpDef[0].dynGrpDef[i].preGrpName,grpDef[0].dynGrpDef[j].name)==0) grpDef[0].dynGrpDef[i].preGrpIdx=j+grpDef[0].nStat;
                }
		if(grpDef[0].dynGrpDef[i].preGrpIdx==-1)
		{
			sprintf(err,"STATIC PRESELECTION GROUP %s FOR DYNAMIC GROUP %s NOT FOUND!\n -> getGrpDefFromJOB\n -> job.c\n",grpDef[0].dynGrpDef[i].preGrpName,grpDef[0].dynGrpDef[i].name);
			fatal(err);
		}
        }

	getIntFromJOB(fn,in,&grpDef[0].nComb,"[combination_groups]",lineCnt);
        allocCombGrpDefs(&grpDef[0].combGrpDef,grpDef[0].nComb);
        for(i=0;i<grpDef[0].nComb;i++)
        {
                getCombGrpFromJOB(fn,in,&grpDef[0].combGrpDef[i],lineCnt);
		for(k=0;k<grpDef[0].combGrpDef[i].nGrpNames;k++)
		{
			grpDef[0].combGrpDef[i].grpIdx[k]=-1;
			for(j=0;j<grpDef[0].nStat;j++)
        		{
				if(strcmp(grpDef[0].combGrpDef[i].grpNames[k],grpDef[0].statGrpDef[j].name)==0) grpDef[0].combGrpDef[i].grpIdx[k]=j;
			}
			for(j=0;j<grpDef[0].nDyn;j++)
			{
				if(strcmp(grpDef[0].combGrpDef[i].grpNames[k],grpDef[0].dynGrpDef[j].name)==0) grpDef[0].combGrpDef[i].grpIdx[k]=j+grpDef[0].nStat;
			}
			if( grpDef[0].combGrpDef[i].grpIdx[k]==-1)
			{
				sprintf(err,"MEMBER GROUP %s FOR COMBINATION GROUP %s NOT FOUND!\n -> getGrpDefFromJOB\n -> job.c\n",grpDef[0].combGrpDef[i].grpNames[k],grpDef[0].combGrpDef[i].name);
				fatal(err);
			}
		}
	}

	return 0;
}

int readJob(char *fn,t_grpDef *grpDef)
{
	FILE *job;
	int lineCnt=0;

	saveOpenRead(&job,fn);

	getGrpDefFromJOB(fn,job,grpDef,&lineCnt);

	fclose(job);

	return 0;
}

