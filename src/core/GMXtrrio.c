
/*
 * $Id: gmx_system_xdr.h,v 1.3 2003/11/17 21:50:40 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */


/*
 * This header file is ONLY used on windows systems, since these do
 * not include the XDR routines present on a unix machine. It will
 * most probably work on other platforms too, but make sure you
 * test that the xtc files produced are ok before using it. 
 *
 * This header file contains Gromacs versions of the definitions for 
 * Sun External Data Representation (XDR) headers and routines.
 *
 * On most UNIX systems this is already present as part of your 
 * system libraries, but since we want to make Gromacs portable to
 * platforms like Microsoft Windows we have created a private version
 * of the necessary routines and distribute them with the Gromacs source.
 * 
 * Although the rest of Gromacs is GPL, you can copy and use the XDR 
 * routines in any way you want as long as you obey Sun's license:
 *
 * Sun RPC is a product of Sun Microsystems, Inc. and is provided for
 * unrestricted use provided that this legend is included on all tape
 * media and as a part of the software program in whole or part.  Users
 * may copy or modify Sun RPC without charge, but are not authorized
 * to license or distribute it to anyone else except as part of a product or
 * program developed by the user.
 *
 * SUN RPC IS PROVIDED AS IS WITH NO WARRANTIES OF ANY KIND INCLUDING THE
 * WARRANTIES OF DESIGN, MERCHANTIBILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE, OR ARISING FROM A COURSE OF DEALING, USAGE OR TRADE PRACTICE.
 *
 * Sun RPC is provided with no support and without any obligation on the
 * part of Sun Microsystems, Inc. to assist in its use, correction,
 * modification or enhancement.
 *
 * SUN MICROSYSTEMS, INC. SHALL HAVE NO LIABILITY WITH RESPECT TO THE
 * INFRINGEMENT OF COPYRIGHTS, TRADE SECRETS OR ANY PATENTS BY SUN RPC
 * OR ANY PART THEREOF.
 *
 * In no event will Sun Microsystems, Inc. be liable for any lost revenue
 * or profits or other special, indirect and consequential damages, even if
 * Sun has been advised of the possibility of such damages.
 *
 * Sun Microsystems, Inc.
 * 2550 Garcia Avenue
 * Mountain View, California  94043
 */ 



#ifndef _gmx_system_xdr_h
#define _gmx_system_xdr_h

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

/*
 * Xdr operations.  XDR_ENCODE causes the type to be encoded into the
 * stream.  XDR_DECODE causes the type to be extracted from the stream.
 * XDR_FREE can be used to release the space allocated by an 
 * XDR_DECODE request.
 */

/* We already have a boolean type in Gromacs, but the XDR library
 * one has a slightly different name (the calls should be identical).
 */
typedef int bool_t;

/* 
 * Aninteger type that is 32 bits wide. Check if int,
 * long or short is 32 bits and die if none of them is :-)
 */
#if (INT_MAX == 2147483647)
   typedef int xdr_int32_t;
   typedef unsigned int xdr_uint32_t;
#elif (LONG_MAX == 2147483647L)
   typedef long xdr_int32_t;
   typedef unsigned long xdr_uint32_t;
#elif (SHRT_MAX == 2147483647)
   typedef short xdr_int32_t;
   typedef unsigned short xdr_uint32_t;
#else
#  error ERROR: No 32 bit wide integer type found!
#endif


enum xdr_op {
  XDR_ENCODE = 0,
  XDR_DECODE = 1,
  XDR_FREE = 2
};

#ifndef FALSE
#      define  FALSE   (0)
#endif
#ifndef TRUE
#      define  TRUE    (1)
#endif


#define BYTES_PER_XDR_UNIT	(4)
/* Macro to round up to units of 4. */
#define XDR_RNDUP(x)  (((x) + BYTES_PER_XDR_UNIT - 1) & ~(BYTES_PER_XDR_UNIT - 1))


/*
 * The XDR handle.
 * Contains operation which is being applied to the stream,
 * an operations vector for the particular implementation (e.g. see xdr_mem.c),
 * and two private fields for the use of the particular implementation.
 */
typedef struct XDR XDR;
struct XDR
  {
    enum xdr_op x_op;		/* operation; fast additional param */
    struct xdr_ops
      {
	bool_t (*x_getlong) (XDR *__xdrs, long *__lp);
	/* get a long from underlying stream */
	bool_t (*x_putlong) (XDR *__xdrs, long *__lp);
	/* put a long to " */
	bool_t (*x_getbytes) (XDR *__xdrs, char *__addr, unsigned int __len);
	/* get some bytes from " */
	bool_t (*x_putbytes) (XDR *__xdrs, char *__addr, unsigned int __len);
	/* put some bytes to " */
	unsigned int (*x_getpostn) (XDR *__xdrs);
	/* returns bytes off from beginning */
	bool_t (*x_setpostn) (XDR *__xdrs, unsigned int __pos);
	/* lets you reposition the stream */
	xdr_int32_t *(*x_inline) (XDR *__xdrs, int __len);
	/* buf quick ptr to buffered data */
	void (*x_destroy) (XDR *__xdrs);
	/* free privates of this xdr_stream */
	bool_t (*x_getint32) (XDR *__xdrs, xdr_int32_t *__ip);
	/* get a int from underlying stream */
	bool_t (*x_putint32) (XDR *__xdrs, xdr_int32_t *__ip);
	/* put a int to " */
      }
     *x_ops;
    char *x_public;		/* users' data */
    char *x_private;		/* pointer to private data */
    char *x_base;		/* private used for position info */
    int x_handy;		/* extra private word */
  };

/*
 * A xdrproc_t exists for each data type which is to be encoded or decoded.
 *
 * The second argument to the xdrproc_t is a pointer to an opaque pointer.
 * The opaque pointer generally points to a structure of the data type
 * to be decoded.  If this pointer is 0, then the type routines should
 * allocate dynamic storage of the appropriate size and return it.
 */

typedef bool_t (*xdrproc_t) (XDR *, void *,...);

/*
 * Operations defined on a XDR handle
 *
 * XDR          *xdrs;
 * xdr_int32_t  *int32p;
 * long         *longp;
 * char         *addr;
 * unsigned int  len;
 * unsigned int  pos;
 */


#define xdr_getint32(xdrs, int32p)                      \
        (*(xdrs)->x_ops->x_getint32)(xdrs, int32p)

#define xdr_putint32(xdrs, int32p)                      \
        (*(xdrs)->x_ops->x_putint32)(xdrs, int32p)

#define xdr_getlong(xdrs, longp)			\
	(*(xdrs)->x_ops->x_getlong)(xdrs, longp)

#define xdr_putlong(xdrs, longp)			\
	(*(xdrs)->x_ops->x_putlong)(xdrs, longp)

#define xdr_getbytes(xdrs, addr, len)			\
	(*(xdrs)->x_ops->x_getbytes)(xdrs, addr, len)

#define xdr_putbytes(xdrs, addr, len)			\
	(*(xdrs)->x_ops->x_putbytes)(xdrs, addr, len)

#define xdr_getpos(xdrs)				\
	(*(xdrs)->x_ops->x_getpostn)(xdrs)

#define xdr_setpos(xdrs, pos)				\
	(*(xdrs)->x_ops->x_setpostn)(xdrs, pos)

#define	xdr_inline(xdrs, len)				\
	(*(xdrs)->x_ops->x_inline)(xdrs, len)

#define	xdr_destroy(xdrs)					\
	do {							\
		if ((xdrs)->x_ops->x_destroy)			\
			(*(xdrs)->x_ops->x_destroy)(xdrs);	\
	} while (0)


extern bool_t xdr_u_short (XDR *__xdrs, unsigned short *__usp);
extern bool_t xdr_int (XDR *__xdrs, int *__ip);
extern bool_t xdr_long (XDR *__xdrs, long *__lp); 
extern bool_t xdr_u_long (XDR *__xdrs, unsigned long *__ulp); 
extern bool_t xdr_bool (XDR *__xdrs, int *__bp);
extern bool_t xdr_opaque (XDR *__xdrs, char *__cp, unsigned int __cnt);
extern bool_t xdr_string (XDR *__xdrs, char **__cpp, unsigned int __maxsize);
extern bool_t xdr_char (XDR *__xdrs, char *__cp);
extern bool_t xdr_u_char (XDR *__xdrs, unsigned char *__cp);
extern bool_t xdr_vector (XDR *__xdrs, char *__basep, unsigned int __nelem,
			  unsigned int __elemsize, xdrproc_t __xdr_elem);
extern bool_t xdr_float (XDR *__xdrs, float *__fp);
extern bool_t xdr_double (XDR *__xdrs, double *__dp);
extern void xdrstdio_create (XDR *__xdrs, FILE *__file, enum xdr_op __xop);

/* free memory buffers for xdr */
extern void xdr_free (xdrproc_t __proc, char *__objp);


#endif /* _gmx_system_xdr_h */

int 
xdropen(XDR *xdrs, const char *filename, const char *type);


int 
xdrclose(XDR *xdrs);

/*
 * for unit alignment
 */
static char xdr_zero[BYTES_PER_XDR_UNIT] = {0, 0, 0, 0};

static xdr_uint32_t xdr_swapbytes(xdr_uint32_t x)
{
  xdr_uint32_t y;
  int i;
  char *px=(char *)&x;
  char *py=(char *)&y;
  
  for(i=0;i<4;i++)
    py[i]=px[3-i];
  
  return y;
}

static xdr_uint32_t xdr_htonl(xdr_uint32_t x)
{
  short s=0x0F00;
  if( *((char *)&s)==(char)0x0F) {
    /* bigendian, do nothing */
    return x;
  } else {
    /* smallendian,swap bytes */
    return xdr_swapbytes(x);
  }
}

static xdr_uint32_t xdr_ntohl(xdr_uint32_t x)
{
  short s=0x0F00;
  if( *((char *)&s)==(char)0x0F) {
    /* bigendian, do nothing */
    return x;
  } else {
    /* smallendian, swap bytes */
    return xdr_swapbytes(x);
  }
}


/*
 * Free a data structure using XDR
 * Not a filter, but a convenient utility nonetheless
 */
void
xdr_free (xdrproc_t proc, char *objp)
{
  XDR x;

  x.x_op = XDR_FREE;
  (*proc) (&x, objp);
}

/*
 * XDR nothing
 */
bool_t
xdr_void (void)
{
  return TRUE;
}

/*
 * XDR integers
 */
bool_t
xdr_int (XDR *xdrs, int *ip)
{

#if INT_MAX < LONG_MAX
  long l;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      l = (long) *ip;
      return xdr_putlong (xdrs, &l);

    case XDR_DECODE:
      if (!xdr_getlong (xdrs, &l))
	{
	  return FALSE;
	}
      *ip = (int) l;
    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
#elif INT_MAX == LONG_MAX
  return xdr_long (xdrs, (long *) ip);
#elif INT_MAX == SHRT_MAX
  return xdr_short (xdrs, (short *) ip);
#else
#error unexpected integer sizes in_xdr_int()
#endif
}


/*
 * XDR unsigned integers
 */
bool_t
xdr_u_int (XDR *xdrs, unsigned int *up)
{
#if UINT_MAX < ULONG_MAX
  unsigned long l;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      l = (unsigned long) * up;
      return xdr_putlong (xdrs, (long *)&l);

    case XDR_DECODE:
      if (!xdr_getlong (xdrs, (long *)&l))
	{
	  return FALSE;
	}
      *up = (unsigned int) l;
    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
#elif UINT_MAX == ULONG_MAX
  return xdr_u_long (xdrs, (unsigned long *) up);
#elif UINT_MAX == USHRT_MAX
  return xdr_short (xdrs, (short *) up);
#else
#error unexpected integer sizes in_xdr_u_int()
#endif
}


/*
 * XDR long integers
 * The definition of xdr_long() is kept for backward
 * compatibility. Instead xdr_int() should be used.
 */
bool_t
xdr_long (XDR *xdrs, long *lp)
{

  if (xdrs->x_op == XDR_ENCODE
      && (sizeof (xdr_int32_t) == sizeof (long)
	  || (xdr_int32_t) *lp == *lp))
    return xdr_putlong (xdrs, lp);

  if (xdrs->x_op == XDR_DECODE)
    return xdr_getlong (xdrs, lp);

  if (xdrs->x_op == XDR_FREE)
    return TRUE;

  return FALSE;
}


/*
 * XDR unsigned long integers
 * The definition of xdr_u_long() is kept for backward
 * compatibility. Instead xdr_u_int() should be used.
 */
bool_t
xdr_u_long (XDR *xdrs, unsigned long *ulp)
{
  switch (xdrs->x_op)
    {
    case XDR_DECODE:
      {
	long int tmp;

	if (xdr_getlong (xdrs, &tmp) == FALSE)
	  return FALSE;

	*ulp = (xdr_uint32_t) tmp;
	return TRUE;
      }

    case XDR_ENCODE:
      if (sizeof (xdr_uint32_t) != sizeof (unsigned long)
	  && (xdr_uint32_t) *ulp != *ulp)
	return FALSE;

      return xdr_putlong (xdrs, (long *) ulp);

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}


/*
 * XDR short integers
 */
bool_t
xdr_short (XDR *xdrs, short *sp)
{
  long l;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      l = (long) *sp;
      return xdr_putlong (xdrs, &l);

    case XDR_DECODE:
      if (!xdr_getlong (xdrs, &l))
	{
	  return FALSE;
	}
      *sp = (short) l;
      return TRUE;

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}


/*
 * XDR unsigned short integers
 */
bool_t
xdr_u_short (XDR *xdrs, unsigned short *usp)
{
  unsigned long l;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      l = (unsigned long) * usp;
      return xdr_putlong (xdrs, (long*)&l);

    case XDR_DECODE:
      if (!xdr_getlong (xdrs, (long*)&l))
	{
	  return FALSE;
	}
      *usp = (unsigned short) l;
      return TRUE;

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}


/*
 * XDR a char
 */
bool_t
xdr_char (XDR *xdrs, char *cp)
{
  int i;

  i = (*cp);
  if (!xdr_int (xdrs, &i))
    {
      return FALSE;
    }
  *cp = i;
  return TRUE;
}

/*
 * XDR an unsigned char
 */
bool_t
xdr_u_char (XDR *xdrs, unsigned char *cp)
{
  unsigned int u;

  u = (*cp);
  if (!xdr_u_int (xdrs, &u))
    {
      return FALSE;
    }
  *cp = u;
  return TRUE;
}

/*
 * XDR booleans
 */
bool_t
xdr_bool (XDR *xdrs, int *bp)
{
#define XDR_FALSE	((long) 0)
#define XDR_TRUE	((long) 1)
  long lb;

  switch (xdrs->x_op)
    {
    case XDR_ENCODE:
      lb = *bp ? XDR_TRUE : XDR_FALSE;
      return xdr_putlong (xdrs, &lb);

    case XDR_DECODE:
      if (!xdr_getlong (xdrs, &lb))
	{
	  return FALSE;
	}
      *bp = (lb == XDR_FALSE) ? FALSE : TRUE;
      return TRUE;

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
#undef XDR_FALSE
#undef XDR_TRUE
}



/*
 * XDR opaque data
 * Allows the specification of a fixed size sequence of opaque bytes.
 * cp points to the opaque object and cnt gives the byte length.
 */
bool_t
xdr_opaque (XDR *xdrs, char *cp, unsigned int cnt)
{
  unsigned int rndup;
  static char crud[BYTES_PER_XDR_UNIT];

  /*
   * if no data we are done
   */
  if (cnt == 0)
    return TRUE;

  /*
   * round byte count to full xdr units
   */
  rndup = cnt % BYTES_PER_XDR_UNIT;
  if (rndup > 0)
    rndup = BYTES_PER_XDR_UNIT - rndup;

  switch (xdrs->x_op)
    {
    case XDR_DECODE:
      if (!xdr_getbytes (xdrs, cp, cnt))
	{
	  return FALSE;
	}
      if (rndup == 0)
	return TRUE;
      return xdr_getbytes (xdrs, (char *)crud, rndup);

    case XDR_ENCODE:
      if (!xdr_putbytes (xdrs, cp, cnt))
	{
	  return FALSE;
	}
      if (rndup == 0)
	return TRUE;
      return xdr_putbytes (xdrs, xdr_zero, rndup);

    case XDR_FREE:
      return TRUE;
    }
  return FALSE;
}


/*
 * XDR null terminated ASCII strings
 * xdr_string deals with "C strings" - arrays of bytes that are
 * terminated by a NULL character.  The parameter cpp references a
 * pointer to storage; If the pointer is null, then the necessary
 * storage is allocated.  The last parameter is the max allowed length
 * of the string as specified by a protocol.
 */
bool_t
xdr_string (xdrs, cpp, maxsize)
     XDR *xdrs;
     char **cpp;
     unsigned int maxsize;
{
  char *sp = *cpp;	/* sp is the actual string pointer */
  unsigned int size;
  unsigned int nodesize;

  /*
   * first deal with the length since xdr strings are counted-strings
   */
  switch (xdrs->x_op)
    {
    case XDR_FREE:
      if (sp == NULL)
	{
	  return TRUE;		/* already free */
	}
      /* fall through... */
    case XDR_ENCODE:
      if (sp == NULL)
	return FALSE;
      size = strlen (sp);
      break;
    case XDR_DECODE:
      break;
    }
  if (!xdr_u_int (xdrs, &size))
    {
      return FALSE;
    }
  if (size > maxsize)
    {
      return FALSE;
    }
  nodesize = size + 1;

  /*
   * now deal with the actual bytes
   */
  switch (xdrs->x_op)
    {
    case XDR_DECODE:
      if (nodesize == 0)
	{
	  return TRUE;
	}
      if (sp == NULL)
	*cpp = sp = (char *) malloc (nodesize);
      if (sp == NULL)
	{
	  (void) fputs ("xdr_string: out of memory\n", stderr);
	  return FALSE;
	}
      sp[size] = 0;
      /* fall into ... */

    case XDR_ENCODE:
      return xdr_opaque (xdrs, sp, size);

    case XDR_FREE:
      free (sp);
      *cpp = NULL;
      return TRUE;
    }
  return FALSE;
}



/* Floating-point stuff */

bool_t
xdr_float(xdrs, fp)
     XDR *xdrs;
     float *fp;
{
	switch (xdrs->x_op) {

	case XDR_ENCODE:
		if (sizeof(float) == sizeof(long))
			return (xdr_putlong(xdrs, (long *)fp));
		else if (sizeof(float) == sizeof(int)) {
			long tmp = *(int *)fp;
			return (xdr_putlong(xdrs, &tmp));
		}
		break;

	case XDR_DECODE:
		if (sizeof(float) == sizeof(long))
			return (xdr_getlong(xdrs, (long *)fp));
		else if (sizeof(float) == sizeof(int)) {
			long tmp;
			if (xdr_getlong(xdrs, &tmp)) {
				*(int *)fp = tmp;
				return (TRUE);
			}
		}
		break;

	case XDR_FREE:
		return (TRUE);
	}
	return (FALSE);
}


bool_t
xdr_double(xdrs, dp)
     XDR *xdrs;
     double *dp;
{

  /* Windows and some other systems dont define double-precision
   * word order in the header files, so unfortunately we have
   * to calculate it!
   */
  static int LSW=-1; /* Least significant fp word */
  
  if(LSW<0) {
    double x=0.987654321; /* Just a number */

    /* Possible representations in IEEE double precision: 
     * (S=small endian, B=big endian)
     * 
     * Byte order, Word order, Hex
     *     S           S       b8 56 0e 3c dd 9a ef 3f    
     *     B           S       3c 0e 56 b8 3f ef 9a dd
     *     S           B       dd 9a ef 3f b8 56 0e 3c
     *     B           B       3f ef 9a dd 3c 0e 56 b8
     */ 
    
    unsigned char ix = *((char *)&x);
    
    if(ix==0xdd || ix==0x3f)
      LSW=1;  /* Big endian word order */
    else if(ix==0xb8 || ix==0x3c)
      LSW=0;  /* Small endian word order */
    else { /* Catch strange errors */
      printf("Error when detecting floating-point word order.\n"
	     "Do you have a non-IEEE system?\n"
	     "If possible, use the XDR libraries provided with your system,\n"
	     "instead of the Gromacs fallback XDR source.\n");
      exit(0);
    }
  }  
  
  switch (xdrs->x_op) {
    
  case XDR_ENCODE:
    if (2*sizeof(long) == sizeof(double)) {
      long *lp = (long *)dp;
      return (xdr_putlong(xdrs, lp+!LSW) &&
	      xdr_putlong(xdrs, lp+LSW));
    } else if (2*sizeof(int) == sizeof(double)) {
      int *ip = (int *)dp;
      long tmp[2];
      tmp[0] = ip[!LSW];
      tmp[1] = ip[LSW];
      return (xdr_putlong(xdrs, tmp) &&
	      xdr_putlong(xdrs, tmp+1));
    }
    break;
    
  case XDR_DECODE:
    if (2*sizeof(long) == sizeof(double)) {
      long *lp = (long *)dp;
      return (xdr_getlong(xdrs, lp+!LSW) &&
	      xdr_getlong(xdrs, lp+LSW));
    } else if (2*sizeof(int) == sizeof(double)) {
      int *ip = (int *)dp;
      long tmp[2];
      if (xdr_getlong(xdrs, tmp+!LSW) &&
	  xdr_getlong(xdrs, tmp+LSW)) {
	ip[0] = tmp[0];
	ip[1] = tmp[1];
	return (TRUE);
      }
    }
    break;
    
  case XDR_FREE:
    return (TRUE);
  }
  return (FALSE);
}


/* Array routines */

/*
 * xdr_vector():
 *
 * XDR a fixed length array. Unlike variable-length arrays,
 * the storage of fixed length arrays is static and unfreeable.
 * > basep: base of the array
 * > size: size of the array
 * > elemsize: size of each element
 * > xdr_elem: routine to XDR each element
 */
bool_t
xdr_vector (xdrs, basep, nelem, elemsize, xdr_elem)
     XDR *xdrs;
     char *basep;
     unsigned int nelem;
     unsigned int elemsize;
     xdrproc_t xdr_elem;
{
#define LASTUNSIGNED	((unsigned int)0-1)
  unsigned int i;
  char *elptr;

  elptr = basep;
  for (i = 0; i < nelem; i++)
    {
      if (!(*xdr_elem) (xdrs, elptr, LASTUNSIGNED))
	{
	  return FALSE;
	}
      elptr += elemsize;
    }
  return TRUE;
#undef LASTUNSIGNED
}



static bool_t xdrstdio_getlong (XDR *, long *);
static bool_t xdrstdio_putlong (XDR *, long *);
static bool_t xdrstdio_getbytes (XDR *, char *, unsigned int);
static bool_t xdrstdio_putbytes (XDR *, char *, unsigned int);
static unsigned int xdrstdio_getpos (XDR *);
static bool_t xdrstdio_setpos (XDR *, unsigned int);
static xdr_int32_t *xdrstdio_inline (XDR *, int);
static void xdrstdio_destroy (XDR *);
static bool_t xdrstdio_getint32 (XDR *, xdr_int32_t *);
static bool_t xdrstdio_putint32 (XDR *, xdr_int32_t *);

/*
 * Ops vector for stdio type XDR
 */
static const struct xdr_ops xdrstdio_ops =
{
  xdrstdio_getlong,		/* deserialize a long int */
  xdrstdio_putlong,		/* serialize a long int */
  xdrstdio_getbytes,       	/* deserialize counted bytes */
  xdrstdio_putbytes,     	/* serialize counted bytes */
  xdrstdio_getpos,		/* get offset in the stream */
  xdrstdio_setpos,		/* set offset in the stream */
  xdrstdio_inline,		/* prime stream for inline macros */
  xdrstdio_destroy,		/* destroy stream */
  xdrstdio_getint32,	/* deserialize a int */
  xdrstdio_putint32		/* serialize a int */
};

/*
 * Initialize a stdio xdr stream.
 * Sets the xdr stream handle xdrs for use on the stream file.
 * Operation flag is set to op.
 */
void
xdrstdio_create (XDR *xdrs, FILE *file, enum xdr_op op)
{
  xdrs->x_op = op;
  /* We have to add the const since the `struct xdr_ops' in `struct XDR'
     is not `const'.  */
  xdrs->x_ops = (struct xdr_ops *) &xdrstdio_ops;
  xdrs->x_private = (char *) file;
  xdrs->x_handy = 0;
  xdrs->x_base = 0;
}

/*
 * Destroy a stdio xdr stream.
 * Cleans up the xdr stream handle xdrs previously set up by xdrstdio_create.
 */
static void
xdrstdio_destroy (XDR *xdrs)
{
  (void) fflush ((FILE *) xdrs->x_private);
  /* xx should we close the file ?? */
}

static bool_t
xdrstdio_getlong (XDR *xdrs, long *lp)
{
  xdr_int32_t mycopy;

  if (fread ((char *) & mycopy, 4, 1, (FILE *) xdrs->x_private) != 1)
    return FALSE;
  *lp = (xdr_int32_t) xdr_ntohl (mycopy);
  return TRUE;
}

static bool_t
xdrstdio_putlong (XDR *xdrs, long *lp)
{
  long mycopy = xdr_htonl (*lp);
  lp = &mycopy;
  if (fwrite ((char *) lp, 4, 1, (FILE *) xdrs->x_private) != 1)
    return FALSE;
  return TRUE;
}

static bool_t
xdrstdio_getbytes (XDR *xdrs, char *addr, unsigned int len)
{
  if ((len != 0) && (fread (addr, (int) len, 1,
			    (FILE *) xdrs->x_private) != 1))
    return FALSE;
  return TRUE;
}

static bool_t
xdrstdio_putbytes (XDR *xdrs, char *addr, unsigned int len)
{
  if ((len != 0) && (fwrite (addr, (int) len, 1,
			     (FILE *) xdrs->x_private) != 1))
    return FALSE;
  return TRUE;
}

static unsigned int
xdrstdio_getpos (XDR *xdrs)
{
  return (unsigned int) ftell ((FILE *) xdrs->x_private);
}

static bool_t
xdrstdio_setpos (XDR *xdrs, unsigned int pos)
{
  return fseek ((FILE *) xdrs->x_private, (long) pos, 0) < 0 ? FALSE : TRUE;
}

static xdr_int32_t *
xdrstdio_inline (XDR *xdrs, int len)
{
  /*
   * Must do some work to implement this: must insure
   * enough data in the underlying stdio buffer,
   * that the buffer is aligned so that we can indirect through a
   * long *, and stuff this pointer in xdrs->x_buf.  Doing
   * a fread or fwrite to a scratch buffer would defeat
   * most of the gains to be had here and require storage
   * management on this buffer, so we don't do this.
   */
  return NULL;
}

static bool_t
xdrstdio_getint32 (XDR *xdrs, xdr_int32_t *ip)
{
  xdr_int32_t mycopy;

  if (fread ((char *) &mycopy, 4, 1, (FILE *) xdrs->x_private) != 1)
    return FALSE;
  *ip = xdr_ntohl (mycopy);
  return TRUE;
}

static bool_t
xdrstdio_putint32 (XDR *xdrs, xdr_int32_t *ip)
{
  xdr_int32_t mycopy = xdr_htonl (*ip);

  ip = &mycopy;
  if (fwrite ((char *) ip, 4, 1, (FILE *) xdrs->x_private) != 1)
    return FALSE;
  return TRUE;
}


/* Unimplemented routines - check the sunrpc still in glibc source,
 * and copy them here if you need them. Make sure that long long and
 * similar datatypes are portable, though, and that they produce the
 * right size on the system you intend to use it on.
 *
 * xdr_hyper, xdr_u_hyper, xdr_longlong_t, xdr_u_longlong_t
 * xdr_enum, xdr_bytes, xdr_union, xdr_netobj, xdr_wrapstring,
 * xdr_array.  / Erik Lindahl 2002-03-25
 */



/*
 * $Id: libxdrf.c,v 1.20.2.2 2005/10/10 19:43:46 lindahl Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#ifdef HAVE_FSEEKO
#define gmx_fseek(A,B,C) fseeko(A,B,C)
#define gmx_ftell(A) ftello(A)
#define gmx_off_t off_t
#else
#define gmx_fseek(A,B,C) fseek(A,B,C)
#define gmx_ftell(A) ftell(A)
#define gmx_off_t int
#endif


#define MAXID 256
static FILE *xdrfiles[MAXID];
static XDR *xdridptr[MAXID];
static char xdrmodes[MAXID];
static unsigned int cnt;


/*___________________________________________________________________________
 |
 | what follows are the C routines for opening, closing xdr streams
 | and the routine to read/write compressed coordinates together
 | with some routines to assist in this task (those are marked
 | static and cannot be called from user programs)
*/
#define MAXABS INT_MAX-2

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x):(y))
#endif
#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
static int magicints[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645,
    812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501,
    8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536,
    82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042,
    8388607, 10568983, 13316085, 16777216 };

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))


/*__________________________________________________________________________
 |
 | xdropen - open xdr file
 |
 | This versions differs from xdrstdio_create, because I need to know
 | the state of the file (read or write) so I can use xdr3dfcoord
 | in eigther read or write mode, and the file descriptor
 | so I can close the file (something xdr_destroy doesn't do).
 |
*/

int xdropen(XDR *xdrs, const char *filename, const char *type) {
    static int init_done = 0;
    enum xdr_op lmode;
    int xdrid;
    char newtype[5];

    if (init_done == 0) {
	for (xdrid = 1; xdrid < MAXID; xdrid++) {
	    xdridptr[xdrid] = NULL;
	}
	init_done = 1;
    }
    xdrid = 1;
    while (xdrid < MAXID && xdridptr[xdrid] != NULL) {
	xdrid++;
    }
    if (xdrid == MAXID) {
	return 0;
    }
    if (*type == 'w' || *type == 'W') {
            strcpy(newtype,"wb+");
	    lmode = XDR_ENCODE;
    } else if (*type == 'a' || *type == 'A') {
            strcpy(newtype,"ab+");
            lmode = XDR_ENCODE;
    } else {
            strcpy(newtype,"rb");
	    lmode = XDR_DECODE;
    }
    xdrfiles[xdrid] = fopen(filename, newtype);
    if (xdrfiles[xdrid] == NULL) {
	xdrs = NULL;
	return 0;
    }
    xdrmodes[xdrid] = *type;
    /* next test isn't usefull in the case of C language
     * but is used for the Fortran interface
     * (C users are expected to pass the address of an already allocated
     * XDR staructure)
     */
    if (xdrs == NULL) {
	xdridptr[xdrid] = (XDR *) malloc((size_t)sizeof(XDR));
	xdrstdio_create(xdridptr[xdrid], xdrfiles[xdrid], lmode);
    } else {
	xdridptr[xdrid] = xdrs;
	xdrstdio_create(xdrs, xdrfiles[xdrid], lmode);
    }
    return xdrid;
}

/*_________________________________________________________________________
 |
 | xdrclose - close a xdr file
 |
 | This will flush the xdr buffers, and destroy the xdr stream.
 | It also closes the associated file descriptor (this is *not*
 | done by xdr_destroy).
 |
*/
 
int xdrclose(XDR *xdrs) {
    int xdrid;
    
    if (xdrs == NULL) {
	fprintf(stderr, "xdrclose: passed a NULL pointer\n");
	exit(1);
    }
    for (xdrid = 1; xdrid < MAXID; xdrid++) {
	if (xdridptr[xdrid] == xdrs) {
	    
	    xdr_destroy(xdrs);
	    fclose(xdrfiles[xdrid]);
	    xdridptr[xdrid] = NULL;
	    return 1;
	}
    } 
    fprintf(stderr, "xdrclose: no such open xdr file\n");
    exit(1);
    
    /* to make some compilers happy: */
    return 0;    
}


