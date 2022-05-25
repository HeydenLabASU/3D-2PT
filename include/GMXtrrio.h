
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



