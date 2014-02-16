
#include <stdio.h>
#include <math.h> 

/* Define float32/int32 as in HDF-4 for portability.  */

#if !defined(int32)
#    define  int32    int
#endif

#if !defined(uint32)
#    define  uint32   unsigned int
#endif

#if !defined(float32)
#    define  float32  float
#endif

#define MAXBITS 20
#define isZERO(I) (!((I) & 0x7fffffff))
#define isINF(I)  (!((I) & 0x7f800000))
#define isNAN(I)  (((I) & 0x7f800000)==0x7f800000)

#ifdef FAST_ISUNDEF
#  define isUndef(A)  ((A) == undef) /* cheap but not robust */
#else
#  define isUndef(A)  ((float32)fabs(undef-(A))<tol) /* robust but expensive */
#endif

#define SKIP(A,I)  (isZERO(I) || isUndef(A) || isNAN(I) || isINF(I)) 

#include "ShaveMantissa.h"   /* protype */

//========================================

float32 SetOffset(float32 minv, float32 maxv)
{
  float32   midv, mnabs,  range;

  range  = (maxv-minv);
  midv   = (maxv+minv)*0.5;
  mnabs  = abs(maxv)>abs(minv) ? abs(minv) : abs(maxv);

  return (range<mnabs) ? midv : midv*(mnabs/range);
}

//========================================

/*
//---------------------------------------------------------------------------
//BOP
//
// !ROUTINE: ShaveMantissa32 - Degrades Precison of 32-bit Float Array 
//
// !INTERFACE:
*/
int ShaveMantissa32 ( a, ain, len, xbits, has_undef, undef, chunksize )

/*
// !INPUT PARAMETERS:
*/

int32   len;        /* size of a[] */
int     xbits;      /* number of bits to excludes from mantissa */
int     has_undef;  /* whether field has missing (undef) values */ 
int32   chunksize;  /* find mid range over chunksize chunks     */ 
float32 undef;      /* missing (undefined) value */
float32 ain[];      /* input array */

/*
// !OUTPUT PARAMETERS:
*/

float32 a[];    // output "shaved" array; can share storage with ain[]

/*
// !DESCRIPTION:  
//
//  This routine returns a lower precision version of the input array
//  {\tt a}, reducing the precision of the mantissa by {\tt xbits}.
//  (The number of bits retained is {\tt nbits = 24 - xbits}, given that
//  32-bit floats in IEEE representation reserves only 24 bits for the
//  mantissa. The purpose of this precision degradation is to promote
//  internal GZIP compression by HDF-4.
//
//  This algorithm produces very similar results as the standard GRIB
//  encoding with fixed number of bits ({\tt nbits = 24 - xbits}),
//  and power of 2 binary scaling.
//
// !REVISION HISTORY:
//
//  08Dec2006  Suarez    First version.
//  09Dec2006  da Silva  Minor revisions for IA64 build, prologue.
//  11Dec2006  da Silva  Merged with Max's newest version handling Inf, NAN
//                       and bug fixes. 
//  18Dec2006  Suarez    Eliminated test for overflow, which I verified was
//                       unnecessary. Treat zeros like undef. Eliminated a
//                       leftover conversion of nan to undef. Corrected macro
//                       for inf.  MAJOR correction to keep code from hanging
//                       when it encountered an invalid value.
//  26Dec2006  Suarez    Added selection of offset based on range and zero
//                       offset. Restored treatment of fields beginning with
//                       invalids and of all constant fields. Fixed bug that
//                       was not cpying input to output arrays.
//EOP
//---------------------------------------------------------------------------
*/

{
  float32   maxv, minv,  offset, *b, *begnxt, *last, tol;
  uint32    round, mask, *i;

  /* sanity checks */

  if ( len < 1 || xbits < 1 ) {
    fprintf(stderr,
	    "ShaveMantissa32: Bad length of mask bits: len= %d, xbits= %d\n", 
	    len, xbits );
    return 1;
  }

  if ( xbits > MAXBITS ) {
    fprintf(stderr,
	    "ShaveMantissa32: Shaving too many bits: %d; maximum allowed is %d\n",
	    xbits, MAXBITS );
    return 2;
  }

  /* if user has not chosen an undef, pick one */

  if ( !has_undef ) undef = (float32) HUGE_VAL;

  /* miscelaneous static combinations */

  tol   = 0.0001*undef;
  mask  = 0xFFFFFFFF<<xbits--;
  round = 1         <<xbits  ;
  last  = &a[len-1];

  // Check input array

  b = a;
  if(ain!=a) {
    if(abs(ain-a)<len) {
      fprintf(stderr,"ShaveMantissa32: Overlapping arrays");
      return 3;
    }

    while(a<=last) *a++=*ain++;
  }

  // Loop over chunks

  while(b<=last) {

    // The beginning of the chunk after the current one

    begnxt = b + chunksize;
    if(begnxt>last) begnxt = last+1;

    // Move to first valid value

    a = b-1;
    while(++a < begnxt) {
      i = (uint32 *)a;
      if(!SKIP(*a,*i)) {maxv=*a; minv=*a; break;}
    }

    // Empty field

    if(a==begnxt) {b=begnxt; continue;}

    // Find man and max valid values

    while(++a<begnxt) {
      i = (uint32 *)a;
      if(!SKIP(*a,*i)) {
	if(*a<minv) minv=*a;
	if(*a>maxv) maxv=*a;
      }
    }

    // Constant field

    if(minv==maxv) {b=begnxt; continue;}

    // Find optimum offset

    offset = SetOffset(minv,maxv);

    // Shave

    a = b-1;
    while(++a<begnxt) {
      i = (uint32 *)a;
      if(!SKIP(*a,*i)) {
	*a -= offset;
	*i = ((*i+round) & mask);
	*a += offset;
      }
    }

    // Prepare for next chunk

    b = begnxt;

  } // End chunk loop

  return 0;
}

//========================================

//    Simple hook for FORTRAN interface.

int SHAVEMANTISSA32 (float32 *a, float32 *ain, int32 *len, int *xbits, 
                     int *has_undef, float32 *undef, int32 *chunksize)
{return (int)ShaveMantissa32(a,ain,*len,*xbits,*has_undef,*undef,*chunksize);}

int SHAVEMANTISSA32_ (float32 *a, float32 *ain, int32 *len, int *xbits, 
                      int *has_undef, float32 *undef, int32 *chunksize)
{return (int)ShaveMantissa32(a,ain,*len,*xbits,*has_undef,*undef,*chunksize);}

int shavemantissa32 (float32 *a, float32 *ain, int32 *len, int *xbits, 
                     int *has_undef, float32 *undef, int32 *chunksize)
{return (int)ShaveMantissa32(a,ain,*len,*xbits,*has_undef,*undef,*chunksize);}

int shavemantissa32_ (float32 *a, float32 *ain, int32 *len, int *xbits, 
                      int *has_undef, float32 *undef, int32 *chunksize)
{return (int)ShaveMantissa32(a,ain,*len,*xbits,*has_undef,*undef,*chunksize);}





