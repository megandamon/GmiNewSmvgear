#ifndef sysAIX
#include <malloc.h>

void QSswap(long long a[], int b[], int i, int j, int m, int n)
{
  long long T;
  int k, S;
  T = a[i]; 
  a[i] = a[j];
  a[j] = T;
  for(k=0;k<m*n;k+=m) {
    S = b[i+k]; 
    b[i+k] = b[j+k];
    b[j+k] = S;
  }
}


   void QuickSort(long long a[], int b[], int l, int r, int m, int n)
   {
        int len;
        int i;
        int j;
        long long v;

        len = r-l+1;
	if (len<=1) return;
	if (len==2) {if (a[l]>a[r]) QSswap(a,b,l,r,m,n); return;}
	if (len==3) {
         if (a[l  ]>a[l+1]) QSswap(a,b,l  ,l+1,m,n);
         if (a[l+1]>a[l+2]) QSswap(a,b,l+1,l+2,m,n);
         if (a[l  ]>a[l+1]) QSswap(a,b,l  ,l+1,m,n);
         return;
        }

        j = r;
        i = l-1;
        v = a[r];
        for(;;)  {
         while(a[++i]<v);
         while(a[--j]>v);
         if (j<i) break;
         QSswap(a,b,i,j,m,n);
        }
        QSswap(a,b,i,r,m,n);
        QuickSort(a,b,l,j,m,n);
        QuickSort(a,b,i+1,r,m,n);
    }


    void QSswap0(long long a[], int i, int j)
    {
            long long T;
            T = a[i]; 
            a[i] = a[j];
            a[j] = T;
    }


   void QuickSort0(long long a[], int l, int r)
   {
        int len;
        int i;
        int j;
        long long v;

        len = r-l+1;
	if (len<=1) return;
	if (len==2) {if (a[l]>a[r]) QSswap0(a,l,r); return;}
	if (len==3) {
         if (a[l  ]>a[l+1]) QSswap0(a,l  ,l+1);
         if (a[l+1]>a[l+2]) QSswap0(a,l+1,l+2);
         if (a[l  ]>a[l+1]) QSswap0(a,l  ,l+1);
         return;
        }

        j = r;
        i = l-1;
        v = a[r];
        for(;;)  {
         while(a[++i]<v);
         while(a[--j]>v);
         if (j<i) break;
         QSswap0(a,i,j);
        }
        QSswap0(a,i,r);
        QuickSort0(a,l,j);
        QuickSort0(a,i+1,r);
    }

// FORTRAN INTERFACES

void QSORT0 (long long a[],          int *r        ) { (void)QuickSort0(a,  0,*r-1      ); }
void QSORT2 (long long a[], int b[], int *r, int *n) { (void)QuickSort (a,b,0,*r-1,*r,*n); }


void QSORT0S (int a[],          int *r        ) 
{
  long long *c;
  int i;
  c=(long long *)malloc(*r*sizeof(long long));
  for(i=0;i<*r;i++) c[i]=(long long)(a[i]);
  (void)QuickSort0(c,0,*r-1);
}


void QSORT2S (int a[], int b[], int *r, int *n) 
{
  long long *c;
  int i;
  c=(long long *)malloc(*r*sizeof(long long));
  for(i=0;i<*r;i++) c[i]=(long long)(a[i]);
  (void)QuickSort (c,b,0,*r-1,*r,*n);
}

//  EXTRA ALIASES FOR OTHER LOADERS

void QSORT0_(long long a[],          int *r        ) { (void)QSORT0 (a,  r  ); }
void QSORT2_(long long a[], int b[], int *r, int *n) { (void)QSORT2 (a,b,r,n); }
void QSORT0S_(int      a[],          int *r        ) { (void)QSORT0S(a,  r  ); }
void QSORT2S_(int      a[], int b[], int *r, int *n) { (void)QSORT2S(a,b,r,n); }

void qsort0_(long long a[],          int *r        ) { (void)QSORT0 (a,  r  ); }
void qsort2_(long long a[], int b[], int *r, int *n) { (void)QSORT2 (a,b,r,n); }
void qsort0s_(int      a[],          int *r        ) { (void)QSORT0S(a,  r  ); }
void qsort2s_(int      a[], int b[], int *r, int *n) { (void)QSORT2S(a,b,r,n); }


void qsort0 (long long a[],          int *r        ) { (void)QSORT0 (a,  r  ); }
void qsort2 (long long a[], int b[], int *r, int *n) { (void)QSORT2 (a,b,r,n); }
void qsort0s (int      a[],          int *r        ) { (void)QSORT0S(a,  r  ); }
void qsort2s (int      a[], int b[], int *r, int *n) { (void)QSORT2S(a,b,r,n); }

#endif
