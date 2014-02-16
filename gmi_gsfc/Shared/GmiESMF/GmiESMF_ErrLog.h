

! The error logging may eventually evolve into a module based
! on the ESMF logger.  For now these macros provide simple
! traceback capability. 

#ifndef GMI_ErrLog_DONE


#define GMI_ErrLog_DONE

#ifdef RETURN_
#undef RETURN_
#endif
#ifdef VERIFY_
#undef VERIFY_
#endif
#ifdef ASSERT_
#undef ASSERT_
#endif

#ifdef IGNORE_
#undef IGNORE_
#endif
 
#define IGNORE_(a) continue

#ifdef I_AM_MAIN

#define VERIFY_(A) if(GMI_VRFY(A,Iam,__LINE__))call GmiESMF_Abort
#define ASSERT_(A) if(GMI_ASRT(A,Iam,__LINE__))call GmiESMF_Abort

#else

#ifdef ANSI_CPP

#define RETURN_(...)   if(GMI_RTRN(__VA_ARGS__,Iam,__LINE__,RC))return
#define VERIFY_(...)   if(GMI_VRFY(__VA_ARGS__,Iam,__LINE__,RC))return
#define ASSERT_(...)   if(GMI_ASRT(__VA_ARGS__,Iam,__LINE__,RC))return

#else

#define RETURN_(A)     if(GMI_RTRN(A,Iam,__LINE__,RC))return
#define VERIFY_(A)     if(GMI_VRFY(A,Iam,__LINE__,RC))return
#define ASSERT_(A)     if(GMI_ASRT(A,Iam,__LINE__,RC))return

#endif
#endif


#endif
