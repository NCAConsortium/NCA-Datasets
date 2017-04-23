      SUBROUTINE PK(ICALL,IDEF,THETA,IREV,EVTREC,NVNT,INDXS,IRGG,GG,          
     X NETAS)
      IMPLICIT DOUBLE PRECISION (A-Z)                                         
      REAL EVTREC                                                             
      SAVE
      INTEGER ICALL,IDEF,IREV,NVNT,INDXS,IRGG,NETAS                           
      DIMENSION IDEF(7,*),THETA(*),EVTREC(IREV,*),INDXS(*),GG(IRGG,31,*)
      COMMON/PROCM1/NEWIND                                                    
      INTEGER NEWIND                                                          
      COMMON/NMPRD1/IERPRD,NETEXT                                             
      COMMON/NMPRD2/ETEXT(3)                                                  
      INTEGER IERPRD,NETEXT                                                   
      CHARACTER*132 ETEXT                                                     
      COMMON/NMPRD4/CRCL,TVCL,GPCL,PPVCL,CL,TVV,V,K10,THALF,MRT,AREA
      COMMON/NMPRD4/AUMC,VSS,IPRE,Y,IRES,IWRE,A00033,A00034,A00035
      COMMON/NMPRD4/A00036,A00037,A00041,A00042,A00046,A00047,A00051
      COMMON/NMPRD4/A00053,A00058,D00001,D00002,D00035,D00036,C00032
      COMMON/NMPRD4/D00037,C00033,D00038,D00039,D00041,D00050
      COMMON/NMPRD4/BBBBBB(00960)
      COMMON/NMPRD7/ETA(30),EPS(30)                                           
      COMMON/ROCM12/MSEC,MFIRST
      INTEGER MSEC,MFIRST
      COMMON/PRCM00/MC0000(6),ME0000(6),MG0000(6),MT0000(6)                   
      INTEGER MC0000,ME0000,MG0000,MT0000                                     
      IF (ICALL.LE.1) THEN                                                    
      MC0000(1)=30
      ME0000(1)=30
      MG0000(1)=080
      MT0000(1)=40
      IDEF(1,001)= -9
      IDEF(1,002)= -1
      IDEF(1,003)=  0
      IDEF(1,004)=  0
      IDEF(2,003)=  0
      IDEF(2,004)=  0
      CALL GETETA(ETA)                                                        
      RETURN                                                                  
      ENDIF                                                                   
      IF (NEWIND.NE.2) THEN
       IF (ICALL.EQ.4) THEN
        CALL SIMETA(ETA)
       ELSE
        CALL GETETA(ETA)
       ENDIF
      ENDIF
      DOSE=EVTREC(NVNT,05)
      AGE=EVTREC(NVNT,13)
      WTKG=EVTREC(NVNT,14)
      SCR=EVTREC(NVNT,15)
      SEX=EVTREC(NVNT,16)
      B00001=140.D0-AGE 
      B00002=SCR*72.D0 
      B00003=0.15D0*SEX 
      B00004=B00001*WTKG/B00002 
      B00005=1.D0-B00003 
      CRCL=B00004*B00005 
      TVCL=THETA(01) 
      B00006=CRCL/80.D0 
      IF(B00006.EQ.0.D0.AND.THETA(02).LE.0.D0)THEN 
      IERPRD=1
      ETEXT(1)='PK SUBROUTINE: ERROR IN COMPUTATION'
      ETEXT(2)='ATTEMPT TO COMPUTE 0**POWER WITH POWER<=0.'
      RETURN
      ENDIF 
      IF(B00006.LT.0.D0)THEN 
      IERPRD=1
      ETEXT(1)='PK SUBROUTINE: ERROR IN COMPUTATION'
      ETEXT(2)='ATTEMPT TO COMPUTE BASE**POWER WITH BASE<0.'
      RETURN
      ENDIF 
      B00007=0.D0 
      IF(B00006.EQ.0.D0)THEN 
      B00007=1.D0 
      ENDIF 
      B00008=1.D0-B00007 
      B00009=B00006+B00007 
      B00010=B00009**THETA(02) 
      B00011=B00008*B00010 
      GPCL=TVCL*B00011 
      PPVCL=DEXP(ETA(01)) 
      CL=GPCL*PPVCL 
C                      A00033 = DERIVATIVE OF CL W.R.T. ETA(01)
      A00033=GPCL*PPVCL 
      TVV=THETA(03) 
      V=TVV 
      K10=CL/V 
      B00013=1.D0/V 
C                      A00035 = DERIVATIVE OF K10 W.R.T. ETA(01)
      A00035=B00013*A00033 
      B00014=DLOG(2.D0) 
      THALF=B00014/K10 
      B00015=-B00014/K10/K10 
C                      A00037 = DERIVATIVE OF THALF W.R.T. ETA(01)
      A00037=B00015*A00035 
      MRT=1.D0/K10 
      B00019=-1.D0/K10/K10 
C                      A00042 = DERIVATIVE OF MRT W.R.T. ETA(01)
      A00042=B00019*A00035 
      AREA=DOSE/CL 
      B00023=-DOSE/CL/CL 
C                      A00047 = DERIVATIVE OF AREA W.R.T. ETA(01)
      A00047=B00023*A00033 
      AUMC=AREA*MRT 
C                      A00052 = DERIVATIVE OF AUMC W.R.T. ETA(01)
      A00052=MRT*A00047 
C                      A00053 = DERIVATIVE OF AUMC W.R.T. ETA(01)
      A00053=AREA*A00042+A00052 
      VSS=V 
      GG(01,1,1)=CL    
      GG(01,02,1)=A00033
      GG(02,1,1)=V     
      IF (MSEC.EQ.1) THEN
C                      A00034 = DERIVATIVE OF A00033 W.R.T. ETA(01)
      A00034=GPCL*PPVCL 
C                      A00036 = DERIVATIVE OF A00035 W.R.T. ETA(01)
      A00036=B00013*A00034 
      B00016=B00014/K10/K10/K10 
C                      A00038 = DERIVATIVE OF B00015 W.R.T. ETA(01)
      A00038=B00016*A00035 
      B00017=B00014/K10/K10/K10 
C                      A00039 = DERIVATIVE OF B00015 W.R.T. ETA(01)
      A00039=B00017*A00035+A00038 
C                      A00040 = DERIVATIVE OF A00037 W.R.T. ETA(01)
      A00040=A00035*A00039 
C                      A00041 = DERIVATIVE OF A00037 W.R.T. ETA(01)
      A00041=B00015*A00036+A00040 
      B00020=1.D0/K10/K10/K10 
C                      A00043 = DERIVATIVE OF B00019 W.R.T. ETA(01)
      A00043=B00020*A00035 
      B00021=1.D0/K10/K10/K10 
C                      A00044 = DERIVATIVE OF B00019 W.R.T. ETA(01)
      A00044=B00021*A00035+A00043 
C                      A00045 = DERIVATIVE OF A00042 W.R.T. ETA(01)
      A00045=A00035*A00044 
C                      A00046 = DERIVATIVE OF A00042 W.R.T. ETA(01)
      A00046=B00019*A00036+A00045 
      B00024=DOSE/CL/CL/CL 
C                      A00048 = DERIVATIVE OF B00023 W.R.T. ETA(01)
      A00048=B00024*A00033 
      B00025=DOSE/CL/CL/CL 
C                      A00049 = DERIVATIVE OF B00023 W.R.T. ETA(01)
      A00049=B00025*A00033+A00048 
C                      A00050 = DERIVATIVE OF A00047 W.R.T. ETA(01)
      A00050=A00033*A00049 
C                      A00051 = DERIVATIVE OF A00047 W.R.T. ETA(01)
      A00051=B00023*A00034+A00050 
C                      A00054 = DERIVATIVE OF A00052 W.R.T. ETA(01)
      A00054=A00047*A00042 
C                      A00055 = DERIVATIVE OF A00052 W.R.T. ETA(01)
      A00055=MRT*A00051+A00054 
C                      A00056 = DERIVATIVE OF A00053 W.R.T. ETA(01)
      A00056=A00042*A00047 
C                      A00057 = DERIVATIVE OF A00053 W.R.T. ETA(01)
      A00057=AREA*A00046+A00056 
C                      A00058 = DERIVATIVE OF A00053 W.R.T. ETA(01)
      A00058=A00055+A00057 
      GG(01,02,02)=A00034
      ENDIF
      RETURN
      END
      SUBROUTINE ERROR (ICALL,IDEF,THETA,IREV,EVTREC,NVNT,INDXS,F,G,HH )      
      IMPLICIT DOUBLE PRECISION (A-Z)                                         
      REAL EVTREC                                                             
      SAVE
      INTEGER ICALL,IDEF,IREV,NVNT,INDXS                                      
      DIMENSION IDEF(*),THETA(*),EVTREC(IREV,*),INDXS(*),G(30,*)
      DIMENSION HH(30,*)
      COMMON/PROCM1/NEWIND                                                    
      INTEGER NEWIND                                                          
      COMMON/NMPRD4/CRCL,TVCL,GPCL,PPVCL,CL,TVV,V,K10,THALF,MRT,AREA
      COMMON/NMPRD4/AUMC,VSS,IPRE,Y,IRES,IWRE,A00033,A00034,A00035
      COMMON/NMPRD4/A00036,A00037,A00041,A00042,A00046,A00047,A00051
      COMMON/NMPRD4/A00053,A00058,D00001,D00002,D00035,D00036,C00032
      COMMON/NMPRD4/D00037,C00033,D00038,D00039,D00041,D00050
      COMMON/NMPRD4/BBBBBB(00960)
      COMMON/NMPRD7/ETA(30),EPS(30)                                           
      COMMON/ROCM12/MSEC,MFIRST
      INTEGER MSEC,MFIRST
      COMMON/ROCM17/NEWL2
      INTEGER NEWL2
      COMMON/PRCM00/MC0000(6),ME0000(6),MG0000(6),MT0000(6)                   
      INTEGER MC0000,ME0000,MG0000,MT0000                                     
      IF (ICALL.LE.1) THEN                                                    
      MC0000(2)=30
      ME0000(2)=30
      MG0000(2)=080
      MT0000(2)=40
      IDEF(2)=-1
      IDEF(3)=00
      RETURN
      ENDIF
      IF (ICALL.EQ.4) THEN
      IF (NEWL2.EQ.1) CALL SIMEPS(EPS)
      ENDIF
      D00001=G(01,1)
      DV=EVTREC(NVNT,09)
      IPRE=F 
      B00001=1.D0+EPS(01) 
      Y=F*B00001+EPS(02) 
C                      D00035 = DERIVATIVE OF Y W.R.T. ETA(01)
      D00035=B00001*D00001 
C                      C00032 = DERIVATIVE OF Y W.R.T. EPS(01)
      C00032=F 
C                      C00033 = DERIVATIVE OF Y W.R.T. EPS(02)
      C00033=1.D0 
C                      D00037 = DERIVATIVE OF C00032 W.R.T. ETA(01)
      D00037=D00001 
      IRES=DV-IPRE 
C                      D00038 = DERIVATIVE OF IRES W.R.T. ETA(01)
      D00038=-D00001 
      IWRE=IRES/IPRE 
      B00002=1.D0/IPRE 
C                      D00040 = DERIVATIVE OF IWRE W.R.T. ETA(01)
      D00040=B00002*D00038 
      B00003=-IRES/IPRE/IPRE 
C                      D00041 = DERIVATIVE OF IWRE W.R.T. ETA(01)
      D00041=B00003*D00001+D00040 
      G(01,1)=D00035  
      HH(01,1)=C00032 
      HH(02,1)=C00033 
      HH(01,02)=D00037
      IF (MSEC.EQ.1) THEN
      D00002=G(01,02)
C                      D00036 = DERIVATIVE OF D00035 W.R.T. ETA(01)
      D00036=B00001*D00002 
C                      D00039 = DERIVATIVE OF D00038 W.R.T. ETA(01)
      D00039=-D00002 
      B00004=-1.D0/IPRE/IPRE 
C                      D00042 = DERIVATIVE OF B00002 W.R.T. ETA(01)
      D00042=B00004*D00001 
C                      D00043 = DERIVATIVE OF D00040 W.R.T. ETA(01)
      D00043=D00038*D00042 
C                      D00044 = DERIVATIVE OF D00040 W.R.T. ETA(01)
      D00044=B00002*D00039+D00043 
      B00006=-1.D0/IPRE/IPRE 
C                      D00045 = DERIVATIVE OF B00003 W.R.T. ETA(01)
      D00045=B00006*D00038 
      B00007=IRES/IPRE/IPRE/IPRE 
C                      D00046 = DERIVATIVE OF B00003 W.R.T. ETA(01)
      D00046=B00007*D00001+D00045 
      B00008=IRES/IPRE/IPRE/IPRE 
C                      D00047 = DERIVATIVE OF B00003 W.R.T. ETA(01)
      D00047=B00008*D00001+D00046 
C                      D00048 = DERIVATIVE OF D00041 W.R.T. ETA(01)
      D00048=D00001*D00047 
C                      D00049 = DERIVATIVE OF D00041 W.R.T. ETA(01)
      D00049=B00003*D00002+D00048 
C                      D00050 = DERIVATIVE OF D00041 W.R.T. ETA(01)
      D00050=D00044+D00049 
      G(01,02)=D00036 
      ENDIF
      F=Y
      RETURN
      END
