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
      COMMON/NMPRD4/TVCL,GPCL,PPVCL,CL,TVV1,GPV1,V1,TVQ,PPVQ,Q,TVV2,V2
      COMMON/NMPRD4/K10,K12,K21,A0,A1,RT1,RT2,BETA,ALPHA,THALF,MRT
      COMMON/NMPRD4/AREA,AUMC,VSS,IPRE,Y,IRES,IWRE,A00033,A00034
      COMMON/NMPRD4/A00037,A00038,A00039,A00040,A00041,A00042,A00043
      COMMON/NMPRD4/A00044,A00045,A00048,A00047,A00046,A00049,A00052
      COMMON/NMPRD4/A00056,A00069,A00105,A00106,A00070,A00107,A00120
      COMMON/NMPRD4/A00156,A00157,A00121,A00158,A00226,A00236,A00235
      COMMON/NMPRD4/A00225,A00232,A00238,A00248,A00247,A00237,A00244
      COMMON/NMPRD4/A00249,A00253,A00255,A00263,A00264,A00256,A00265
      COMMON/NMPRD4/D00001,D00002,D00003,D00004,D00005,D00042,D00045
      COMMON/NMPRD4/D00044,D00041,D00043,C00032,D00047,D00046,C00033
      COMMON/NMPRD4/D00048,D00050,D00051,D00049,D00052,D00055,D00075
      COMMON/NMPRD4/D00076,D00056,D00079,BBBBBB(00902)
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
      CR60=EVTREC(NVNT,13)
      BSA=EVTREC(NVNT,14)
      ALB=EVTREC(NVNT,15)
      TVCL=THETA(01) 
      B00001=1.D0+THETA(02)*CR60 
      GPCL=TVCL*B00001 
      PPVCL=DEXP(ETA(01)) 
      CL=GPCL*PPVCL 
C                      A00033 = DERIVATIVE OF CL W.R.T. ETA(01)
      A00033=GPCL*PPVCL 
      TVV1=THETA(03) 
      B00003=ALB/34.D0 
      IF(B00003.EQ.0.D0.AND.THETA(04).LE.0.D0)THEN 
      IERPRD=1
      ETEXT(1)='PK SUBROUTINE: ERROR IN COMPUTATION'
      ETEXT(2)='ATTEMPT TO COMPUTE 0**POWER WITH POWER<=0.'
      RETURN
      ENDIF 
      IF(B00003.LT.0.D0)THEN 
      IERPRD=1
      ETEXT(1)='PK SUBROUTINE: ERROR IN COMPUTATION'
      ETEXT(2)='ATTEMPT TO COMPUTE BASE**POWER WITH BASE<0.'
      RETURN
      ENDIF 
      B00005=0.D0 
      IF(B00003.EQ.0.D0)THEN 
      B00005=1.D0 
      ENDIF 
      B00006=1.D0-B00005 
      B00007=B00003+B00005 
      B00008=B00007**THETA(04) 
      B00009=B00006*B00008 
      B00004=BSA*B00009 
      GPV1=TVV1*B00004 
      V1=GPV1 
      TVQ=THETA(05) 
      PPVQ=DEXP(ETA(02)) 
      Q=TVQ*PPVQ 
C                      A00037 = DERIVATIVE OF Q W.R.T. ETA(02)
      A00037=TVQ*PPVQ 
      TVV2=THETA(06) 
      V2=TVV2 
      K10=CL/V1 
      B00011=1.D0/V1 
C                      A00039 = DERIVATIVE OF K10 W.R.T. ETA(01)
      A00039=B00011*A00033 
      K12=Q/V1 
      B00012=1.D0/V1 
C                      A00041 = DERIVATIVE OF K12 W.R.T. ETA(02)
      A00041=B00012*A00037 
      K21=Q/V2 
      B00013=1.D0/V2 
C                      A00043 = DERIVATIVE OF K21 W.R.T. ETA(02)
      A00043=B00013*A00037 
      A0=K10*K21 
C                      A00045 = DERIVATIVE OF A0 W.R.T. ETA(01)
      A00045=K21*A00039 
C                      A00046 = DERIVATIVE OF A0 W.R.T. ETA(02)
      A00046=K10*A00043 
      A1=K10+K12+K21 
C                      A00052 = DERIVATIVE OF A1 W.R.T. ETA(02)
      A00052=A00043+A00041 
      B00016=A1*A1-4.D0*A0 
      B00018=DSQRT(B00016) 
      B00017=A1+B00018 
      RT1=B00017/2.D0 
C                      A00057 = DERIVATIVE OF B00016 W.R.T. ETA(02)
      A00057=A1*A00052 
C                      A00058 = DERIVATIVE OF B00016 W.R.T. ETA(01)
      A00058=A1*A00039 
C                      A00059 = DERIVATIVE OF B00016 W.R.T. ETA(02)
      A00059=A1*A00052+A00057 
C                      A00060 = DERIVATIVE OF B00016 W.R.T. ETA(01)
      A00060=A1*A00039+A00058 
C                      A00061 = DERIVATIVE OF B00016 W.R.T. ETA(02)
      A00061=-4.D0*A00046+A00059 
C                      A00062 = DERIVATIVE OF B00016 W.R.T. ETA(01)
      A00062=-4.D0*A00045+A00060 
      B00019=.5D0/B00018 
C                      A00063 = DERIVATIVE OF B00018 W.R.T. ETA(01)
      A00063=B00019*A00062 
C                      A00064 = DERIVATIVE OF B00018 W.R.T. ETA(02)
      A00064=B00019*A00061 
C                      A00067 = DERIVATIVE OF B00017 W.R.T. ETA(02)
      A00067=A00064+A00052 
C                      A00068 = DERIVATIVE OF B00017 W.R.T. ETA(01)
      A00068=A00063+A00039 
      B00020=1.D0/2.D0 
C                      A00069 = DERIVATIVE OF RT1 W.R.T. ETA(01)
      A00069=B00020*A00068 
C                      A00070 = DERIVATIVE OF RT1 W.R.T. ETA(02)
      A00070=B00020*A00067 
      B00026=A1*A1-4.D0*A0 
      B00028=DSQRT(B00026) 
      B00027=A1-B00028 
      RT2=B00027/2.D0 
C                      A00108 = DERIVATIVE OF B00026 W.R.T. ETA(02)
      A00108=A1*A00052 
C                      A00109 = DERIVATIVE OF B00026 W.R.T. ETA(01)
      A00109=A1*A00039 
C                      A00110 = DERIVATIVE OF B00026 W.R.T. ETA(02)
      A00110=A1*A00052+A00108 
C                      A00111 = DERIVATIVE OF B00026 W.R.T. ETA(01)
      A00111=A1*A00039+A00109 
C                      A00112 = DERIVATIVE OF B00026 W.R.T. ETA(02)
      A00112=-4.D0*A00046+A00110 
C                      A00113 = DERIVATIVE OF B00026 W.R.T. ETA(01)
      A00113=-4.D0*A00045+A00111 
      B00029=.5D0/B00028 
C                      A00114 = DERIVATIVE OF B00028 W.R.T. ETA(01)
      A00114=B00029*A00113 
C                      A00115 = DERIVATIVE OF B00028 W.R.T. ETA(02)
      A00115=B00029*A00112 
C                      A00118 = DERIVATIVE OF B00027 W.R.T. ETA(02)
      A00118=-A00115+A00052 
C                      A00119 = DERIVATIVE OF B00027 W.R.T. ETA(01)
      A00119=-A00114+A00039 
      B00030=1.D0/2.D0 
C                      A00120 = DERIVATIVE OF RT2 W.R.T. ETA(01)
      A00120=B00030*A00119 
C                      A00121 = DERIVATIVE OF RT2 W.R.T. ETA(02)
      A00121=B00030*A00118 
      Q00000=0.D0 
      IF(RT1.LT.RT2)THEN 
      Q00000=1.D0 
      ENDIF 
      BETA=Q00000*RT1 
      ALPHA=Q00000*RT2 
C                      A00169 = DERIVATIVE OF BETA W.R.T. ETA(01)
      A00169=Q00000*A00069 
C                      A00170 = DERIVATIVE OF BETA W.R.T. ETA(02)
      A00170=Q00000*A00070 
C                      A00171 = DERIVATIVE OF ALPHA W.R.T. ETA(01)
      A00171=Q00000*A00120 
C                      A00172 = DERIVATIVE OF ALPHA W.R.T. ETA(02)
      A00172=Q00000*A00121 
      Q00002=0.D0 
      Q00003=1.D0 
      IF(RT2.LT.RT1)THEN 
      Q00002=1.D0 
      Q00003=0.D0 
      ENDIF 
      B00044=Q00002*RT2+Q00003*BETA 
      B00045=Q00002*RT1+Q00003*ALPHA 
C                      A00189 = DERIVATIVE OF B00044 W.R.T. ETA(01)
      A00189=Q00002*A00120 
C                      A00190 = DERIVATIVE OF B00044 W.R.T. ETA(02)
      A00190=Q00002*A00121 
C                      A00191 = DERIVATIVE OF B00044 W.R.T. ETA(02)
      A00191=Q00003*A00170+A00190 
C                      A00192 = DERIVATIVE OF B00044 W.R.T. ETA(01)
      A00192=Q00003*A00169+A00189 
C                      A00193 = DERIVATIVE OF B00045 W.R.T. ETA(01)
      A00193=Q00002*A00069 
C                      A00194 = DERIVATIVE OF B00045 W.R.T. ETA(02)
      A00194=Q00002*A00070 
C                      A00195 = DERIVATIVE OF B00045 W.R.T. ETA(02)
      A00195=Q00003*A00172+A00194 
C                      A00196 = DERIVATIVE OF B00045 W.R.T. ETA(01)
      A00196=Q00003*A00171+A00193 
      BETA=B00044 
      ALPHA=B00045 
      B00046=DLOG(2.D0) 
      THALF=B00046/BETA 
      B00047=-B00046/BETA/BETA 
C                      A00225 = DERIVATIVE OF THALF W.R.T. ETA(02)
      A00225=B00047*A00191 
C                      A00226 = DERIVATIVE OF THALF W.R.T. ETA(01)
      A00226=B00047*A00192 
      B00050=BETA 
      MRT=1.D0/BETA 
      B00051=-1.D0/BETA/BETA 
C                      A00237 = DERIVATIVE OF MRT W.R.T. ETA(02)
      A00237=B00051*A00191 
C                      A00238 = DERIVATIVE OF MRT W.R.T. ETA(01)
      A00238=B00051*A00192 
      B00054=BETA 
      AREA=DOSE/CL 
      B00055=-DOSE/CL/CL 
C                      A00249 = DERIVATIVE OF AREA W.R.T. ETA(01)
      A00249=B00055*A00033 
      AUMC=AREA*MRT 
C                      A00254 = DERIVATIVE OF AUMC W.R.T. ETA(01)
      A00254=MRT*A00249 
C                      A00255 = DERIVATIVE OF AUMC W.R.T. ETA(01)
      A00255=AREA*A00238+A00254 
C                      A00256 = DERIVATIVE OF AUMC W.R.T. ETA(02)
      A00256=AREA*A00237 
      VSS=V1+V2 
      GG(01,1,1)=CL    
      GG(01,02,1)=A00033
      GG(02,1,1)=V1    
      GG(03,1,1)=Q     
      GG(03,03,1)=A00037
      GG(04,1,1)=V2    
      IF (MSEC.EQ.1) THEN
C                      A00034 = DERIVATIVE OF A00033 W.R.T. ETA(01)
      A00034=GPCL*PPVCL 
C                      A00038 = DERIVATIVE OF A00037 W.R.T. ETA(02)
      A00038=TVQ*PPVQ 
C                      A00040 = DERIVATIVE OF A00039 W.R.T. ETA(01)
      A00040=B00011*A00034 
C                      A00042 = DERIVATIVE OF A00041 W.R.T. ETA(02)
      A00042=B00012*A00038 
C                      A00044 = DERIVATIVE OF A00043 W.R.T. ETA(02)
      A00044=B00013*A00038 
C                      A00047 = DERIVATIVE OF A00045 W.R.T. ETA(02)
      A00047=A00039*A00043 
C                      A00048 = DERIVATIVE OF A00045 W.R.T. ETA(01)
      A00048=K21*A00040 
C                      A00049 = DERIVATIVE OF A00046 W.R.T. ETA(02)
      A00049=K10*A00044 
C                      A00056 = DERIVATIVE OF A00052 W.R.T. ETA(02)
      A00056=A00042+A00044 
C                      A00071 = DERIVATIVE OF A00057 W.R.T. ETA(02)
      A00071=A00052*A00052 
C                      A00072 = DERIVATIVE OF A00057 W.R.T. ETA(02)
      A00072=A1*A00056+A00071 
C                      A00073 = DERIVATIVE OF A00058 W.R.T. ETA(02)
      A00073=A00039*A00052 
C                      A00074 = DERIVATIVE OF A00058 W.R.T. ETA(01)
      A00074=A00039*A00039 
C                      A00075 = DERIVATIVE OF A00058 W.R.T. ETA(01)
      A00075=A1*A00040+A00074 
C                      A00076 = DERIVATIVE OF A00059 W.R.T. ETA(02)
      A00076=A00052*A00052 
C                      A00077 = DERIVATIVE OF A00059 W.R.T. ETA(02)
      A00077=A1*A00056+A00076 
C                      A00078 = DERIVATIVE OF A00059 W.R.T. ETA(02)
      A00078=A00072+A00077 
C                      A00079 = DERIVATIVE OF A00060 W.R.T. ETA(02)
      A00079=A00039*A00052 
C                      A00080 = DERIVATIVE OF A00060 W.R.T. ETA(01)
      A00080=A00039*A00039 
C                      A00081 = DERIVATIVE OF A00060 W.R.T. ETA(01)
      A00081=A1*A00040+A00080 
C                      A00082 = DERIVATIVE OF A00060 W.R.T. ETA(01)
      A00082=A00075+A00081 
C                      A00083 = DERIVATIVE OF A00060 W.R.T. ETA(02)
      A00083=A00073+A00079 
C                      A00084 = DERIVATIVE OF A00061 W.R.T. ETA(02)
      A00084=-4.D0*A00049 
C                      A00085 = DERIVATIVE OF A00061 W.R.T. ETA(02)
      A00085=A00078+A00084 
C                      A00086 = DERIVATIVE OF A00062 W.R.T. ETA(01)
      A00086=-4.D0*A00048 
C                      A00087 = DERIVATIVE OF A00062 W.R.T. ETA(02)
      A00087=-4.D0*A00047 
C                      A00088 = DERIVATIVE OF A00062 W.R.T. ETA(02)
      A00088=A00083+A00087 
C                      A00089 = DERIVATIVE OF A00062 W.R.T. ETA(01)
      A00089=A00082+A00086 
      B00025=-.5D0/B00018/B00018 
C                      A00090 = DERIVATIVE OF B00019 W.R.T. ETA(02)
      A00090=B00025*A00064 
C                      A00091 = DERIVATIVE OF B00019 W.R.T. ETA(01)
      A00091=B00025*A00063 
C                      A00092 = DERIVATIVE OF A00063 W.R.T. ETA(01)
      A00092=A00062*A00091 
C                      A00093 = DERIVATIVE OF A00063 W.R.T. ETA(02)
      A00093=A00062*A00090 
C                      A00094 = DERIVATIVE OF A00063 W.R.T. ETA(01)
      A00094=B00019*A00089+A00092 
C                      A00095 = DERIVATIVE OF A00063 W.R.T. ETA(02)
      A00095=B00019*A00088+A00093 
C                      A00096 = DERIVATIVE OF A00064 W.R.T. ETA(02)
      A00096=A00061*A00090 
C                      A00097 = DERIVATIVE OF A00064 W.R.T. ETA(02)
      A00097=B00019*A00085+A00096 
C                      A00101 = DERIVATIVE OF A00067 W.R.T. ETA(02)
      A00101=A00056+A00097 
C                      A00104 = DERIVATIVE OF A00068 W.R.T. ETA(01)
      A00104=A00040+A00094 
C                      A00105 = DERIVATIVE OF A00069 W.R.T. ETA(01)
      A00105=B00020*A00104 
C                      A00106 = DERIVATIVE OF A00069 W.R.T. ETA(02)
      A00106=B00020*A00095 
C                      A00107 = DERIVATIVE OF A00070 W.R.T. ETA(02)
      A00107=B00020*A00101 
C                      A00122 = DERIVATIVE OF A00108 W.R.T. ETA(02)
      A00122=A00052*A00052 
C                      A00123 = DERIVATIVE OF A00108 W.R.T. ETA(02)
      A00123=A1*A00056+A00122 
C                      A00124 = DERIVATIVE OF A00109 W.R.T. ETA(02)
      A00124=A00039*A00052 
C                      A00125 = DERIVATIVE OF A00109 W.R.T. ETA(01)
      A00125=A00039*A00039 
C                      A00126 = DERIVATIVE OF A00109 W.R.T. ETA(01)
      A00126=A1*A00040+A00125 
C                      A00127 = DERIVATIVE OF A00110 W.R.T. ETA(02)
      A00127=A00052*A00052 
C                      A00128 = DERIVATIVE OF A00110 W.R.T. ETA(02)
      A00128=A1*A00056+A00127 
C                      A00129 = DERIVATIVE OF A00110 W.R.T. ETA(02)
      A00129=A00123+A00128 
C                      A00130 = DERIVATIVE OF A00111 W.R.T. ETA(02)
      A00130=A00039*A00052 
C                      A00131 = DERIVATIVE OF A00111 W.R.T. ETA(01)
      A00131=A00039*A00039 
C                      A00132 = DERIVATIVE OF A00111 W.R.T. ETA(01)
      A00132=A1*A00040+A00131 
C                      A00133 = DERIVATIVE OF A00111 W.R.T. ETA(01)
      A00133=A00126+A00132 
C                      A00134 = DERIVATIVE OF A00111 W.R.T. ETA(02)
      A00134=A00124+A00130 
C                      A00135 = DERIVATIVE OF A00112 W.R.T. ETA(02)
      A00135=-4.D0*A00049 
C                      A00136 = DERIVATIVE OF A00112 W.R.T. ETA(02)
      A00136=A00129+A00135 
C                      A00137 = DERIVATIVE OF A00113 W.R.T. ETA(01)
      A00137=-4.D0*A00048 
C                      A00138 = DERIVATIVE OF A00113 W.R.T. ETA(02)
      A00138=-4.D0*A00047 
C                      A00139 = DERIVATIVE OF A00113 W.R.T. ETA(02)
      A00139=A00134+A00138 
C                      A00140 = DERIVATIVE OF A00113 W.R.T. ETA(01)
      A00140=A00133+A00137 
      B00035=-.5D0/B00028/B00028 
C                      A00141 = DERIVATIVE OF B00029 W.R.T. ETA(02)
      A00141=B00035*A00115 
C                      A00142 = DERIVATIVE OF B00029 W.R.T. ETA(01)
      A00142=B00035*A00114 
C                      A00143 = DERIVATIVE OF A00114 W.R.T. ETA(01)
      A00143=A00113*A00142 
C                      A00144 = DERIVATIVE OF A00114 W.R.T. ETA(02)
      A00144=A00113*A00141 
C                      A00145 = DERIVATIVE OF A00114 W.R.T. ETA(01)
      A00145=B00029*A00140+A00143 
C                      A00146 = DERIVATIVE OF A00114 W.R.T. ETA(02)
      A00146=B00029*A00139+A00144 
C                      A00147 = DERIVATIVE OF A00115 W.R.T. ETA(02)
      A00147=A00112*A00141 
C                      A00148 = DERIVATIVE OF A00115 W.R.T. ETA(02)
      A00148=B00029*A00136+A00147 
C                      A00151 = DERIVATIVE OF A00118 W.R.T. ETA(02)
      A00151=-A00148 
C                      A00152 = DERIVATIVE OF A00118 W.R.T. ETA(02)
      A00152=A00056+A00151 
C                      A00153 = DERIVATIVE OF A00119 W.R.T. ETA(02)
      A00153=-A00146 
C                      A00154 = DERIVATIVE OF A00119 W.R.T. ETA(01)
      A00154=-A00145 
C                      A00155 = DERIVATIVE OF A00119 W.R.T. ETA(01)
      A00155=A00040+A00154 
C                      A00156 = DERIVATIVE OF A00120 W.R.T. ETA(01)
      A00156=B00030*A00155 
C                      A00157 = DERIVATIVE OF A00120 W.R.T. ETA(02)
      A00157=B00030*A00153 
C                      A00158 = DERIVATIVE OF A00121 W.R.T. ETA(02)
      A00158=B00030*A00152 
C                      A00173 = DERIVATIVE OF A00169 W.R.T. ETA(01)
      A00173=Q00000*A00105 
C                      A00174 = DERIVATIVE OF A00169 W.R.T. ETA(02)
      A00174=Q00000*A00106 
C                      A00175 = DERIVATIVE OF A00170 W.R.T. ETA(02)
      A00175=Q00000*A00107 
C                      A00176 = DERIVATIVE OF A00171 W.R.T. ETA(01)
      A00176=Q00000*A00156 
C                      A00177 = DERIVATIVE OF A00171 W.R.T. ETA(02)
      A00177=Q00000*A00157 
C                      A00178 = DERIVATIVE OF A00172 W.R.T. ETA(02)
      A00178=Q00000*A00158 
C                      A00197 = DERIVATIVE OF A00189 W.R.T. ETA(01)
      A00197=Q00002*A00156 
C                      A00198 = DERIVATIVE OF A00189 W.R.T. ETA(02)
      A00198=Q00002*A00157 
C                      A00199 = DERIVATIVE OF A00190 W.R.T. ETA(02)
      A00199=Q00002*A00158 
C                      A00200 = DERIVATIVE OF A00191 W.R.T. ETA(02)
      A00200=Q00003*A00175 
C                      A00201 = DERIVATIVE OF A00191 W.R.T. ETA(02)
      A00201=A00199+A00200 
C                      A00202 = DERIVATIVE OF A00192 W.R.T. ETA(02)
      A00202=Q00003*A00174 
C                      A00203 = DERIVATIVE OF A00192 W.R.T. ETA(01)
      A00203=Q00003*A00173 
C                      A00204 = DERIVATIVE OF A00192 W.R.T. ETA(02)
      A00204=A00198+A00202 
C                      A00205 = DERIVATIVE OF A00192 W.R.T. ETA(01)
      A00205=A00197+A00203 
C                      A00206 = DERIVATIVE OF A00193 W.R.T. ETA(01)
      A00206=Q00002*A00105 
C                      A00207 = DERIVATIVE OF A00193 W.R.T. ETA(02)
      A00207=Q00002*A00106 
C                      A00208 = DERIVATIVE OF A00194 W.R.T. ETA(02)
      A00208=Q00002*A00107 
C                      A00209 = DERIVATIVE OF A00195 W.R.T. ETA(02)
      A00209=Q00003*A00178 
C                      A00210 = DERIVATIVE OF A00195 W.R.T. ETA(02)
      A00210=A00208+A00209 
C                      A00211 = DERIVATIVE OF A00196 W.R.T. ETA(02)
      A00211=Q00003*A00177 
C                      A00212 = DERIVATIVE OF A00196 W.R.T. ETA(01)
      A00212=Q00003*A00176 
C                      A00213 = DERIVATIVE OF A00196 W.R.T. ETA(02)
      A00213=A00207+A00211 
C                      A00214 = DERIVATIVE OF A00196 W.R.T. ETA(01)
      A00214=A00206+A00212 
      B00048=B00046/B00050/B00050/B00050 
C                      A00227 = DERIVATIVE OF B00047 W.R.T. ETA(02)
      A00227=B00048*A00191 
C                      A00228 = DERIVATIVE OF B00047 W.R.T. ETA(01)
      A00228=B00048*A00192 
      B00049=B00046/B00050/B00050/B00050 
C                      A00229 = DERIVATIVE OF B00047 W.R.T. ETA(02)
      A00229=B00049*A00191+A00227 
C                      A00230 = DERIVATIVE OF B00047 W.R.T. ETA(01)
      A00230=B00049*A00192+A00228 
C                      A00231 = DERIVATIVE OF A00225 W.R.T. ETA(02)
      A00231=A00191*A00229 
C                      A00232 = DERIVATIVE OF A00225 W.R.T. ETA(02)
      A00232=B00047*A00201+A00231 
C                      A00233 = DERIVATIVE OF A00226 W.R.T. ETA(01)
      A00233=A00192*A00230 
C                      A00234 = DERIVATIVE OF A00226 W.R.T. ETA(02)
      A00234=A00192*A00229 
C                      A00235 = DERIVATIVE OF A00226 W.R.T. ETA(02)
      A00235=B00047*A00204+A00234 
C                      A00236 = DERIVATIVE OF A00226 W.R.T. ETA(01)
      A00236=B00047*A00205+A00233 
      B00052=1.D0/B00054/B00054/B00054 
C                      A00239 = DERIVATIVE OF B00051 W.R.T. ETA(02)
      A00239=B00052*A00191 
C                      A00240 = DERIVATIVE OF B00051 W.R.T. ETA(01)
      A00240=B00052*A00192 
      B00053=1.D0/B00054/B00054/B00054 
C                      A00241 = DERIVATIVE OF B00051 W.R.T. ETA(02)
      A00241=B00053*A00191+A00239 
C                      A00242 = DERIVATIVE OF B00051 W.R.T. ETA(01)
      A00242=B00053*A00192+A00240 
C                      A00243 = DERIVATIVE OF A00237 W.R.T. ETA(02)
      A00243=A00191*A00241 
C                      A00244 = DERIVATIVE OF A00237 W.R.T. ETA(02)
      A00244=B00051*A00201+A00243 
C                      A00245 = DERIVATIVE OF A00238 W.R.T. ETA(01)
      A00245=A00192*A00242 
C                      A00246 = DERIVATIVE OF A00238 W.R.T. ETA(02)
      A00246=A00192*A00241 
C                      A00247 = DERIVATIVE OF A00238 W.R.T. ETA(02)
      A00247=B00051*A00204+A00246 
C                      A00248 = DERIVATIVE OF A00238 W.R.T. ETA(01)
      A00248=B00051*A00205+A00245 
      B00056=DOSE/CL/CL/CL 
C                      A00250 = DERIVATIVE OF B00055 W.R.T. ETA(01)
      A00250=B00056*A00033 
      B00057=DOSE/CL/CL/CL 
C                      A00251 = DERIVATIVE OF B00055 W.R.T. ETA(01)
      A00251=B00057*A00033+A00250 
C                      A00252 = DERIVATIVE OF A00249 W.R.T. ETA(01)
      A00252=A00033*A00251 
C                      A00253 = DERIVATIVE OF A00249 W.R.T. ETA(01)
      A00253=B00055*A00034+A00252 
C                      A00257 = DERIVATIVE OF A00254 W.R.T. ETA(01)
      A00257=A00249*A00238 
C                      A00258 = DERIVATIVE OF A00254 W.R.T. ETA(02)
      A00258=A00249*A00237 
C                      A00259 = DERIVATIVE OF A00254 W.R.T. ETA(01)
      A00259=MRT*A00253+A00257 
C                      A00260 = DERIVATIVE OF A00255 W.R.T. ETA(01)
      A00260=A00238*A00249 
C                      A00261 = DERIVATIVE OF A00255 W.R.T. ETA(01)
      A00261=AREA*A00248+A00260 
C                      A00262 = DERIVATIVE OF A00255 W.R.T. ETA(02)
      A00262=AREA*A00247 
C                      A00263 = DERIVATIVE OF A00255 W.R.T. ETA(01)
      A00263=A00259+A00261 
C                      A00264 = DERIVATIVE OF A00255 W.R.T. ETA(02)
      A00264=A00258+A00262 
C                      A00265 = DERIVATIVE OF A00256 W.R.T. ETA(02)
      A00265=AREA*A00244 
      GG(01,02,02)=A00034
      GG(03,03,03)=A00038
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
      COMMON/NMPRD4/TVCL,GPCL,PPVCL,CL,TVV1,GPV1,V1,TVQ,PPVQ,Q,TVV2,V2
      COMMON/NMPRD4/K10,K12,K21,A0,A1,RT1,RT2,BETA,ALPHA,THALF,MRT
      COMMON/NMPRD4/AREA,AUMC,VSS,IPRE,Y,IRES,IWRE,A00033,A00034
      COMMON/NMPRD4/A00037,A00038,A00039,A00040,A00041,A00042,A00043
      COMMON/NMPRD4/A00044,A00045,A00048,A00047,A00046,A00049,A00052
      COMMON/NMPRD4/A00056,A00069,A00105,A00106,A00070,A00107,A00120
      COMMON/NMPRD4/A00156,A00157,A00121,A00158,A00226,A00236,A00235
      COMMON/NMPRD4/A00225,A00232,A00238,A00248,A00247,A00237,A00244
      COMMON/NMPRD4/A00249,A00253,A00255,A00263,A00264,A00256,A00265
      COMMON/NMPRD4/D00001,D00002,D00003,D00004,D00005,D00042,D00045
      COMMON/NMPRD4/D00044,D00041,D00043,C00032,D00047,D00046,C00033
      COMMON/NMPRD4/D00048,D00050,D00051,D00049,D00052,D00055,D00075
      COMMON/NMPRD4/D00076,D00056,D00079,BBBBBB(00902)
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
      D00002=G(02,1)
      DV=EVTREC(NVNT,09)
      IPRE=F 
      B00001=DEXP(EPS(01)) 
      Y=F*B00001+EPS(02) 
C                      D00041 = DERIVATIVE OF Y W.R.T. ETA(02)
      D00041=B00001*D00002 
C                      D00042 = DERIVATIVE OF Y W.R.T. ETA(01)
      D00042=B00001*D00001 
C                      C00032 = DERIVATIVE OF Y W.R.T. EPS(01)
      C00032=F*B00001 
C                      C00033 = DERIVATIVE OF Y W.R.T. EPS(02)
      C00033=1.D0 
C                      D00046 = DERIVATIVE OF C00032 W.R.T. ETA(02)
      D00046=B00001*D00002 
C                      D00047 = DERIVATIVE OF C00032 W.R.T. ETA(01)
      D00047=B00001*D00001 
      IRES=DV-IPRE 
C                      D00048 = DERIVATIVE OF IRES W.R.T. ETA(01)
      D00048=-D00001 
C                      D00049 = DERIVATIVE OF IRES W.R.T. ETA(02)
      D00049=-D00002 
      IWRE=IRES/IPRE 
      B00002=1.D0/IPRE 
C                      D00053 = DERIVATIVE OF IWRE W.R.T. ETA(02)
      D00053=B00002*D00049 
C                      D00054 = DERIVATIVE OF IWRE W.R.T. ETA(01)
      D00054=B00002*D00048 
      B00003=-IRES/IPRE/IPRE 
C                      D00055 = DERIVATIVE OF IWRE W.R.T. ETA(01)
      D00055=B00003*D00001+D00054 
C                      D00056 = DERIVATIVE OF IWRE W.R.T. ETA(02)
      D00056=B00003*D00002+D00053 
      G(01,1)=D00042  
      G(02,1)=D00041  
      HH(01,1)=C00032 
      HH(02,1)=C00033 
      HH(01,02)=D00047
      HH(01,03)=D00046
      IF (MSEC.EQ.1) THEN
      D00003=G(01,02)
      D00004=G(02,02)
      D00005=G(02,03)
C                      D00043 = DERIVATIVE OF D00041 W.R.T. ETA(02)
      D00043=B00001*D00005 
C                      D00044 = DERIVATIVE OF D00042 W.R.T. ETA(02)
      D00044=B00001*D00004 
C                      D00045 = DERIVATIVE OF D00042 W.R.T. ETA(01)
      D00045=B00001*D00003 
C                      D00050 = DERIVATIVE OF D00048 W.R.T. ETA(01)
      D00050=-D00003 
C                      D00051 = DERIVATIVE OF D00048 W.R.T. ETA(02)
      D00051=-D00004 
C                      D00052 = DERIVATIVE OF D00049 W.R.T. ETA(02)
      D00052=-D00005 
      B00004=-1.D0/IPRE/IPRE 
C                      D00057 = DERIVATIVE OF B00002 W.R.T. ETA(01)
      D00057=B00004*D00001 
C                      D00058 = DERIVATIVE OF B00002 W.R.T. ETA(02)
      D00058=B00004*D00002 
C                      D00059 = DERIVATIVE OF D00053 W.R.T. ETA(02)
      D00059=D00049*D00058 
C                      D00060 = DERIVATIVE OF D00053 W.R.T. ETA(02)
      D00060=B00002*D00052+D00059 
C                      D00061 = DERIVATIVE OF D00054 W.R.T. ETA(02)
      D00061=D00048*D00058 
C                      D00062 = DERIVATIVE OF D00054 W.R.T. ETA(01)
      D00062=D00048*D00057 
C                      D00063 = DERIVATIVE OF D00054 W.R.T. ETA(02)
      D00063=B00002*D00051+D00061 
C                      D00064 = DERIVATIVE OF D00054 W.R.T. ETA(01)
      D00064=B00002*D00050+D00062 
      B00006=-1.D0/IPRE/IPRE 
C                      D00065 = DERIVATIVE OF B00003 W.R.T. ETA(02)
      D00065=B00006*D00049 
C                      D00066 = DERIVATIVE OF B00003 W.R.T. ETA(01)
      D00066=B00006*D00048 
      B00007=IRES/IPRE/IPRE/IPRE 
C                      D00067 = DERIVATIVE OF B00003 W.R.T. ETA(01)
      D00067=B00007*D00001+D00066 
C                      D00068 = DERIVATIVE OF B00003 W.R.T. ETA(02)
      D00068=B00007*D00002+D00065 
      B00008=IRES/IPRE/IPRE/IPRE 
C                      D00069 = DERIVATIVE OF B00003 W.R.T. ETA(01)
      D00069=B00008*D00001+D00067 
C                      D00070 = DERIVATIVE OF B00003 W.R.T. ETA(02)
      D00070=B00008*D00002+D00068 
C                      D00071 = DERIVATIVE OF D00055 W.R.T. ETA(02)
      D00071=D00001*D00070 
C                      D00072 = DERIVATIVE OF D00055 W.R.T. ETA(01)
      D00072=D00001*D00069 
C                      D00073 = DERIVATIVE OF D00055 W.R.T. ETA(01)
      D00073=B00003*D00003+D00072 
C                      D00074 = DERIVATIVE OF D00055 W.R.T. ETA(02)
      D00074=B00003*D00004+D00071 
C                      D00075 = DERIVATIVE OF D00055 W.R.T. ETA(01)
      D00075=D00064+D00073 
C                      D00076 = DERIVATIVE OF D00055 W.R.T. ETA(02)
      D00076=D00063+D00074 
C                      D00077 = DERIVATIVE OF D00056 W.R.T. ETA(02)
      D00077=D00002*D00070 
C                      D00078 = DERIVATIVE OF D00056 W.R.T. ETA(02)
      D00078=B00003*D00005+D00077 
C                      D00079 = DERIVATIVE OF D00056 W.R.T. ETA(02)
      D00079=D00060+D00078 
      G(01,02)=D00045 
      G(02,02)=D00044 
      G(02,03)=D00043 
      ENDIF
      F=Y
      RETURN
      END