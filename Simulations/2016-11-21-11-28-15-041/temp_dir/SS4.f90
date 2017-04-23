!*********************************COPYRIGHT******************************************
!                                                                                   !
!       THE NONMEM SYSTEM MAY BE DISTRIBUTED ONLY BY ICON DEVELOPMENT               !
!       SOLUTIONS.                                                                  !
!                                                                                   !
!       COPYRIGHT BY ICON DEVELOPMENT SOLUTIONS                                     !
!       2009-2013 ALL RIGHTS RESERVED.                                              !
!                                                                                   !
!       DO NOT ATTEMPT TO MODIFY CODE WITHOUT FIRST CONSULTING WITH                 !
!       ICON DEVELOPMENT SOLUTIONS.                                                 !
!                                                                                   !
!************************************************************************************
!
!-----------------------------HISTORY------------------------------------------------
! VERSION     : NONMEM VII
! AUTHOR      : ALISON J. BOECKMANN
! CREATED ON  : AUG/1984
! LANGUAGE    : FORTRAN 90/95
! LAST UPDATE : JUL/2008 - COMMON BLOCKS REPLACED WITH MODULES
!               NOV/2008 - INTRODUCED HEADER INFORMATIONS AND RESTRUCTURED AS PER
!                          THE NONMEM STANDARDS
!               FEB/2010 - CHANGED SIZES TO PRSIZES
!               FEB/2011 - INTEGRATED 7.2BETA5.8B MODIFICATIONS
!
!----------------------------- SS4.F90 ----------------------------------------------
!
! SUBROUTINE SS
!
! DESCRIPTION : Bi-exponential model with absorption and urine compartments
!
! ARGUMENTS   : NONE
!
! CALLED BY   : SSS - Supervisor of steady state
!
! CALLS       : NONE
!
! ALGORITHM   : - Save first derivatives before getting second derivatives
!               - Second derivatives
!               - If SSID /= 0, then Normal entry; Else
!                 - Set SSID=4 and RETURN
!               - Normal entry
!                 - Check for errors in basic PK parameters
!                 - Branch if constant infusion rather than multiple doses
!                 - Compute common terms and expressions
!                 - Set up coefficients FA,GC,GP,HC,HP,FC,FP for the equations
!                 - Loop on ALPHA,BETA: RT(I)=ALPHA, RT(J)=BETA within loop
!                   - Compute denominators for fractions
!                 - Set X0,X1,X2 according to dose type and compartment
!                 - Compute multiple bolus doses
!                 - Else, compute multiple infusions
!                 - Reduce system of three equations to two equations
!                 - Solve the simultaneous linear equations
!                 - Compute constant or background infusion
!                 - Estimate derivatives
!                   - Compute common terms and expressions
!                   - Loop on ALPHA, BETA
!                   - Find multiple bolus doses
!                   - Else, find multiple infusions
!                   - Set up the simultaneous linear equations
!                   - Compute constant or background infusion
!                   - Save first derivatives for ETA(K)
!                 - Loop for second derivatives
!                   - Compute multiple bolus doses
!                   - Compute multiple infusions
!                   - Solve simultaneous equations
!                   - Get constant or background infusion
!
! MODULES USED: PRSIZES,NMPRD_INT,PRCOM_INT,PRMOD_INT,PROCM_INT,NMPRD_REAL,PRCM_REAL,
!               PRCOM_REAL,PRMOD_REAL,PROCM_REAL,NMPRD_CHAR,PRCOM_LOG   
!
! CONTAINS    : NONE
!
! LOCAL'S     : A10,A10E,A12,A12E,A21,A21E,A21R,A21RE,AA,AAE,D10,D10E,DA,DADIF,DADIFE,
!               DADS,DADSE,DAE,DDIF,DDIFA,DDIFAE,DDIFE,DEDIF,DEDIFE,DEN,DPRD,DPROD,
!               DPRODE,DR,DRE,EX,EXA,EXADL,EXADLE,EXAE,EXDEL,EXDELE,EXE,FA,FAE,FC,
!               FCA,FCAE,FCE,FCR,FCRE,FP,FPA,FPAE,FPE,FPR,FPRE,G,G2,GC,GCE,GE,GG,GP,
!               GPE,GR,GRE,HC,HCE,HP,HPE,I,J,K,L,MDA,MDE,MDI,NEWA,NEWAE,NEWC,NEWP,P,
!               PE,PPE,RT,RTE,S,SE,T0,T0E,T1,T1E,T2,T2E,TEMP,TMA,TMB,TMC,TMD,X0,X0E,
!               X1,X1E,X2,X2E,Y,Y1,Y1E,Y2,Y2E,YE,Z,Z1,Z1E,Z2,Z2E,ZE
!
! Save first derivatives before getting second derivatives
!               A101,A121,A211,A21R1,AA1,D101,DA1,DADIF1,DADS1,DDIF1,DDIFA1,DEDIF1,
!               DENE1,DPRD1,DPROD1,DR1,EX1,EXA1,EXADL1,EXDEL1,FA1,FC1,FP1,FPA1,FPR1,
!               GE1,GR1,MDA1,MDE1,MDI1,NEWA1,NEWC1,NEWP1,RT1,S1,T01,T11,T21,TEMP1,
!               TMA1,TMB1,TMC1,TMD1,X11,X21,Y11,Y21,YE1,Z11,Z21,ZE1
!
! Second derivatives
!               A102,A122,A212,A21R2,AA2,D,D1,D102,D2,DA2,DADIF2,DADS2,DDIF2,DDIFA2,DEDIF2,
!               DENE2,DPRD2,DPROD2,DR2,EX2,EXA2,EXADL2,EXDEL2,FA2,FC2,FCA2,FCR2,FP2,FPA2,
!               FPR2,GC2,GE2,GP2,GR2,HC2,HP2,MDA2,MDE2,MDI2,NEWA2,NEWC2,NEWP2,P2,RT2,S2,
!               SSRT,SSRT1,SSRT2,T02,T12,T22,TEMP2,TMA2,TMB2,TMC2,TMD2,X02,X12,X22,Y12,Y22,
!               YE2,Z12,Z22,ZE2
!
!---------------------------- END OF HEADER -----------------------------------------
!
      SUBROUTINE SS
!      
      USE PRSIZES,    ONLY: ISIZE,DPSIZE,PE
! INTEGER
      USE NMPRD_INT,  ONLY: IERPRD   
      USE PRCOM_INT,  ONLY: ITSC,NPETAS,SSID,SSC,SV
      USE PRMOD_INT,  ONLY: IKE,IK12,IK21,KP,KC,IKA,KA
      USE PROCM_INT,  ONLY: IDXETA
! REAL
      USE NMPRD_REAL, ONLY: SQINFN
      USE PRCM_REAL,  ONLY: SSAE2,SSRE2,SDELE2,RHOE2     
      USE PRCOM_REAL, ONLY: D2DTE,DT,DTE,RHO,RHOE,SDEL,SDELE,SSA,SSAE,SSR,SSRE,G3, &
                            ZERO,ONE
      USE PRMOD_REAL, ONLY: TWO,FOUR
      USE PROCM_REAL, ONLY: AMNT,DAETA,D2AETA      
! CHARACTER
      USE NMPRD_CHAR, ONLY: ETEXT
! LOGICAL          
      USE PRCOM_LOG,  ONLY: NOETAS,SECOND
!
      IMPLICIT NONE
!      
      SAVE
!           
!------------------------------------------------------------------------------------
!     COMMON /CM73/ REPS,TOL,INFNTY,SQINFN,EPS
!     DOUBLE PRECISION REPS,TOL,INFNTY,SQINFN
!     REAL EPS
!     COMMON /NMPRD1/ IERPRD,NETEXT
!     COMMON /NMPRD2/ ETEXT(3)
!     INTEGER IERPRD,NETEXT
!     CHARACTER*132 ETEXT
!     COMMON /PRCOM0/ NP,NBP,YFORM
!     COMMON /PRCOM0/ MAXKF,IFORM
!     COMMON /PRCOM0/ IDC,IDO,MAXIC,ISV,IINST,ITURN
!     COMMON /PRCOM0/ JTIME,JCONT,JEVENT,JAMT,JRATE,JSS,JDELTA
!     COMMON /PRCOM0/ JCOMPT,JCOMPF,JERROR,SSC,KREC,JMORE,JDUM
!     COMMON /PRCOM1/ NOETAS,SECOND
!     COMMON /PRCOM2/ IBF,IRR,IS,ID,ITSC,IFR,ILAG
!     COMMON /PRCOM3/ ITRANS,IRGG,IREV,NPETAS,NPEPS
!     COMMON /PRCOM4/ G3,HH,DELTA,DT,DTE
!     COMMON /PRCOM4/ YMULT,ZERO,ONE,XR,XD,TSTART,DTSTAR
!     COMMON /PRCOM4/ DDELTA,D2DELT,ADTE,D2ADTE
!     COMMON /PRCOM5/ ISPEC,DCTR,BETA,DD
!     COMMON /PRCOM5/ IP,IPOOL,IHEAD,INEXT,IBACK,SV
!     COMMON /PRCOM6/ IA,IAA,IAEA,IRA,IREA,IDA,IDEA,R,RE
!     COMMON /PRCOM6/ RHO,RHOE,SDEL,SDELE,SSA,SSAE,SSR,SSRE
!     COMMON /PRCOM6/ SAMT,SDEL1
!     COMMON /PRCOM6/ I2AEA,I2REA,I2DEA,R2E,D2DTE,D2TSTA
!     COMMON /PRCM6A/ SSAE2,SSRE2,SDELE2,RHOE2
!     COMMON /PRCOM7/ ADVID,SSID
!     DOUBLE PRECISION DELTA,G3,HH
!     DOUBLE PRECISION SDELE,RHOE,SSAE,SSRE,YMULT,ZERO,XR,XD
!     DOUBLE PRECISION SSAE2,SSRE2,SDELE2,RHOE2
!     DOUBLE PRECISION ONE,TSTART,DTSTAR(PE)
!     DOUBLE PRECISION DDELTA(PE),D2DELT(PE,PE),ADTE(PE),D2ADTE(PE,PE)
!     DOUBLE PRECISION IA(90),IAA(90),IAEA(90,PE),IRA(90),IDA(90)
!     DOUBLE PRECISION IREA(90,PE),IDEA(90,PE),R(PC),RE(PC,PE)
!     DOUBLE PRECISION I2REA(90,PE,PE),I2DEA(90,PE,PE),I2AEA(90,PE,PE)
!     DOUBLE PRECISION R2E(PC,PE,PE),D2DTE(PE,PE)
!     DOUBLE PRECISION D2TSTA(PE,PE)
!     DOUBLE PRECISION DT,DTE(PE),RHO,SDEL,SSA,SSR
!     DOUBLE PRECISION SAMT,SDEL1
!     DIMENSION SDELE(PE),RHOE(PE),SSAE(PE),SSRE(PE)
!     DIMENSION SSAE2(PE,PE),SSRE2(PE,PE),SDELE2(PE,PE),RHOE2(PE,PE)
!     DIMENSION G3(PG+1,PE+1,PE+1),HH(PE,PE)
!     INTEGER IRGG,IREV,ITRANS,NPETAS,NPEPS
!     INTEGER JCONT,JTIME,JEVENT,JAMT,JRATE,JSS,JDELTA
!     INTEGER JCOMPT,JCOMPF,JERROR
!     INTEGER IDC,IDO,NP,NBP,SSC,KREC,JMORE,JDUM
!     INTEGER ISV(PC),IBF(PC),IRR(PC),SV(PC)
!     INTEGER IINST(PC),ITURN(PC),ITSC,IFR,ILAG(PC),IS(PC),ID(PC)
!     INTEGER ADVID,SSID,MAXKF,IFORM(PG+1),YFORM,MAXIC
!     INTEGER BETA(90),IPOOL(90),IP,IHEAD,INEXT(90),IBACK(90)
!     INTEGER ISPEC,DD(90),DCTR
!     LOGICAL NOETAS,SECOND
!     COMMON /PROCM4/ A,DAETA,D2AETA
!     DOUBLE PRECISION A,DAETA,D2AETA
!     DIMENSION A(PC),DAETA(PC,PE),D2AETA(PC,PE,PE)
!     COMMON /PRBI/ TWO,FOUR,IKE,IK12,IK21,IV,KP,KC
!     COMMON /PRBI/ IKA,KA
! FOR LVOUT FEATURE
!     COMMON /PROCM5/ NACTIV,M(0:PE)
!     INTEGER NACTIV,M
! LOCAL VARIABLES START HERE
!     INTEGER KP,KC,IKE,IK12,IK21,IV
!     INTEGER KA,IKA
!     DOUBLE PRECISION TWO,FOUR
!------------------------------------------------------------------------------------
!
! Local Variables 
!
      INTEGER(KIND=ISIZE) :: I,J,K,L
! 
      REAL(KIND=DPSIZE)   :: A10,A10E,A12,A12E,A21,A21E,A21R(2),A21RE(2),AA,AAE,D10,     &
                             D10E,DA,DADIF(2),DADIFE(2),DADS(2),DADSE(2),DAE,DDIF(2),    &
                             DDIFA(2),DDIFAE(2),DDIFE(2),DEDIF(2),DEDIFE(2),DEN,DPRD,    &
                             DPROD,DPRODE,DR(2),DRE(2),EX(2),EXA,EXADL,EXADLE,EXAE,      &
                             EXDEL(2),EXDELE(2),EXE(2),FA,FAE,FC,FCA,FCAE,FCE,FCR(2),    &
                             FCRE(2),FP,FPA,FPAE,FPE,FPR(2),FPRE(2),G,G2,GC(2),GCE(2),   &
                             GE,GG,GP(2),GPE(2),GR(2),GRE(2),HC(2),HCE(2),HP(2),HPE(2),  &
                             MDA(2),MDE(2),MDI(2),NEWA,NEWAE,NEWC,NEWP,P,PPE,RT(2),      &
                             RTE(2),S,SE,T0,T0E,T1,T1E,T2,T2E,TEMP,TMA(2),TMB(2),TMC(2), &
                             TMD(2),X0,X0E,X1,X1E,X2,X2E,Y,Y1,Y1E,Y2,Y2E,YE,Z,Z1,Z1E,    &
                             Z2,Z2E,ZE
! Save first derivatives before getting second derivatives
      REAL(KIND=DPSIZE)   :: A101(PE),A121(PE),A211(PE),A21R1(2,PE),AA1(PE),D101(PE),    &
                             DA1(PE),DADIF1(2,PE),DADS1(2,PE),DDIF1(2,PE),DDIFA1(2,PE),  &
                             DEDIF1(2,PE),DENE1(PE),DPRD1(PE),DPROD1(PE),DR1(2,PE),      &
                             EX1(2,PE),EXA1(PE),EXADL1(PE),EXDEL1(2,PE),FA1(PE),FC1(PE), &
                             FP1(PE),FPA1(PE),FPR1(2,PE),GE1(PE),GR1(2,PE),MDA1(2,PE),   &
                             MDE1(2,PE),MDI1(2,PE),NEWA1(PE),NEWC1(PE),NEWP1(PE),        &
                             RT1(2,PE),S1(PE),T01(PE),T11(PE),T21(PE),TEMP1(PE),         &
                             TMA1(2,PE),TMB1(2,PE),TMC1(2,PE),TMD1(2,PE),X11(PE),X21(PE),&
                             Y11(PE),Y21(PE),YE1(PE),Z11(PE),Z21(PE),ZE1(PE)             
! Second derivatives
      REAL(KIND=DPSIZE)   :: A102,A122,A212,A21R2(2),AA2,D,D1(PE),D102,D2,DA2,DADIF2(2), &
                             DADS2(2),DDIF2(2),DDIFA2(2),DEDIF2(2),DENE2,DPRD2,DPROD2,   &
                             DR2(2),EX2(2),EXA2,EXADL2,EXDEL2(2),FA2,FC2,FCA2,FCR2(2),   &
                             FP2,FPA2,FPR2(2),GC2(2),GE2,GP2(2),GR2(2),HC2(2),HP2(2),    &
                             MDA2(2),MDE2(2),MDI2(2),NEWA2,NEWC2,NEWP2,P2,RT2(2),S2,SSRT,&
                             SSRT1(PE),SSRT2,T02,T12,T22,TEMP2,TMA2(2),TMB2(2),TMC2(2),  &
                             TMD2(2),X02,X12,X22,Y12,Y22,YE2,Z12,Z22,ZE2        
!     
      G2(I,J,K)=G3(I,IDXETA(J)+1,IDXETA(K)+1) ! Statement function
      GG(I,J)=G3(I,IDXETA(J-1)+1,1)
!
! Steady state variables:
! From PRED: 
! SSC     = Compartment number for dose
! DT      = Capital Delta (Dosing interval)
! SDEL    = Small Delta (Duration of infusion)
! SSA     = Dose amount
! SSR     = Rate of constant infusion
! RHO     = Rate of multiple infusion
! RHOE    = Partial derivative of RHO wrt ETA(K)
! DTE(K)  = Partial derivative of DT wrt ETA(K)
! SELE(K) = Partial derivative of SDEL wrt ETA(K)
! SSAE(K) = Partial derivative of SSA wrt ETA(K)
! SSRE(K) = Partial derivative of SSR wrt ETA(K)
!
      IF (SSID == 0) THEN
        SSID=4; GO TO 999  
      END IF
      
! Normal entry
      A10=GG(ITSC,1)*GG(IKE,1)
      A12=GG(ITSC,1)*GG(IK12,1)
      A21=GG(ITSC,1)*GG(IK21,1)
      AA=ONE
!      
      IF (SV(KA) /= 0) AA=GG(ITSC,1)*GG(IKA,1)
!      
! Check for errors in basic PK parameters
      IF (GG(ITSC,1) <= ZERO)  THEN
        ETEXT(2)='PK PARAMETER FOR TIME SCALE IS NON-POSITIVE'
        IERPRD=1; GO TO 999  
      ELSE IF (A10 <= ZERO) THEN
        ETEXT(2)='PK PARAMETER FOR K IS NON-POSITIVE'
        IERPRD=1; GO TO 999  
      ELSE IF (A12 <= ZERO) THEN
        ETEXT(2)='PK PARAMETER FOR K23 IS NON-POSITIVE'
        IERPRD=1; GO TO 999  
      ELSE IF (A21 <= ZERO) THEN
        ETEXT(2)='PK PARAMETER FOR K32 IS NON-POSITIVE'
        IERPRD=1; GO TO 999   
      ELSE IF (SV(KA) /= 0) THEN
        IF (AA <= ZERO) THEN
          ETEXT(2)='PK PARAMETER FOR KA IS NON-POSITIVE'
          IERPRD=1; GO TO 999   
        END IF
      END IF
      
      DA=ONE/AA
      D10=ONE/A10
      G=A12+A10
      Z=A12/A21
      Y=G/A21
      
! Branch if constant infusion rather than multiple doses
      IF (DT /= ZERO) THEN
        S=A10+A12+A21
        P=A21*A10
        IF (S >= SQINFN) THEN
          ETEXT(1)='K+K23+K32 IS TOO LARGE. THE CHARACTERISTIC EQUATION'
          ETEXT(2)='CANNOT BE SOLVED.'
          IERPRD=1; GO TO 999   
        END IF
        D=S*S-FOUR*P
        TEMP=SQRT(D)
        RT(1)=(S+TEMP)/TWO
        RT(2)=TWO*P/(S+TEMP)
        IF (RT(2) == ZERO) THEN
          ETEXT(1)='A ROOT OF THE CHARACTERISTIC EQUATION IS ZERO BECAUSE'
          ETEXT(2)='K*K32 IS MUCH SMALLER THAN (K+K23+K32)**2.'
          ETEXT(3)='PERHAPS K OR K32 IS VERY SMALL, OR K23 IS VERY LARGE.'
          IERPRD=1; GO TO 999  
        END IF
! Common terms and expressions
        DO I=1,2
          DR(I)=ONE/RT(I)
          EX(I)=EXP(-(DT+SDEL)*RT(I))
          A21R(I)=A21-RT(I)
          GR(I)=G-RT(I)
          IF (RHO == ZERO) CYCLE
          EXDEL(I)=EXP(-DT*RT(I))
        END DO
! Set up coefficients FA,GC,GP,HC,HP,FC,FP for the equations:
!       A(KA)=A(KA)*FA+CA
!       A(KC)=A(KC)*GC+A(KP)*HC+A(KA)*FC+CC
!       A(KP)=A(KP)*GP+A(KP)*HP+A(KP)*FP+CP
! Constant terms CA,CC,CP will be developed later as part of X'S
! Loop on ALPHA,BETA: RT(I)=ALPHA, RT(J)=BETA within loop
        J=3
        DO I=1,2
          J=J-1
          DDIF(I)=ONE/(RT(J)-RT(I))    ! Denominators for fractions
          IF (RHO /= ZERO) THEN
            DEDIF(I)=D10*DDIF(I)
            DADIF(I)=DR(I)*DDIF(I)
            DADS(I)=DADIF(I)*DDIF(I)
          END IF
          GC(I)=A21R(I)*DDIF(I)*EX(I)
          HC(I)=A21*DDIF(I)*EX(I)
          GP(I)=A12*DDIF(I)*EX(I)
          HP(I)=GR(I)*DDIF(I)*EX(I)
        END DO
!       
        IF (SV(KA) /= 0) THEN
          EXA=EXP(-(DT+SDEL)*AA)
          IF (RHO /= ZERO) EXADL=EXP(-DT*AA)         
          DDIFA(1:2)=ONE/(AA-RT(1:2))
          FCR(1:2)=AA*A21R(1:2)*DDIFA(1:2)*DDIF(1:2)*EX(1:2)
          FPR(1:2)=DDIFA(1:2)*DDIF(1:2)*EX(1:2)
!          
          DPROD=DDIFA(1)*DDIFA(2)
          FA=EXA
          FCA=AA*(A21-AA)*DPROD*EXA
          FPA=DPROD*EXA
          FC=FCA+FCR(1)+FCR(2)
          FP=(FPA+FPR(1)+FPR(2))*AA*A12
        END IF
        
! Set X0,X1,X2 according to dose type and compartment
        X0=ZERO; X1=ZERO; X2=ZERO
!            
        IF (RHO == ZERO) THEN   ! Multiple bolus doses
          IF (SSC == KA) X0=-SSA
          IF (SSC == KC) X1=-SSA
          IF (SSC == KP) X2=-SSA
        ELSE                    ! Multiple infusions
          MDE(1)=DEDIF(1)*EXDEL(1)
          MDE(2)=DEDIF(2)*EXDEL(2)
          MDA(1)=DADS(1)*EX(1)
          MDA(2)=DADS(2)*EX(2)
          IF (SSC /= KA) THEN
            IF (SSC == KP) THEN 
              T1 = (G+A21R(1))*MDE(1)                   &
                  +(G+A21R(2))*MDE(2)                   &
                  -(A21R(1)*A21+GR(1)*A21)*MDA(1)       &
                  -(A21R(2)*A21+GR(2)*A21)*MDA(2)
              X1 =-RHO*T1
              T2 = (A12+Y*GR(1))*MDE(1)                 &
                  +(A12+Y*GR(2))*MDE(2)                 &
                  -(GR(1)*GR(1)+A12*A21)*MDA(1)         &
                  -(GR(2)*GR(2)+A12*A21)*MDA(2)
              X2 =-RHO*T2
            ELSE
              T1 = (A12+A21R(1))*MDE(1)                 &
                  +(A12+A21R(2))*MDE(2)                 &
                  -(A21R(1)*A21R(1)+A12*A21)*MDA(1)     &
                  -(A21R(2)*A21R(2)+A12*A21)*MDA(2)
              X1 =-RHO*T1
              T2 = (A12+Z*GR(1))*MDE(1)                 &
                  +(A12+Z*GR(2))*MDE(2)                 &
                  -(GR(1)*A12+A21R(1)*A12)*MDA(1)       &
                  -(GR(2)*A12+A21R(2)*A12)*MDA(2)
              X2 =-RHO*T2
            END IF 
          ELSE
            T0 = DA*(EXADL-EXA)
            X0 =-RHO*T0
            MDI(1) = DDIF(1)*EXDEL(1)
            MDI(2) = DDIF(2)*EXDEL(2)
            DPRD   = DPROD*(EXADL-EXA)
            TMA(1) = (A21R(1)+A12)*D10+A21R(1)*DDIFA(1)
            TMA(2) = (A21R(2)+A12)*D10+A21R(2)*DDIFA(2)
            TMB(1) = AA*A21R(1)*A21R(1)+A12*AA*A21
            TMB(2) = AA*A21R(2)*A21R(2)+A12*AA*A21
            T1 = (A21-AA)*DPRD                          &
                +TMA(1)*MDI(1)                          &
                +TMA(2)*MDI(2)                          &
                -TMB(1)*DDIFA(1)*MDA(1)                 &
                -TMB(2)*DDIFA(2)*MDA(2)
            X1 =-RHO*T1
            TMC(1) = (A12+Z*GR(1))*D10+A12*DDIFA(1)
            TMC(2) = (A12+Z*GR(2))*D10+A12*DDIFA(2)
            TMD(1) = AA*A12*(A21R(1)+GR(1))
            TMD(2) = AA*A12*(A21R(2)+GR(2))
            T2 = A12*DPRD                               &
                +TMC(1)*MDI(1)                          &
                +TMC(2)*MDI(2)                          &
                -TMD(1)*DDIFA(1)*MDA(1)                 &
                -TMD(2)*DDIFA(2)*MDA(2)
            X2 =-RHO*T2
! Reduce system of three equations to two equations
!           X1 = A(KC)*Y1+A(KP)*Z2
!           X2 = A(KC)*Y2+A(KP)*Z2
          END IF
        END IF
!    
        IF (SV(KA) /= 0) THEN
          IF (FA-ONE == ZERO) THEN
            ETEXT(1)='ELIMINATION FROM SOME COMPARTMENT IS NEGLIGIBLE.'
            ETEXT(2)='STEADY-STATE IS NOT ACHIEVABLE.'
            IERPRD=1; GO TO 999   
          END IF
          NEWA=X0/(FA-ONE)
          AMNT(KA)=AMNT(KA)+NEWA
          X1=X1-FC*NEWA
          X2=X2-FP*NEWA
        END IF
!        
! Solve the simultaneous linear equations
        Y1=GC(1)+GC(2)-ONE
        Y2=GP(1)+GP(2)
        Z1=HC(1)+HC(2)
        Z2=HP(1)+HP(2)-ONE
        DEN=ONE/(Y1*Z2-Y2*Z1)
        NEWC=(X1*Z2-X2*Z1)*DEN
        NEWP=(Y1*X2-Y2*X1)*DEN
        AMNT(KC)=AMNT(KC)+NEWC
        AMNT(KP)=AMNT(KP)+NEWP
      END IF
      
! Constant or background infusion
      IF (SSR /= ZERO) THEN
        IF (SSC /= KP) THEN 
          IF (SSC == KA) AMNT(KA)=AMNT(KA)+SSR*DA
          SSRT=SSR*D10
          AMNT(KC)=AMNT(KC)+SSRT
          AMNT(KP)=AMNT(KP)+SSRT*Z
        ELSE
          AMNT(KC)=AMNT(KC)+SSR*D10
          AMNT(KP)=AMNT(KP)+SSR*Y*D10
        END IF  
      END IF
!      
      IF (NOETAS) GO TO 999
!      
! Derivatives
      DO K=1,NPETAS
        A10E = GG(ITSC,1)*GG(IKE,K+1)+GG(ITSC,K+1)*GG(IKE,1)
        A12E = GG(ITSC,1)*GG(IK12,K+1)+GG(ITSC,K+1)*GG(IK12,1)
        A21E = GG(ITSC,1)*GG(IK21,K+1)+GG(ITSC,K+1)*GG(IK21,1)
        AAE  = GG(ITSC,1)*GG(IKA,K+1)+GG(ITSC,K+1)*GG(IKA,1)
        DAE  =-DA*AAE/AA
        D10E =-D10*A10E/A10
        GE   = A12E+A10E
        ZE   = (A12E-Z*A21E)/A21
        YE   = (GE-Y*A21E)/A21
        IF (DT /= ZERO) THEN
          SE    = A10E+A12E+A21E
          PPE   = A21*A10E+A21E*A10
          D1(K) = TWO*S*SE-FOUR*PPE
          TEMP1(K) = D1(K)/(TWO*TEMP)
          RTE(1) = (SE+TEMP1(K))/TWO
          RTE(2) = (TWO*PPE-RT(2)*(SE+TEMP1(K)))/(S+TEMP)
! Common terms and expressions
          DO I=1,2
            DRE(I)  =-DR(I)/RT(I)*RTE(I)
            EXE(I)  =-EX(I)*((DT+SDEL)*RTE(I)+(DTE(K)+SDELE(K))*RT(I))
            A21RE(I)= A21E-RTE(I)
            GRE(I)  = GE-RTE(I)
            IF (RHO == ZERO) CYCLE
            EXDELE(I) =-EXDEL(I)*(DT*RTE(I)+DTE(K)*RT(I))
          END DO
!    
          J=3
          DO I=1,2   ! Loop on ALPHA, BETA
            J=J-1
            DDIFE(I) =-(RTE(J)-RTE(I))*DDIF(I)*DDIF(I)
            IF (RHO /= ZERO) THEN
              DEDIFE(I) = D10*DDIFE(I)+D10E*DDIF(I)
              DADIFE(I) = DR(I)*DDIFE(I)+DRE(I)*DDIF(I)
              DADSE(I)  = DADIF(I)*DDIFE(I)+DADIFE(I)*DDIF(I)
            END IF
            GCE(I) = A21R(I)*DDIF(I)*EXE(I)+A21R(I)*DDIFE(I)*EX(I)+A21RE(I)*DDIF(I)*EX(I)
            HCE(I) = A21*DDIF(I)*EXE(I)+A21*DDIFE(I)*EX(I)+A21E*DDIF(I)*EX(I)
            GPE(I) = A12*DDIF(I)*EXE(I)+A12*DDIFE(I)*EX(I)+A12E*DDIF(I)*EX(I)
            HPE(I) = GR(I)*DDIF(I)*EXE(I)+GR(I)*DDIFE(I)*EX(I)+GRE(I)*DDIF(I)*EX(I)
          END DO
!    
          IF (SV(KA) /= 0) THEN
            EXAE = -EXA*((DT+SDEL)*AAE+(DTE(K)+SDELE(K))*AA)
            IF (RHO /= ZERO) EXADLE =-EXADL*(DT*AAE+DTE(K)*AA)
            DDIFAE(1:2)=-DDIFA(1:2)*DDIFA(1:2)*(AAE-RTE(1:2))
            FPRE(1:2)  = DDIFA(1:2)*DDIF(1:2)*EXE(1:2)+DDIFA(1:2)*DDIFE(1:2)*EX(1:2)    &
                        +DDIFAE(1:2)*DDIF(1:2)*EX(1:2)
            FCRE(1:2)  = FPR(1:2)*(AA*A21RE(1:2)+AAE*A21R(1:2))+FPRE(1:2)*AA*A21R(1:2)
!           
            DPRODE = DDIFA(1)*DDIFAE(2)+DDIFAE(1)*DDIFA(2)
            FAE    = EXAE
            FCAE   = AA*(A21-AA)*DPROD*EXAE+AA*(A21-AA)*DPRODE*EXA                      &
                    +AA*(A21E-AAE)*DPROD*EXA+AAE*(A21-AA)*DPROD*EXA
            FPAE   = DPROD*EXAE+DPRODE*EXA
            FCE    = FCAE+FCRE(1)+FCRE(2)
            FPE    = (FPA+FPR(1)+FPR(2))*(AA*A12E+AAE*A12)+(FPAE+FPRE(1)+FPRE(2))*AA*A12
          END IF
!              
          X0E=ZERO; X1E=ZERO; X2E=ZERO
!   
          IF (RHO == ZERO) THEN   ! Multiple bolus doses
            IF (SSC == KA) X0E=-SSAE(K)
            IF (SSC == KC) X1E=-SSAE(K)
            IF (SSC == KP) X2E=-SSAE(K)
          ELSE                    ! Multiple infusions
            MDE1(1,K) = DEDIF(1)*EXDELE(1)+DEDIFE(1)*EXDEL(1)
            MDE1(2,K) = DEDIF(2)*EXDELE(2)+DEDIFE(2)*EXDEL(2)
            MDA1(1,K) = DADS(1)*EXE(1)+DADSE(1)*EX(1)
            MDA1(2,K) = DADS(2)*EXE(2)+DADSE(2)*EX(2)
            IF (SSC /= KA) THEN
              IF (SSC == KP) THEN 
                T1E = (G+A21R(1))*MDE1(1,K)+(GE+A21RE(1))*MDE(1)                &
                     +(G+A21R(2))*MDE1(2,K)+(GE+A21RE(2))*MDE(2)                &
                     -(A21R(1)*A21+GR(1)*A21)*MDA1(1,K)                         &
                     -(A21RE(1)*A21+GRE(1)*A21+A21R(1)*A21E+GR(1)*A21E)*MDA(1)  &
                     -(A21R(2)*A21+GR(2)*A21)*MDA1(2,K)                         &
                     -(A21RE(2)*A21+GRE(2)*A21+A21R(2)*A21E+GR(2)*A21E)*MDA(2)
                X1E =-RHO*T1E-RHOE(K)*T1
                T2E = (A12+Y*GR(1))*MDE1(1,K)+(A12E+Y*GRE(1)+YE*GR(1))*MDE(1)   &
                     +(A12+Y*GR(2))*MDE1(2,K)+(A12E+Y*GRE(2)+YE*GR(2))*MDE(2)   &
                     -(GR(1)*GR(1)+A12*A21)*MDA1(1,K)                           &
                     -(TWO*GR(1)*GRE(1)+A12*A21E+A12E*A21)*MDA(1)               &
                     -(GR(2)*GR(2)+A12*A21)*MDA1(2,K)                           &
                     -(TWO*GR(2)*GRE(2)+A12*A21E+A12E*A21)*MDA(2)
                X2E =-RHO*T2E-RHOE(K)*T2
              ELSE
                T1E = (A12+A21R(1))*MDE1(1,K)+(A12E+A21RE(1))*MDE(1)            &
                     +(A12+A21R(2))*MDE1(2,K)+(A12E+A21RE(2))*MDE(2)            &
                     -(A21R(1)*A21R(1)+A12*A21)*MDA1(1,K)                       &
                     -(TWO*A21R(1)*A21RE(1)+A12*A21E+A12E*A21)*MDA(1)           &
                     -(A21R(2)*A21R(2)+A12*A21)*MDA1(2,K)                       &
                     -(TWO*A21R(2)*A21RE(2)+A12*A21E+A12E*A21)*MDA(2)
                X1E =-RHO*T1E-RHOE(K)*T1
                T2E = (A12+Z*GR(1))*MDE1(1,K)+(A12E+Z*GRE(1)+ZE*GR(1))*MDE(1)   &
                     +(A12+Z*GR(2))*MDE1(2,K)+(A12E+Z*GRE(2)+ZE*GR(2))*MDE(2)   &
                     -(GR(1)+A21R(1))*A12*MDA1(1,K)                             &
                     -((GR(1)+A21R(1))*A12E+(GRE(1)+A21RE(1))*A12)*MDA(1)       &
                     -(GR(2)+A21R(2))*A12*MDA1(2,K)                             &
                     -((GR(2)+A21R(2))*A12E+(GRE(2)+A21RE(2))*A12)*MDA(2)
                X2E =-RHO*T2E-RHOE(K)*T2
              END IF  
            ELSE
              T0E = DA*(EXADLE-EXAE)+DAE*(EXADL-EXA)
              X0E =-RHO*T0E-RHOE(K)*T0
              MDI1(1,K) = DDIF(1)*EXDELE(1)+DDIFE(1)*EXDEL(1)
              MDI1(2,K) = DDIF(2)*EXDELE(2)+DDIFE(2)*EXDEL(2)
              DPRD1(K)  = DPRODE*(EXADL-EXA)+DPROD*(EXADLE-EXAE)
              TMA1(1,K) = (A21RE(1)+A12E)*D10+A21RE(1)*DDIFA(1)                 &
                         +(A21R(1)+A12)*D10E+A21R(1)*DDIFAE(1)
              TMA1(2,K) = (A21RE(2)+A12E)*D10+A21RE(2)*DDIFA(2)                 &
                         +(A21R(2)+A12)*D10E+A21R(2)*DDIFAE(2)
              TMB1(1,K) = AAE*A21R(1)*A21R(1)+A12E*AA*A21                       &       
                         +AA*A21RE(1)*A21R(1)+A12*AAE*A21                       &
                         +AA*A21R(1)*A21RE(1)+A12*AA*A21E
              TMB1(2,K) = AAE*A21R(2)*A21R(2)+A12E*AA*A21                       &
                         +AA*A21RE(2)*A21R(2)+A12*AAE*A21                       &
                         +AA*A21R(2)*A21RE(2)+A12*AA*A21E
              T1E = (A21E-AAE)*DPRD+(A21-AA)*DPRD1(K)                           &
                   +TMA1(1,K)*MDI(1)+TMA(1)*MDI1(1,K)                           &
                   +TMA1(2,K)*MDI(2)+TMA(2)*MDI1(2,K)                           &
                   -TMB1(1,K)*DDIFA(1)*MDA(1)                                   &
                   -TMB(1)*DDIFAE(1)*MDA(1)-TMB(1)*DDIFA(1)*MDA1(1,K)           &
                   -TMB1(2,K)*DDIFA(2)*MDA(2)                                   &    
                   -TMB(2)*DDIFAE(2)*MDA(2)-TMB(2)*DDIFA(2)*MDA1(2,K)
              X1E =-RHO*T1E-RHOE(K)*T1
              TMC1(1,K) = (A12E+ZE*GR(1)+Z*GRE(1))*D10+A12E*DDIFA(1)            &
                         +(A12+Z*GR(1))*D10E+A12*DDIFAE(1)
              TMC1(2,K) = (A12E+ZE*GR(2)+Z*GRE(2))*D10+A12E*DDIFA(2)            &
                         +(A12+Z*GR(2))*D10E+A12*DDIFAE(2)
              TMD1(1,K) = AAE*A12*(A21R(1)+GR(1))                               &
                         +AA*A12E*(A21R(1)+GR(1))                               &
                         +AA*A12*(A21RE(1)+GRE(1))
              TMD1(2,K) = AAE*A12*(A21R(2)+GR(2))                               &
                         +AA*A12E*(A21R(2)+GR(2))                               &
                         +AA*A12*(A21RE(2)+GRE(2))
              T2E = A12E*DPRD+A12*DPRD1(K)                                      &
                   +TMC1(1,K)*MDI(1)+TMC(1)*MDI1(1,K)                           &
                   +TMC1(2,K)*MDI(2)+TMC(2)*MDI1(2,K)                           &
                   -TMD1(1,K)*DDIFA(1)*MDA(1)                                   &
                   -TMD(1)*DDIFAE(1)*MDA(1)-TMD(1)*DDIFA(1)*MDA1(1,K)           &
                   -TMD1(2,K)*DDIFA(2)*MDA(2)                                   &
                   -TMD(2)*DDIFAE(2)*MDA(2)-TMD(2)*DDIFA(2)*MDA1(2,K)
              X2E =-RHO*T2E-RHOE(K)*T2
            END IF 
          END IF 
! Set up the simultaneous linear eqns
          IF (SV(KA) /= 0) THEN
            NEWAE =(X0E-NEWA*FAE)/(FA-ONE)
            DAETA(KA,K) = DAETA(KA,K)+NEWAE
            X1E = X1E-FC*NEWAE-FCE*NEWA
            X2E = X2E-FP*NEWAE-FPE*NEWA
          END IF     
          Y1E = GCE(1)+GCE(2)
          Y2E = GPE(1)+GPE(2)
          Z1E = HCE(1)+HCE(2)
          Z2E = HPE(1)+HPE(2)
          DENE1(K) =-DEN*DEN*(Y1E*Z2+Y1*Z2E-Y2E*Z1-Y2*Z1E)
          NEWC1(K) = (X1*Z2-X2*Z1)*DENE1(K)+(X1E*Z2+X1*Z2E-X2E*Z1-X2*Z1E)*DEN
          NEWP1(K) = (Y1*X2-Y2*X1)*DENE1(K)+(Y1E*X2+Y1*X2E-Y2E*X1-Y2*X1E)*DEN
          DAETA(KC,K) = DAETA(KC,K)+NEWC1(K)
          DAETA(KP,K) = DAETA(KP,K)+NEWP1(K)
        END IF   
!          
! Constant or background infusion
        IF (SSR /= ZERO) THEN
          IF (SSC /= KP) THEN
            IF (SSC == KA) DAETA(KA,K)=DAETA(KA,K)+SSR*DAE+SSRE(K)*DA
            SSRT1(K) = SSR*D10E+SSRE(K)*D10
            DAETA(KC,K) = DAETA(KC,K)+SSRT1(K)
            DAETA(KP,K) = DAETA(KP,K)+SSRT*ZE+SSRT1(K)*Z
          ELSE
            DAETA(KC,K) = DAETA(KC,K)+SSR*D10E+SSRE(K)*D10
            DAETA(KP,K) = DAETA(KP,K)+SSR*Y*D10E+SSR*YE*D10+SSRE(K)*Y*D10
          END IF  
        END IF
!        
        IF (.NOT.SECOND) CYCLE
!        
! Save first derivatives for ETA(K)
        A101(K)=A10E
        A121(K)=A12E
        A211(K)=A21E
        AA1(K)=AAE
        DA1(K)=DAE
        D101(K)=D10E
        GE1(K)=GE
        ZE1(K)=ZE
        YE1(K)=YE
        IF (DT == ZERO) CYCLE
        S1(K)=SE
        RT1(1,K)=RTE(1)
        RT1(2,K)=RTE(2)
        DO I=1,2
          DR1(I,K)=DRE(I)
          EX1(I,K)=EXE(I)
          A21R1(I,K)=A21RE(I)
          GR1(I,K)=GRE(I)
          IF (RHO /= ZERO) EXDEL1(I,K)=EXDELE(I)
          DDIF1(I,K)=DDIFE(I)
          IF (RHO /= ZERO) THEN
            DEDIF1(I,K)=DEDIFE(I)
            DADIF1(I,K)=DADIFE(I)
            DADS1(I,K)=DADSE(I)
          END IF
        END DO
        IF (SV(KA) /= 0) THEN
          EXA1(K)=EXAE
          IF (RHO /= ZERO) EXADL1(K)=EXADLE
          DDIFA1(1,K)=DDIFAE(1)
          FPR1(1,K)=FPRE(1)
          DDIFA1(2,K)=DDIFAE(2)
          FPR1(2,K)=FPRE(2)
          DPROD1(K)=DPRODE
          FA1(K)=FAE
          FPA1(K)=FPAE
          FC1(K)=FCE
          FP1(K)=FPE
          NEWA1(K)=NEWAE
        END IF
        X11(K)=X1E
        X21(K)=X2E
        IF (RHO /= ZERO) THEN
          IF (SSC == KA) T01(K)=T0E
          T11(K)=T1E
          T21(K)=T2E
        END IF
        Y11(K)=Y1E
        Y21(K)=Y2E
        Z11(K)=Z1E
        Z21(K)=Z2E
      END DO
!      
      IF (.NOT.SECOND) GO TO 999
!      
      DO K=1,NPETAS    ! Loop for second derivatives
        DO J=K,NPETAS
          A102 = GG(ITSC,J+1)*GG(IKE,K+1)+G2(ITSC,J,K)*GG(IKE,1)+GG(ITSC,1)*G2(IKE,J,K)      &
                +GG(ITSC,K+1)*GG(IKE,J+1)
          A122 = GG(ITSC,J+1)*GG(IK12,K+1)+G2(ITSC,J,K)*GG(IK12,1)+GG(ITSC,1)*G2(IK12,J,K)   &
                +GG(ITSC,K+1)*GG(IK12,J+1)
          A212 = GG(ITSC,J+1)*GG(IK21,K+1)+G2(ITSC,J,K)*GG(IK21,1)+GG(ITSC,1)*G2(IK21,J,K)   &
                +GG(ITSC,K+1)*GG(IK21,J+1)
          AA2  = GG(ITSC,J+1)*GG(IKA,K+1)+G2(ITSC,J,K)*GG(IKA,1)+GG(ITSC,1)*G2(IKA,J,K)      &
                +GG(ITSC,K+1)*GG(IKA,J+1)
          DA2  = (-DA1(J)*AA1(K)-DA*AA2-DA1(K)*AA1(J))/AA
          D102 = (-D101(J)*A101(K)-D10*A102-D101(K)*A101(J))/A10
          GE2  = A122+A102
          ZE2  = (A122-Z*A212-ZE1(K)*A211(J)-ZE1(J)*A211(K))/A21
          YE2  = (GE2-Y*A212-YE1(K)*A211(J)-YE1(J)*A211(K))/A21 
          IF (DT /= ZERO) THEN
            S2 = A102+A122+A212
            P2 = A211(J)*A101(K)+A212*A10+A21*A102+A211(K)*A101(J)
            D2 = TWO*S*S2+TWO*S1(J)*S1(K)-FOUR*P2
            TEMP2  = (TEMP*D2-D1(K)*TEMP1(J))/(TWO*D)
            RT2(1) = (S2+TEMP2)/TWO
            RT2(2) = (S2-TEMP2)/TWO
            DO I=1,2
              DR2(I) = (-DR(I)*RT2(I)-DR1(I,J)*RT1(I,K)-DR1(I,K)*RT1(I,J))/RT(I)
              EX2(I) =-EX1(I,J)*((DT+SDEL)*RT1(I,K)+(DTE(K)+SDELE(K))*RT(I))                 &
                      -EX(I)*((DTE(J)+SDELE(J))*RT1(I,K)+(D2DTE(J,K)+SDELE2(J,K))*RT(I))     &
                      -EX(I)*((DT+SDEL)*RT2(I)+(DTE(K)+SDELE(K))*RT1(I,J))
              A21R2(I) = A212-RT2(I)
              GR2(I) = GE2-RT2(I)
              IF (RHO == ZERO) CYCLE
              EXDEL2(I) =-EXDEL1(I,J)*(DTE(K)*RT(I)+DT*RT1(I,K))-EXDEL(I)*(DTE(K)*RT1(I,J)   &
                         +D2DTE(J,K)*RT(I)+DTE(J)*RT1(I,K)+DT*RT2(I))
            END DO
            L=3
            DO I=1,2
              L=L-1
              DDIF2(I) =-(RT2(L)-RT2(I))*DDIF(I)*DDIF(I)-(RT1(L,K)                           &
                        -RT1(I,K))*TWO*DDIF(I)*DDIF1(I,J)
              IF (RHO /= ZERO) THEN 
                DEDIF2(I) = D101(J)*DDIF1(I,K)+D102*DDIF(I)+D10*DDIF2(I)                     &
                           +D101(K)*DDIF1(I,J)  
                DADIF2(I) = DR1(I,J)*DDIF1(I,K)+DR2(I)*DDIF(I)+DR(I)*DDIF2(I)                &
                           +DR1(I,K)*DDIF1(I,J)
                DADS2(I)  = DADIF1(I,J)*DDIF1(I,K)+DADIF2(I)*DDIF(I)                         &
                           +DADIF(I)*DDIF2(I)+DADIF1(I,K)*DDIF1(I,J)
              END IF
!              
              GC2(I) = A21R1(I,J)*DDIF(I)*EX1(I,K)+A21R(I)*DDIF1(I,J)*EX1(I,K)               &
                      +A21R(I)*DDIF(I)*EX2(I)+A21R1(I,J)*DDIF1(I,K)*EX(I)                    &
                      +A21R(I)*DDIF2(I)*EX(I)+A21R(I)*DDIF1(I,K)*EX1(I,J)                    &
                      +A21R2(I)*DDIF(I)*EX(I)+A21R1(I,K)*DDIF1(I,J)*EX(I)                    &
                      +A21R1(I,K)*DDIF(I)*EX1(I,J)
              HC2(I) = A211(J)*DDIF(I)*EX1(I,K)+A211(J)*DDIF1(I,K)*EX(I)+A212*DDIF(I)*EX(I)  &
                      +A21*DDIF1(I,J)*EX1(I,K)+A21*DDIF2(I)*EX(I)+A211(K)*DDIF1(I,J)*EX(I)   &
                      +A21*DDIF(I)*EX2(I)+A21*DDIF1(I,K)*EX1(I,J)+A211(K)*DDIF(I)*EX1(I,J)
              GP2(I) = A121(J)*DDIF(I)*EX1(I,K)+A121(J)*DDIF1(I,K)*EX(I)+A122*DDIF(I)*EX(I)  &
                      +A12*DDIF1(I,J)*EX1(I,K)+A12*DDIF2(I)*EX(I)+A121(K)*DDIF1(I,J)*EX(I)   &
                      +A12*DDIF(I)*EX2(I)+A12*DDIF1(I,K)*EX1(I,J)+A121(K)*DDIF(I)*EX1(I,J)
              HP2(I) = GR1(I,J)*DDIF(I)*EX1(I,K)+GR1(I,J)*DDIF1(I,K)*EX(I)                   &
                      +GR2(I)*DDIF(I)*EX(I)+GR(I)*DDIF1(I,J)*EX1(I,K)+GR(I)*DDIF2(I)*EX(I)   &
                      +GR1(I,K)*DDIF1(I,J)*EX(I)+GR(I)*DDIF(I)*EX2(I)                        &
                      +GR(I)*DDIF1(I,K)*EX1(I,J)+GR1(I,K)*DDIF(I)*EX1(I,J)
            END DO
!            
            IF (SV(KA) /= 0) THEN
              EXA2 =-EXA1(J)*((DT+SDEL)*AA1(K)+(DTE(K)+SDELE(K))*AA)                         &
                    -EXA*((DTE(J)+SDELE(J))*AA1(K)+(D2DTE(J,K)+SDELE2(J,K))*AA               &
                    +(DT+SDEL)*AA2+(DTE(K)+SDELE(K))*AA1(J))   
!                     
              IF (RHO /= ZERO) EXADL2 =-EXADL1(J)*(DT*AA1(K)+DTE(K)*AA)-EXADL*(DTE(J)*AA1(K) &
                                        +DT*AA2+D2DTE(J,K)*AA+DTE(K)*AA1(J))
              DO I=1,2
                DDIFA2(I) = -TWO*DDIFA(I)*DDIFA1(I,J)*(AA1(K)-RT1(I,K))                      &
                            -DDIFA(I)*DDIFA(I)*(AA2-RT2(I))
                FPR2(I) = DDIFA1(I,J)*DDIF(I)*EX1(I,K)+DDIFA1(I,J)*DDIF1(I,K)*EX(I)          &
                         +DDIFA2(I)*DDIF(I)*EX(I)+DDIFA(I)*DDIF1(I,J)*EX1(I,K)               &
                         +DDIFA(I)*DDIF2(I)*EX(I)+DDIFA1(I,K)*DDIF1(I,J)*EX(I)               &
                         +DDIFA(I)*DDIF(I)*EX2(I)+DDIFA(I)*DDIF1(I,K)*EX1(I,J)               &
                         +DDIFA1(I,K)*DDIF(I)*EX1(I,J)
                FCR2(I) = FPR1(I,J)*(AA*A21R1(I,K)+AA1(K)*A21R(I))                           &
                         +FPR(I)*(AA1(J)*A21R1(I,K)+AA2*A21R(I))                             &
                         +FPR(I)*(AA*A21R2(I)+AA1(K)*A21R1(I,J))                             &
                         +FPR2(I)*AA*A21R(I)+FPR1(I,K)*AA1(J)*A21R(I)                        &    
                         +FPR1(I,K)*AA*A21R1(I,J)
              END DO
              DPROD2 = DDIFA1(1,J)*DDIFA1(2,K)+DDIFA2(1)*DDIFA(2)+DDIFA(1)*DDIFA2(2)         &
                      +DDIFA1(1,K)*DDIFA1(2,J)
              FA2  = EXA2
              FCA2 = AA1(J)*(A21-AA)*DPROD*EXA1(K)+AA1(J)*(A21-AA)*DPROD1(K)*EXA             &
                    +AA1(J)*(A211(K)-AA1(K))*DPROD*EXA+AA2*(A21-AA)*DPROD*EXA                &
                    +AA*(A211(J)-AA1(J))*DPROD*EXA1(K)+AA*(A211(J)-AA1(J))*DPROD1(K)*EXA     &
                    +AA*(A212-AA2)*DPROD*EXA+AA1(K)*(A211(J)-AA1(J))*DPROD*EXA               &
                    +AA*(A21-AA)*DPROD1(J)*EXA1(K)+AA*(A21-AA)*DPROD2*EXA                    &
                    +AA*(A211(K)-AA1(K))*DPROD1(J)*EXA+AA1(K)*(A21-AA)*DPROD1(J)*EXA         &
                    +AA*(A21-AA)*DPROD*EXA2+AA*(A21-AA)*DPROD1(K)*EXA1(J)                    &
                    +AA*(A211(K)-AA1(K))*DPROD*EXA1(J)+AA1(K)*(A21-AA)*DPROD*EXA1(J)
              FPA2 = DPROD1(J)*EXA1(K)+DPROD2*EXA+DPROD*EXA2+DPROD1(K)*EXA1(J)
              FC2  = FCA2+FCR2(1)+FCR2(2)
              FP2  = (FPA1(J)+FPR1(1,J)+FPR1(2,J))*(AA*A121(K)+AA1(K)*A12)                   &
                    +(FPA+FPR(1)+FPR(2))*(AA1(J)*A121(K)+AA*A122+AA2*A12+AA1(K)*A121(J))     &
                    +(FPA2+FPR2(1)+FPR2(2))*AA*A12+(FPA1(K)+FPR1(1,K)                        &
                    +FPR1(2,K))*(AA1(J)*A12+AA*A121(J))
            END IF
!            
            X02=ZERO; X12=ZERO; X22=ZERO
!            
            IF (RHO == ZERO) THEN  ! Multiple bolus doses
              IF (SSC == KA) X02=-SSAE2(J,K)
              IF (SSC == KC) X12=-SSAE2(J,K)
              IF (SSC == KP) X22=-SSAE2(J,K)
            ELSE                   ! Multiple infusions
              MDE2(1) = DEDIF1(1,J)*EXDEL1(1,K)+DEDIF2(1)*EXDEL(1)                      &
                       +DEDIF(1)*EXDEL2(1)+DEDIF1(1,K)*EXDEL1(1,J)      
              MDE2(2) = DEDIF1(2,J)*EXDEL1(2,K)+DEDIF2(2)*EXDEL(2)                      &
                       +DEDIF(2)*EXDEL2(2)+DEDIF1(2,K)*EXDEL1(2,J)
              MDA2(1) = DADS1(1,J)*EX1(1,K)+DADS2(1)*EX(1)                              &
                       +DADS(1)*EX2(1)+DADS1(1,K)*EX1(1,J)   
              MDA2(2) = DADS1(2,J)*EX1(2,K)+DADS2(2)*EX(2)                              &
                       +DADS(2)*EX2(2)+DADS1(2,K)*EX1(2,J)
              IF (SSC /= KA) THEN 
                IF (SSC == KP) THEN 
                  T12 = (GE1(J)+A21R1(1,J))*MDE1(1,K)+(G+A21R(1))*MDE2(1)               &
                       +(GE2+A21R2(1))*MDE(1)+(GE1(K)+A21R1(1,K))*MDE1(1,J)             &
                       +(GE1(J)+A21R1(2,J))*MDE1(2,K)+(G+A21R(2))*MDE2(2)               &
                       +(GE2+A21R2(2))*MDE(2)+(GE1(K)+A21R1(2,K))*MDE1(2,J)             &
                       -(A21R1(1,J)*A21+GR1(1,J)*A21)*MDA1(1,K)                         &  
                       -(A21R(1)*A211(J)+GR(1)*A211(J))*MDA1(1,K)                       &
                       -(A21R(1)*A21+GR(1)*A21)*MDA2(1)                                 &
                       -((A21R2(1)*A21+GR2(1)*A21+A21R1(1,J)*A211(K)+GR1(1,J)*A211(K)   &  
                       +A21R1(1,K)*A211(J)+GR1(1,K)*A211(J)+A21R(1)*A212+GR(1)*A212))   &
                       *MDA(1)-(A21R1(1,K)*A21+GR1(1,K)*A21+A21R(1)*A211(K)             &
                       +GR(1)*A211(K))*MDA1(1,J)                                        &
                       -(A21R1(2,J)*A21+GR1(2,J)*A21)*MDA1(2,K)                         &
                       -(A21R(2)*A211(J)+GR(2)*A211(J))*MDA1(2,K)                       &
                       -(A21R(2)*A21+GR(2)*A21)*MDA2(2)                                 &  
                       -((A21R2(2)*A21+GR2(2)*A21+A21R1(2,J)*A211(K)+GR1(2,J)*A211(K)   &
                       +A21R1(2,K)*A211(J)+GR1(2,K)*A211(J)+A21R(2)*A212+GR(2)*A212))   &
                       *MDA(2)-(A21R1(2,K)*A21+GR1(2,K)*A21+A21R(2)*A211(K)             &
                       +GR(2)*A211(K))*MDA1(2,J)
                  X12 =-RHOE(J)*T11(K)-RHOE2(J,K)*T1-RHO*T12-RHOE(K)*T11(J)
                  T22 = (A121(J)+YE1(J)*GR(1)+Y*GR1(1,J))*MDE1(1,K)                     &   
                       +(A12+Y*GR(1))*MDE2(1)                                           &
                       +(A122+YE1(J)*GR1(1,K)+YE2*GR(1)+Y*GR2(1)+YE1(K)*GR1(1,J))*MDE(1)&
                       +(A121(K)+Y*GR1(1,K)+YE1(K)*GR(1))*MDE1(1,J)                     &    
                       +(A121(J)+YE1(J)*GR(2)+Y*GR1(2,J))*MDE1(2,K)                     &
                       +(A12+Y*GR(2))*MDE2(2)                                           &
                       +(A122+YE1(J)*GR1(2,K)+YE2*GR(2)                                 &
                       +Y*GR2(2)+YE1(K)*GR1(2,J))*MDE(2)                                &
                       +(A121(K)+Y*GR1(2,K)+YE1(K)*GR(2))*MDE1(2,J)                     &    
                       -(GR1(1,J)*GR(1)+A121(J)*A21+GR(1)*GR1(1,J)+A12*A211(J))*MDA1(1,K)&
                       -(GR(1)*GR(1)+A12*A21)*MDA2(1)                                   &
                       -(TWO*GR1(1,J)*GR1(1,K)+A121(J)*A211(K)+A122*A21                 &
                       +TWO*GR(1)*GR2(1)+A12*A212+A121(K)*A211(J))*MDA(1)               &
                       -(TWO*GR(1)*GR1(1,K)+A12*A211(K)+A121(K)*A21)*MDA1(1,J)          &
                       -(TWO*GR1(2,J)*GR(2)+A121(J)*A21+A12*A211(J))*MDA1(2,K)          &
                       -(GR(2)*GR(2)+A12*A21)*MDA2(2)                                   &
                       -(TWO*GR1(2,J)*GR1(2,K)+A121(J)*A211(K)+A122*A21                 &
                       +TWO*GR(2)*GR2(2)+A12*A212+A121(K)*A211(J))*MDA(2)               &
                       -(TWO*GR(2)*GR1(2,K)+A12*A211(K)+A121(K)*A21)*MDA1(2,J)
                  X22 =-RHOE(J)*T21(K)-RHOE2(J,K)*T2-RHO*T22-RHOE(K)*T21(J)
                ELSE 
                  T12 = (A121(J)+A21R1(1,J))*MDE1(1,K)+(A12+A21R(1))*MDE2(1)            &
                       +(A122+A21R2(1))*MDE(1)+(A121(K)+A21R1(1,K))*MDE1(1,J)           &
                       +(A121(J)+A21R1(2,J))*MDE1(2,K)                                  &
                       +(A12+A21R(2))*MDE2(2)                                           &
                       +(A122+A21R2(2))*MDE(2)                                          &
                       +(A121(K)+A21R1(2,K))*MDE1(2,J)                                  &
                       -(TWO*A21R(1)*A21R1(1,J)+A121(J)*A21+A12*A211(J))*MDA1(1,K)      &
                       -(A21R(1)*A21R(1)+A12*A21)*MDA2(1)                               &
                       -(TWO*A21R1(1,J)*A21R1(1,K)+A121(J)*A211(K)+A122*A21             &
                       +TWO*A21R(1)*A21R2(1)+A12*A212+A121(K)*A211(J))*MDA(1)           &
                       -(TWO*A21R(1)*A21R1(1,K)+A12*A211(K)+A121(K)*A21)*MDA1(1,J)      &
                       -(TWO*A21R(2)*A21R1(2,J)+A12*A211(J)+A121(J)*A21)*MDA1(2,K)      &
                       -(A21R(2)*A21R(2)+A12*A21)*MDA2(2)                               &
                       -(TWO*A21R1(2,J)*A21R1(2,K)+A121(J)*A211(K)+A122*A21             &
                       +TWO*A21R(2)*A21R2(2)+A12*A212+A121(K)*A211(J))*MDA(2)           &
                       -(TWO*A21R(2)*A21R1(2,K)+A12*A211(K)+A121(K)*A21)*MDA1(2,J)
                  X12 =-RHOE(J)*T11(K)-RHOE2(J,K)*T1-RHO*T12-RHOE(K)*T11(J)
                  T22 = (A121(J)+ZE1(J)*GR(1)+Z*GR1(1,J))*MDE1(1,K)                     &
                       +(A12+Z*GR(1))*MDE2(1)                                           &
                       +(A122+Z*GR2(1)+ZE1(J)*GR1(1,K)+ZE2*GR(1)+ZE1(K)*GR1(1,J))*MDE(1)&
                       +(A121(K)+Z*GR1(1,K)+ZE1(K)*GR(1))*MDE1(1,J)                     &    
                       +(A121(J)+ZE1(J)*GR(2)+Z*GR1(2,J))*MDE1(2,K)                     &
                       +(A12+Z*GR(2))*MDE2(2)                                           &
                       +(A122+Z*GR2(2)+ZE1(J)*GR1(2,K)+ZE2*GR(2)+ZE1(K)*GR1(2,J))*MDE(2)&
                       +(A121(K)+Z*GR1(2,K)+ZE1(K)*GR(2))*MDE1(2,J)                     &    
                       -(GR1(1,J)+A21R1(1,J))*A12*MDA1(1,K)                             &
                       -(GR(1)+A21R(1))*A121(J)*MDA1(1,K)-(GR(1)+A21R(1))*A12*MDA2(1)   &
                       -((GR1(1,J)+A21R1(1,J))*A121(K)+(GR2(1)+A21R2(1))*A12)*MDA(1)    &
                       -((GR(1)+A21R(1))*A122+(GR1(1,K)+A21R1(1,K))*A121(J))*MDA(1)     &
                       -((GR(1)+A21R(1))*A121(K)+(GR1(1,K)+A21R1(1,K))*A12)*MDA1(1,J)   &
                       -(GR1(2,J)+A21R1(2,J))*A12*MDA1(2,K)                             &
                       -(GR(2)+A21R(2))*A121(J)*MDA1(2,K)-(GR(2)+A21R(2))*A12*MDA2(2)   &
                       -((GR1(2,J)+A21R1(2,J))*A121(K)+(GR2(2)+A21R2(2))*A12)*MDA(2)    &
                       -((GR(2)+A21R(2))*A122+(GR1(2,K)+A21R1(2,K))*A121(J))*MDA(2)     &
                       -((GR(2)+A21R(2))*A121(K)+(GR1(2,K)+A21R1(2,K))*A12)*MDA1(2,J)   
                  X22 =-RHOE(J)*T21(K)-RHOE2(J,K)*T2-RHO*T22-RHOE(K)*T21(J)
                END IF
              ELSE
                T02 = DA1(J)*(EXADL1(K)-EXA1(K))+DA2*(EXADL-EXA)                        &
                     +DA*(EXADL2-EXA2)+DA1(K)*(EXADL1(J)-EXA1(J))
                X02 =-RHOE(J)*T01(K)-RHOE2(J,K)*T0-RHO*T02-RHOE(K)*T01(J)
                MDI2(1) = DDIF1(1,J)*EXDEL1(1,K)+DDIF2(1)*EXDEL(1)                      &
                         +DDIF(1)*EXDEL2(1)+DDIF1(1,K)*EXDEL1(1,J)
                MDI2(2) = DDIF1(2,J)*EXDEL1(2,K)+DDIF2(2)*EXDEL(2)                      &
                         +DDIF(2)*EXDEL2(2)+DDIF1(2,K)*EXDEL1(2,J)
                DPRD2   = DPROD2*(EXADL-EXA)+DPROD1(J)*(EXADL1(K)-EXA1(K))              &
                         +DPROD1(K)*(EXADL1(J)-EXA1(J))+DPROD*(EXADL2-EXA2)
                TMA2(1) = (A21R2(1)+A122)*D10+A21R2(1)*DDIFA(1)                         &
                         +(A21R1(1,K)+A121(K))*D101(J)+A21R1(1,K)*DDIFA1(1,J)           &
                         +(A21R1(1,J)+A121(J))*D101(K)+A21R1(1,J)*DDIFA1(1,K)           &
                         +(A21R(1)+A12)*D102+A21R(1)*DDIFA2(1)
                TMA2(2) = (A21R2(2)+A122)*D10+A21R2(2)*DDIFA(2)                         &
                         +(A21R1(2,K)+A121(K))*D101(J)+A21R1(2,K)*DDIFA1(2,J)           &
                         +(A21R1(2,J)+A121(J))*D101(K)+A21R1(2,J)*DDIFA1(2,K)           &
                         +(A21R(2)+A12)*D102+A21R(2)*DDIFA2(2)
                TMB2(1) = AA2*A21R(1)*A21R(1)+AA1(K)*A21R1(1,J)*A21R(1)                 &
                         +AA1(K)*A21R(1)*A21R1(1,J)+A122*AA*A21+A121(K)*AA1(J)*A21      &
                         +A121(K)*AA*A211(J)+AA1(J)*A21R1(1,K)*A21R(1)                  &
                         +AA*A21R2(1)*A21R(1)+AA*A21R1(1,K)*A21R1(1,J)                  &
                         +A121(J)*AA1(K)*A21+A12*AA2*A21+A12*AA1(K)*A211(J)             &
                         +AA1(J)*A21R(1)*A21R1(1,K)+AA*A21R1(1,J)*A21R1(1,K)            &
                         +AA*A21R(1)*A21R2(1)+A121(J)*AA*A211(K)                        &
                         +A12*AA1(J)*A211(K)+A12*AA*A212
                TMB2(2) = AA2*A21R(2)*A21R(2)+AA1(K)*A21R1(2,J)*A21R(2)                 &
                         +AA1(K)*A21R(2)*A21R1(2,J)+A122*AA*A21+A121(K)*AA1(J)*A21      &
                         +A121(K)*AA*A211(J)+AA1(J)*A21R1(2,K)*A21R(2)                  &
                         +AA*A21R2(2)*A21R(2)+AA*A21R1(2,K)*A21R1(2,J)                  &
                         +A121(J)*AA1(K)*A21+A12*AA2*A21+A12*AA1(K)*A211(J)             &
                         +AA1(J)*A21R(2)*A21R1(2,K)+AA*A21R1(2,J)*A21R1(2,K)            &
                         +AA*A21R(2)*A21R2(2)+A121(J)*AA*A211(K)                        &
                         +A12*AA1(J)*A211(K)+A12*AA*A212
                T12 = (A212-AA2)*DPRD+(A211(K)-AA1(K))*DPRD1(J)                         &
                     +(A211(J)-AA1(J))*DPRD1(K)+(A21-AA)*DPRD2                          &
                     +TMA2(1)*MDI(1)+TMA1(1,J)*MDI1(1,K)                                &
                     +TMA1(1,K)*MDI1(1,J)+TMA(1)*MDI2(1)                                &
                     +TMA2(2)*MDI(2)+TMA1(2,J)*MDI1(2,K)                                &
                     +TMA1(2,K)*MDI1(2,J)+TMA(2)*MDI2(2)                                &
                     -TMB2(1)*DDIFA(1)*MDA(1)-TMB1(1,K)*DDIFA1(1,J)*MDA(1)              &
                     -TMB1(1,K)*DDIFA(1)*MDA1(1,J)-TMB1(1,J)*DDIFA1(1,K)*MDA(1)         &
                     -TMB(1)*DDIFA2(1)*MDA(1)-TMB(1)*DDIFA1(1,K)*MDA1(1,J)              &
                     -TMB1(1,J)*DDIFA(1)*MDA1(1,K)-TMB(1)*DDIFA1(1,J)*MDA1(1,K)         &    
                     -TMB(1)*DDIFA(1)*MDA2(1)-TMB2(2)*DDIFA(2)*MDA(2)                   &
                     -TMB1(2,K)*DDIFA1(2,J)*MDA(2)-TMB1(2,K)*DDIFA(2)*MDA1(2,J)         &
                     -TMB1(2,J)*DDIFA1(2,K)*MDA(2)-TMB(2)*DDIFA2(2)*MDA(2)              &
                     -TMB(2)*DDIFA1(2,K)*MDA1(2,J)-TMB1(2,J)*DDIFA(2)*MDA1(2,K)         &
                     -TMB(2)*DDIFA1(2,J)*MDA1(2,K)-TMB(2)*DDIFA(2)*MDA2(2)
                X12 =-RHOE(J)*T11(K)-RHOE2(J,K)*T1-RHO*T12-RHOE(K)*T11(J)
                TMC2(1) = (A121(K)+ZE1(K)*GR(1)+Z*GR1(1,K))*D101(J)                     &
                         +(A122+ZE2*GR(1)+ZE1(J)*GR1(1,K)+ZE1(K)*GR1(1,J)+Z*GR2(1))*D10 &
                         +(A121(J)+ZE1(J)*GR(1)+Z*GR1(1,J))*D101(K)                     &
                         +(A12+Z*GR(1))*D102                                            &
                         +A122*DDIFA(1)+A121(K)*DDIFA1(1,J)                             &
                         +A121(J)*DDIFA1(1,K)+A12*DDIFA2(1)
                TMC2(2) = (A121(K)+ZE1(K)*GR(2)+Z*GR1(2,K))*D101(J)                     &
                         +(A122+ZE2*GR(2)+ZE1(J)*GR1(2,K)+ZE1(K)*GR1(2,J)+Z*GR2(2))*D10 &
                         +(A121(J)+ZE1(J)*GR(2)+Z*GR1(2,J))*D101(K)                     &    
                         +(A12+Z*GR(2))*D102                                            &
                         +A122*DDIFA(2)+A121(K)*DDIFA1(2,J)                             &
                         +A121(J)*DDIFA1(2,K)+A12*DDIFA2(2)
                TMD2(1) = AA2*A12*(A21R(1)+GR(1))                                       &
                         +AA1(K)*A121(J)*(A21R(1)+GR(1))                                &
                         +AA1(K)*A12*(A21R1(1,J)+GR1(1,J))                              &
                         +AA1(J)*A121(K)*(A21R(1)+GR(1))                                &
                         +AA*A122*(A21R(1)+GR(1))                                       &
                         +AA*A121(K)*(A21R1(1,J)+GR1(1,J))                              &
                         +AA1(J)*A12*(A21R1(1,K)+GR1(1,K))                              &
                         +AA*A121(J)*(A21R1(1,K)+GR1(1,K))                              &
                         +AA*A12*(A21R2(1)+GR2(1))
                TMD2(2) = AA2*A12*(A21R(2)+GR(2))                                       &
                         +AA1(K)*A121(J)*(A21R(2)+GR(2))                                &
                         +AA1(K)*A12*(A21R1(2,J)+GR1(2,J))                              &
                         +AA1(J)*A121(K)*(A21R(2)+GR(2))                                &
                         +AA*A122*(A21R(2)+GR(2))                                       &
                         +AA*A121(K)*(A21R1(2,J)+GR1(2,J))                              &
                         +AA1(J)*A12*(A21R1(2,K)+GR1(2,K))                              &
                         +AA*A121(J)*(A21R1(2,K)+GR1(2,K))                              &
                         +AA*A12*(A21R2(2)+GR2(2))
                T22 = A122*DPRD+A121(J)*DPRD1(K)+A121(K)*DPRD1(J)+A12*DPRD2             &
                     +TMC2(1)*MDI(1)+TMC1(1,J)*MDI1(1,K)                                &
                     +TMC1(1,K)*MDI1(1,J)+TMC(1)*MDI2(1)                                &
                     +TMC2(2)*MDI(2)+TMC1(2,J)*MDI1(2,K)                                &
                     +TMC1(2,K)*MDI1(2,J)+TMC(2)*MDI2(2)                                &
                     -TMD2(1)*DDIFA(1)*MDA(1)-TMD1(1,K)*DDIFA1(1,J)*MDA(1)              &
                     -TMD1(1,K)*DDIFA(1)*MDA1(1,J)                                      &
                     -TMD1(1,J)*DDIFA1(1,K)*MDA(1)-TMD(1)*DDIFA2(1)*MDA(1)              &
                     -TMD(1)*DDIFA1(1,K)*MDA1(1,J)                                      &
                     -TMD1(1,J)*DDIFA(1)*MDA1(1,K)-TMD(1)*DDIFA1(1,J)*MDA1(1,K)         &
                     -TMD(1)*DDIFA(1)*MDA2(1)                                           &        
                     -TMD2(2)*DDIFA(2)*MDA(2)-TMD1(2,K)*DDIFA1(2,J)*MDA(2)              &
                     -TMD1(2,K)*DDIFA(2)*MDA1(2,J)                                      &
                     -TMD1(2,J)*DDIFA1(2,K)*MDA(2)-TMD(2)*DDIFA2(2)*MDA(2)              &
                     -TMD(2)*DDIFA1(2,K)*MDA1(2,J)                                      &
                     -TMD1(2,J)*DDIFA(2)*MDA1(2,K)-TMD(2)*DDIFA1(2,J)*MDA1(2,K)         &
                     -TMD(2)*DDIFA(2)*MDA2(2)   
                X22 =-RHOE(J)*T21(K)-RHOE2(J,K)*T2-RHO*T22-RHOE(K)*T21(J)
              END IF  
            END IF
!                
! Simultaneous equations
            IF (SV(KA) /= 0) THEN
              NEWA2 = (X02-NEWA*FA2-NEWA1(J)*FA1(K)-NEWA1(K)*FA1(J))/(FA-ONE)
              D2AETA(KA,J,K) = D2AETA(KA,J,K)+NEWA2
              X12 = X12-FC1(J)*NEWA1(K)-FC2*NEWA-FC*NEWA2-FC1(K)*NEWA1(J)
              X22 = X22-FP1(J)*NEWA1(K)-FP2*NEWA-FP*NEWA2-FP1(K)*NEWA1(J)
            END IF
            Y12 = GC2(1)+GC2(2)
            Y22 = GP2(1)+GP2(2)
            Z12 = HC2(1)+HC2(2)
            Z22 = HP2(1)+HP2(2)
            DENE2 =-TWO*DEN*DENE1(J)*(Y1*Z21(K)+Y11(K)*Z2-Y2*Z11(K)-Y21(K)*Z1)          &
                   -DEN*DEN*(Y11(J)*Z21(K)+Y12*Z2-Y21(J)*Z11(K)-Y22*Z1                  &
                   +Y1*Z22+Y11(K)*Z21(J)-Y2*Z12-Y21(K)*Z11(J))
            NEWC2 = (X1*Z2-X2*Z1)*DENE2                                                 &
                   +(X11(J)*Z2-X21(J)*Z1+X1*Z21(J)-X2*Z11(J))*DENE1(K)                  &
                   +(X11(J)*Z21(K)+X12*Z2-X21(J)*Z11(K)-X22*Z1                          &
                   +X1*Z22+X11(K)*Z21(J)-X2*Z12-X21(K)*Z11(J))*DEN                      &
                   +(X1*Z21(K)+X11(K)*Z2-X2*Z11(K)-X21(K)*Z1)*DENE1(J)
            NEWP2 = (Y1*X2-Y2*X1)*DENE2                                                 &
                   +(Y11(J)*X2-Y21(J)*X1+Y1*X21(J)-Y2*X11(J))*DENE1(K)                  &
                   +(Y11(J)*X21(K)+Y12*X2-Y21(J)*X11(K)-Y22*X1                          &
                   +Y1*X22+Y11(K)*X21(J)-Y2*X12-Y21(K)*X11(J))*DEN                      &
                   +(Y1*X21(K)+Y11(K)*X2-Y2*X11(K)-Y21(K)*X1)*DENE1(J)
            D2AETA(KC,J,K) = D2AETA(KC,J,K)+NEWC2
            D2AETA(KP,J,K) = D2AETA(KP,J,K)+NEWP2
          END IF
!          
! Constant or background infusion
          IF (SSR == ZERO) CYCLE
          IF (SSC /= KP) THEN
            IF (SSC == KA) D2AETA(KA,J,K) = D2AETA(KA,J,K)+SSRE(J)*DA1(K)+SSRE2(J,K)*DA &
                                           +SSR*DA2+SSRE(K)*DA1(J)
            SSRT2 = SSRE(J)*D101(K)+SSRE2(J,K)*D10+SSR*D102+SSRE(K)*D101(J)
            D2AETA(KC,J,K) = D2AETA(KC,J,K)+SSRT2
            D2AETA(KP,J,K) = D2AETA(KP,J,K)+SSRT1(J)*ZE1(K)+SSRT*ZE2+SSRT2*Z            &
                            +SSRT1(K)*ZE1(J)
            CYCLE
          END IF
            D2AETA(KC,J,K) = D2AETA(KC,J,K)+SSRE(J)*D101(K)+SSRE2(J,K)*D10              &
                            +SSR*D102+SSRE(K)*D101(J)
            D2AETA(KP,J,K) = D2AETA(KP,J,K)                                             &
                             +SSRE(J)*Y*D101(K)+SSRE(J)*YE1(K)*D10+SSRE2(J,K)*Y*D10     &
                             +SSR*YE1(J)*D101(K)+SSR*YE2*D10+SSRE(K)*YE1(J)*D10         &
                             +SSR*Y*D102+SSR*YE1(K)*D101(J)+SSRE(K)*Y*D101(J)
        END DO
      END DO 
!      
  999 RETURN 
! 
      END SUBROUTINE SS
