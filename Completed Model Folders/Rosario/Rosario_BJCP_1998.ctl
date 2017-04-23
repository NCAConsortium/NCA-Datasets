;Model Desc: 2 CMT IVB LINEAR
;Project Name: NCAC SIMS
;Author: BM
;QC: JK
;QC2: BM

$PROB RUN# 1

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT NTIM NTAD TIME TAD DOSE 
       AMT SS II DV MDV EVID 
       CMT CR60 BSA ALB

$DATA NMDATA IGNORE=@ ; NMDATA = RosarioRun1.csv, RosarioRun2.csv

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Intravenous Dosage Amount Data Item (mg)
; DV:   Dependent Variable (mg/L);
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing);
; CR60: Creatinine Clearance (with 60 umol/L cutoff)
; BSA:  Body Surface Area (m**2)
; ALB:  Albumin (g/L)

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************
$SUBROUTINE ADVAN3 TRANS4

$PK

    TVCL   = THETA(1)                                             ; POPULATION CLEARACE
    GPCL   = TVCL * (1 + THETA(2) * CR60)                         ; GROUP CLEARANCE
    PPVCL  = EXP(ETA(1))                                          ; BSV CLEARANCE
    CL     = GPCL * PPVCL                                         ; INDIVIDUAL CLEARANCE

    TVV1   = THETA(3)                                             ; POPULATION CENTRAL VOLUME
    GPV1   = TVV1 * (BSA * (ALB/34) ** THETA(4))                  ; GROUP CENTRAL VOLUME 
    V1     = GPV1                                                 ; INDIVIDUAL CENTRAL VOLUME

    TVQ    = THETA(5)                                             ; POPULATION DISTRIBUTIONAL CLEARANCE
    PPVQ   = EXP(ETA(2))                                          ; BSV DISTRIBUTIONAL CLEARANCE
    Q      = TVQ * PPVQ                                           ; INDIVIDUAL DISTRIBUTIONAL CLEARANCE

    TVV2   = THETA(6)                                             ; POPULATION INTERCOMPARTMENTAL CLEARANCE 2 
    V2     = TVV2                                                 ; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE 2
      
;------------------------------------------------------------
;----           TIME AFTER DOSE CALCULATION              ----
;------------------------------------------------------------



;------------------------------------------------------------
;----           SECONDARY MODEL PARAMETERS               ----
;------------------------------------------------------------
    K10   = CL/V1
    K12   = Q/V1
    K21   = Q/V2
    
    A0 = K10 * K21
    A1 = K10 + K12 + K21

    RT1 = (A1 + SQRT(A1*A1 - 4* A0)) / 2
    RT2 = (A1 - SQRT(A1*A1 - 4* A0)) / 2

    IF(RT1.LT.RT2) THEN
      BETA  = RT1
      ALPHA = RT2
    ENDIF
    IF(RT2.LT.RT1) THEN
      BETA  = RT2
      ALPHA = RT1
    ENDIF

    THALF  = LOG(2)/BETA                                ; HALF-LIFE(h)
    MRT  = 1/BETA                                       ; MEAN RESIDENT TIME (h) 
    AREA = DOSE/CL                                      ; MEAN RESIDENT TIME (mg*hr/L) 
    AUMC = AREA*MRT                                     ; AREA UNDER THE CURVE (mg*hr**2/L)
    VSS  = V1 + V2                                      ; VOLUME STEADY-STATE (L)

;------------------------------------------------------------
;----               RESIDUAL ERROR MODEL                 ----
;------------------------------------------------------------
$ERROR

IPRE = F
Y    = F * EXP(EPS(1)) + EPS(2)
IRES = DV-IPRE
IWRE = IRES/IPRE

;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************
$THETA
0.88                                 ; TV CLEARANCE (L/H) CL
0.043                                ; COV CR60 ~ CL
8.59                                 ; TV V1 (L)
-0.39                                ; COV ALB ~ V1
1.30                                 ; TV Q (L/H)
9.79                                 ; TV V2 (L)
  
$OMEGA
0.0336                               ; [E] BETWEEN SUBJECT VARIABILITY (CL)
0.0765                               ; [E] BETWEEN SUBJECT VARIABILITY (Q)

$SIGMA
0.025                                ; [P] PROPORTIONAL RESIDUAL ERROR
0.107                                ; [P] ADDITIVE RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************
$SIM ONLYSIM (8675309)

$TABLE                                         ; STANDARD TABLE
       ID NTIM NTAD TIME TAD DOSE AMT EVID MDV CMT 
       CR60 BSA ALB
       CL Q V1 V2 ETA1 ETA2
       THALF MRT AREA AUMC VSS
       IRES IPRE                               ;DV PRED RES WRES - APPENDED
       NOPRINT ONEHEADER FILE=sdtab

