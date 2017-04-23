;Model Desc: 2 CMT INFUSION LINEAR
;Project Name: NCAC SIMS
;Author: BM
;QC: JK
;QC2: BM

$PROB RUN# 1

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT NTIM NTAD TIME TAD DOSE AMT DUR RATE=DROP SS II DV MDV 
EVID CMT WTKG CRCL

$DATA NMDATA IGNORE=@ ; NMDATA = JonssonRun1.csv, JonssonRun2.csv

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Intravenous Dosage Amount Data Item (mg)
; DUR:  Infusion Duration (hr)
; DV:   Dependent Variable (mg/L);
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing);
; CRCL: Creatinine Clearance (mL/min)
; WTKG: Body Weight (kg)

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************
$SUBROUTINE ADVAN3 TRANS4

$PK

    TVCL   = THETA(1)                                             ; POPULATION CLEARACE
    IF (CRCL.GT.40) THEN
      GPCL = TVCL * (1 + THETA(2) * (CRCL-70))                    ; GROUP CLEARANCE CRCL>40
    END IF
    IF (CRCL.LE.40) THEN
      GPCL = THETA(3)                                             ; GROUP CLEARANCE CRCL</=40
    ENDIF
    PPVCL  = EXP(ETA(1))                                          ; BSV CLEARANCE
    CL     = GPCL * PPVCL                                         ; INDIVIDUAL CLEARANCE

    TVV1   = THETA(4)                                             ; POPULATION CENTRAL VOLUME
    GPV1   = TVV1 * (1 + THETA(5) * (WTKG-75))                    ; GROUP CENTRAL VOLUME
    PPVV1  = EXP(ETA(2))                                          ; BSV VOLUME
    V1     = GPV1 * PPVV1                                         ; INDIVIDUAL CENTRAL VOLUME

    TVQ    = THETA(6)                                             ; POPULATION DISTRIBUTIONAL CLEARANCE
    Q      = TVQ                                                  ; INDIVIDUAL DISTRIBUTIONAL CLEARANCE

    TVV2   = THETA(7)                                             ; POPULATION INTERCOMPARTMENTAL CLEARANCE 2 
    V2     = TVV2                                                 ; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE 2

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

IPRE=F
Y    = F * (1 + EPS(1))
IRES = DV-IPRE
IWRE = IRES/IPRE

;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************
$THETA
4.54                                 ; TV CLEARANCE (L/H) CL
0.012                                ; COV CRCL ~ CL
2.91                                 ; GPCL FOR CRCL<40
7.76                                 ; TV VOLUME (L) V1
0.0198                               ; COV WTKG ~ V1
13.1                                 ; TV Q (L/H)
7.17                                 ; TV V2 (L)
  
$OMEGA BLOCK (2)
0.0515                               ; [E] BETWEEN SUBJECT VARIABILITY (CL)
0.0236                               ; [F] CORRELATION CL~V1
0.1484                               ; [E] BETWEEN SUBJECT VARIABILITY (V1)

$SIGMA
0.02722                              ; [P] PROPORTIONAL RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************
$SIM ONLYSIM (8675309)

$TABLE                                       ; STANDARD TABLE
       ID NTIM NTAD TIME TAD DOSE AMT DUR EVID MDV CMT 
       CRCL WTKG
       CL Q V1 V2 ETA1 ETA2
       THALF MRT AREA AUMC VSS
       IRES IPRE                          ;DV PRED RES WRES - APPENDED
       NOPRINT ONEHEADER FILE=sdtab


