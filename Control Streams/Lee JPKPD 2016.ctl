;Model Desc: 3 CMT IVB LINEAR
;Project Name: NCAC SIMS
;Author: BM
;QC: JK
;QC2: BM
;Notes: Scale needed for central compartment with dose in mg
;       and DV in ng/mL. Also, authors used ADVAN6 but shouldn't
;       lead to any issues in sim.
$PROB RUN# 1

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT TIME TAD DOSE AMT SS II DV MDV EVID CMT AGE WTKG

$DATA NMDATA IGNORE=@ ; NMDATA = LeeRun1REP1.csv

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Intravenous or subcutaneous Dosage Amount Data Item (mg)
; DV:   Dependent Variable (ng/mL);
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing);
; AGE:  Age (yr)
; WTKG: Weight (kg)

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************

$SUBROUTINE ADVAN11 TRANS4

$PK

    TVCL  = THETA(1)                                              ; POPULATION CLEARANCE
    GPCL  = TVCL * (WTKG/60)**0.75                                ; ALLOMETRIC WTKG ~ CL
    IF (AGE.GE.THETA(2)) THEN
    GP2CL = GPCL * EXP(THETA(3) * (AGE-THETA(2)))                 ; COV AGE >= INF ~ CL
    ELSE
    GP2CL = GPCL                                                  ; COV AGE < INF ~ CL
    ENDIF
    PPVCL = EXP(ETA(1))
    CL    = GP2CL * PPVCL                                         ; INDIVIDUAL CLEARANCE

    TVV1   = THETA(4)                                             ; POPULATION CENTRAL VOLUME
    GPV1   = TVV1 * (WTKG/60)                                       ; ALLOMETRIC WTKG ~ V1
    V1     = GPV1                                                 ; INDIVIDUAL CENTRAL VOLUME

    TVV2   = THETA(5)                                             ; POPULATION VOLUME 2
    GPV2   = TVV2 * (WTKG/60)                                       ; ALLOMETRIC WTKG ~ V3
    V2     = GPV2                                                 ; INDIVIDUAL VOLUME 2

    TVV3   = THETA(6)                                             ; POPULATION VOLUME 3
    GPV3   = TVV3 * (WTKG/60)                                       ; ALLOMETRIC WTKG ~ V3
    PPVV3  = EXP(ETA(2))                                          ; BSV VOLUME 3 WTKG
    V3     = GPV3 * PPVV3                                         ; INDIVIDUAL VOLUME 3

    TVQ2   = THETA(7)                                             ; POPULATION INTERCOMPARTMENTAL CLEARANCE 2 
    GPQ2   = TVQ2 * (WTKG/60)**0.75                                 ; GROUP INTERCOMPARTMENTAL CLEARANCE 2 WTKG
    PPVQ2  = EXP(ETA(3))                                          ; BSV INTERCOMPARTMENTAL CLEARANCE 2 WTKG
    Q2     = GPQ2 * PPVQ2                                         ; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE 2

    TVQ3   = THETA(8)                                             ; POPULATION INTERCOMPARTMENTAL CLEARANCE 3
    GPQ3   = TVQ3 * (WTKG/60)**0.75                                 ; GROUP INTERCOMPARTMENTAL CLEARANCE 3 WTKG
    PPVQ3  = EXP(ETA(4))                                          ; BSV INTERCOMPARTMENTAL CLEARANCE 3 WTKG
    Q3     = GPQ3 * PPVQ3                                         ; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE 3

    S1     = V1/1000                                              ; SCALING FACTOR - AMT in mg, DV in ng/mL (ug/L)

;------------------------------------------------------------
;----           SECONDARY MODEL PARAMETERS               ----
;------------------------------------------------------------

    K10 = CL/V1
    K12 = Q2/V1
    K21 = Q2/V2
    K13 = Q3/V1
    K31 = Q3/V3

    A0 = K10 * K21 * K31
    A1 = K10 * K31 + K21 * K31 + K21 * K13 + K10 * K21 + K31 * K12
    A2 = K10 + K12 + K13 + K21 + K31
   
    PPAR = A1 - (A2 * A2 / 3)
    QPAR = (2 * A2 * A2 * A2 / 27) - (A1 * A2 / 3) + A0
    G1 = SQRT(-(PPAR * PPAR * PPAR) / 27)
    PHI = ACOS((-QPAR/2) / G1)/3
    G2 = 2 * EXP(LOG(G1) / 3)
    PI = 4 * ATAN(1)

    RT1 = -(COS(PHI) * G2 - A2 / 3)
    RT2 = -(COS(PHI + 2 * PI/3) * G2 - A2 / 3)
    RT3 = -(COS(PHI + 4 * PI/3) * G2 - A2 / 3)
   
    IF(RT1.LT.RT2.AND.RT1.LT.RT3) L3 = RT1
    IF(RT2.LT.RT1.AND.RT2.LT.RT3) L3 = RT2
    IF(RT3.LT.RT1.AND.RT3.LT.RT2) L3 = RT3

    THALF  = LOG(2)/L3                                  ; HALF-LIFE(h)
    MRT  = 1/L3                                         ; MEAN RESIDENT TIME (h) 
    AREA = DOSE/CL                                      ; MEAN RESIDENT TIME (mg*hr/L) 
    AUMC = AREA*MRT                                     ; AREA UNDER THE CURVE (mg*hr**2/L)
    VSS  = V1 + V2 + V3                                 ; VOLUME STEADY-STATE (L)

;------------------------------------------------------------
;----               RESIDUAL ERROR MODEL                 ----
;------------------------------------------------------------
$ERROR

IPRE=F
Y    = F + F * EPS(1)
IRES = DV-IPRE
IWRE = IRES/IPRE

;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************

$THETA
0.19                                ; TV CLEARANCE (L/HR) CL
57                                  ; COV AGE INF ~ CL
0.029                               ; COV AGE ~ CL
4.95                                ; TV V1 (L)
31.9                                ; TV V2 (L)
88.3                                ; TV V3 (L)
1.69                                ; TV Q2 (L/HR)
0.895                               ; TV Q3 (L/HR)
  
$OMEGA
0.1355                              ; [E] BETWEEN SUBJECT VARIABILITY (CL)
0.2786                              ; [E] BETWEEN SUBJECT VARIABILITY (V3)
0.0987                              ; [E] BETWEEN SUBJECT VARIABILITY (Q2)
0.3004                              ; [E] BETWEEN SUBJECT VARIABILITY (Q3)

$SIGMA
0.0151                              ; [P] PROPORTIONAL RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************

$SIM ONLYSIM (8675309)

$TABLE                                        ; STANDARD TABLE
       ID TIME AMT DV EVID MDV AGE WTKG
       CL V1 V2 V3 Q2 Q3 ETA1 ETA2 ETA3 ETA4
       THALF MRT AREA AUMC VSS
       IWRE IRES IPRE                          ;DV PRED RES WRES - APPENDED
       NOPRINT ONEHEADER FILE=sdtab


