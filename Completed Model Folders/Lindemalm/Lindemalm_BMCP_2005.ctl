;Model Desc: 3 CMT EV LINEAR  
;Project Name: NCAC SIMS
;Author: JK
;QC: 
;Compound: Cladribine

$PROB RUN# SIM Lindemalm SC Route

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT TIME TAD DOSE AMT SS II DV MDV EVID CMT

$DATA NMDATA IGNORE=@ ; NMDATA = LINDEMALMRUN1REP1

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Subcutaneous Dosage Amount Data Item (mg)
; DV:   Dependent Variable (mg/L)
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing)
; ADDL: Additional number of doses given after the first dosing record
; II:   Dosing Interval
; F:    Bioavailabilty = 1 for SC Route

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************

$SUBROUTINE ADVAN12 TRANS4

$PK

TVCL   = THETA(1)                       ; POPULATION CLEARANCE
GPCL   = TVCL                           ; GROUP CLEARANCE
PPVCL  = EXP(ETA(1))                    ; BSV CLEARANCE
CL     = GPCL * PPVCL                   ; INDIVIDUAL CLEARANCE

TVV2   = THETA(2)                       ; POPULATION CENTRAL VOLUME
GPV2   = TVV2                           ; GROUP CENTRAL VOLUME WT
PPVV2  = EXP(ETA(2))
V2     = GPV2 * PPVV2                   ; INDIVIDUAL CENTRAL VOLUME

TVQ3   = THETA(3)                       ; POPULATION INTERCOMPARTMENTAL CLEARANCE 2 
GPQ3   = TVQ3                           ; GROUP INTERCOMPARTMENTAL CLEARANCE 2 WT
PPVQ3  = EXP(ETA(3))                    ; BSV INTERCOMPARTMENTAL CLEARANCE 2 WT
Q3     = GPQ3 * PPVQ3                   ; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE 2

TVV3   = THETA(4)                       ; POPULATION VOLUME 2
GPV3   = TVV3                           ; GROUP VOLUME 2 WT
PPVV3  = EXP(ETA(4))
V3     = GPV3 * PPVV3                   ; INDIVIDUAL VOLUME 2

TVQ4   = THETA(5)                       ; POPULATION INTERCOMPARTMENTAL CLEARANCE 3
GPQ4   = TVQ4                           ; GROUP INTERCOMPARTMENTAL CLEARANCE 3 WT
PPVQ4  = EXP(ETA(5))                    ; BSV INTERCOMPARTMENTAL CLEARANCE 3 WT
Q4     = GPQ4 * PPVQ4                   ; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE 3	

TVV4   = THETA(6)                       ; POPULATION VOLUME 3
GPV4   = TVV4                           ; GROUP VOLUME 3 WT
PPVV4  = EXP(ETA(6))                    ; BSV VOLUME 3 WT
V4     = GPV4 * PPVV4                   ; INDIVIDUAL VOLUME 3

TVKA   = THETA(7)                       ; POPULATION SC KA (ABSORPTION RATE CONSTANT)
GPKA   = TVKA                           ; GROUP SC KA
KA     = GPKA                           ; INDIVIDUAL SC KA

      
;------------------------------------------------------------
;----           TIME AFTER DOSE CALCULATION              ----
;------------------------------------------------------------



;------------------------------------------------------------
;----           SECONDARY MODEL PARAMETERS               ----
;------------------------------------------------------------

K10 = CL/V2
K12 = Q3/V2
K21 = Q3/V3
K13 = Q4/V2
K31 = Q4/V4

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

     T12 = 0.693/L3
     MRT = 1/L3
     AREA = DOSE/CL
     AUMC = AREA*MRT
     VSS = V2 + V3 + V4
    
;------------------------------------------------------------
;----               RESIDUAL ERROR MODEL                 ----
;------------------------------------------------------------
$ERROR

IPRE = F
Y    = F + F*EPS(1)
IRES = DV-IPRE
IWRE = IRES/IPRE

;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************

$THETA
39.3                    ; CLEARANCE (L/HR)
71.7                    ; CENTRAL VOLUME (L) 
51.1                    ; INTERCOMPARTMENTAL CLEARANCE 2 (L/HR)
475                     ; PERIPHERAL VOLUME 2 (L)
105                     ; INTERCOMPARTMENTAL CLEARANCE 3 (L/HR)
73.6                    ; PERIPHERAL VOLUME 3 (L)
2.48                    ; KA SC ABSORPTION RATE CONSTANT (1/HR)

$OMEGA
0.226                   ; [E] BETWEEN SUBJECT VARIABILITY (CL)
0.109                   ; [E] BETWEEN SUBJECT VARIABILITY (V2)
0.316                   ; [E] BETWEEN SUBJECT VARIABILITY (Q3)
0.399                   ; [E] BETWEEN SUBJECT VARIABILITY (V3)
0.316                   ; [E] BETWEEN SUBJECT VARIABILITY (Q4)
0.316                   ; [E] BETWEEN SUBJECT VARIABILITY (V4)

$SIGMA
0.162                   ; [P] PROPORTIONAL RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************

$SIM ONLYSIM (5471207)

$TABLE                                       ; STANDARD TABLE
       ID TIME II AMT DOSE DV CMT EVID PRED MDV IWRE IRES
       IPRE CL V2 Q3 V3 Q4 V4 KA T12 MRT AUMC VSS AREA
       ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 
       NOPRINT ONEHEADER FILE=sdtab
