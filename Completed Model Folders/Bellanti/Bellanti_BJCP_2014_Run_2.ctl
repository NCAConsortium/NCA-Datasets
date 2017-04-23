;Model Desc: 1 CMT EV LINEAR W/ LAG 
;Project Name: NCAC SIMS
;Author: JK
;QC: 

$PROB RUN# SIM Bellanti

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT TIME TAD DOSE AMT SS II DV MDV EVID CMT SEX

$DATA NMDATA IGNORE=@ ; BELLANTIRUN1REP1

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Oral Dosage Amount Data Item (mg)
; DV:   Dependent Variable (mg/L)
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing)
; ADDL: Additional number of doses given after the first dosing record
; II:   Dosing Interval
; SEX:  GENDER (0 = FEMALE, 1 = MALE)

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************

$SUBROUTINE ADVAN2 TRANS2

$PK


TVKA  = THETA(1)            ; POPULATION ABSORPTION CONSTANT
GPKA  = TVKA
PPVKA = EXP(ETA(3))
KA    = GPKA * PPVKA        ; INDIVIDUAL ABSORPTION CONSTANT

TVCL  = THETA(2)            ; POPULATION CLEARANCE
GPCL  = TVCL 
PPVCL = EXP(ETA(1))
CL    = GPCL * PPVCL        ; INDIVIDUAL CLEARANCE

IF(SEX.EQ.0)TVV = THETA(3)  ; FEMALE CENTRAL VOLUME
IF(SEX.EQ.1)TVV = THETA(4)  ; MALE CENTRAL VOLUME 
GPV   = TVV 
PPVV  = EXP(ETA(2))
V     = GPV * PPVV          ; INDIVIDUAL CENTRAL VOLUME

ALAG1 = THETA(5)            ; LAG TIME INTO ABSORPTION COMPARTMENT

      
;------------------------------------------------------------
;----           TIME AFTER DOSE CALCULATION              ----
;------------------------------------------------------------



;------------------------------------------------------------
;----           SECONDARY MODEL PARAMETERS               ----
;------------------------------------------------------------

    K10  = CL/V
    T12  = 0.693/K10
    MRT  = 1/K10
    AREAF = DOSE/CL ; Doesn't account for bioavailabilty
    AUMC = AREAF*MRT
    VSSF  = V

;------------------------------------------------------------
;----               RESIDUAL ERROR MODEL                 ----
;------------------------------------------------------------
$ERROR

W = 2.4
IPRE=F
Y    = F + W*EPS(1)
IRES = DV-IPRE
IWRE = IRES/IPRE

;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************

$THETA
8.2           ; TV KA (1/HR)
30.8          ; CLEARANCE (L/HR) CL
65.3          ; FEMALE VOLUME (L)
78.4          ; MALE VOLUME (L)
0.146         ; LAG TIME (1/HR)

$OMEGA BLOCK (2)
0.0554        ; [E] BETWEEN SUBJECT VARIABILITY (CL/F)
0.0339 0.245  ; [E] BETWEEN SUBJECT VARIABILITY (V/F)
 
$OMEGA 
0.689         ; [E] BETWEEN SUBJECT VARIABILITY (KA)


$SIGMA
0.00566       ; [P] PROPORTIONAL RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************

$SIM ONLYSIM (3773315)

$TABLE ID TIME SEX II AMT DOSE DV CMT EVID PRED MDV IWRE IRES
IPRE CL V T12 MRT AREAF AUMC VSSF ETA1 ETA2 ETA3  
NOPRINT ONEHEADER FILE=sdtab


