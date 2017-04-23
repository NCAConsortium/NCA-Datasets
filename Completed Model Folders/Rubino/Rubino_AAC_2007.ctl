;Model Desc: 1 CMT EV LINEAR 
;Project Name: NCAC SIMS
;Author: JK
;QC: 

$PROB RUN# SIM Rubino

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT TIME TAD DOSE AMT SS II DV MDV EVID CMT BSA WT

$DATA NMDATA IGNORE=@ ; RUBINORUN1REP1

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Oral Dosage Amount Data Item (mg)
; DV:   Dependent Variable (ug/mL)
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing)
; ADDL: Additional number of doses given after the first dosing record
; II:   Dosing Interval
; BSA:  Body surface area (m**2)
; WT:   Weight (kg)

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************

$SUBROUTINE ADVAN2 TRANS2

$PK


TVKA  = THETA(1)       ; POPULATION ABSORPTION CONSTANT
GPKA  = TVKA
PPVKA = EXP(ETA(3))
KA    = GPKA * PPVKA   ; INDIVIDUAL ABSORPTION CONSTANT

TVCL  = THETA(2)       ; POPULATION CLEARANCE
GPCL  = TVCL * BSA     ; GROUP CLEARANCE
PPVCL = EXP(ETA(1))
CL    = GPCL * PPVCL   ; INDIVIDUAL CLEARANCE

TVV   = THETA(3)       ; POPULATION CENTRAL VOLUME
GPV   = TVV * WT       ; GROUP CENTRAL VOLUME
PPVV  = EXP(ETA(2))
V     = GPV * PPVV     ; INDIVIDUAL CENTRAL VOLUME

      
;------------------------------------------------------------
;----           TIME AFTER DOSE CALCULATION              ----
;------------------------------------------------------------



;------------------------------------------------------------
;----           SECONDARY MODEL PARAMETERS               ----
;------------------------------------------------------------

    K10 = CL / V
    T12 = 0.693 / K10
    MRT = 1/K10
    AREAF = DOSE/CL
    AUMC = AREAF*MRT
    VSSF = V
    
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
0.936        ; TV KA (1/HR)
8.46         ; CLEARANCE/F (L/HR/M**2) CL
2.15         ; VOLUME/F (L/KG)

$OMEGA BLOCK(2)
0.0845       ; [E] BETWEEN SUBJECT VARIABILITY (CL/F)
0.0578 0.058 ; [E] BETWEEN SUBJECT VARIABILITY (V/F)
 
$OMEGA 
0.534        ; [E] BETWEEN SUBJECT VARIABILITY (KA)


$SIGMA
0.063        ; [P] PROPORTIONAL RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************

$SIM ONLYSIM (8960714)

$TABLE ID TIME BSA WT II DOSE AMT DV CMT EVID MDV PRED IWRE 
IRES IPRE CL V KA T12 MRT AREAF AUMC VSSF ETA1 ETA2 ETA3  
NOPRINT ONEHEADER FILE=sdtab
