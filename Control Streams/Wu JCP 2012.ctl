;Model Desc: 1 CMT INFUSION NONLINEAR
;Project Name: NCAC SIMS
;Author: BM
;QC: JK
;QC2: BM
;Notes: Scaling factor for central compartment needed for 
;       mg into ug. Need a concentration in differential to
;       have the units match up. C1 = A(1)/V, then include
;       C1 into the differential where A(1) is. Prop error
;       should be 0.020736 based on article.

$PROB RUN# 1

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT TIME TAD DOSE AMT DUR RATE SS II DV MDV EVID CMT BSA LT

$DATA NMDATA IGNORE=C ; NMDATA = analysis.data.set_28JUN2016.csv

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Intravenous or subcutaneous Dosage Amount Data Item (mg)
; RATE: Infusion Rate (mg/hr)
; DV:   Dependent Variable (ug/L);
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing)
; BSA:  Body Surface Area (m**2)
; LT:   Presence of Liver Tumor (0 - Absent, 1 - Present)

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************

$SUBROUTINE ADVAN6 TRANS1 TOL=4

$MODEL COMP=(CENTRAL,DEFDOSE)           ;cmt=1

$PK

    TVV   = THETA(1)                                             ; POPULATION CLEARACE
    PPVV  = EXP(ETA(1))
    V     = TVV * PPVV

    IF (LT.EQ.0) THEN
    TVVMAX = THETA(2)
    END IF
    IF (LT.EQ.1) THEN
    TVVMAX = THETA(3)
    END IF
    PPVVMA  = EXP(ETA(2))                                        ; BSV CLEARANCE
    VMAX     = TVVMAX * PPVVMA                                   ; INDIVIDUAL CLEARANCE

    TVKM   = THETA(4)                                             ; POPULATION CENTRAL VOLUME
    KM     = TVKM                                                 ; INDIVIDUAL CENTRAL VOLUME

    S1     = V/1000                                               ; SCALING FACTOR ug/L to mg/L

$DES
    C1 = A(1)/V

    DADT(1) = -((VMAX*C1)/(KM+C1))

;------------------------------------------------------------
;----           SECONDARY MODEL PARAMETERS               ----
;------------------------------------------------------------

VSS = V
    
;------------------------------------------------------------
;----               RESIDUAL ERROR MODEL                 ----
;------------------------------------------------------------
$ERROR

IPRE=F
Y    = F * (1 + EPS(1)) + EPS(2)
IRES = DV-IPRE
IWRE = IRES/IPRE

;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************

$THETA
3.63                                 ; TV VOLUME (L) V
97.4                                 ; TV VMAX (ug/h) NO LIVER TUMOR
150                                  ; TV VMAX (ug/h) LIVER TUMOR
992                                  ; KM (ug/L)

  
$OMEGA
0.4177                               ; [E] BETWEEN SUBJECT VARIABILITY (V)
1.649                                ; [E] BETWEEN SUBJECT VARIABILITY (VMAX)

$SIGMA
0.02074                              ; [P] PROPORTIONAL RESIDUAL ERROR
10.9                                 ; [P] ADDITIVE RESIDUAL ERROR


;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************

$SIM ONLYSIM (8675309) 

$TABLE                                       ; STANDARD TABLE
       ID TIME AMT DV EVID MDV BSA LT 
       V VMAX KM ETA1 ETA2
       VSS
       IWRE IRES IPRE                          ;DV PRED RES WRES - APPENDED
       NOPRINT ONEHEADER FILE=sdtab



