;Model Desc: 2 CMT IVB LINEAR
;Project Name: NCAC SIMS
;Author: BM
;QC:

$PROB RUN# 1

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TIME AMT DV EVID MDV CR60 BSA ALB

$DATA NMDATA IGNORE=C ; NMDATA = analysis.data.set_28JUN2016.csv

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Intravenous or subcutaneous Dosage Amount Data Item (IU)
; DV:   Dependent Variable (IU/dL);
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing);
; CR60: Creatinine Clearance (with 60 umol/L cutoff)
; BSA:  Body Surface Area (m**2)
; ALB:  Albumin (g/L)

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************

$SUBROUTINE ADVAN6 TRANS1

$PK

    TVCL   = THETA(1)            								  ; POPULATION CLEARACE
  	GPCL   = TVCL * (1 + THETA(2) * CR60)						  ; GROUP CLEARANCE
    PPVCL  = EXP(ETA(1))										  ; BSV CLEARANCE
    CL     = GPCL * PPVCL        								  ; INDIVIDUAL CLEARANCE

    TVV1   = THETA(3)            							      ; POPULATION CENTRAL VOLUME
	GPV1   = TVV1 * (BSA * (ALB/34) ** THETA(4)		     	      ; GROUP CENTRAL VOLUME 
    V1     = GPV1   										      ; INDIVIDUAL CENTRAL VOLUME

    TVQ    = THETA(5)                                             ; POPULATION DISTRIBUTIONAL CLEARANCE
    PPVQ   = EXP(ETA(2))                                          ; BSV DISTRIBUTIONAL CLEARANCE
    Q      = TVQ * PPVQ                                           ; INDIVIDUAL DISTRIBUTIONAL CLEARANCE

    TVV2   = THETA(6)                                             ; POPULATION INTERCOMPARTMENTAL CLEARANCE 2 
    V2     = TVV2                                                 ; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE 2

    S1     = V                                                    ; SCALING FACTOR
    S2     = V2                                                   ; SCALING FACTOR
      
;------------------------------------------------------------
;----           TIME AFTER DOSE CALCULATION              ----
;------------------------------------------------------------



;------------------------------------------------------------
;----           SECONDARY MODEL PARAMETERS               ----
;------------------------------------------------------------


    
;------------------------------------------------------------
;----               RESIDUAL ERROR MODEL                 ----
;------------------------------------------------------------
$ERROR

IPRE=F
Y    = F * EXP(EPS(1)) + EPS(2)
IRES = DV-IPRE
IWRE = IRES/IPRE

;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************

$THETA
0.88				            	 ; TV CLEARANCE (L/H) CL
0.043								 ; COV CR60 ~ CL
8.59							     ; TV V1 (L)
-0.39        				     	 ; COV ALB ~ V1
1.30        				     	 ; TV Q (L/H)
9.79        			  	   	     ; TV V2 (L)
  
$OMEGA
0.0342                     ; [E] BETWEEN SUBJECT VARIABILITY (CL)
0.0795                     ; [E] BETWEEN SUBJECT VARIABILITY (Q)

$SIGMA
0.025                     ; [P] PROPORTIONAL RESIDUAL ERROR
0.107                     ; [P] ADDITIVE RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************

$SIM ONLYSIM (8675309) SUBPROBLEMS=1000

$TABLE                                       ; STANDARD TABLE
       ID TIME AMT DV EVID MDV CR60 BSA ALB DV PRED MDV IWRE IRES
       IPRE CL Q V1 V2 ETA1 ETA2
       NOPRINT ONEHEADER FILE=sdtab


