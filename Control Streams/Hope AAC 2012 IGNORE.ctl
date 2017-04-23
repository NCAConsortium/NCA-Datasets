;Model Desc: 2 CMT IV/EV NONLINEAR W/ LAG 
;Project Name: NCAC SIMS
;Author: JK
;QC: 

$PROB RUN# SIM HOPE

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TIME AMT DV EVID MDV ADDL II WT

$DATA NMDATA LRECL=??? IGNORE=C ; NMDATA = XXX.csv

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Intravenous or subcutaneous Dosage Amount Data Item (IU)
; DV:   Dependent Variable (IU/dL)
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing)
; ADDL: Additional number of doses given after the first dosing record
; II:	Dosing Interval
; WT: Weight (kg)

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************

$SUBROUTINE ADVAN11 TRANS4

$PK

    TVCL  = THETA(1)           																	; POPULATION CLEARANCE
    GPCL  = TVCL *(WT/70)**0.75
	PPVCL = EXP(ETA(1))
    CL    = GPCL * PPVCL        																	; INDIVIDUAL CLEARANCE

    TVV1  = THETA(4)  																	; POPULATION CENTRAL VOLUME
   GPV1   = TVV1 * (WT/70)
   PPVV1  = EXP(ETA(2))
   V1     = GPV1 * PPVV1     																        ; INDIVIDUAL CENTRAL VOLUME

    TVV2   = THETA(5)                                           ; POPULATION VOLUME 2
	GPV2   = TVV2 * (WT/70)
	PPVV2  = EXP(ETA(3))
    V2     = GPV2 * PPVV2                                               ; INDIVIDUAL VOLUME 2

    TVV3   = THETA(6)                                             ; POPULATION VOLUME 3
    GPV3   = TVV3 * (WT/70)                                       ; GROUP VOLUME 3 WT
    PPVV3  = EXP(ETA(4))                                          ; BSV VOLUME 3 WT
    V3     = GPV3 * PPVV3                                         ; INDIVIDUAL VOLUME 3

    TVQ2   = THETA(7)                                             ; POPULATION INTERCOMPARTMENTAL CLEARANCE 2 
    GPQ2   = TVQ2 * (WT/70)**0.75                                 ; GROUP INTERCOMPARTMENTAL CLEARANCE 2 WT
    PPVQ2  = EXP(ETA(5))                                          ; BSV INTERCOMPARTMENTAL CLEARANCE 2 WT
    Q2     = GPQ2 * PPVQ2                                         ; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE 2

    TVQ3   = THETA(8)                                             ; POPULATION INTERCOMPARTMENTAL CLEARANCE 3
    GPQ3   = TVQ3 * (WT/70)**0.75                                 ; GROUP INTERCOMPARTMENTAL CLEARANCE 3 WT
    PPVQ3  = EXP(ETA(6))                                          ; BSV INTERCOMPARTMENTAL CLEARANCE 3 WT
    Q3     = GPQ3 * PPVQ3                                         ; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE 3

    S1     = V                                                    ; SCALING FACTOR
    S2     = V2                                                   ; SCALING FACTOR
    S3     = V3                                                   ; SCALING FACTOR
      
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
Y    = F + F*EPS(1)
IRES = DV-IPRE
IWRE = IRES/IPRE

;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************

$THETA
0.19				            	 ; TV CLEARANCE (L/HR) CL
57									  		 ; COV AGE INF ~ CL
0.029									     ; COV AGE ~ CL
4.95        				     	 ; TV V1
31.9        				     	 ; TV V2
88.3        			  	   	 ; TV V3
1.69        		  		   	 ; TV Q2
0.895        				   	   ; TV Q3
  
$OMEGA
0.1452                     ; [E] BETWEEN SUBJECT VARIABILITY (CL)
0.3215                     ; [E] BETWEEN SUBJECT VARIABILITY (V3)
0.1037                     ; [E] BETWEEN SUBJECT VARIABILITY (Q2)
0.3505                     ; [E] BETWEEN SUBJECT VARIABILITY (Q3)




$SIGMA
0.0151                     ; [P] PROPORTIONAL RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************

$SIM ONLYSIM (8675309) SUBPROBLEMS=1000

$TABLE                                       ; STANDARD TABLE
       ID TIME AMT DV EVID MDV WTKG AGE DV PRED MDV IWRE IRES
       IPRE CL V1 V2 V3 Q2 Q3 ETA1 ETA2 ETA3 ETA4 
       NOPRINT ONEHEADER FILE=sdtab


