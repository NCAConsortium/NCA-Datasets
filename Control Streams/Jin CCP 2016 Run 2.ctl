;Model Desc: 2 CMT ORAL LINEAR W/ LAG
;Project Name: NCAC SIMS Jin
;Author: JK
;QC: BM

$PROB RUN# SIM JIN

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT NTIM NTAD TIME TAD DOSE AMT SS II DV MDV EVID 
CMT WT PAT  

$DATA NMDATA IGNORE=@ ; NMDATA = JINRUN1REP1.csv
  
; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Oral Dosage Amount Data Item (mg)
; DV:   Dependent Variable (ng/mL)
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing)
; ADDL: Number of additional doses given after the first
; II:   Dosing Interval 
; WT:   Weight (kg)
; PAT:  Patient Status (1=patient or 2=healthy volunteer)
; DOSE: Dose given to the subject

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************  

$SUBROUTINE ADVAN4 TRANS4

$PK
 CLT=THETA(1) ; POPULATION CLEARANCE
 IF (PAT.EQ.2) CLT=THETA(8)
 QT=THETA(3)  ; POPULATION INTERCOMPARTMENTAL CLEARANCE
 IF (PAT.EQ.2) QT=THETA(9)

 FWT=WT
; IF (ID.EQ.161) FWT=75
; IF (ID.EQ.211) FWT=75
; IF (ID.EQ.421) FWT=75
; IF (ID.EQ.434) FWT=75
; IF (ID.EQ.480) FWT=75
 SWT=THETA(10)*LOG(FWT/75)

 ; PK parameters
 CL=EXP(CLT+SWT+ETA(1))         ; INDIVIDUAL CLEARANCE
 V2=EXP(THETA(2)+ETA(2))        ; INDIVIDUAL VOLUME
 Q=EXP(QT+ETA(3))               ; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE
 V3=EXP(THETA(4)+ETA(4))        ; INDIVIDUAL PERIPHERAL VOLUME
 KA=EXP(THETA(5)+ETA(5))        ; INDIVIDUAL ABSORPTION RATE CONSTANT
 ALAG1=EXP(THETA(6)+ETA(6))     ; INDIVIDUAL LAG TIME
 F1=EXP(THETA(7)*LOG(DOSE/150)) ; INDIVIDUAL INFLUENCE OF DOSE ON BIOAVAILABILITY

 S2=V2/1000      ; SCALING FACTOR
 
;------------------------------------------------------------
;----           SECONDARY MODEL PARAMETERS               ----
;------------------------------------------------------------
 
 K20=CL/V2
 K23=Q/V2
 K32=Q/V3
  
A0 = K20 * K32
A1 = K20 + K23 + K32

RT1 = (A1 + SQRT(A1*A1 - 4* A0)) / 2
RT2 = (A1 - SQRT(A1*A1 - 4* A0)) / 2

IF(RT1.LT.RT2) THEN
    BETA = RT1
    ALPHA = RT2
ENDIF
IF(RT2.LT.RT1) THEN
  BETA = RT2
  ALPHA = RT1
ENDIF

    T12 = 0.693 / BETA
    MRT = 1/BETA
    AREA = (F1*DOSE/CL)*1000 ; (ng*h/mL) and bioavailability
    AUMC = AREA*MRT
    VSS = V2 + V3 
 
;------------------------------------------------------------
;----               RESIDUAL ERROR MODEL                 ----
;------------------------------------------------------------
$ERROR
  IPRED=F
  IF(F.GT.0) THEN
  IPRED=LOG(F)
  ELSE
  IPRED=0 
  ENDIF
  Y=IPRED+EPS(1)
  IRES=DV-IPRED
  IWRES=IRES/IPRED
  
;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************  

$THETA 
 (2.7)      ; CL (L/hr)
 (3.12)     ; V2 (L)
 (2.47)     ; Q  (L/hr)
 (4.29)     ; V3 (L)
 (-.729)    ; KA (1/hr)
 (-1.398)   ; LAG (hr)
 (-.262)    ; DOSE ON F1
 (2.98)     ; CL HV
 (2.06)     ; Q HV
 (0.245)    ; BW on CL

$OMEGA BLOCK(2) 
.146               ;IIV ON CL
.112 .725          ;IIV ON V2
$OMEGA BLOCK(2) 
.151               ;IIV ON Q
.231 .538          ;IIV ON V3
$OMEGA .147 .207   ;IIV ON KA AND LAG

$SIGMA .286        ; RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************

$SIM ONLYSIM (4621978)

$TABLE ID TIME AMT DOSE CMT DV EVID MDV WT PAT II PRED IRES IWRES 
IPRED CL V2 Q V3 KA ALAG1 F1 T12 MRT AREA AUMC VSS 
NOPRINT ONEHEADER FILE=sdtab
