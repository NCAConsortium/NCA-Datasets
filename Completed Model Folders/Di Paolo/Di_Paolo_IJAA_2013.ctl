;Model Desc: 1 CMT INFUSION LINEAR
;Project Name: NCAC SIMS
;Author: BM
;QC: JK
;QC2: BM

$PROB RUN# 1

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT NTIM NTAD TIME TAD DOSE AMT DUR RATE=DROP 
       SS=DROP II=DROP DV MDV EVID CMT AGE WTKG SCR SEX

$DATA NMDATA IGNORE=@ ; NMDATA = DiPaoloRun1.csv, DiPaoloRun2.csv

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Intravenous Dosage Amount Data Item (mg)
; DUR:  Infusion Duration (hr)
; DV:   Dependent Variable (mg/L)
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing)
; AGE:  Age (yr)
; WTKG: Body Weight (kg)
; SCR:  Serum Creatinine (mg/dL)
; SEX:  Gender (0 - Male, 1 - Female)

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************
$SUBROUTINE ADVAN1 TRANS2

$PK

    CRCL   = ((140-AGE)*WTKG/(SCR*72))*(1-(0.15*SEX))    ; CALCULATE CRCL

    TVCL   = THETA(1)                                    ; POPULATION CLEARACE
    GPCL   = TVCL * (CRCL/80) ** THETA(2)                ; GROUP CLEARANCE
    PPVCL  = EXP(ETA(1))                                 ; BSV CLEARANCE
    CL     = GPCL * PPVCL                                ; INDIVIDUAL CLEARANCE

    TVV    = THETA(3)                                    ; POPULATION CENTRAL VOLUME
    V      = TVV                                         ; INDIVIDUAL CENTRAL VOLUME

;------------------------------------------------------------
;----           SECONDARY MODEL PARAMETERS               ----
;------------------------------------------------------------
    K10 = CL/V                                           ; ELIMINATION RATE CONSTANT (1/h)

    THALF  = LOG(2)/K10                                  ; HALF-LIFE(h)
    MRT    = 1/K10                                       ; MEAN RESIDENT TIME (h) 
    AREA   = DOSE/CL                                     ; MEAN RESIDENT TIME (mg*hr/L) 
    AUMC   = AREA*MRT                                    ; AREA UNDER THE CURVE (mg*hr**2/L)
    VSS    = V                                           ; VOLUME STEADY-STATE (L)

;------------------------------------------------------------
;----               RESIDUAL ERROR MODEL                 ----
;------------------------------------------------------------
$ERROR

IPRE = F
Y    = F *(1 + EPS(1)) + EPS(2)
IRES = DV-IPRE
IWRE = IRES/IPRE

;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************
$THETA
0.8016                               ; TV CLEARANCE (L/H) CL
0.2026                               ; COV CRCL ~ CL
12.29                                ; TV V (L)
  
$OMEGA
0.0421                               ; [E] BETWEEN SUBJECT VARIABILITY (CL)

$SIGMA
0.1316                               ; [P] PROPORTIONAL RESIDUAL ERROR
1.422                                ; [P] ADDITIVE RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************
$SIM ONLYSIM (8675309)

$TABLE                                         ; STANDARD TABLE
       ID NTIM NTAD TIME TAD DOSE AMT DUR EVID MDV CMT
       AGE WTKG SCR SEX CRCL
       CL V ETA1 
       THALF MRT AREA AUMC VSS
       IRES IPRE                          ;DV PRED RES WRES - APPENDED
       NOPRINT ONEHEADER FILE=sdtab

