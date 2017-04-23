;Model Desc: 1 CMT EV NONLINEAR  
;Project Name: NCAC SIMS
;Author: JK
;QC: BM

$PROB RUN# SIM Ding

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT NTIM NTAD TIME TAD DOSE AMT SS II DV MDV EVID 
CMT BW TDD CBZ

$DATA NMDATA IGNORE=@ ; NMDATA = DingRun2REP1

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Oral Dosage Amount Data Item (mg)
; DV:   Dependent Variable (mg/L)
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing)
; ADDL: Additional number of doses given after the first dosing record
; II:   Dosing Interval
; BW:   Weight (kg)
; TDD:  Total Daily Dose (mg/kg)
; CBZ:  Carbamazepine coadmin = 1, otherwise = 0

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************

$SUBROUTINE ADVAN6 TOL=3

$MODEL NCOMPARTMENTS=2
NPARAMETERS=8
COMP (DEPOT, DEFDOSE)
COMP (CENTRAL, DEFOBS)

$PK

TVEMAX= THETA(4)                                     ; POPULATION EMAX 
EM    = TVEMAX                                       ; INDIVIDUAL EMAX

TVGAM = THETA(5)                                     ; POP FIXED HILL COEFFICIENT
GAM   = TVGAM                                        ; IND FIXED HILL COEFFICIENT

TDD50 = THETA(6)                                     ; POP TDD50

TVCL  = THETA(1)                                     ; POPULATION CLEARANCE
GPCL1 = TVCL*((BW/70)**THETA(7))*(THETA(8)**CBZ)     ; GRP CL
GPCL2 = GPCL1*(1+(EM*TDD**GAM)/(TDD50**GAM+TDD**GAM)); GRP CL
PPVCL = EXP(ETA(1))                                  ; BSV CLEARANCE
CL    = GPCL2 * PPVCL                                ; INDIVIDUAL CLEARANCE

TVV1  = THETA(2)                                     ; POPULATION CENTRAL VOLUME
GPV1  = TVV1*(BW/70)                                 ; GROUP CENTRAL VOLUME WT
V1    = GPV1                                         ; INDIVIDUAL CENTRAL VOLUME

TVKA  = THETA(3)                                     ; POPULATION PO TAB KA (ABSORPTION RATE)
GPKA  = TVKA                                         ; GROUP PO TAB KA
KA    = GPKA                                         ; INDIVIDUAL PO TAB KA





$DES

C1=A(1)
C2=A(2)/V1
DADT(1)=-KA*C1
DADT(2)=(KA*C1)-(CL*C2)


      
;------------------------------------------------------------
;----           TIME AFTER DOSE CALCULATION              ----
;------------------------------------------------------------



;------------------------------------------------------------
;----           SECONDARY MODEL PARAMETERS               ----
;------------------------------------------------------------

;K = CL/V1

;    T12 = 0.693 / K
;    MRT = 1/K
;    AREA = DOSE/CL
;    AUMC = AREA*MRT
     VSSF = V1

 
;------------------------------------------------------------
;----               RESIDUAL ERROR MODEL                 ----
;------------------------------------------------------------
$ERROR

IPRE = F
Y    = F + EPS(1)
IRES = DV-IPRE
IWRE = IRES/IPRE

;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************

$THETA
0.300   ; CLEARANCE (L/HR)
22.2    ; CENTRAL VOLUME (L) 
1.57    ; KA PO TAB ABSORPTION RATE CONSTANT (1/HR)
2.8     ; EMAX
1.68    ; HILL COEFFICIENT
37.4    ; TDD50 (MG/KG) (TDD WHEN EMAX IS INCREASED 50%)
0.667   ; EXPONENT ON WEIGHT FOR CLEARANCE 
1.43    ; EFFECT OF CBZ ON CLEARANCE

$OMEGA
0.178   ; [E] BETWEEN SUBJECT VARIABILITY (CL)

$SIGMA
13.3    ; [A] ADDITIVE RESIDUAL ERROR (MG/L)

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************

$SIM ONLYSIM (3776054)

$TABLE ID TIME II DOSE AMT DV CMT EVID PRED MDV 
IWRE IRES IPRE CL V1 KA EM GAM TDD50 ETA1 BW CBZ TDD VSSF
NOPRINT ONEHEADER FILE=sdtab
