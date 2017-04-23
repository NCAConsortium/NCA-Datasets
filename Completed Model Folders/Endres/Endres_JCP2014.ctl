;Model Desc: 2 CMT EV NONLINEAR  
;Project Name: NCAC SIMS
;Author: JK
;QC: BM
;Notes

$PROB RUN# SIM ENDRES

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT TIME TAD DOSE AMT SS II DV MDV EVID CMT TBW AGE DIAG

$DATA NMDATA IGNORE=@ ; EndresRun1REP1

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Subcutaneous Dosage Amount Data Item (mg)
; DV:   Dependent Variable (ug/mL)
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing)
; ADDL: Additional number of doses given after the first dosing record
; II:   Dosing Interval
; TBW:  Total Body Weight (kg)
; AGE:  Age(yr)
; DIAG: Disease diagnosis (0 = Healthy Volunteer, 1 = Psoriasis)

;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************

$SUBROUTINE ADVAN6 TOL=3

$MODEL NCOMPARTMENTS=3
NPARAMETERS=14
COMP (DEPOT,DEFDOSE)
COMP (CENTRAL, DEFOBS, NOOFF)
COMP (PERIF1 NOOFF NODOSE)

$PK

TVCL    = THETA(1)                   ; POPULATION CLEARANCE
GPCL=TVCL*((TBW/90)**THETA(8))*((AGE/40)**THETA(9))*(THETA(10)**DIAG) ; GRP CL
PPVCL   = EXP(ETA(1))                ; BSV CLEARANCE
CL      = GPCL * PPVCL               ; INDIVIDUAL CLEARANCE

TVV1    = THETA(2)                   ; POPULATION CENTRAL VOLUME
GPV1    = TVV1*((TBW/90)**THETA(11)) ; GROUP CENTRAL VOLUME WT
PPVV1   = EXP(ETA(2))
V1      = GPV1 * PPVV1               ; INDIVIDUAL CENTRAL VOLUME

TVQ     = THETA(3)                   ; POPULATION INTERCOMPARTMENTAL CLEARANCE 2 
GPQ     = TVQ*((TBW/90)**THETA(8))   ; GROUP INTERCOMPARTMENTAL CLEARANCE 2 WT
PPVQ    = EXP(ETA(3))                ; BSV INTERCOMPARTMENTAL CLEARANCE 2 WT
Q       = GPQ * PPVQ                 ; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE 2

TVV2    = THETA(4)                   ; POPULATION VOLUME 2
GPV2    = TVV2*((TBW/90)**THETA(11)) ; GROUP VOLUME 2 
PPVV2   = EXP(ETA(4))
V2      = GPV2 * PPVV2               ; INDIVIDUAL VOLUME 2

TVKA    = THETA(5)                   ; POPULATION SC KA (ABSORPTION RATE)
GPKA    = TVKA                       ; GROUP SC KA
PPVKA   = EXP(ETA(5))                ; BSV KA
KA      = GPKA * PPVKA               ; INDIVIDUAL SC KA

TVVMAX  = THETA(6)                   ; POPULATION VMAX 
BWVMAX  = (TBW/90)**THETA(12)
GPVMAX  = TVVMAX*BWVMAX*((AGE/40)**THETA(13))*(THETA(14)**DIAG)
PPVMAX  = EXP(ETA(6))
VMAX    = GPVMAX * PPVMAX            ; INDIVIDUAL VMAX

KM      = THETA(7)                   ; IND FIXED KM

$DES

C1=A(1)
C2=A(2)/V1
C3=A(3)/V2
DADT(1)=-KA*C1
DADT(2)=(KA*C1)+(Q*(C3-C2))-(CL*C2)-(C2*VMAX/(KM+C2))
DADT(3)=Q*(C2-C3)

      
;------------------------------------------------------------
;----           TIME AFTER DOSE CALCULATION              ----
;------------------------------------------------------------



;------------------------------------------------------------
;----           SECONDARY MODEL PARAMETERS               ----
;------------------------------------------------------------

;    K = CL/V1

;    T12 = 0.693 / K
;    MRT = 1/K
;    AREA = DOSE/CL
;    AUMC = AREA*MRT
     VSSF = V1 + V2
    
;------------------------------------------------------------
;----               RESIDUAL ERROR MODEL                 ----
;------------------------------------------------------------
$ERROR

IPRE=F
Y    = F + F*EPS(1) + EPS(2)
IRES = DV-IPRE
IWRE = IRES/IPRE

;************************************************************
;****           INITIAL PARAMETER ESTIMATES              ****
;************************************************************

$THETA
0.223                   ; CLEARANCE (L/DAY)
4.62                    ; CENTRAL VOLUME (L) 
0.697                   ; INTERCOMPARTMENTAL CLEARANCE  (L/DAY)
1.84                    ; PERIPHERAL VOLUME 2 (L)
0.236                   ; KA SC ABSORPTION RATE CONSTANT (1/DAY)
5.40                    ; VMAX (MG/DAY)
0.02                    ; KM (UG/ML) 
0.598                   ; TBW ON CL & Q
0.468                   ; AGE ON CL
0.927                   ; DIAG ON CL
0.849                   ; TBW ON V1 & V2
1.12                    ; TBW ON VMAX
-0.320                  ; AGE ON VMAX
1.04                    ; DIAG ON VMAX

$OMEGA
0.391                   ; [E] BETWEEN SUBJECT VARIABILITY (CL)
0.395                   ; [E] BETWEEN SUBJECT VARIABILITY (V1)
0.0223                  ; [E] BETWEEN SUBJECT VARIABILITY (Q)
0.55                    ; [E] BETWEEN SUBJECT VARIABILITY (V2)
0.289                   ; [E] BETWEEN SUBJECT VARIABILITY (KA)
0.0649                  ; [E] BETWEEN SUBJECT VARIABILITY (VMAX)

$SIGMA
0.0142                  ; [P] PROPORTIONAL RESIDUAL ERROR
0.384                   ; [A] ADDITIVE RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************

$SIM ONLYSIM (7284605)

$TABLE ID TIME II AMT DOSE DV CMT EVID PRED MDV IWRE IRES DIAG
IPRE CL V1 Q V2 KA VMAX KM ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 
TBW AGE VSSF
NOPRINT ONEHEADER FILE=sdtab
