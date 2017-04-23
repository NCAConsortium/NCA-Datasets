;Model Desc: 3 CMT IV NONLINEAR  
;Project Name: NCAC SIMS
;Author: JK
;QC: BM
;--> I think you'll either need to include a scaling factor or convert to ug for 
;    the dose in your simulation dataset
;--> Not sure if the non-linear equations are right for compartment two - for the 
;    movement of drug from V2 into V1, it looks like C1 is in the denominator. Should 
;    it be C2? Likewise, should C2 be in the demoninator for drug leaving V2? May
;    just have to write the equations out in your differential
; JK Scaling added for V1
;	 Made 2 distributional clearances (Q1A and Q1B) for drug to and from V1 to V2

$PROB RUN# SIM CHEN

;************************************************************
;****                   DATASET                          ****
;************************************************************
$INPUT C ID TRT TIME TAD DOSE AMT DUR RATE SS II DV MDV EVID CMT ALB

$DATA  ChenJCPRUN1REP1.csv IGNORE=@

; C:    Comment (Any row starting with a 'C' will be ignored)
; ID:   NONMEM Subject Identification Number 
; TIME: Actual Time Field
; AMT:  Intravenous Dosage Amount Data Item (mg)
; DV:   Dependent Variable (ug/L)
; EVID: Event Identification (0 = Obs, 1 = Dosing)
; MDV:  Missing Dependent Variable (0 = Value, 1 = Missing)
; ADDL: Additional number of doses given after the first dosing record
; II:	Dosing Interval
; DUR:  Duration of infustion (h)
; RATE: Rate of infusion (mg/h) 0.5 hour infusion (Rate = 2x Dose)
; ALB:  Albumin level (g/dL)


;************************************************************
;****         PHARMACOKINETIC POPULATION MODEL           ****
;************************************************************

$SUBROUTINE ADVAN6 TOL=3

$MODEL 	NCOMPARTMENTS=3
		NPARAMETERS=14
		COMP (CENTRAL, DEFDOSE, DEFOBS)
		COMP (PERIF1 NOOFF NODOSE)
		COMP (PERIF2 NOOFF NODOSE)
		
$PK

	TVVMEL	= THETA(1)            				; POPULATION CENTRAL MAX ELIMINATION
  	GPVMEL	= TVVMEL * ((ALB/3.9)**THETA(9))	; GRP VMEL
    PPVVMEL = EXP(ETA(1))						; BSV VMEL
    VMEL    = GPVMEL * PPVVMEL        			; INDIVIDUAL VMEL
	
	TVKMEL	= THETA(2)            				; POPULATION CONC @ 50% OF VMEL
  	GPKMEL	= TVKMEL							; GRP KMEL
    KMEL    = GPKMEL 			       			; INDIVIDUAL KMEL

    TVV1   = THETA(3)            				; POPULATION CENTRAL VOLUME
	GPV1   = TVV1								; GROUP CENTRAL VOLUME WT
	PPVV1  = EXP(ETA(2))
    V1     = GPV1 * PPVV1        				; INDIVIDUAL CENTRAL VOLUME

	TVVMTR	= THETA(4)            				; POPULATION INTERCMT MAX DISTRIB
  	GPVMTR	= TVVMTR							; GRP VMTR
    PPVVMTR = EXP(ETA(3))						; BSV VMTR
    VMTR    = GPVMTR * PPVVMTR        			; INDIVIDUAL VMTR
	
	TVKMTR	= THETA(5)            				; POPULATION CONC @ 50% OF VMTR
  	GPKMTR	= TVKMTR							; GRP KMTR
    KMTR    = GPKMTR 			       			; INDIVIDUAL KMTR	
	
    TVV2   = THETA(6)                       	; POPULATION VOLUME 2
    GPV2   = TVV2						        ; GROUP VOLUME 2 
	PPVV2  = EXP(ETA(4))
    V2     = GPV2 * PPVV2                   	; INDIVIDUAL VOLUME 2
	
	TVQ2   = THETA(7)                       	; POPULATION INTERCOMPARTMENTAL CLEARANCE 2 
    GPQ2   = TVQ2						        ; GROUP INTERCOMPARTMENTAL CLEARANCE 2 WT
    Q2     = GPQ2			                   	; INDIVIDUAL INTERCOMPARTMENTAL CLEARANCE 2
	
    TVV3   = THETA(8)                       	; POPULATION VOLUME 2
    GPV3   = TVV3						        ; GROUP VOLUME 2 
    V3     = GPV3			                   	; INDIVIDUAL VOLUME 2
	
	S1 = V1/1000
	
$DES
	
	C1=A(1)/V1
	C2=A(2)/V2
	C3=A(3)/V3
	CL=(VMEL/(KMEL+C1))
	Q1A=(VMTR/(HMTR+C1)) ;DIST CL FROM CENTRAL TO V2
	Q1B=(VMTR/(HMTR+C2)) ;DIST CL FROM V2 TO CENTRAL
	DADT(1)=Q1*(C2-C1) + Q2*(C3-C1)-(CL*C1)
	DADT(2)=(Q1A*C1) - (Q1B*C2)
	DADT(3)=Q2*(C1-C3)

      
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
8070				            ; VMEL (UG/H)
40.2							; KMEL (UG/L) 
15.8							; CENTRAL VOLUME  (L)
325000							; VMTR (UG/H)
4260 							; KMTR (UG/L)
1650							; PERIPHERAL VOLUME 2 (L)
41.6							; DISCTRIBUTION CLEARANCE (L/H) 					
75.4							; PERIPHERAL VOLUME 3 (L)
0.554							; ALBUMIN ON VMEL


$OMEGA
0.0497                     		; [E] BETWEEN SUBJECT VARIABILITY (VMEL)
0.2181                   		; [E] BETWEEN SUBJECT VARIABILITY (V1)
0.0713                   	    ; [E] BETWEEN SUBJECT VARIABILITY (VMTR)
0.065                     		; [E] BETWEEN SUBJECT VARIABILITY (V2)					

$SIGMA
0.06125                     	; [P] PROPORTIONAL RESIDUAL ERROR

;************************************************************
;****        ESTIMATION AND SIMULATION ROUTINES          ****
;************************************************************

$SIM ONLYSIM (4567182)

$TABLE                                       ; STANDARD TABLE
       ID TIME DOSE DUR RATE II AMT DV CMT EVID MDV DV PRED MDV IWRE IRES
       IPRE KMTR KMEL VMTR VMEL V1 V2 Q3 V3 ETA1 ETA2 ETA3 ETA4 ALB
       NOPRINT ONEHEADER FILE=sdtab