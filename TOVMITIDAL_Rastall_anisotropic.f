c ╭────── · · ✧ · · ──────╮
      PROGRAM TOVSOLVER     
c ╰────── · · ✧ · · ──────╯

c           ˚₊  calculates moment of inertia & tidal deformability 
c           ˚₊  EoS: Chaplygin dark star with anisotropy (Horvat et al.)
c           ˚₊  framework: Rastall gravity (modified as effective matter)

c           ˚₊  python version available
c           ˚₊  github.com/jessiechd
c           ˚₊  last updated: 2025
     

      IMPLICIT DOUBLE PRECISION (A-H,K,L,M,O-Z) 
      INTEGER IL
      CHARACTER(LEN=100) :: fn, fn2

c---  !change anisotropy, rastall, PCC parameters here: 
      alp = 0.0                 ! anisotropy factor (isotropic = 0)`T
      rast = 0.05               ! Rastall parameter (GR = 0)
      PCCmax = 500              ! max PCC to stop printing
      PCCprof = 101.D0          ! PCC point to save star profile data 

      WRITE(fn,'(A,F3.1,A,F4.2,A)') 'file_a=',alp,"_b=",rast,'.txt'
      PRINT *, fn
      OPEN (unit=8,STATUS='unknown',file=fn)

      WRITE(fn2,'(A,F3.1,A,F4.2,A)') 'profile_a=',alp,"_b=",rast,'.txt'
      PRINT *, fn2
      OPEN (unit=9,STATUS='unknown',file=fn2)



      ALPHA = alp * 1.D0    
      RAS = rast * 1.D0    

      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15

c---  begin iterations from PCC=1 to PCCmax
      DO 10 IL=1,PCCmax,2     
      PCC=1.0D0*IL                   ! pressure in center
      TUCT=0.1D-8                    ! trail initial condition for matric component nu

c---  call subroutine that computes initial condition
      CALL TOVI(PCC,TUCT,ROS,GMOS,MNURT,ALPHA,RAS)

      RNS=ROS*1.D3                    ! ROS = radius (km)
      DBS=1.D0-2.D0*GS*GMOS*MSS/RNS

      KC=0.5D0*DLOG(DBS)-MNURT
      ROT= TUCT + KC                  ! ROT = correct initial condition for matric component nu

      PMIN=1.0D-9

c---  call subroutine that computes momin & tidal deformability variables
      CALL TOVMI(PCC,ROT,ROS2,GMOS2,MNURT,OMEGA,KAPPA,ALPHA,RAS,YR2,XP)
      
      
c---  process k2, and get lambda parameters for tidal deformability
      MASST = GMOS2
      k2 = k_2(XP,MASST,YR2) 

      C=GS*MASST*MSS/XP
      C2=C*C
      C3=C2*C
      C4=C3*C
      C5=C4*C
      LMBD=2.D0*K2/(3.D0*C5) 

      Const=1.D-36
      XP5=XP*XP*XP*XP*XP*1D10
      GB=6.674D-8
      lambda=2.D0*K2*XP5*Const/(3.D0*GB) 

      RNS2=ROS2*1.D3

c---  get moment of inertia
      MOMIN=KAPPA/(GS*GMOS*MSS*RNS2*RNS2)               ! I/MR2 (dimensionless)
      MI=MOMIN*1.98892D33*1.D10*GMOS*ROS2*ROS2/1.0D45   ! moment of inertia (10^45 g.cm2)

      COMP = GMOS2/ROS2                                 ! compactness (M/R)
      ED = FED(PCC)                                     ! energy density (MeV/fm3)

      GMR = 1-2.D0*GMOS2*1.48D0/ROS2
      ZSUR = 1/SQRT(GMR) - 1                            ! surface gravitational redshift

      ANISO = alpha*PCC*2.D0*GS*GMOS2*MSS/ROS2         ! anisotropy factor  (MeV/fm3)
      PT = PCC + ANISO                                 ! tangential pressure (MeV/fm3)

      Pr_eff = P_eff(PCC,ED,RAS,ANISO)              
      ED_eff = E_eff(PCC,ED,RAS,ANISO)
      Pt_eff = Pr_eff - PCC + PT
         
      DEDPr = DEDP(PCC)
      DEDPt = DEDP(PT)
      V2r = 1.D0/DEDPr                         ! speed of sound (radial)
      V2t = 1.D0/DEDPt                         ! speed of sound (tangential)

      DEDPr_eff = DEDP_eff(RAS,DSDP,DEDPr)
      DEDPt_eff = DEDP_eff(RAS,DSDP,DEDPt)
      V2r_eff = 1.D0/DEDPr_eff                 ! speed of sound (radial); effective
      V2t_eff = 1.D0/DEDPt_eff                 ! speed of sound (tangential); effective

      EC1 = ED_eff - Pr_eff - 2.D0*Pt_eff       ! energy condition: EC1>0
      EC2 = ED_eff + Pr_eff + 2.D0*Pt_eff       ! energy condition: EC2>0

      IF (EC1 <= 0.D0) THEN
            WRITE(*,*) "EC1 out of bounds. Stopping program."
            WRITE(*,*) "ED - Pr - 2Pt > 0 not fulfilled."
      END IF

      IF (EC2 <= 0.D0) THEN
            WRITE(*,*) "EC2 out of bounds. Stopping program."
            WRITE(*,*) "ED + Pr + 2Pt > 0 not fulfilled."
      END IF
      
      WRITE(8,*)IL,PCC,ROS2,GMOS2,PT,ED,LOG(ED),COMP,ZSUR,Pr_eff,Pt_eff,
     &          ED_eff,V2r_eff,V2t_eff,EC1,EC2,MOMIN,MI,k2,LMBD,lambda
            

      WRITE(*,*) "Computing: PCC =", IL, "/", PCCmax, "(MeV/fm3)"
      
      
 10   CONTINUE
       
      STOP
      END 


c--- ₊˚ ‿︵‿︵‿︵‿︵‿︵‿︵୨୧ · · ♡ · · ୨୧‿︵‿︵‿︵‿︵‿︵‿︵ ˚₊ ---

c---  subroutine TOVMI: compute momen inersia & tidal deformability -----
     
      SUBROUTINE TOVMI(PCC,TUCT,ROS,GMOS,MNURT,OMEGA,KAPPA,A,R,YR,XP)
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IN    
      DIMENSION YA(10), EK(5,10), Y(10)

      HC  = 197.327D0
      ALPHA = A 
      RAS = R 

      IM=7              ! number of equations
      IN=IM-1           ! number of 1st order differential equations

      YL=2.D0 

      Y(1)=PCC          ! pressure 
      Y(2)=0.1D-8       ! star mass ˚₊‧꒰ა ☆ ໒꒱ ‧₊˚
      Y(3)=TUCT         ! matric nu 
      Y(4)=0.1D-8       ! omega 
      Y(5)=0.1D-8       ! kappa 
      Y(6)=YL           ! y tidal
      Y(7)=FED(PCC)     ! energy density

      PU=1.0D-1         ! interval of radius for printing
      NS=32             ! number of steps in one print
      XL=30.0D3         ! maximum radius to stop calculation

      H=PU/NS           
      XP=1.0D-3
      HH=H/(2.0D0)

      LI=0              ! line number initialization ⊹₊｡ꕤ˚₊⊹

      IF (PCC == 101.D0) THEN
            WRITE(*,*) "Computing star profile (PCC = 101 MeV/fm3)."
      END IF
 
 28   LI=LI+1           
      
      DO N=1,NS
         XB=XP          ! old radius (k1) 
         XP=XP+H        ! new radius (k4)
         XM=XB+HH       ! midpoint radius (k2, k3)
 
C---     compute K1, L1, M1, N1, O1
         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO

         XA=XB
         CALL FUNCX(EK,J,YA,XA,H,YL,ALPHA,RAS)

C---     compute K2, L2, M2, N2, O2
         J=2
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
         P0=YA(1)
         ED=FED(P0)
         YA(7)=ED
         
         XA=XM
         CALL FUNCX(EK,J,YA,XA,H,YL,ALPHA,RAS)
        
C---     compute K3, L3, M3, N3, O3
         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
         P0=YA(1)
         ED=FED(P0)
         YA(7)=ED    
         
         XA=XM
         CALL FUNCX(EK,J,YA,XA,H,YL,ALPHA,RAS)

C---     compute K4, L4, M4, N4, O4
         J=4
         DO I=1,IN
            YA(I)=Y(I)+EK(3,I)
         END DO

         P0=YA(1)
         ED=FED(P0)
         YA(7)=ED
         
         XA=XP
         CALL FUNCX(EK,J,YA,XA,H,YL,ALPHA,RAS)

C---     4th order RK scheme ✩₊˚.⋆☾⋆⁺₊✧
         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
         P0=Y(1)
         ED=FED(P0) 
         Y(7)=ED
          
      END DO

      PS=Y(1)
      ROS=(XP/1.D3)
      GMOS=Y(2)
      MNURT=Y(3)
      OMET=Y(4)
      KAPAT=Y(5)
      YR=Y(6)

      PMIN=0.7D-8

      IF (PCC == PCCprof) THEN
            WRITE(9,*) PS, ROS, GMOS, MNURT, ED
      END IF

      IF (PS .GT. PMIN) GOTO 28  

      IF (PCC == PCCprof) THEN
            WRITE(*,*) "Star profile saved."
      END IF   

c---  process omega & kappa for moment of inertia ૮꒰ ˶• ༝ •˶꒱ა ♡
      ETA=1.D0/(OMET+2.D0*KAPAT/(XP*XP*XP))
      OMEGA=ETA*OMET
      KAPPA=ETA*KAPAT



      RETURN
      END
c------------------------ end of subroutine TOVMI -----------------------

c--- ₊˚ ‿︵‿︵‿︵‿︵‿︵‿︵୨୧ · · ♡ · · ୨୧‿︵‿︵‿︵‿︵‿︵‿︵ ˚₊ ---

c---  subroutine TOVI: compute TOV initial condition for nu  ------------
      SUBROUTINE TOVI(PCC,TUCT,ROS,GMOS,MNURT,ALPHA,RAS)
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IN    
      DIMENSION YA(10), EK(5,10), Y(10)

      HC  = 197.327D0

      IM=4
      IN=IM-1  

      Y(1)=PCC          ! pressure
      Y(2)=0.1D-8       ! star mass 
      Y(3)=TUCT         ! matric nu
      Y(4)=FED(PCC)     ! energy density
   
   
      PU=1.0D-1
      NS=16
      XL=20.0D3

      H=PU/NS
      XP=1.0D-3
      HH=H/(2.0D0)

      LI=0              ! line number initialization

 28   LI=LI+1

      DO N=1,NS
         XB=XP          ! old radius (k1) 
         XP=XP+H        ! new radius (k4)
         XM=XB+HH       ! midpoint radius (k2, k3)

C---     compute K1, L1, M1
         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO

         XA=XB
         CALL FUNCT(EK,J,YA,XA,H,ALPHA,RAS)

C---     compute K2, L2, M2
         J=2
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
         P0=YA(1)
         ED=FED(P0)
         YA(4)=ED
    
         XA=XM
         CALL FUNCT(EK,J,YA,XA,H,ALPHA,RAS)

C---     compute K3, L3, M3
         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
         P0=YA(1)
         ED=FED(P0)
         YA(4)=ED 
        
         XA=XM
         CALL FUNCT(EK,J,YA,XA,H,ALPHA,RAS)

C---     compute K4, L4, M4
         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO
         P0=YA(1)
         ED=FED(P0)
         YA(4)=ED

         XA=XP
         CALL FUNCT(EK,J,YA,XA,H,ALPHA,RAS)

C---     4th order RK scheme ✩₊˚.⋆☾⋆⁺₊✧
         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
         P0=Y(1)
         ED=FED(P0)
         Y(4)=ED 
      
       END DO

       PS=Y(1)
       PMIN=1.0D-9

      IF (PS .GT. PMIN  ) GOTO 28

      ROS=(XP/1.D3)
      GMOS=Y(2)
      MNURT=Y(3)
      
      RETURN
      END
c------------------------ end of subroutine TOVI ------------------------

c--- ₊˚ ‿︵‿︵‿︵‿︵‿︵‿︵୨୧ · · ♡ · · ୨୧‿︵‿︵‿︵‿︵‿︵‿︵ ˚₊ ---

c---  subroutine FUNCX: differential equations for momin, tidal  --------
      SUBROUTINE FUNCX(EK,J,YA,XA,H,YL,ALPHA,RAS)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(5,10), YA(10)
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15

      PRESS=YA(1)
      MASST=YA(2)
      MNU=YA(3)
      OMGT=YA(4)
      KAPT=YA(5)
      YTIDAL=YA(6)
      EDEN=YA(7) 
      EXPMNU=EXP(MNU)
      EXPMMNU=EXP(-MNU)

      ANI = ALPHA*PRESS*2.D0*GS*MASST*MSS/XA      ! anisotropy σ = Pt-Pr

      PRESS_EFF = P_eff(PRESS,EDEN,RAS,ANI)  
      EDEN_EFF = E_eff(PRESS,EDEN,RAS,ANI)

      DEDP1=DEDP(PRESS)                           ! dE/dP
      DSDP=ALPHA*2.D0*GS*MASST*MSS/XA             ! dσ/dP


      DPDR=-GS*EDEN*MASST*MSS*H/(XA*XA)
     &        *(1.D0+PRESS/EDEN)
     &       *(1.D0+4.D0*PI*XA*XA*XA*PRESS_EFF/(MSS*MASST))
     &       /(1.D0-2.D0*GS*MASST*MSS/XA)
     &       + 2*H*ANI/XA

      Y = 1.D0 - 3.D0*RAS - 2.D0*RAS*DSDP + RAS*DEDP1

      A=ALPHA*2.D0*GS*MASST*MSS/XA + 1.D0         ! A = dPt/dPr

      SMGM=(1.0D0-2.0D0*GS*MASST*MSS/XA)
      SRP=(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MSS*MASST))
      DQYT=SRP*SRP/(SMGM*SMGM)
      FYT=(1-4.D0*PI*XA*XA*GS*(EDEN_EFF - PRESS_EFF))/SMGM

      DEDPEFF = DEDP_eff(RAS,DSDP,DEDP1)
      SUKU=4.D0*EDEN_EFF+8.D0*PRESS_EFF+((DEDPEFF+1.D0)*(EDEN+PRESS)/A)

      QYT=4.D0*PI*XA*XA*GS*SUKU
      QYT=(QYT-YL*(YL+1.D0))/SMGM
      QYT=QYT-4.D0*(GS*MASST*MSS)*(GS*MASST*MSS)*DQYT/(XA*XA)


      EK(J,1) = DPDR/Y                                         ! dP/dr

      EK(J,2)=4.D0*PI*XA*XA*EDEN_EFF*H/MSS                     ! dm/dr
 
      EK(J,3)=GS*MASST*MSS*H/(XA*XA)                           ! dν/dr
     &       *(1.D0+4.D0*PI*XA*XA*XA*PRESS_EFF/(MSS*MASST))
     &       /(1.D0-2.D0*GS*MASST*MSS/XA)

      EK(J,4)=6.D0*EXPMNU*H/DSQRT(1.D0-2.D0*GS*MASST*MSS/XA)   ! dω/dr
     &     *KAPT*(1-ANI/(EDEN+PRESS))/(XA*XA*XA*XA)

      EK(J,5)=GS*8.D0*PI*(XA*XA*XA*XA)*EXPMMNU*EDEN*H/3.D0     ! dκ/dr
     &        *(1.D0+PRESS/EDEN)*OMGT
     &        /DSQRT(1.D0-2.D0*GS*MASST*MSS/XA)

      EK(J,6)=-(YTIDAL*YTIDAL+YTIDAL*FYT+QYT)*H/XA             ! dy/dr
 
      RETURN
      END
c------------------------ end of subroutine FUNCX -----------------------

c--- ₊˚ ‿︵‿︵‿︵‿︵‿︵‿︵୨୧ · · ♡ · · ୨୧‿︵‿︵‿︵‿︵‿︵‿︵ ˚₊ ---

c---  subroutine FUNCT: TOV equations for initial condition for nu  -----
      SUBROUTINE FUNCT(EK,J,YA,XA,H,ALPHA,RAS)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(5,10), YA(10)
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
  
      EDEN=YA(4)  
      PRESS=YA(1)
      MASST=YA(2)
      MNU=YA(3)

      ANI = ALPHA*PRESS*2.D0*GS*MASST*MSS/XA        ! anisotropy σ = Pt-Pr

      PRESS_EFF = P_eff(PRESS,EDEN,RAS,ANI)
      EDEN_EFF = E_eff(PRESS,EDEN,RAS,ANI)

      DEDP1=DEDP(PRESS)                             ! dE/dP
      DSDP=alpha*2.D0*GS*MASST*MSS/XA               ! dσ/dP

      DPDR=-GS*EDEN*MASST*MSS*H/(XA*XA)
     &        *(1.D0+PRESS/EDEN)
     &       *(1.D0+4.D0*PI*XA*XA*XA*PRESS_EFF/(MSS*MASST))
     &       /(1.D0-2.D0*GS*MASST*MSS/XA)
     &       + 2*H*ANI/XA

      Y = 1.D0 - 3.D0*RAS - 2.D0*RAS*DSDP + RAS*DEDP1

      EK(J,1) = DPDR/Y                                         ! dP/dr

      EK(J,2)=4.D0*PI*XA*XA*EDEN_EFF*H/MSS                     ! dm/dr
 
      EK(J,3)=GS*MASST*MSS*H/(XA*XA)                           ! dν/dr
     &       *(1.D0+4.D0*PI*XA*XA*XA*PRESS_EFF/(MSS*MASST))
     &       /(1.D0-2.D0*GS*MASST*MSS/XA)

      RETURN
      END
c------------------------ end of subroutine FUNCT -----------------------

c--- ₊˚ ‿︵‿︵‿︵‿︵‿︵‿︵୨୧ · · ♡ · · ୨୧‿︵‿︵‿︵‿︵‿︵‿︵ ˚₊ ---

c---  function FED: energy density as a function of pressure (DS CG)  ---
      FUNCTION FED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      REAL :: A, B
      conv = 7.5534D11              ! unit conversion to MeV/fm3

c---  !change Chaplygin dark star parameters here
      A = 0.3D0                     ! dark star CG parameter A
      B = 6.0D-20 * (conv**2)       ! dark star CG parameter B
    
      FED = (P0 + SQRT(P0**2 + 4*A*B))/(A*2)

      RETURN
      END
c-------------------------- end of function FED -------------------------

c--- ₊˚ ‿︵‿︵‿︵‿︵‿︵‿︵୨୧ · · ♡ · · ୨୧‿︵‿︵‿︵‿︵‿︵‿︵ ˚₊ ---

c---  function DEDP: derivative/turunan FED terhadap pressure  ----------
      FUNCTION DEDP(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL FED

      h  = 1.D-5                    ! note: xa = P
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h 
      C1=FED(x4)
      C2=FED(x2)
   
      DEDP = (FED(x4)-8.D0*FED(x2)+8.D0*FED(x1)-FED(x3))/(12.D0*h)
      RETURN
      END
c------------------------- end of function DEDP -------------------------

c--- ₊˚ ‿︵‿︵‿︵‿︵‿︵‿︵୨୧ · · ♡ · · ୨୧‿︵‿︵‿︵‿︵‿︵‿︵ ˚₊ ---

c---  function DEDP_eff: turunan FED_eff terhadap pressure_eff  ---------
      FUNCTION DEDP_eff(RAS,S,DEDP1)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      
      Y=3.D0 + 2.D0*S - DEDP1
      num = DEDP1 + RAS*Y
      deno = 1.D0 - RAS*Y

      DEDP_eff = num/deno
      RETURN
      END
c----------------------- end of function DEDP_eff -----------------------

c--- ₊˚ ‿︵‿︵‿︵‿︵‿︵‿︵୨୧ · · ♡ · · ୨୧‿︵‿︵‿︵‿︵‿︵‿︵ ˚₊ ---

c---  function P_eff: = Pressure with Ras. parameter & anisotropy  ------
      FUNCTION P_eff(PRESS,EDEN,RAS,ANI)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)

      T = 3.D0*PRESS + 2.D0*ANI - EDEN
      P_eff = PRESS - RAS*T

      RETURN
      END
c------------------------- end of function P_eff ------------------------

c--- ₊˚ ‿︵‿︵‿︵‿︵‿︵‿︵୨୧ · · ♡ · · ୨୧‿︵‿︵‿︵‿︵‿︵‿︵ ˚₊ ---

c---  function E_eff: = energy den. with Ras. parameter & anisotropy  ---
      FUNCTION E_eff(PRESS,EDEN,RAS,ANI)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)

      T = 3.D0*PRESS + 2.D0*ANI - EDEN
      E_eff = EDEN + RAS*T

      RETURN
      END
c------------------------- end of function E_eff ------------------------

c--- ₊˚ ‿︵‿︵‿︵‿︵‿︵‿︵୨୧ · · ♡ · · ୨୧‿︵‿︵‿︵‿︵‿︵‿︵ ˚₊ ---

C      ================================================================
C      ELECTRIC TIDAL LOVE NUMBER
C      ===============================================================


      FUNCTION k_2(XP,MASST,YR2) 
C     *********************************************************
C     DEFINES love number K2 and surcial love number H2
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 

      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15

      C=GS*MASST*MSS/XP
      C2=C*C
      C3=C2*C
      C4=C3*C
      C5=C4*C
      
      FKUR1=(4.D0*(YR2+1.D0)*C4 +(6.D0*YR2-4.D0)*C3+(26.D0-22.D0*YR2)*C2
     &     + 3.D0*(5.D0*YR2-8.D0)*C -3.D0*YR2+6.D0)*2.D0*C
      
      FKUR2=-3.D0*(1.D0-2.D0*C)*(1.D0-2.D0*C)*(2.D0*C*(YR2-1.D0)-YR2+2)
     &      *DLOG(1.D0/(1.D0-2.D0*C))

      FKUR=FKUR1+FKUR2

      K2 = (8.D0/5.D0)*(1.D0-2.D0*C)*(1.D0-2.D0*C)*C5
      K2 = K2*(2.D0*C*(YR2-1.D0)-YR2+2.D0)/FKUR
      k_2 = K2 * 1.D0
      RETURN
      END