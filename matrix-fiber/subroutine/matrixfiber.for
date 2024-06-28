C ======================================================================
C User subroutine VUEL for Abaqus:
C designed for simulating the fracture of composite materials by a double-phased field
C Type 1: C3D16 phase-field tetrahedral element
C Type 2: C3D32 phase-field brick element
C Type 3: C3D8R reduced 
C the freedom degree in thermal-displacement is x,y,z,T.
C for the hexahedron elements,  the 2 × 2 × 2 integral or 1 × 1 × 1 integral scheme is used here.
C Jing Lee, 2024
C ======================================================================
C PROPS(1) = Density (rho)
C PROPS(2) = Longitudinal Young’s modulus (E1) 
C PROPS(3) = Transverse Young’s modulus (E2,E3)
C PROPS(4) = Shear modulus (G12,G13)
C PROPS(5) = Shear modulus (G23)
C PROPS(6) = Poisson’s ratio (V12,V13)
C PROPS(7) = Poisson’s ratio (V23)
C PROPS(8) = Transverse tensile strength (ST1) 
C PROPS(9) = Transverse shear strength (ST2)
C PROPS(10) = In-plane shear strength (STL) 
C PROPS(11) = Longitudinal tensile strength(YL)
C PROPS(12) = Transverse tensile critical energy release rate(GcT1)
C PROPS(13) = Transverse shear critical energy release rate (GcT2)
C PROPS(14) = In-plane shear critical energy release rate (GcTL)
C PROPS(15) = Longitudinal tensile critical energy release rate(GcL)
C PROPS(16) = matrix phase field internal length scales (l0m)
C PROPS(17) = fiber phase field internal length scales(l0L)
C PROPS(18) = artificial viscosity parameters for matrix(etam)
C PROPS(19) = artificial viscosity parameters for fiber(etaL)
C PROPS(20) = Material direction vector component (nx)
C PROPS(21) = Material direction vector component (ny)
C PROPS(22) = Material direction vector component (nz)
C
C ---- Used variables ---------------------
C N_ELEM - number of elements used in the model divided
C            by 3 - (N_stress_couple_phase_matrix+N_VUMATHT+N_phase_fiber)/3 (to be changed for each model)           
C NSTV,NSTV1,NSTV2 - overall solution dependent variables (displacements,strains,stresses,energies, phase,history)

C ======================================================================
      SUBROUTINE VUEL(nblock,rhs,amass,dtimeStable,svars,nsvars,
     1                energy,
     2                nnode,ndofel,props,nprops,jprops,njprops,
     3                coords,ncrd,u,du,v,a,
     4                jtype,jElem,
     5                time,period,dtimeCur,dtimePrev,kstep,kinc,
     6                lflags,
     7                dMassScaleFactor,
     8                predef,npredef,
     9                jdltyp, adlmag)
C The following is the standard definition of the abaqus interface. don't make any changes!
      include 'vaba_param.inc'
    
      

C     operational code keys
      parameter ( jMassCalc            = 1,
     *            jIntForceAndDtStable = 2,
     *            jExternForce         = 3)

C     flag indices
      parameter (iProcedure = 1,
     *           iNlgeom    = 2,
     *           iOpCode    = 3,
     *           nFlags     = 3)

C     energy array indices
      parameter ( iElPd = 1,
     *            iElCd = 2,
     *            iElIe = 3,
     *            iElTs = 4,
     *            iElDd = 5,
     *            iElBv = 6,
     *            iElDe = 7,
     *            iElHe = 8,
     *            iElKe = 9,
     *            iElTh = 10,
     *            iElDmd = 11,
     *            iElDc = 12,
     *            nElEnergy = 12)

C     predefined variables indices
      parameter ( iPredValueNew = 1,
     *            iPredValueOld = 2,
     *            nPred         = 2)

C     time indices
      parameter (iStepTime  = 1,
     *           iTotalTime = 2,
     *           nTime      = 2)
c     procedure flags
      parameter ( jDynExplicit = 17 )
      parameter ( jCoupTS = 74 )

      dimension rhs(nblock,ndofel),amass(nblock,ndofel,ndofel),
     1          dtimeStable(nblock),
     2          svars(nblock,nsvars),energy(nblock,nElEnergy),
     3          props(nprops),jprops(njprops),
     4          jelem(nblock),time(nTime),lflags(nFlags),
     5          coords(nblock,nnode,ncrd),
     6          u(nblock,ndofel), du(nblock,ndofel),
     7          v(nblock,ndofel), a(nblock, ndofel),
     8          dMassScaleFactor(nblock),  
     9          predef(nblock,nnode,npredef,nPred),
     *          adlmag(nblock)
C   Declare custom variables here

C     ==================================================================
      parameter(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,THREE=3.D0,
     1   TOLER=1.0D-8,FOUR=4.D0,EIGHT=8.d0,RP25=0.25D0,HALF=0.5D0,
     2   SIX=6.D0,N_ELEM=1,NSTV1=24,NSTV2=2,NSTV=NSTV1+NSTV2)
C    c0 = integral(@(d) 4* sqrt(2*d-d.*d), 0, 1); normalization constant
      real(8),parameter:: a2=-0.5D0,p=2.D0,c0=3.14D0,beita=35
      
      INTEGER I,J,L,K,K1,K2,K3,K4,IX,IY,IZ,INPT,Innode,INODE,Index

      REAL*8 AINTW(8),XII(8,3),XI(3),dNdxi(nnode,3),
     1   VJACOB(3,3),dNdx(nnode,3),VJABOBINV(3,3),AN(nnode),
     2   BP(3,nnode),DP(3),sdv(NSTV),BB(6,nnode*3),VNI(3,nnode*3),
     3   ULOC(3),strain(6),stress(6),myenergy(6),dir_f(3)

      REAL*8 DTM
      REAL*8 rho,E1,E2,G12,G23,V12,V23,ST1,
     1     ST2,STL,SYL,GcT1,GcT2,GcTL,GcL,l0m,
     2     l0L,etaM,etaL,nx,ny,nz,nnorm,ProjDp

      REAL*8 wT1d,wT2d,wTLd,WLd,wT1d_dd,wT2d_dd,wTLd_dd,wLd_dd
      REAL*8 PhiT1,PhiT2,PhiTL,PhiL,PhiT1_H,PhiT2_H,PhiTL_H,PhiL_H
      REAL*8 phase,phase_m,phase_l
      REAL*8 C11,C22,C44,C55,C12,C23
      REAL*8 a1_T1,a1_T2,a1_TL
      COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
      logical :: mydebug = .true.
      logical :: enable_c3d8r = .true.
      real*8 :: print_time=1E-5,next_print_time=0.d0
C      ======================for multithread---------------------------
      integer numThreads,NUMPROCESSES,KPROCESSNUM,myThreadID     
      CALL VGETNUMCPUS( NUMPROCESSES )
      CALL VGETRANK( KPROCESSNUM )
      numThreads = GETNUMTHREADS()
      myThreadID = get_thread_id()
c     to be completed
C     =============================

c     real*8 USRVAR_ARRAY(*)
c      pointer(ptr_USRVAR,USRVAR_ARRAY)
c      logical :: allocated = .false.
c      save ptr_USRVAR, allocated

c      if (.not. allocated) then
c        ptr_USRVAR = SMAFloatArrayCreate(1, N_ELEM*NSTV*8, zero)
c        allocated = .true.
c      endif
       call MutexInit( 1 )
c      logical :: firstrun = .false.
c      integer tempvar
c      if (firstrun) then
c        write(*,*) "write a number"
c        read(*,*) tempvar
c        firstrun = .false.
c      endif
C     ==================================================================
C     Material parameters
C     ==================================================================
      rho=props(1)
      E1=props(2)
      E2=props(3)
      G12=props(4)
      G23=props(5)
      V12=props(6)
      V23=props(7)
      ST1=props(8)
      ST2=props(9)
      STL=props(10)
      SYL=props(11)
      GcT1=props(12)
      GcT2=props(13)
      GcTL=props(14)
      GcL=props(15)
      l0m=props(16)
      l0L=props(17)
      etaM=props(18)
      etaL=props(19)  
      nx=props(20)
      ny=props(21)
      nZ=props(22)
      nnorm = sqrt((nx**two+ny**two+nz**two))
      nx = nx/nnorm
      ny = ny/nnorm
      nz = nz/nnorm
      dir_f(1) = nx
      dir_f(2) = ny
      dir_f(3) = nz
      call Cmatrix(C11,C22,C44,C55,C12,C23,E1,E2,G12,G23,V12,V23)
      do kblock = 1, nblock
        dtimeStable(kblock) = 1e-8;
      enddo
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
      if (jtype.EQ.ONE) then
        XII(1,1) = ONE/FOUR
        XII(1,2) = ONE/FOUR
        XII(1,3) = ONE/FOUR
        Innode=1
        AINTW(1) = ONE/SIX
      elseif (jtype.EQ.TWO) then
        if (enable_c3d8r) then
         XII(1,1) = zero
         XII(1,2) = zero
         XII(1,3) = zero
         Innode = 1.0
         do I=1,Innode
          AINTW(I) = EIGHT
         enddo
        else 
         XII(1,1) = MONE/THREE**HALF
         XII(1,2) = MONE/THREE**HALF
         XII(1,3) = MONE/THREE**HALF
         XII(2,1) = ONE/THREE**HALF
         XII(2,2) = MONE/THREE**HALF
         XII(2,3) = MONE/THREE**HALF
         XII(3,1) = ONE/THREE**HALF
         XII(3,2) = ONE/THREE**HALF
         XII(3,3) = MONE/THREE**HALF
         XII(4,1) = MONE/THREE**HALF
         XII(4,2) = ONE/THREE**HALF
         XII(4,3) = MONE/THREE**HALF
         XII(5,1) = MONE/THREE**HALF
         XII(5,2) = MONE/THREE**HALF
         XII(5,3) = ONE/THREE**HALF
         XII(6,1) = ONE/THREE**HALF
         XII(6,2) = MONE/THREE**HALF
         XII(6,3) = ONE/THREE**HALF
         XII(7,1) = ONE/THREE**HALF
         XII(7,2) = ONE/THREE**HALF
         XII(7,3) = ONE/THREE**HALF
         XII(8,1) = MONE/THREE**HALF
         XII(8,2) = ONE/THREE**HALF
         XII(8,3) = ONE/THREE**HALF
         Innode = 8.0
         do I=1,Innode
            AINTW(I) = ONE
         enddo
        endif
       endif
C     ==================================================================
C     Calculating properties at each element
C     ==================================================================    
C     ==================================================================
C     Calculating for the displacement and phase_matrix
C     ==================================================================   
      do kblock = 1, nblock
      if (jelem(kblock).le.N_ELEM) then
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
        do INPT=1,Innode
C     Initializing solution dependent variables (phase,history)
            do I=1,NSTV1
             sdv(I)=svars(kblock,NSTV1*(INPT-1)+I)
            enddo
C     Local coordinates of the integration point
            XI(1) = XII(INPT,1)
            XI(2) = XII(INPT,2) 
            XI(3) = XII(INPT,3) 
            
C     Shape functions and local derivatives
            if (jtype.EQ.ONE) then
             call SHAPEFUNT(AN,dNdxi,XI)
            elseif (jtype.EQ.TWO) then             
             call SHAPEFUN(AN,dNdxi,XI)
            endif
C   Shape functions for vector
            IZ=ZERO
            DO I = 1,NNODE
             IX=IZ+1
             IY=IX+1
             IZ=IY+1
             VNI(1,IX)=AN(I)
             VNI(2,IX)=ZERO
             VNI(3,IX)=ZERO
             VNI(1,IY)=ZERO
             VNI(2,IY)=AN(I)
             VNI(3,IY)=ZERO
             VNI(1,IZ)=ZERO
             VNI(2,IZ)=ZERO
             VNI(3,IZ)=AN(I)
            END DO
C     Jacobian
C     Note that ncrd and nnode are in a different order than uel
            do I = 1,3
             do J = 1,3
              VJACOB(I,J) = ZERO
              do K = 1,nnode
               VJACOB(I,J) = VJACOB(I,J)+coords(kblock,K,I)*dNdxi(K,J) 
     1                      
              enddo
             enddo
            enddo       
            DTM = ZERO
            DTM = VJACOB(1,1)*VJACOB(2,2)*VJACOB(3,3)+VJACOB(1,2)*
     1       VJACOB(2,3)*VJACOB(3,1)+VJACOB(1,3)*VJACOB(2,1)*
     2       VJACOB(3,2)-VJACOB(3,1)*VJACOB(2,2)*VJACOB(1,3)-
     3       VJACOB(3,2)*VJACOB(2,3)*VJACOB(1,1)-VJACOB(3,3)*
     4       VJACOB(2,1)*VJACOB(1,2)
C     
            if (DTM.lt.ZERO) then
             WRITE(7,*) 'Negative Jacobian',DTM
             call XPLB_EXIT
            endif

C     Inverse of Jacobian
            VJABOBINV(1,1)=(VJACOB(2,2)*VJACOB(3,3)-VJACOB(2,3)*
     1       VJACOB(3,2))/DTM
            VJABOBINV(1,2)=-(VJACOB(1,2)*VJACOB(3,3)-VJACOB(3,2)*
     1       VJACOB(1,3))/DTM
            VJABOBINV(1,3)=(VJACOB(1,2)*VJACOB(2,3)-VJACOB(1,3)*
     1       VJACOB(2,2))/DTM
            VJABOBINV(2,1)=-(VJACOB(2,1)*VJACOB(3,3)-VJACOB(2,3)*
     1       VJACOB(3,1))/DTM
            VJABOBINV(2,2)=(VJACOB(1,1)*VJACOB(3,3)-VJACOB(1,3)*
     1       VJACOB(3,1))/DTM
            VJABOBINV(2,3)=-(VJACOB(1,1)*VJACOB(2,3)-VJACOB(1,3)*
     1       VJACOB(2,1))/DTM
            VJABOBINV(3,1)=(VJACOB(2,1)*VJACOB(3,2)-VJACOB(2,2)*
     1       VJACOB(3,1))/DTM
            VJABOBINV(3,2)=-(VJACOB(1,1)*VJACOB(3,2)-VJACOB(1,2)*
     1       VJACOB(3,1))/DTM
            VJABOBINV(3,3)=(VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*
     1       VJACOB(2,1))/DTM       
C     Derivatives of shape functions respect to global ccordinates
            do K = 1,nnode
             do I = 1,3
                dNdx(K,I) = ZERO
              do J = 1,3
               dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)                         
              enddo
             enddo
            enddo
C     Calculating B matrix (B=LN) for phase
            do INODE=1,nnode
             BP(1,INODE)=dNdx(INODE,1)
             BP(2,INODE)=dNdx(INODE,2)
             BP(3,INODE)=dNdx(INODE,3)
            enddo
C     Calculating B matrix (B=LN) for displacement
            IZ=ZERO
            do INODE=1,NNODE
                IX=IZ+1
                IY=IX+1
                IZ=IY+1
                BB(1,IX)= dNdx(INODE,1)
                BB(2,IX)= ZERO
                BB(3,IX)= ZERO
                BB(4,IX)= dNdx(INODE,2)
                BB(5,IX)= dNdx(INODE,3)
                BB(6,IX)= ZERO
                BB(1,IY)= ZERO
                BB(2,IY)= dNdx(INODE,2)
                BB(3,IY)= ZERO
                BB(4,IY)= dNdx(INODE,1)
                BB(5,IY)= ZERO
                BB(6,IY)= dNdx(INODE,3)
                BB(1,IZ)= ZERO
                BB(2,IZ)= ZERO
                BB(3,IZ)= dNdx(INODE,3)
                BB(4,IZ)= ZERO
                BB(5,IZ)= dNdx(INODE,1)
                BB(6,IZ)= dNdx(INODE,2)
            enddo
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
            phase=ZERO
            do I=1,nnode
                phase=phase+AN(I)*U(kblock,4*I)
            enddo
            if (phase>one) then
              phase=one
            else if (phase<zero) then
              phase=zero
            endif
  
            do I=1,3
                DP(I)=ZERO
            enddo
            do I=1,3
                do J=1,nnode
                    DP(I)=DP(I)+BP(I,J)*U(kblock,4*J)
                enddo
            enddo
            ProjDp = ZERO
            do I=1,3
              ProjDp = ProjDp + DP(I)*dir_f(I)
            enddo
C     ==================================================================
C     Nodal displacements
C     ==================================================================
            do J=1,3
             ULOC(J)=ZERO
            enddo
            do J=1,3
             do I=1,nnode
              do K=1,3                  
               ULOC(J)=ULOC(J)+VNI(J,3*(I-1)+K)*U(kblock,4*(I-1)+K)               
              enddo
             enddo
            enddo
C update the phased and displacement state variable
            sdv(21) = phase
            do J=1,3
                sdv(J)=ULOC(J)
            enddo
   
            call DegraFun_matrix(wT1d,wT2d,wTLd,wT1d_dd,wT2d_dd,
     1       wTLd_dd,a2,p,E1,E2,G12,G23,V12,V23,GcT1,GcT2,GcTL,
     2       c0,l0m,ST1,ST2,STL,phase)
C     ==================================================================
C     Calculated crack driving force
C     ==================================================================  
            call MutexLock( 1 )  
            if (time(itotalTime).le.zero) then
              phiT1_H= HALF*(ST1**2)/(C22+C23)
              phiT2_H= HALF*(ST2**2)/C44
              phiTL_H= HALF*(STL**2)/C55
            else            
              PhiT1=USRVAR(jelem(kblock),17,INPT)
              PhiT1_H=USRVAR(jelem(kblock),22,INPT)
              if (PhiT1 > PhiT1_H) then
              phiT1_H = PhiT1              
              endif

              PhiT2=USRVAR(jelem(kblock),18,INPT)
              PhiT2_H=USRVAR(jelem(kblock),23,INPT)
              if (PhiT2 > PhiT2_H) then
                  phiT2_H = PhiT2                 
              endif

              PhiTL=USRVAR(jelem(kblock),19,INPT)
              PhiTL_H=USRVAR(jelem(kblock),24,INPT)
              if (PhiTL > PhiTL_H) then
                  phiTL_H = PhiTL                 
              endif
            endif
            sdv(22) = phiT1_H
            sdv(23) = phiT2_H
            sdv(24) = phiTL_H
            drive = -wT1d_dd*PhiT1_H/GcT1-wT2d_dd*PhiT2_H/
     1          GcT2-wTLd_dd*PhiTL_H/GcTL   
C     ==================================================================
C     Calculating strain
C     ==================================================================
            do J=1,6
              strain(J)=ZERO
            enddo
            do I=1,6
             do J=1,nnode
              do K=1,3
                strain(I)=strain(I)+BB(I,3*(J-1)+K)*U(kblock,4*(J-1)+K)              
              enddo  
             enddo
            enddo
            do J=1,6
              sdv(J+3)=strain(J)
            enddo
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
            phase_m = USRVAR(jelem(kblock),21,INPT)
            phase_l = USRVAR(jelem(kblock),25,INPT)
C     ==================================================================
C     Calculating stresses and energy
C     ==================================================================
            call DegraFun_matrix(wT1d,wT2d,wTLd,wT1d_dd,wT2d_dd,
     1       wTLd_dd,a2,p,E1,E2,G12,G23,V12,V23,GcT1,GcT2,GcTL,
     2       c0,l0m,ST1,ST2,STL,phase_m)
            call DegraFun_fiber(WLd,wLd_dd,a2,p,E1,E2,G12,G23,V12,
     1       V23,GcL,c0,l0L,SYL,phase_l)

            call compute_stress_energy(stress,myenergy,strain,
     1       E1,E2,G12,G23,V12,V23,wT1d,wT2d,wTLd,WLd,nx,ny,nz) 

            do J=1,6
              sdv(J+9)=stress(j)
            enddo
            do J=1,5
              sdv(J+15)=myenergy(j)       
            enddo      
            
             
            if ( lflags(iOpCode).eq.jMassCalc ) then         
C     ==================================================================
C     lumped mass and capacitance matrix
C     ==================================================================
             do I=1,nnode
              amass(kblock,4*I,4*I)=amass(kblock,4*I,4*I)+etam*AN(I)*
     1         AN(I)*AINTW(INPT)*DTM
              do J=1,3
               Index = 4*(I-1)+J
               amass(kblock,Index,Index)=amass(kblock,Index,Index)+
     1          AN(I)*rho*AN(I)*AINTW(INPT)*DTM
              enddo
             enddo
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
           elseif ( lflags(iOpCode).eq.jIntForceAndDtStable) then             
             do I=1,nnode
              do J=1,3
               rhs(kblock,4*I)=rhs(kblock,4*I)+two/c0*l0m*BP(J,I)*
     1          (DP(J)+beita*dir_f(J)*two*ProjDp)*AINTW(INPT)*DTM
              enddo
              rhs(kblock,4*I)=rhs(kblock,4*I)-(drive-
     1         two/c0*(one-phase)/l0m)*AN(I)*AINTW(INPT)*DTM
             enddo

             

             do K1=1,nnode
              do I=1,3
               Index=4*(K1-1)+I
               do K4=1,6           
                rhs(kblock,Index)=rhs(kblock,Index)+AINTW(INPT)*
     1               BB(K4,3*(K1-1)+I)*STRESS(K4)*DTM 
               enddo
              enddo
            enddo
            endif
            
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================    
            do I=1,NSTV1
             svars(kblock,NSTV1*(INPT-1)+I)=sdv(I)
             USRVAR(jelem(kblock),I,INPT)=svars(kblock,NSTV1*(INPT-1)+I)                         
            enddo
            call MutexUnlock( 1 ) 
        enddo
      else
C     ==================================================================
C     Calculating for phase_fiber
C     ==================================================================  
        do INPT=1,Innode
C     Initializing solution dependent variables (phase,history)
            do I=1,NSTV2
             sdv(I)=svars(kblock,NSTV2*(INPT-1)+I)
            enddo
C     Local coordinates of the integration point
            XI(1) = XII(INPT,1)
            XI(2) = XII(INPT,2) 
            XI(3) = XII(INPT,3) 
C     Shape functions and local derivatives
            if (jtype.EQ.ONE) then
             call SHAPEFUNT(AN,dNdxi,XI)
            elseif (jtype.EQ.TWO) then
             call SHAPEFUN(AN,dNdxi,XI)
            endif
C     Jacobian
C     Note that ncrd and nnode are in a different order than uel
            do I = 1,3
             do J = 1,3
              VJACOB(I,J) = ZERO
              do K = 1,nnode
               VJACOB(I,J) = VJACOB(I,J)+coords(kblock,K,I)*dNdxi(K,J) 
     1                      
              enddo
             enddo
            enddo       
            DTM = ZERO
            DTM = VJACOB(1,1)*VJACOB(2,2)*VJACOB(3,3)+VJACOB(1,2)*
     1       VJACOB(2,3)*VJACOB(3,1)+VJACOB(1,3)*VJACOB(2,1)*
     2       VJACOB(3,2)-VJACOB(3,1)*VJACOB(2,2)*VJACOB(1,3)-
     3       VJACOB(3,2)*VJACOB(2,3)*VJACOB(1,1)-VJACOB(3,3)*
     4       VJACOB(2,1)*VJACOB(1,2)
C     
            if (DTM.lt.ZERO) then
             WRITE(7,*) 'Negative Jacobian',DTM
             call XPLB_EXIT
            endif

C     Inverse of Jacobian
            VJABOBINV(1,1)=(VJACOB(2,2)*VJACOB(3,3)-VJACOB(2,3)*
     1       VJACOB(3,2))/DTM
            VJABOBINV(1,2)=-(VJACOB(1,2)*VJACOB(3,3)-VJACOB(3,2)*
     1       VJACOB(1,3))/DTM
            VJABOBINV(1,3)=(VJACOB(1,2)*VJACOB(2,3)-VJACOB(1,3)*
     1       VJACOB(2,2))/DTM
            VJABOBINV(2,1)=-(VJACOB(2,1)*VJACOB(3,3)-VJACOB(2,3)*
     1       VJACOB(3,1))/DTM
            VJABOBINV(2,2)=(VJACOB(1,1)*VJACOB(3,3)-VJACOB(1,3)*
     1       VJACOB(3,1))/DTM
            VJABOBINV(2,3)=-(VJACOB(1,1)*VJACOB(2,3)-VJACOB(1,3)*
     1       VJACOB(2,1))/DTM
            VJABOBINV(3,1)=(VJACOB(2,1)*VJACOB(3,2)-VJACOB(2,2)*
     1       VJACOB(3,1))/DTM
            VJABOBINV(3,2)=-(VJACOB(1,1)*VJACOB(3,2)-VJACOB(1,2)*
     1       VJACOB(3,1))/DTM
            VJABOBINV(3,3)=(VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*
     1       VJACOB(2,1))/DTM       
C     Derivatives of shape functions respect to global ccordinates
            do K = 1,nnode
             do I = 1,3
                dNdx(K,I) = ZERO
              do J = 1,3
               dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)                         
              enddo
             enddo
            enddo
C     Calculating B matrix (B=LN) for phase
            do INODE=1,nnode
             BP(1,INODE)=dNdx(INODE,1)
             BP(2,INODE)=dNdx(INODE,2)
             BP(3,INODE)=dNdx(INODE,3)
            enddo
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
            phase=ZERO
            do I=1,nnode
                phase=phase+AN(I)*U(kblock,4*I)
            enddo 
            if (phase>one) then
              phase=one
            elseif (phase<zero) then
              phase=zero
            endif      
            do I=1,3
                DP(I)=ZERO
            enddo
            do I=1,3
                do J=1,nnode
                    DP(I)=DP(I)+BP(I,J)*U(kblock,4*J)
                enddo
            enddo
C update the phase state variable
            sdv(1) = phase
            call DegraFun_fiber(WLd,wLd_dd,a2,p,E1,E2,G12,G23,V12,V23,
     1       GcL,c0,l0L,SYL,phase)  
C     ==================================================================
C     Calculated crack driving force
C     ================================================================== 
           call MutexLock( 1 )          
            if (time(itotalTime).le.zero) then
              phiL_H = HALF*(SYL**2)/C11 
            else    
              PhiL=USRVAR(jelem(kblock)-N_ELEM*2,20,INPT)
              PhiL_H=USRVAR(jelem(kblock)-N_ELEM*2,26,INPT)
              if (PhiL > PhiL_H) then
                PhiL_H = PhiL                
              endif
            endif
            sdv(2) = PhiL_H
            drive = -wLd_dd*PhiL_H/GcL

            if ( lflags(iOpCode).eq.jMassCalc ) then         
C     ==================================================================
C     lumped mass and capacitance matrix
C     ==================================================================
             do I=1,nnode
              amass(kblock,4*I,4*I)=amass(kblock,4*I,4*I)+AN(I)*etaL*
     1         AN(I)*AINTW(INPT)*DTM
              do J=1,3
               Index = 4*(I-1)+J
C     ==================================================================
C     dummy displacement mass
C     ==================================================================     
               amass(kblock,Index,Index)=amass(kblock,Index,Index)+
     1          AN(I)*rho*AN(I)*AINTW(INPT)*DTM
              enddo
             enddo
       

            elseif ( lflags(iOpCode).eq.jIntForceAndDtStable) then             
             do I=1,nnode
              do J=1,3
               rhs(kblock,4*I)=rhs(kblock,4*I)+two/c0*l0L*BP(J,I)*DP(J)*
     1          AINTW(INPT)*DTM
              enddo
              rhs(kblock,4*I)=rhs(kblock,4*I)-(drive-
     1         two/c0*(one-phase)/l0L)*AN(I)*AINTW(INPT)*DTM
             enddo
            endif
                
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================    
            do I=1,NSTV2
             svars(kblock,NSTV2*(INPT-1)+I)=sdv(I)
             USRVAR(jelem(kblock)-N_ELEM*2,NSTV1+I,INPT)=svars(kblock,
     1        NSTV2*(INPT-1)+I)                         
            enddo
            call MutexUnlock( 1 ) 
        enddo
      endif
C     enforcing the Irreversible condition


      enddo    
      if (mydebug.and.time(itotalTime)>next_print_time) then
c       write(7,*) "time",time(itotalTime)
c       write(7,*) ,"kblock,jelem(kblock)",kblock,jelem(kblock)
c       write(7,*) " numThreads", numThreads
       next_print_time =  next_print_time + print_time
      endif  
      
      RETURN
      end
C --------------------------------------------------------------      
C Shape functions for tetrahedral elements
C --------------------------------------------------------------      
      subroutine SHAPEFUNT(AN,dNdxi,XI)
        include 'vaba_param.inc'
        Real*8 AN(4),dNdxi(4,3)
        Real*8 XI(3)
        parameter(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

C     Values of shape functions as a function of local coord.
        AN(1) = ONE-XI(1)-XI(2)-XI(3)
        AN(2) = XI(1)
        AN(3) = XI(2)
        AN(4) = XI(3)
C
C     Derivatives of shape functions respect to local coordinates
        do I=1,4
        do J=1,3
            dNdxi(I,J) =  ZERO
        enddo
        enddo
        dNdxi(1,1) =  MONE
        dNdxi(1,2) =  MONE
        dNdxi(1,3) =  MONE
C
        dNdxi(2,1) =  ONE
        dNdxi(2,2) =  ZERO
        dNdxi(2,3) =  ZERO
C
        dNdxi(3,1) =  ZERO
        dNdxi(3,2) =  ONE
        dNdxi(3,3) =  ZERO
C
        dNdxi(4,1) =  ZERO
        dNdxi(4,2) =  ZERO
        dNdxi(4,3) =  ONE
C
        RETURN
        end
C --------------------------------------------------------------      
C Shape functions for brick elements
C --------------------------------------------------------------      
      subroutine SHAPEFUN(AN,dNdxi,XI)
        include 'vaba_param.inc'
        REAL*8 AN(8),dNdxi(8,3)
        REAL*8 XI(3)
        parameter(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0,EIGHT=8.D0)

C     Values of shape functions as a function of local coord.
        AN(1) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE-XI(3))
        AN(2) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE-XI(3))
        AN(3) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE-XI(3))
        AN(4) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE-XI(3))
        AN(5) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE+XI(3))
        AN(6) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE+XI(3))
        AN(7) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE+XI(3))
        AN(8) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE+XI(3))
        
C     Derivatives of shape functions respect to local coordinates
        do I=1,8
        do J=1,3
            dNdxi(I,J) =  ZERO
        enddo
        enddo
        dNdxi(1,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
        dNdxi(1,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
        dNdxi(1,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
        dNdxi(2,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
        dNdxi(2,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
        dNdxi(2,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
        dNdxi(3,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
        dNdxi(3,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
        dNdxi(3,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
        dNdxi(4,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
        dNdxi(4,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
        dNdxi(4,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
        dNdxi(5,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
        dNdxi(5,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
        dNdxi(5,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
        dNdxi(6,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
        dNdxi(6,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
        dNdxi(6,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
        dNdxi(7,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
        dNdxi(7,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
        dNdxi(7,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
        dNdxi(8,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
        dNdxi(8,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
        dNdxi(8,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
        
        RETURN
        end
C --------------------------------------------------------------      
C Elastic matrix transformation
C --------------------------------------------------------------     
      subroutine Cmatrix(C11,C22,C44,C55,C12,C23,E1,E2,G12,G23,V12,V23)
      include 'vaba_param.inc'
      real*8 C11,C22,C44,C55,C12,C23,E1,E2,G12,G23,V12,V23
      parameter (ONE=1.D0,TWO=2.D0)
      C11 = -(E1*(E1 - E1*V23))/(TWO*E2*V12**TWO - E1 + E1*V23)
      C22 = -(E2*(- E2*V12**TWO + E1))/(TWO*E2*V12**TWO*V23 +  
     1  TWO*E2*V12**TWO + E1*V23**TWO - E1)
      C44 = G12
      C12 = -(E1*E2*V12)/(TWO*E2*V12**TWO - E1 + E1*V23)
      C23 = -(E2**TWO*V12**TWO + E1*V23*E2)/(TWO*E2*V12**TWO*V23 + 
     1  TWO*E2*V12**TWO + E1*V23**TWO- E1)
      C55 = ONE/TWO*(C22-C23)
      end
C --------------------------------------------------------------      
C Matrix and fiber degradation function
C -------------------------------------------------------------- 
      subroutine DegraFun_matrix(wT1d,wT2d,wTLd,wT1d_dd,wT2d_dd,
     1   wTLd_dd,a2,p,E1,E2,G12,G23,V12,V23,GcT1,GcT2,GcTL,
     2   c0,l0m,ST1,ST2,STL,phase)     
      include 'vaba_param.inc'
      parameter (ONE=1.D0,TWO=2.D0)
      real*8 wT1d,wT2d,wTLd,wT1d_dd,wT2d_dd,
     1  wTLd_dd,a2,p,E1,E2,G12,G23,V12,V23,GcT1,GcT2,GcTL,
     2  c0,l0m,ST1,ST2,STL,phase
      real*8 C11,C22,C44,C55,C12,C23,a1_T1,a1_T2,a1_TL

      call Cmatrix(C11,C22,C44,C55,C12,C23,E1,E2,G12,G23,V12,V23)
C degradation function for transverse tensile failure
      call Co_a1(a1_T1,TWO*(C22+C23),GcT1,ST1,c0,l0m)
      call DegraFun(wT1d,wT1d_dd,phase,a1_T1,a2,p)
C degradation function for transverse shear failure
      call Co_a1(a1_T2,C55,GcT2,ST2,c0,l0m)
      call DegraFun(wT2d,wT2d_dd,phase,a1_T2,a2,p)
C degradation function for transverse longitudinal failure
      call Co_a1(a1_TL,C44,GcTL,STL,c0,l0m)
      call DegraFun(wTLd,wTLd_dd,phase,a1_TL,a2,p)
                

      end

      subroutine DegraFun_fiber(WLd,wLd_dd,a2,p,E1,E2,G12,G23,V12,V23,
     1   GcL,c0,l0L,SYL,phase)  
      include 'vaba_param.inc'
            real*8 WLd,wLd_dd,a2,p,E1,E2,G12,G23,V12,V23,
     1   GcL,c0,l0L,SYL,phase
            real*8 C11,C22,C44,C55,C12,C23,a1_L
      
            call Cmatrix(C11,C22,C44,C55,C12,C23,E1,E2,G12,G23,V12,V23)
C degradation function for longitudinal failure
            call Co_a1(a1_L,C11,GcL,SYL,c0,l0L)
            call DegraFun(wLd,wLd_dd,phase,a1_L,a2,p)
      end
C --------------------------------------------------------------      
C degradation function
C -------------------------------------------------------------- 
      subroutine DegraFun(wd,wd_dd,phase,a1,a2,p)
        include 'vaba_param.inc'
        parameter (ONE=1.d0,TWO=2.d0)
        real*8 wd,wd_dd,phase,a1,a2,p
        wd = (ONE-phase)**p/((ONE-phase)**p+a1*phase+a1*a2*phase*phase)
        wd_dd = - ((ONE-phase)**p * (a1 - p * (ONE-phase)**(p - ONE)+ 
     1      TWO * a1 * a2 * phase))/(a1 * phase + (ONE - phase)**p + 
     2      a1 * a2 * phase**TWO)**TWO - (p * (ONE - phase)**(p - ONE))/
     3       (a1 * phase + (ONE - phase)**p + a1 * a2 * phase**TWO)         
      end
C
C --------------------------------------------------------------      
C Calculate the degenerate function coefficient a1
C --------------------------------------------------------------     
      subroutine Co_a1(a1,E,Gc,ft,c0,l0)
        include 'vaba_param.inc'
        real*8 a1,E,Gc,ft,c0,l0
        parameter (FOUR=4.D0)
        a1 = FOUR*E*Gc/c0/l0/ft/ft
      end

      subroutine compute_stress_energy(stress,energy,strain,E1,E2,G12,
     1    G23,V12,V23,wT1d,wT2d,wTLd,WLd,nx,ny,nz)
        include 'vaba_param.inc'
        real*8 E1,E2,G12,G23,V12,V23,
     1      wT1d,wT2d,wTLd,WLd,nx,ny,nz
        
        real*8 stress(6),strain(6),energy(5),dphi_dstrainilon(4),
     1      dede(4,6),rotate_strain(6),rotate_stress(6)
        parameter (ONE=1.D0,TWO=2.D0,HALF=0.5D0)
        real*8 C11,C22,C44,C55,C12,C23,eL,eT1,eT2,eTL
        call Cmatrix(C11,C22,C44,C55,C12,C23,E1,E2,G12,G23,V12,V23)
        call rotate_components(strain,nx,ny,nz,rotate_strain)
        call character_strain(eL,eT1,eT2,eTL,rotate_strain)

        dphi_dstrainilon(1)=wLd*C11*(abs(eL)+eL)*HALF +
     1    C11*(-abs(eL)+eL)*HALF+wLd*wT1d*TWO*C12*eT1    
        dphi_dstrainilon(2)=TWO*wT1d*(C22+C23)*(abs(eT1)+eT1)*HALF +
     1     TWO*(C22+C23)*(-abs(eT1)+eT1)*HALF+wLd*wT1d*TWO*C12*eL
        dphi_dstrainilon(3)=wTLd*C44*eTL
        dphi_dstrainilon(4)=wT2d*C55*eT2
        do I=1,4
         do J=1,6
          dede(I,J)=0
         enddo
        enddo 
        dede(1,1)=ONE
        dede(2,2)=HALF
        dede(2,3)=HALF
        if (eTL>0) then
          dede(3,4)=rotate_strain(4)/eTL
          dede(3,6)=rotate_strain(6)/eTL
        endif
        if (eT2>0) then
          dede(4,2)=(rotate_strain(2)-rotate_strain(3))/eT2
          dede(4,3)=-dede(4,2)
          dede(4,5)=rotate_strain(5)/eT2
        endif
        do J = 1, 6
            stress(J) = 0
            do I = 1, 4
               stress(J) = stress(J) + dphi_dstrainilon(I)*dede(I, J)
            enddo
        enddo 
        call inv_rotate_components(stress,nx,ny,nz,rotate_stress)
        stress=rotate_stress
C 1:total_energy 2:ET1 3:ET2 4:ETL 5:EL
         energy(1) = HALF*C11*eL**TWO + TWO*C12*eL*eT1 + (C22+C23)*  
     1    eT1**TWO + half*C44*eTL**TWO + HALF*C55*eT2
         energy(2) = (C22+C23)*((abs(eT1)+eT1)*HALF)**TWO
         energy(3) = HALF*C55*eT2**TWO
         energy(4) = HALF*C44*eTL**TWO
         energy(5) = HALF*C11*((abs(eL)+eL)*HALF)**TWO
      end

      subroutine character_strain(eL,eT1,eT2,eTL,strain)
      include 'vaba_param.inc'
      real*8 eL,eT1,eT2,eTL
      real*8 strain(6)
      parameter (TWO=2.D0)
        eL = strain(1)
        eT1 = (strain(2)+strain(3))/TWO
        eT2 = sqrt((strain(2)-strain(3))**TWO+strain(5)**TWO)
        eTL = sqrt(strain(4)**TWO+strain(6)**TWO)
      end

      subroutine rotate_components(vector, nx, ny, nz,
     1     rotate_vector)
          include 'vaba_param.inc'
          real*8, intent(in) :: vector(6)
          real*8, intent(in) :: nx, ny, nz           
          real*8, intent(out) :: rotate_vector(6)
  
          real*8  full_matrix(3,3),full_rotate_matrix(3,3),rotate(3,3)
          parameter (two=2.d0)
          call vector_to_matrix(vector, full_matrix)
          rotate = reshape([nx,ny,nz,-ny,nx,0.d0,nx*nz,-ny*nz,
     1   (nx**two+ny**two)], [3, 3])
          full_rotate_matrix=matmul(transpose(rotate),
     1      matmul(full_matrix,rotate))
          call matrix_to_vector(full_rotate_matrix,rotate_vector)
      
        end subroutine rotate_components

      subroutine inv_rotate_components(vector, nx, ny, nz,
     1     inv_rotate_vector)
        include 'vaba_param.inc'

        real*8, intent(in) :: vector(6)
        real*8, intent(in) :: nx, ny, nz               
        real*8, intent(out) :: inv_rotate_vector(6)
        parameter (two=2.d0)
        real*8  full_matrix(3,3),full_rotate_matrix(3,3),rotate(3,3)

        call vector_to_matrix(vector, full_matrix)
        rotate = reshape([nx,ny,nz,-ny,nx,0.d0,nx*nz,-ny*nz,
     1   (nx**two+ny**two)], [3, 3])
        full_rotate_matrix=matmul(rotate,
     1      matmul(full_matrix,transpose(rotate)))
        call matrix_to_vector(full_rotate_matrix,inv_rotate_vector)
       
        end subroutine inv_rotate_components

      subroutine matrix_to_vector(matrix, vector)
        include 'vaba_param.inc'
        real*8, intent(in) ::  matrix(3,3) 
        real*8, intent(out) :: vector(6)  
        
        vector(1) = matrix(1,1)
        vector(2) = matrix(2,2)
        vector(3) = matrix(3,3)
        vector(4) = matrix(1,2)
        vector(5) = matrix(2,3)
        vector(6) = matrix(1,3)
      end subroutine matrix_to_vector

      subroutine vector_to_matrix(vector,matrix)
        include 'vaba_param.inc'
        real*8, intent(in) ::  vector(6)  
        real*8, intent(out) ::  matrix(3,3)
        
        matrix(1,1) = vector(1) 
        matrix(2,2) = vector(2)
        matrix(3,3) = vector(3)
        matrix(1,2) = vector(4)
        matrix(2,3) = vector(5)
        matrix(1,3) = vector(6)
        matrix(2,1) = matrix(1,2)
        matrix(3,2) = matrix(2,3)
        matrix(3,1) = matrix(1,3)
      end subroutine vector_to_matrix
      
      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, jInfoArray,
     2  stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
      parameter (i_info_AnnealFlag = 1, 
     *     i_info_Intpt    = 2, ! Integration station number
     *     i_info_layer  = 3, ! Layer number
     *     i_info_kspt   = 4, ! Section point number in current layer
     *     i_info_effModDefn = 5, ! =1 if Bulk/ShearMod need to be defined
     *     i_info_ElemNumStartLoc   = 6) ! Start loc of user element number
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), dtArray(2*(nblock)+1), strainInc(nblock,
     2  ndir+nshr),relSpinInc(nblock,nshr), tempOld(nblock), 
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock), jInfoArray(*)
C
      parameter(N_ELEM=1,NSTV=26)
      parameter( zero = 0., one = 1., two = 2., three = 3.,
     1  third = one/three, half = .5, twoThirds = two/three,
     2  threeHalfs = 1.5 )
      COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
      character*80 cmname
C
      pointer (ptrjElemNum, jElemNum)
      dimension jElemNum(nblock)
      dimension DDSDDE(ndir+nshr,ndir+nshr)
C
      lAnneal = jInfoArray(i_info_AnnealFlag) 
      iLayer = jInfoArray(i_info_layer)
      kspt   = jInfoArray(i_info_kspt)
      intPt  = jInfoArray(i_info_Intpt)
      iUpdateEffMod = jInfoArray(i_info_effModDefn)
      iElemNumStartLoc = jInfoArray(i_info_ElemNumStartLoc)
      ptrjElemNum = loc(jInfoArray(iElemNumStartLoc))

      EMOD=PROPS(1)
      ENU=PROPS(2)
      EG=EMOD/(TWO*(ONE+ENU))
      EG2=EG*TWO
      ELAM=EG2*ENU/(ONE-TWO*ENU)
      NTENS = ndir+nshr
      
C
C   Stiffness tensor
C
        DO K1=1, NTENS
            DO K2=1, NTENS
            DDSDDE(K2, K1)=0.0
            enddo
        enddo
C
        DO K1=1, ndir
            DO K2=1, ndir
            DDSDDE(K2, K1)=ELAM
            enddo
            DDSDDE(K1, K1)=EG2+ELAM
        enddo
C
        DO K1=ndir+1, NTENS
            DDSDDE(K1, K1)=EG
        enddo

        do i = 1, nblock
        do K1 = 1,NTENS
            stressNew(i,K1) = stressOld(i,K1)
        enddo
C Trial stress
        do k1=1,NTENS
         do k2=1,NTENS
            stressNew(i,K2) = stressNew(i,K2)+DDSDDE(K2,K1)*
     1       strainInc(i,K1)
         enddo
        enddo
      enddo

      do kblock = 1, nblock
        IF (intPt.EQ.3) THEN
         intPt=4
        ELSEIF (intPt.EQ.4) THEN
         intPt=3
        endIF
        IF (intPt.EQ.7) THEN
         intPt=8
        ELSEIF (intPt.EQ.8) THEN
         intPt=7
        endIF
           
        do I=1,nstatev
         stateNew(kblock,I)=USRVAR(jElemNum(kblock)-N_ELEM,I,intPt)
        enddo
      enddo
      return
      end
    
    



     