!
! Copyright (c) Stanford University, The Regents of the University of
!               California, and others.
!
! All Rights Reserved.
!
! See Copyright-SimVascular.txt for additional details.
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject
! to the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------------------------------------------------------------------
!
!     This is for solving fluid transport equation solving Navier-Stokes
!     equations. Dirichlet boundary conditions are either treated
!     strongly or weakly.
!
!--------------------------------------------------------------------


      SUBROUTINE IFEM_INIT(Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) iM

      DO iM=1, nMsh
         msh(iM)%iGC = 0
         CALL GETNSTENCIL(msh(iM))   
      END DO   

      DO iM=2, nMsh
         CALL IFEM_FINDCLOSEST(msh(1), msh(iM), Dg)
         IF( intrp .EQ. ifemIntrp_MLS) THEN 
            CALL IFEM_FINDMLSW(msh(1), msh(iM), Dg)
         END IF
      END DO

      IF( nMsh .EQ. 2) THEN
         CALL IFEM_FINDHIDDEN_ND(msh(1), msh(2), Dg)
         CALL IFEM_FINDHIDDEN_ELM_FLUID(msh(1), msh(2), Dg)
      ELSE IF( nMsh .EQ. 3) THEN

!        If we need the hidden nodes from the solid
!         CALL IFEM_FINDHIDDEN_SOL(msh(1), msh(3), Dg)
!        To use lstDirNd
C        CALL IFEM_FINDHIDDEN_ND(msh(1), msh(3), Dg)
         CALL IFEM_FINDHIDDEN_ND(msh(1), msh(2), Dg)
!        Compute hidden elm and lstDirNd         
         CALL IFEM_FINDHIDDEN_ELM_FS(msh(1), Dg) 


      END IF

      IF( intrp .EQ. ifemIntrp_FEM) THEN 
         CALL IFEM_FIND_Nd2Elm_BG(msh(1), msh(2), Dg)
         CALL IFEM_FIND_Nd2Elm_FG(msh(1), msh(2), Dg)
      END IF

!     Allocation of arrays to store unknown in Fg and Bg meshes 
      ALLOCATE(msh(1)%YgBG(nsd,msh(1)%nNo))
      ALLOCATE(msh(2)%YgFG(nsd+1,msh(2)%nNo))
      IF( nMsh .EQ. 3) ALLOCATE(msh(3)%YgFG(nsd+1,msh(3)%nNo))

      RETURN
      END SUBROUTINE IFEM_INIT


!####################################################################


      SUBROUTINE IFEM_UPDATE(Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) iM

      DO iM=1, nMsh
         msh(iM)%iGC = 0
         CALL GETNSTENCIL(msh(iM))   
      END DO   

      DO iM=2, nMsh
         CALL IFEM_FINDCLOSEST(msh(1), msh(iM), Dg)
         IF( intrp .EQ. ifemIntrp_MLS) THEN 
            CALL IFEM_FINDMLSW(msh(1), msh(iM), Dg)
         END IF
      END DO

      IF( nMsh .EQ. 2) THEN
         CALL IFEM_FINDHIDDEN_ND(msh(1), msh(2), Dg)
         CALL IFEM_FINDHIDDEN_ELM_FLUID(msh(1), msh(2), Dg)
      ELSE IF( nMsh .EQ. 3) THEN

!        If we need the hidden nodes from the solid
!         CALL IFEM_FINDHIDDEN_SOL(msh(1), msh(3), Dg)
         CALL IFEM_FINDHIDDEN_ND(msh(1), msh(2), Dg)
         CALL IFEM_FINDHIDDEN_ELM_FS(msh(1), Dg)

      END IF

      IF( intrp .EQ. ifemIntrp_FEM) THEN 
         CALL IFEM_FIND_Nd2Elm_BG(msh(1), msh(2), Dg)
         CALL IFEM_FIND_Nd2Elm_FG(msh(1), msh(2), Dg)
      END IF

!     Allocation of arrays to store unknown in Fg and Bg meshes 
      IF(.NOT.ALLOCATED(msh(1)%YgBG)) THEN 
         ALLOCATE(msh(1)%YgBG(nsd,msh(1)%nNo))
      END IF
      IF(.NOT.ALLOCATED(msh(2)%YgFG)) THEN 
         ALLOCATE(msh(2)%YgFG(nsd+1,msh(2)%nNo))
      END IF
      IF( nMsh .EQ. 3) THEN 
         IF( .NOT.ALLOCATED(msh(3)%YgFG) ) THEN 
            ALLOCATE(msh(3)%YgFG(nsd+1,msh(3)%nNo))
         END IF 
      END IF

      RETURN
      END SUBROUTINE IFEM_UPDATE

!####################################################################
!####################################################################

      SUBROUTINE CONSTRUCT_FLUID_MM(lM, Ag, Yg, iM, iter)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(IN) :: iM, iter

      LOGICAL :: vmsStab
      INTEGER(KIND=IKIND) a, e, g, l, Ac, eNoN, cPhys, AcLc
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)
      TYPE(fsType) :: fs(2)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   bfl(:,:), lR(:,:), lK(:,:,:), ylFg(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xwl(:,:), xql(:,:), Nwx(:,:),
     2   Nwxx(:,:), Nqx(:,:)

      eNoN = lM%eNoN

      IF (lM%nFs .EQ. 1) THEN
         vmsStab = .TRUE.
      ELSE
         vmsStab = .FALSE.
      END IF

!     l = 3, if nsd==2 ; else 6;
      l = nsymd

!     FLUID: dof = nsd+1
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   bfl(nsd,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))

      IF( iM .EQ. 2) ALLOCATE(ylFg(nsd+1,eNoN)) 

!     Loop over all elements of mesh
      DO e=1, lM%nEl

!        Jump alloc if elem is hidden for bh mesh 
         IF ( iter .GE. 3) THEN 
            IF((iM .EQ. 1) .AND. (lM%lstHdnElm(e) .EQ. 1)) CYCLE
         END IF

!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_fluid) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Create local copies
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)   = Ac
            xl(:,a)  = x(:,Ac)
            al(:,a)  = Ag(:,Ac)
            yl(:,a)  = Yg(:,Ac)
            bfl(:,a) = Bf(:,Ac)
         END DO

!        Create local copies if we are in the FG mesh 
         IF( iM .GE. 2) THEN 
            DO a=1, eNoN
               Ac = lM%IEN(a,e)
               AcLc = lM%lN(Ac)

C                ylFg(1:nsd,a) = Yg(1:nsd,Ac)
C                ylFg(nsd+1,a) = lM%YgFG(nsd+1,AcLc)

               ylFg(:,a) = lM%YgFG(:,AcLc)
            END DO
         END IF

!        Initialize residue and tangents
         lR = 0._RKIND
         lK = 0._RKIND

!        Set function spaces for velocity and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 1)

!        Define element coordinates appropriate for function spaces
         ALLOCATE(xwl(nsd,fs(1)%eNoN), Nwx(nsd,fs(1)%eNoN),
     2      Nwxx(l,fs(1)%eNoN))
         ALLOCATE(xql(nsd,fs(2)%eNoN), Nqx(nsd,fs(2)%eNoN))
         xwl(:,:) = xl(:,:)
         xql(:,:) = xl(:,1:fs(2)%eNoN)
         Nwx      = 0._RKIND
         Nqx      = 0._RKIND
         Nwxx     = 0._RKIND

!        Gauss integration 1
         DO g=1, fs(1)%nG
            IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
               CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF

            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

               CALL GNNxx(l, fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g),
     2            fs(1)%Nxx(:,:,g), xwl, Nwx, Nwxx)
            END IF
            w = fs(1)%w(g) * Jac

            IF (nsd .EQ. 3) THEN
               CALL FLUID3D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl,
     3            bfl, lR, lK)
               
               IF( iM .GE. 2) THEN 
                  CALL IFEM_2DRES(fs(1)%eNoN, fs(2)%eNoN, w,
     2               fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, yl, ylFg, lR)
               END IF

             ELSE IF (nsd .EQ. 2) THEN
               CALL FLUID2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl,
     3            bfl, lR, lK)
               IF( iM .GE. 2) THEN 

C                   CALL IFEM_2DRES(fs(1)%eNoN, fs(2)%eNoN, w,
C      2               fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, yl, ylFg, lR)
               END IF
            END IF
         END DO ! g: loop

!        Set function spaces for velocity and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 2)

!        Gauss integration 2
         DO g=1, fs(2)%nG
            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF

            IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
               CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = fs(2)%w(g) * Jac

            IF (nsd .EQ. 3) THEN
               CALL FLUID3D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl,
     3            bfl, lR, lK)

            ELSE IF (nsd .EQ. 2) THEN
               CALL FLUID2D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl,
     3            bfl, lR, lK)
            END IF
         END DO ! g: loop

         DEALLOCATE(xwl, xql, Nwx, Nwxx, Nqx)

!        Assembly
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
         END IF
#endif
      END DO ! e: loop

      DEALLOCATE(ptr, xl, al, yl, bfl, lR, lK)
      IF(ALLOCATED(ylFg)) DEALLOCATE(ylFg)  

      CALL DESTROY(fs(1))
      CALL DESTROY(fs(2))


      RETURN
      END SUBROUTINE CONSTRUCT_FLUID_MM
!####################################################################
      SUBROUTINE IFEM_2DRES(eNoNw, eNoNq, w, Nw, Nq, Nwx,
     2   Nqx, yl, ylFg, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Nw(eNoNw), Nq(eNoNq),
     2   Nwx(2,eNoNw), Nqx(2,eNoNq), yl(nsd+1,eNoNw), ylFg(nsd+1,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw) 

      INTEGER(KIND=IKIND) a, b, k
      REAL(KIND=RKIND) amd, wl, wr, rho, 
     2   gam, mu, mu_s, mu_g, p, u(2), ux(2,2), es(2,2), rM(2,2), T1


      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)

      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      amd  = eq(cEq)%am/T1
      wl   = w*T1
      wr   = w*rho

!     Note that indices are not selected based on the equation because
!     fluid equation always come first
!     Velocity and its gradients, inertia (acceleration & body force)
      u   = 0._RKIND
      ux  = 0._RKIND
      DO a=1, eNoNw

         u(1)    = u(1) + Nw(a)*ylFg(1,a)
         u(2)    = u(2) + Nw(a)*ylFg(2,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*ylFg(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*ylFg(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*ylFg(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*ylFg(2,a)

      END DO

!     Pressure 
      p  = 0._RKIND
      DO a=1, eNoNq
         p  = p + Nq(a)*ylFg(3,a)
      END DO

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)
      es(2,1) = ux(2,1) + ux(1,2)
      es(1,2) = es(2,1)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(2,1)*es(2,1)
     2    + es(1,2)*es(1,2) + es(2,2)*es(2,2)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
!     The returned mu_g := (d\mu / d\gamma)
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, mu, mu_s, mu_g)
  
      rM(1,1) = mu*es(1,1) - p
      rM(2,1) = mu*es(2,1) 

      rM(1,2) = mu*es(1,2) 
      rM(2,2) = mu*es(2,2) - p


!     Local residue
      DO a=1, eNoNw
         lR(1,a) = lR(1,a) - w*(Nwx(1,a)*rM(1,1)
     2      + Nwx(2,a)*rM(2,1))

         lR(2,a) = lR(2,a) - w*(Nwx(1,a)*rM(1,2)
     2      + Nwx(2,a)*rM(2,2))
      END DO

      RETURN
      END SUBROUTINE IFEM_2DRES

!####################################################################
!####################################################################
!     Add contribution from background fluid stresses on the boundary
      SUBROUTINE SET_BG_BCNEU_TO_FG(Yg, Dg, iM)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(IN) :: iM

      INTEGER(KIND=IKIND) iFa, jM

      IF( iM .EQ. 1 ) RETURN
      
      jM = 2

      DO iFa = 1, msh(jM)%nFa
         CALL CONSTRUCT_NEUL(jM, msh(jM)%fa(iFa), Yg, Dg)
      END DO

      RETURN
      END SUBROUTINE SET_BG_BCNEU_TO_FG
!####################################################################
!     Construct Neumann BCs
      SUBROUTINE CONSTRUCT_NEUL(iM, lFa, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iM
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, Ec, cPhys, eNoN, AcLc, b, 
     2             eNoNb, bRob 
      REAL(KIND=RKIND) w, h, nV(nsd), y(tDof), Jac, ksix(nsd,nsd), 
     2             st(nsd), es(nsd,nsd), u(nsd), p, ux(nsd,nsd), vis,
     3             esN(nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:), ptrBd(:), 
     2             maplocBdElm(:)
      REAL(KIND=RKIND), ALLOCATABLE :: N(:), hl(:), yl(:,:), lR(:,:),
     2             lK(:,:,:), ylFg(:,:), xl(:,:), Nx(:,:), qpBd(:,:), 
     3             xRefBnd(:,:), RefQp(:,:)


      eNoNb = lFa%eNoN
      eNoN  = msh(iM)%eNoN
      bRob = 0.5_RKIND

      DO e=1, lFa%nEl
C          write(*,*)" face elem local ", e
!        Get the id element in the global ordering of the element owing this face  
         Ec    = lFa%gE(e)
         cDmn  = DOMAIN(msh(iM), cEq, Ec)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         vis   = eq(cEq)%dmn(cDmn)%visc%mu_i

C          write(*,*)" Ec ", Ec
C          write(*,*)" vis ", vis

         ALLOCATE( ptr(eNoN), ptrBd(eNoNb), N(eNoN), yl(tDof,eNoN),
     2      lR(dof,eNoN), lK(dof*dof,eNoN,eNoN), ylFg(nsd+1,eNoN), 
     3      maplocBdElm(eNoNb), Nx(nsd,eNoN), qpBd(nsd,lFa%nG), 
     4      xRefBnd(nsd,eNoNb), RefQp(nsd-1,lFa%nG), xl(nsd,eNoN) )
         lK = 0._RKIND
         lR = 0._RKIND
         Nx = 0._RKIND
         N  = 0._RKIND

!        Get coordinates elements, local solution, local bg vel 
!        and define local map from face to elem 
         DO a=1, eNoN
            Ac = msh(iM)%IEN(a,Ec)

            xl(:,a) = x(:,Ac)
            ptr(a)  = Ac
            yl(:,a) = Yg(:,Ac)
            IF (mvMsh) xl(:,a) = xl(:,a) + Do(nsd+2:2*nsd+1,Ac)

!           We store in ylFg the local vel and pressure from the BG mesh 
            AcLc = msh(iM)%lN(Ac)
            ylFg(:,a) = msh(iM)%YgFG(:,AcLc)
         END DO

         DO a=1, eNoNb
            Ac      = lFa%IEN(a,e)
            ptrBd(a)  = Ac
         END DO

!        Building map from local face id to local id in the bulm element 
         DO a=1, eNoNb
            Ac = ptrBd(a)
            DO b=1, eNoN
               IF( ptr(b) .EQ. Ac ) maplocBdElm(a) = b
            END DO
         END DO


!        Map quad point into parent bulk element
!        Fisrt, get coordinates of bnd elem in the ref bulk elem 
         CALL GETCFace(msh(iM)%eType, eNoN, eNoNb, maplocBdElm, xRefBnd)
!        Get quad point of the ref bnd elem 
         CALL GETGIP_RefBDELM(nsd-1, lFa%eType, eNoNb, lFa%nG, RefQp)         
!        Now, we map the quad point from fre bnd elm to ref bnd face (bulk elm)
         CALL GETP_BDELM(nsd-1, lFa%eType, eNoNb, lFa%nG, xRefBnd, 
     2       RefQp, qpBd)


C          write(*,*)" before mapping "
C             write(*,*)" Ac = ", Ac
C             write(*,*)" xl = ", xl
C             write(*,*)" ptr = ", ptr
C             write(*,*)" yl = ", yl
C             write(*,*)" ylFg = ", ylFg
C             write(*,*)" ptrBd = ", ptrBd
C             write(*,*)" maplocBdElm = ", maplocBdElm

C             write(*,*)" xRefBnd = ", xRefBnd
C             write(*,*)" RefQp = ", RefQp
C             write(*,*)" qpBd = ", qpBd

         DO g=1, lFa%nG

!           Evaluate gradient at this quad point 
!           The grad is const, we can evaluate it in any quad point
            CALL GNN(eNoN, nsd, msh(iM)%fs(1)%Nx(:,:,g), xl, Nx, Jac,
     2            ksix)

!           Evaluate test function at this quad point 
            CALL GETN(nsd, msh(iM)%eType, eNoN, qpBd(:,g), N)

            CALL GNNB(lFa, e, g, nsd-1, eNoNb, lFa%Nx(:,:,g), nV)
            Jac = SQRT(NORM(nV))
            nV  = nV/Jac
            w   = lFa%w(g)*Jac

C             write(*,*)" nV = ", nV
C             write(*,*)" w = ", w

!           Compute sigma(u,p) n 
            u = 0._RKIND
            p = 0._RKIND
            ux = 0._RKIND
            es = 0._RKIND
            esN = 0._RKIND
            st = 0._RKIND

            DO a=1, eNoN
               u(1) = u(1) + N(a)*ylFg(1,a)
               u(2) = u(2) + N(a)*ylFg(2,a)

               p = p + N(a)*ylFg(3,a)

               ux(1,1) = ux(1,1) + Nx(1,a)*ylFg(1,a)
               ux(2,1) = ux(2,1) + Nx(1,a)*ylFg(2,a)
               ux(1,2) = ux(1,2) + Nx(2,a)*ylFg(1,a)
               ux(2,2) = ux(2,2) + Nx(2,a)*ylFg(2,a)
            END DO

!           Strain rate tensor 2*e_ij := (u_ij + u_ji)
            es(1,1) = ux(1,1) + ux(1,1)
            es(2,1) = ux(2,1) + ux(1,2)
            es(2,2) = ux(2,2) + ux(2,2)
            es(1,2) = es(2,1)

            esN(1) = vis*(es(1,1)*nV(1) + es(1,2)*nV(2))
            esN(2) = vis*(es(2,1)*nV(1) + es(2,2)*nV(2))

            st(1) = esN(1) - p*nV(1) 
            st(2) = esN(2) - p*nV(2)


            IF (nsd .EQ. 2) THEN
               DO a=1, eNoN
                  lR(1,a) = lR(1,a) - w*N(a)*st(1) 
                  lR(2,a) = lR(2,a) - w*N(a)*st(2) 
               END DO
            ELSE
               DO a=1, eNoN
                  lR(1,a) = lR(1,a) - w*N(a)*st(1) 
                  lR(2,a) = lR(2,a) - w*N(a)*st(2) 
                  lR(3,a) = lR(3,a) - w*N(a)*st(3)
               END DO
            END IF

!           For Robin type BC:  
C             IF (nsd .EQ. 2) THEN
C                DO a=1, eNoN
C                   lR(1,a) = lR(1,a) + 
C      2               bRob * w * N(a) * u(1) *( u(1)*nV(1) + u(2)*nV(2) )
C                   lR(2,a) = lR(2,a) + 
C      2               bRob * w * N(a) * u(2) *( u(1)*nV(1) + u(2)*nV(2) )             
C                END DO
C             ELSE
C                DO a=1, eNoN
C                   lR(1,a) = lR(1,a) + 
C      2               bRob * w * N(a) * u(1) *( u(1)*nV(1) + u(2)*nV(2) )
C                   lR(2,a) = lR(2,a) + 
C      2               bRob * w * N(a) * u(2) *( u(1)*nV(1) + u(2)*nV(2) )   
C                   lR(3,a) = lR(3,a) + 
C      2               bRob * w * N(a) * u(3) *( u(1)*nV(1) + u(2)*nV(2) )     
C                END DO
C             END IF           


         END DO

!        Now doing the assembly part
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK,  lR)
#ifdef WITH_TRILINOS
         END IF
#endif
         DEALLOCATE(ptr, ptrBd, N, yl, lR, lK, ylFg, maplocBdElm, 
     2              qpBd, xRefBnd, RefQp, xl, Nx)

      END DO

      RETURN
      END SUBROUTINE CONSTRUCT_NEUL

!####################################################################
!####################################################################
!     Add contribution from background fluid stresses on the boundary
      SUBROUTINE SET_BCNEU_TO_FG(Yg, Dg, iM)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(IN) :: iM

      INTEGER(KIND=IKIND) iFa, jM

      IF( iM .EQ. 1 ) RETURN
      IF( iM .EQ. 3 ) RETURN
      
      jM = 2

      DO iFa = 1, msh(jM)%nFa
         IF( iFa .LE. 5 ) CYCLE 
         CALL CONSTRUCT_NEUL(jM, msh(jM)%fa(iFa), Yg, Dg)
      END DO

      RETURN
      END SUBROUTINE SET_BCNEU_TO_FG
!####################################################################
!     Construct Neumann BCs
      SUBROUTINE CONSTRUCT_RK_NEUL(iM, lFa, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iM
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, Ec, cPhys, eNoN, AcLc, b, 
     2             eNoNb
      REAL(KIND=RKIND) w, h, nV(nsd), y(tDof), Jac, ksix(nsd,nsd), 
     2             st(nsd), es(nsd,nsd), u(nsd), p, ux(nsd,nsd), vis,
     3             esN(nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:), ptrBd(:), 
     2             maplocBdElm(:)
      REAL(KIND=RKIND), ALLOCATABLE :: N(:), hl(:), yl(:,:), lR(:,:),
     2             lK(:,:,:), ylFg(:,:), xl(:,:), Nx(:,:), qpBd(:,:), 
     3             xRefBnd(:,:), RefQp(:,:)


      eNoNb = lFa%eNoN
      eNoN  = msh(iM)%eNoN

      DO e=1, lFa%nEl
!        Get the id element in the global ordering of the element owing this face  
         Ec    = lFa%gE(e)
         cDmn  = DOMAIN(msh(iM), cEq, Ec)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         vis   = eq(cEq)%dmn(cDmn)%visc%mu_i

         ALLOCATE( ptr(eNoN), ptrBd(eNoNb), N(eNoN), yl(tDof,eNoN),
     2      lR(dof,eNoN), lK(dof*dof,eNoN,eNoN), ylFg(nsd+1,eNoN), 
     3      maplocBdElm(eNoNb), Nx(nsd,eNoN), qpBd(nsd,lFa%nG), 
     4      xRefBnd(nsd,eNoNb), RefQp(nsd-1,lFa%nG), xl(nsd,eNoN) )
         lK = 0._RKIND
         lR = 0._RKIND
         Nx = 0._RKIND
         N  = 0._RKIND

!        Get coordinates elements, local solution, local bg vel 
!        and define local map from face to elem 
         DO a=1, eNoN
            Ac = msh(iM)%IEN(a,Ec)

            xl(:,a) = x(:,Ac)
            ptr(a)  = Ac
            yl(:,a) = Yg(:,Ac)
            IF (mvMsh) xl(:,a) = xl(:,a) + Do(nsd+2:2*nsd+1,Ac)

!           We store in ylFg the local vel and pressure from the BG mesh 
            AcLc = msh(iM)%lN(Ac)
            ylFg(:,a) = msh(iM)%YgFG(:,AcLc)
         END DO

         DO a=1, eNoNb
            Ac      = lFa%IEN(a,e)
            ptrBd(a)  = Ac
         END DO

!        Building map from local face id to local id in the bulm element 
         DO a=1, eNoNb
            Ac = ptrBd(a)
            DO b=1, eNoN
               IF( ptr(b) .EQ. Ac ) maplocBdElm(a) = b
            END DO
         END DO


!        Map quad point into parent bulk element
!        Fisrt, get coordinates of bnd elem in the ref bulk elem 
         CALL GETCFace(msh(iM)%eType, eNoN, eNoNb, maplocBdElm, xRefBnd)
!        Get quad point of the ref bnd elem 
         CALL GETGIP_RefBDELM(nsd-1, lFa%eType, eNoNb, lFa%nG, RefQp)         
!        Now, we map the quad point from fre bnd elm to ref bnd face (bulk elm)
         CALL GETP_BDELM(nsd-1, lFa%eType, eNoNb, lFa%nG, xRefBnd, 
     2       RefQp, qpBd)

         DO g=1, lFa%nG

!           Evaluate gradient at this quad point 
!           The grad is const, we can evaluate it in any quad point
            CALL GNN(eNoN, nsd, msh(iM)%fs(1)%Nx(:,:,g), xl, Nx, Jac,
     2            ksix)

!           Evaluate test function at this quad point 
            CALL GETN(nsd, msh(iM)%eType, eNoN, qpBd(:,g), N)


            CALL GNNB(lFa, e, g, nsd-1, eNoNb, lFa%Nx(:,:,g), nV)
            Jac = SQRT(NORM(nV))
            nV  = nV/Jac
            w   = lFa%w(g)*Jac

!--         Foreground part
!           Compute sigma(u,p) n 
            u = 0._RKIND
            p = 0._RKIND
            ux = 0._RKIND
            es = 0._RKIND
            esN = 0._RKIND
            st = 0._RKIND

            DO a=1, eNoN
               u(1) = u(1) + N(a)*yl(1,a)
               u(2) = u(2) + N(a)*yl(2,a)

               p = p + N(a)*ylFg(3,a)

               ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a)
               ux(2,1) = ux(2,1) + Nx(1,a)*yl(2,a)
               ux(1,2) = ux(1,2) + Nx(2,a)*yl(1,a)
               ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a)
            END DO

!           Strain rate tensor 2*e_ij := (u_ij + u_ji)
            es(1,1) = ux(1,1) + ux(1,1)
            es(2,1) = ux(2,1) + ux(1,2)
            es(2,2) = ux(2,2) + ux(2,2)
            es(1,2) = es(2,1)

            esN(1) = vis*(es(1,1)*nV(1) + es(1,2)*nV(2))
            esN(2) = vis*(es(2,1)*nV(1) + es(2,2)*nV(2))

            st(1) = esN(1) - p*nV(1) 
            st(2) = esN(2) - p*nV(2)

C !--         Background part
C             uBG = 0._RKIND
C             pBG = 0._RKIND
C             stBG = 0._RKIND

C             DO a = 1, eNoN
C                uBG(1) = uBG(1) + N(a)*yl(1,a)
C                uBG(2) = uBG(2) + N(a)*yl(2,a)

C                pBG = pBG + N(a)*yl(3,a)

C                ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a)
C                ux(2,1) = ux(2,1) + Nx(1,a)*yl(2,a)
C                ux(1,2) = ux(1,2) + Nx(2,a)*yl(1,a)
C                ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a)
C             END DO

C !           Strain rate tensor 2*e_ij := (u_ij + u_ji)
C             es(1,1) = ux(1,1) + ux(1,1)
C             es(2,1) = ux(2,1) + ux(1,2)
C             es(2,2) = ux(2,2) + ux(2,2)
C             es(1,2) = es(2,1)

C             esN(1) = vis*(es(1,1)*nV(1) + es(1,2)*nV(2))
C             esN(2) = vis*(es(2,1)*nV(1) + es(2,2)*nV(2))

C             st(1) = esN(1) - p*nV(1) 
C             st(2) = esN(2) - p*nV(2)



!--         Assembly 
            IF (nsd .EQ. 2) THEN
               DO a=1, eNoN
                  lR(1,a) = lR(1,a) - w*N(a)*st(1) 
                  lR(2,a) = lR(2,a) - w*N(a)*st(2) 
               END DO
            ELSE
               DO a=1, eNoN
                  lR(1,a) = lR(1,a) - w*N(a)*st(1) 
                  lR(2,a) = lR(2,a) - w*N(a)*st(2) 
                  lR(3,a) = lR(3,a) - w*N(a)*st(3)
               END DO
            END IF


         END DO

!        Now doing the assembly part
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK,  lR)
#ifdef WITH_TRILINOS
         END IF
#endif
         DEALLOCATE(ptr, ptrBd, N, yl, lR, lK, ylFg, maplocBdElm, 
     2              qpBd, xRefBnd, RefQp, xl, Nx)

      END DO

      RETURN
      END SUBROUTINE CONSTRUCT_RK_NEUL


!####################################################################
!     Returns Gauss integration points in current coordinates from 
!     reference face (dimFace) to current face (dimFace)
      SUBROUTINE GETGIP_SElmBDELM(insd, eType, eNoNF, nG, crdFace, xi)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: insd, eType, nG, eNoNF
      REAL(KIND=RKIND), INTENT(IN) :: crdface(insd,eNoNF)
      REAL(KIND=RKIND), INTENT(OUT) :: xi(insd, nG)

      REAL(KIND=RKIND) :: s, t, in, jq, xiRef(insd,nG), N(eNoNF,nG)

      IF (eType .EQ. eType_NRB) RETURN

C       write(*,*)" inside GETGIP_SELM "

!     3D elements
      SELECT CASE(eType)
      CASE(eType_TRI3)
         s = 2._RKIND/3._RKIND
         t = 1._RKIND/6._RKIND
         xiRef(1,1) = t; xiRef(2,1) = t
         xiRef(1,2) = s; xiRef(2,2) = t
         xiRef(1,3) = t; xiRef(2,3) = s

!     2D elements
      CASE(eType_LIN1)
         s = 1._RKIND/SQRT(3._RKIND)
         xiRef(1,1) = -s
         xiRef(1,2) =  s

      END SELECT

C       write(*,*)" sub face ref elemem = ", crdFace

      DO in = 1, eNoNF
!        Get nodal shape function in the reference face element 
         CALL GETN(insd, eType, eNoNF, xiRef(:,in), N(:,in))
      END DO

!     N_ij is the weight that we will use to define the quad point 
!     in the new subElement in the reference bulk element with 
!     coordinates = crdFace
      xi = 0._RKIND

!     Loop over quad point          
      DO jq = 1, nG
!        Loop over element nodes 
         DO in = 1, eNoNF
            xi(:, jq) = xi(:, jq) + crdFace(:,in)*N(in,jq)
         END DO
      END DO  

C       write(*,*)" quad point sub elemem ", xi

      RETURN
      END SUBROUTINE GETGIP_SElmBDELM
!--------------------------------------------------------------------
      SUBROUTINE GETGIP_RefBDELM(insd, eType, eNoNF, nG, xiRef)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: insd, eType, nG, eNoNF
      REAL(KIND=RKIND), INTENT(OUT) :: xiRef(insd, nG)

      REAL(KIND=RKIND) :: s, t

      IF (eType .EQ. eType_NRB) RETURN

!     3D elements
      SELECT CASE(eType)
      CASE(eType_TRI3)
         s = 2._RKIND/3._RKIND
         t = 1._RKIND/6._RKIND
         xiRef(1,1) = t; xiRef(2,1) = t
         xiRef(1,2) = s; xiRef(2,2) = t
         xiRef(1,3) = t; xiRef(2,3) = s

!     2D elements
      CASE(eType_LIN1)
         s = 1._RKIND/SQRT(3._RKIND)
         xiRef(1,1) = -s
         xiRef(1,2) =  s

      END SELECT

      RETURN
      END SUBROUTINE GETGIP_RefBDELM
!--------------------------------------------------------------------
!     Returns shape functions at given point (reference coords)
      SUBROUTINE GETN(insd, eType, eNoN, xi, N)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: insd, eType, eNoN
      REAL(KIND=RKIND), INTENT(IN)  :: xi(insd)
      REAL(KIND=RKIND), INTENT(OUT) :: N(eNoN)

      IF (eType .EQ. eType_NRB) RETURN

!     3D face elements
      SELECT CASE(eType)
      CASE(eType_TRI3)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = 1._RKIND - xi(1) - xi(2)

!     2D face elements
      CASE(eType_LIN1)
         N(1) = (1._RKIND - xi(1))*0.5_RKIND
         N(2) = (1._RKIND + xi(1))*0.5_RKIND

      CASE(eType_TET4)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = xi(3)
         N(4) = 1._RKIND - xi(1) - xi(2) - xi(3)

      END SELECT

      RETURN
      END SUBROUTINE GETN
!--------------------------------------------------------------------
!     Returns coordinates of a face (with locl id locId) in the ref bulk elem
      SUBROUTINE GETCFace(eType, eNoN, eNoNb, locId, xi)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eType, eNoN, eNoNb, 
     2                                   locId(eNoNb)
      REAL(KIND=RKIND), INTENT(OUT)  :: xi(nsd,eNoNb)
      
      REAL(KIND=RKIND) :: xRefElm(nsd,eNoN)

      IF (eType .EQ. eType_NRB) RETURN

!     2D face elements
      SELECT CASE(eType)
      CASE(eType_TRI3)

         xRefElm(:,1) = (/1._RKIND, 0._RKIND/)
         xRefElm(:,2) = (/0._RKIND, 1._RKIND/)
         xRefElm(:,3) = (/0._RKIND, 0._RKIND/)
!     3D face elements
      CASE(eType_TET4)

         xRefElm(:,1) = (/1._RKIND, 0._RKIND, 0._RKIND/)
         xRefElm(:,2) = (/0._RKIND, 1._RKIND, 0._RKIND/)
         xRefElm(:,3) = (/0._RKIND, 0._RKIND, 1._RKIND/)
         xRefElm(:,4) = (/0._RKIND, 0._RKIND, 0._RKIND/)

      END SELECT

      xi(:,:) = xRefElm(:, locId(:))

      RETURN
      END SUBROUTINE GETCFace
!--------------------------------------------------------------------
!     Returns quad points in current coordinates (current face, dimFace) 
!     from points (nbr of point = nbr of quad points) in the reference 
!     face (dimFace-1) 
      SUBROUTINE GETP_BDELM(insd, eType, eNoNF, nG, crdFace, xin, xout)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: insd, eType, nG, eNoNF
      REAL(KIND=RKIND), INTENT(IN) :: crdFace(insd+1, eNoNF), 
     2              xin(insd, nG)
      REAL(KIND=RKIND), INTENT(OUT) :: xout(insd+1, nG)

      REAL(KIND=RKIND) :: in, jq, N(eNoNF,nG)

C       write(*,*)" face sub elemem = ", crdFace

      DO in = 1, eNoNF
!        Get nodal shape function in the reference face element 
         CALL GETN(insd, eType, eNoNF, xin(:,in), N(:,in))
      END DO

!     N_ij is the weight that we will use to define the quad point 
!     in the new subElement in the reference bulk element with 
!     coordinates = crdFace
      xout = 0._RKIND

!     Loop over quad point          
      DO jq = 1, nG
!        Loop over element nodes 
         DO in = 1, eNoNF
            xout(:, jq) = xout(:, jq) + crdFace(:,in)*N(in,jq)
         END DO
      END DO  

C       write(*,*)" quad point sub elemem ", xout

      RETURN
      END SUBROUTINE GETP_BDELM
!####################################################################
!####################################################################

      SUBROUTINE CONSTRUCT_FSI_MM(lM, Ag, Yg, Dg, iM, iter)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(IN) :: iM, iter

      LOGICAL :: vmsStab
      INTEGER(KIND=IKIND) a, e, g, l, Ac, eNoN, cPhys, iFn, nFn, AcLc
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)
      TYPE(fsType) :: fs(2)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   dl(:,:), bfl(:,:), fN(:,:), pS0l(:,:), pSl(:), ya_l(:),
     3   lR(:,:), lK(:,:,:), lKd(:,:,:), ylFg(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xwl(:,:), xql(:,:), Nwx(:,:),
     2   Nwxx(:,:), Nqx(:,:)

      eNoN = lM%eNoN
      nFn  = lM%nFn
      IF (nFn .EQ. 0) nFn = 1

      IF (lM%nFs .EQ. 1) THEN
         vmsStab = .TRUE.
      ELSE
         vmsStab = .FALSE.
      END IF

!     l = 3, if nsd==2 ; else 6;
      l = nsymd

      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), bfl(nsd,eNoN), fN(nsd,nFn), pS0l(nsymd,eNoN),
     3   pSl(nsymd), ya_l(eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN),
     4   lKd(dof*nsd,eNoN,eNoN))

      IF( iM .GE. 2) ALLOCATE(ylFg(nsd+1,eNoN)) 

!     Loop over all elements of mesh
      DO e=1, lM%nEl

!        Jump alloc if elem is hidden for bh mesh 
C          IF ( iter .GE. 2) THEN 
            IF((iM .EQ. 1) .AND. (lM%lstHdnElm(e) .EQ. 1)) CYCLE
C          END IF

         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF ((cPhys .NE. phys_fluid)  .AND.
     2       (cPhys .NE. phys_lElas)  .AND.
     3       (cPhys .NE. phys_struct) .AND.
     4       (cPhys .NE. phys_ustruct)) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Create local copies
         fN   = 0._RKIND
         pS0l = 0._RKIND
         ya_l = 0._RKIND
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)   = Ac
            xl(:,a)  = x(:,Ac)
            al(:,a)  = Ag(:,Ac)
            yl(:,a)  = Yg(:,Ac)
            dl(:,a)  = Dg(:,Ac)
            bfl(:,a) = Bf(:,Ac)
            IF (ALLOCATED(lM%fN)) THEN
               DO iFn=1, nFn
                  fN(:,iFn) = lM%fN((iFn-1)*nsd+1:iFn*nsd,e)
               END DO
            END IF
            IF (ALLOCATED(pS0)) pS0l(:,a) = pS0(:,Ac)
            IF (cem%cpld) ya_l(a) = cem%Ya(Ac)
         END DO

!        Create local copies if we are in the FG mesh 
         IF( iM .GE. 2) THEN 
            DO a=1, eNoN
               Ac = lM%IEN(a,e)
               AcLc = lM%lN(Ac)

C                ylFg(1:nsd,a) = Yg(1:nsd,Ac)
C                ylFg(nsd+1,a) = lM%YgFG(nsd+1,AcLc)

               ylFg(:,a) = lM%YgFG(:,AcLc)
            END DO
         END IF

!        For FSI, fluid domain should be in the current configuration
         IF (cPhys .EQ. phys_fluid) THEN
            xl(:,:) = xl(:,:) + dl(nsd+2:2*nsd+1,:)
         END IF

!        Initialize residue and tangents
         lR  = 0._RKIND
         lK  = 0._RKIND
         lKd = 0._RKIND

!        Set function spaces for velocity and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 1)

!        Define element coordinates appropriate for function spaces
         ALLOCATE(xwl(nsd,fs(1)%eNoN), Nwx(nsd,fs(1)%eNoN),
     2      Nwxx(l,fs(1)%eNoN))
         ALLOCATE(xql(nsd,fs(2)%eNoN), Nqx(nsd,fs(2)%eNoN))
         xwl(:,:) = xl(:,:)
         xql(:,:) = xl(:,1:fs(2)%eNoN)
         Nwx      = 0._RKIND
         Nqx      = 0._RKIND
         Nwxx     = 0._RKIND

!        Gauss integration 1
         DO g=1, fs(1)%nG
            IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
               CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF

            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

               CALL GNNxx(l, fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g),
     2            fs(1)%Nxx(:,:,g), xwl, Nwx, Nwxx)
            END IF
            w = fs(1)%w(g) * Jac

            IF (nsd .EQ. 3) THEN
               SELECT CASE (cPhys)
               CASE (phys_fluid)
                  CALL FLUID3D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               CASE (phys_lElas)
                  CALL LELAS3D(fs(1)%eNoN, w, fs(1)%N(:,g), Nwx, al, dl,
     2               bfl, pS0l, pSl, lR, lK)

               CASE (phys_struct)
                  CALL STRUCT3D(fs(1)%eNoN, nFn, w, fs(1)%N(:,g), Nwx,
     2               al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK)

               CASE (phys_ustruct)
                  CALL USTRUCT3D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, nFn,
     2               w, Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, al, yl,
     3               dl, bfl, fN, ya_l, lR, lK, lKd)

               END SELECT

            ELSE IF (nsd .EQ. 2) THEN
               SELECT CASE (cPhys)
               CASE (phys_fluid)
                  CALL FLUID2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

C                   IF( iM .GE. 2) THEN 
C                      CALL IFEM_2DRES(fs(1)%eNoN, fs(2)%eNoN, w,
C      2               fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, yl, ylFg, lR)
C                   END IF

               CASE (phys_lElas)
                  CALL LELAS2D(fs(1)%eNoN, w, fs(1)%N(:,g), Nwx, al, dl,
     2               bfl, pS0l, pSl, lR, lK)

               CASE (phys_struct)
                  CALL STRUCT2D(fs(1)%eNoN, nFn, w, fs(1)%N(:,g), Nwx,
     2               al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK)

               CASE (phys_ustruct)
                  CALL USTRUCT2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, nFn,
     2               w, Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, al, yl,
     3               dl, bfl, fN, ya_l, lR, lK, lKd)

C                   IF( iM .GE. 2) THEN 
C                      CALL IFEM_2DRES(fs(1)%eNoN, fs(2)%eNoN, w,
C      2               fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, yl, ylFg, lR)
C                   END IF

               END SELECT
            END IF
         END DO ! g: loop

!        Set function spaces for velocity and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 2)

!        Gauss integration 2
         DO g=1, fs(2)%nG
            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF

            IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
               CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = fs(2)%w(g) * Jac

            IF (nsd .EQ. 3) THEN
               SELECT CASE (cPhys)
               CASE (phys_fluid)
                  CALL FLUID3D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               CASE (phys_ustruct)
                  CALL USTRUCT3D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl,
     3               dl, bfl, lR, lK, lKd)
               END SELECT

            ELSE IF (nsd .EQ. 2) THEN
               SELECT CASE (cPhys)
               CASE (phys_fluid)
                  CALL FLUID2D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               CASE (phys_ustruct)
                  CALL USTRUCT2D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl,
     3               dl, bfl, lR, lK, lKd)
               END SELECT
            END IF
         END DO ! g: loop

         DEALLOCATE(xwl, xql, Nwx, Nwxx, Nqx)

!        Assembly
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            IF (cPhys .EQ. phys_ustruct) err = "Cannot assemble "//
     2         "USTRUCT using Trilinos"
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            IF (cPhys .EQ. phys_ustruct) THEN
               CALL USTRUCT_DOASSEM(eNoN, ptr, lKd, lK, lR)
            ELSE
               CALL DOASSEM(eNoN, ptr, lK, lR)
            END IF
#ifdef WITH_TRILINOS
         END IF
#endif
      END DO ! e: loop

      DEALLOCATE(ptr, xl, al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK,
     2   lKd)

      CALL DESTROY(fs(1))
      CALL DESTROY(fs(2))

      RETURN
      END SUBROUTINE CONSTRUCT_FSI_MM


!####################################################################
!####################################################################
!     Find closest lM node for each node of the foreground mesh lMFg, TODO for parallel
      SUBROUTINE IFEM_FINDCLOSEST(lMBg, lMFg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)


      INTEGER(KIND=IKIND) :: a, e, i, j, iM, Ac, Acn, find
      REAL(KIND=RKIND), ALLOCATABLE :: xFgCur(:,:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: elmBg(:), nodeBg(:)
      REAL(KIND=RKIND), ALLOCATABLE :: poly(:,:)
      REAL(KIND=RKIND) :: xp(nsd), xs(nsd), minb(nsd), maxb(nsd)
      REAL(KIND=RKIND) :: diam, maxDist
      LOGICAL :: flag = .FALSE.
      
      find = 0      


!     Update current foreground (from Lag to Eul) and fluid location (in case of ALE)
      ALLOCATE(xFgCur(nsd,lMFg%nNo))
      ALLOCATE(elmBg(lMBg%nNo))
      ALLOCATE(nodeBg(lMBg%nNo))
      ALLOCATE(poly(nsd,lMBg%eNoN))

      xFgCur = 0._RKIND
      elmBg = 0
      nodeBg = 0



      maxDist = TINY(maxDist)
!     compute solid mesh diamter (in the current configuration)
      DO e=1, lMFg%nEl

         Ac = lMFg%IEN(1,e)
         Acn = lMFg%IEN(2,e)
         diam = DIST(x(:,Ac),x(:,Acn))
         IF(diam .GT. maxDist) maxDist = diam

         Ac = lMFg%IEN(2,e)
         Acn = lMFg%IEN(3,e)
         diam = DIST(x(:,Ac),x(:,Acn))
         IF(diam .GT. maxDist) maxDist = diam

         Ac = lMFg%IEN(1,e)
         Acn = lMFg%IEN(3,e)
         diam = DIST(x(:,Ac),x(:,Acn))
         IF(diam .GT. maxDist) maxDist = diam

      END DO


      DO a=1, lMFg%nNo
         Ac = lMFg%gN(a)
         xFgCur(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
      END DO

!     Create a bounding box around of the current solid location 
      minb = HUGE(minb)
      maxb = TINY(maxb)

      DO i=1, nsd
         minb(i) = MINVAL(xFgCur(i,:)) - maxDist
         maxb(i) = MAXVAL(xFgCur(i,:)) + maxDist
      END DO

!     loop over the background mesh to search if any of the 
!     foreground node is inside 
      DO e=1, lMBg%nEl ! fluid elem id
         flag = .FALSE.

         DO a=1, lMBg%eNoN
            Ac = lMBg%IEN(a,e)

            xp(:) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)

!           is the fluid element in the solid BBox? 
            IF((xp(1) .LE. maxb(1)) .AND. (xp(1) .GE. minb(1)) 
     2                 .AND. (xp(2) .LE. maxb(2)) 
     3                 .AND. (xp(2) .GE. minb(2))) THEN 

!               The node is inside the Bounding Box
                flag = .TRUE.
                EXIT
            END IF 
         END DO

         IF (flag) THEN
!               Loop over the solid node, and search if is inside the fluid element 
C                 write(*,*) "Element ", e, " inside BBox, with vertices: "

            DO a=1, lMBg%eNoN
               Ac = lMBg%IEN(a,e)
               poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
C                   write(*,*) poly(:,a)
            END DO  

            DO a=1, lMFg%nNo
               xs = xFgCur(:,a)

               find = IN_POLY(xs,poly)
C                   write(*,*)"find node solid ", a, " is ", find

!                 is the solid node in this fluis element 
               IF (find .EQ. 1) THEN 
                  elmBg(a) = e 
!                    If inside, serach for the closest point   
!                    local (element) id of the closest point
                  nodeBg(a) = CLOSEST_POINT(xs,poly) 
                  nodeBg(a) = lMBg%IEN(nodeBg(a),e) 
C                      write(*,*)"Closes id is ", msh(iM)%IEN(nodeBg(a),e)
C                      EXIT
               END IF
            END DO
         END IF

      END DO

      IF( ALLOCATED(lMFg%clsBgElm) ) DEALLOCATE(lMFg%clsBgElm)
      IF( ALLOCATED(lMFg%clsBgNd) ) DEALLOCATE(lMFg%clsBgNd)

      ALLOCATE(lMFg%clsBgElm(lMFg%nNo))
      ALLOCATE(lMFg%clsBgNd(lMFg%nNo))

      lMFg%clsBgNd = nodeBg
      lMFg%clsBgElm = elmBg

!     ---------------------------
!     Printing 
C       DO a=1, lMFg%nNo
C          write(*,*)"closest foreground mesh id ", a,  
C      2    " is node id ", lMFg%clsBgNd(a), " in elem ", 
C      3    lMFg%clsBgElm(a)
C       END DO

C        write(*,*) "fg nodes are into bg elem: ", lMFg%clsBgElm 
C       write(*,*) "closest local id is : ", lMFg%clsBgNd
!     ---------------------------

      DEALLOCATE(xFgCur)
      DEALLOCATE(elmBg)
      DEALLOCATE(nodeBg)
      DEALLOCATE(poly)

      RETURN
      END SUBROUTINE IFEM_FINDCLOSEST

!####################################################################

!     Find all hidden nodes of background mesh wrt foreground mesh, TODO for parallel
      SUBROUTINE IFEM_FINDHIDDEN_ND(lMBg, lMFg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

      INTEGER(KIND=IKIND) :: a, e, i, Ac, aBg, find, nbrHid, count
      INTEGER(KIND=IKIND) :: FlagBg(lMBg%nNo), FlagBgElm(lMBg%nEl)
      REAL(KIND=RKIND) :: poly(nsd,lMBg%eNoN), xFgCur(nsd,lMFg%nNo)
      REAL(KIND=RKIND) :: xp(nsd), minb(nsd), maxb(nsd)
      REAL(KIND=RKIND) :: diam, maxDist
      LOGICAL :: flag = .FALSE.
      

!     Create a list of flag of the FlagBg
      FlagBg = 0

!     Create a BBox around foreground mesh 
      DO a=1, lMFg%nNo
         Ac = lMFg%gN(a)
         xFgCur(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
      END DO

!     Create a bounding box around of the current solid location 
      minb = HUGE(minb)
      maxb = TINY(maxb)

      CALL GETMESHDIAM(lMFg)

      DO i=1, nsd
         minb(i) = MINVAL(xFgCur(i,:)) - lMFg%diam 
         maxb(i) = MAXVAL(xFgCur(i,:)) + lMFg%diam
      END DO

!     Finf hidden nodes 
!     For each node of the background mesh nMesh

      DO aBg = 1, lMBg%nNo
         Ac = lMBg%gN(aBg)
         xp(:) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)

!        we check if it is inside the augmented bbox, if no, for sure is outside 
!        if it's inside, we serach an element of the foreground mesh 
!        that contains the node 

         IF((xp(1) .LE. maxb(1)) .AND. (xp(1) .GE. minb(1)) 
     2                 .AND. (xp(2) .LE. maxb(2)) 
     3                 .AND. (xp(2) .GE. minb(2))) THEN 
!           The node is inside the Bounding Box
            flag = .TRUE.
         END IF 


         IF (flag) THEN

            DO e=1, lMFg%nEl

                DO a=1, lMFg%eNoN
                    Ac = lMFg%IEN(a,e)
                    poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
                END DO  


                find = IN_POLY(xp,poly)

!               is the solid node in this fluis element 
                IF (find .EQ. 1) THEN 
                    FlagBg(aBg) = 1 
C                     write(*,*)" node BG ", aBg, " is hidden"
                    EXIT
                END IF
            END DO
         END IF
      END DO
C       write(*,*)" FlagBg = ", FlagBg

      nbrHid = SUM(FlagBg)

      IF( ALLOCATED( lMBg%lstHdnNd )) DEALLOCATE(lMBg%lstHdnNd)
      ALLOCATE(lMBg%lstHdnNd(nbrHid))
      count = 0 

      DO aBg = 1, lMBg%nNo
        IF( FlagBg(aBg) .EQ. 1 ) THEN 
            count = count + 1
            lMBg%lstHdnNd(count) = lMBg%gN(aBg)
        END IF
      END DO


      RETURN
      END SUBROUTINE IFEM_FINDHIDDEN_ND

!----------------------------------------------------
!     Find in which element (of the FG mesh) each hidden node of the BG mesh is 
!     Find in which element (local ordering of the FG mesh)
      SUBROUTINE IFEM_FIND_Nd2Elm_BG(lMBg, lMFg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

      INTEGER(KIND=IKIND) :: i, e, a, find, aBg, Ac, idBG, nbrHnd
      REAL(KIND=RKIND) :: xp(nsd), poly(nsd,lMFg%eNoN)

C       nbrHnd = SIZE(lMBg%lstDirNd)
      nbrHnd = SIZE(lMBg%lstHdnNd)

      IF( ALLOCATED(lMBg%lstNdToELm)) DEALLOCATE(lMBg%lstNdToELm)
      ALLOCATE(lMBg%lstNdToELm(nbrHnd))
      lMBg%lstNdToELm = 0

      DO i = 1, nbrHnd

!       Get global id of the hidden node 
C         idBG = lMBg%lstDirNd(i) 
        idBG = lMBg%lstHdnNd(i) 
        xp(:) = x(:,idBG) + lD(nsd+2:2*nsd+1,idBG)

!       Loop over FG mesh to find the belonging element 
        DO e = 1, lMFg%nEl

            DO a=1, lMFg%eNoN
                Ac = lMFg%IEN(a,e)
                poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
            END DO  

            find = IN_POLY(xp,poly)

!           is the node in this fluid element?
            IF (find .EQ. 1) THEN 
                lMBg%lstNdToELm(i) = e
                EXIT
            END IF
        END DO

      END DO

      RETURN
      END SUBROUTINE IFEM_FIND_Nd2Elm_BG

!----------------------------------------------------
!     Find in which element (of the BG mesh) each node of the FG mesh is 
      SUBROUTINE IFEM_FIND_Nd2Elm_FG(lMBg, lMFg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

      INTEGER(KIND=IKIND) :: i, e, a, find, aBg, Ac, idFG
      REAL(KIND=RKIND) :: xp(nsd), poly(nsd,lMBg%eNoN)

      IF( ALLOCATED(lMFg%lstNdToELm)) DEALLOCATE(lMFg%lstNdToELm)
      ALLOCATE(lMFg%lstNdToELm(lMFg%nNo))
      lMFg%lstNdToELm = 0

      DO i = 1, lMFg%nNo

!       Get global id of the hidden node 
        idFG = lMFg%gN(i) 
        xp(:) = x(:,idFG) + lD(nsd+2:2*nsd+1,idFG)

!       Loop over BG mesh to find the belonging element 
        DO e = 1, lMBg%nEl

            DO a=1, lMBg%eNoN
                Ac = lMBg%IEN(a,e)
                poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
            END DO  

            find = IN_POLY(xp,poly)

!           is the node in this fluid element?
            IF (find .EQ. 1) THEN 
                lMFg%lstNdToELm(i) = e
                EXIT
            END IF
        END DO

      END DO

      RETURN
      END SUBROUTINE IFEM_FIND_Nd2Elm_FG

!----------------------------------------------------
!     Find hidden nodes of background mesh wrt foreground mesh, TODO for parallel
      SUBROUTINE IFEM_FINDHIDDEN_ELM_FLUID(lMBg, lMFg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

      INTEGER(KIND=IKIND) :: a, e, i, Ac, nbrHid, count

      nbrHid = SIZE(lMBg%lstHdnNd)
!     Find hidden element 

      IF( ALLOCATED( lMBg%lstHdnElm )) DEALLOCATE(lMBg%lstHdnElm)
      ALLOCATE( lMBg%lstHdnElm(lMBg%nEl) )
      lMBg%lstHdnElm = 0

!     For each element of the background mesh nMesh
      DO e = 1, lMBg%nEl

!        Loop on the nodes and check if they are all hidden, if yes, that's an hidden elm 
         count = 0 

         DO a=1, lMBg%eNoN
            Ac = lMBg%IEN(a,e)
            DO i = 1, nbrHid
               IF( lMBg%lstHdnNd(i) .EQ. Ac) count = count + 1
            END DO 
         END DO

         IF( count .EQ. lMBg%eNoN ) lMBg%lstHdnElm(e) = 1
C          IF( count .EQ. lMBg%eNoN ) THEN 
C             write(*,*)" lMBg%lstHdnElm(", e, ") = 1"
C          END IF

      END DO

      RETURN
      END SUBROUTINE IFEM_FINDHIDDEN_ELM_FLUID

!----------------------------------------------------
!     Find hidden nodes of background mesh wrt foreground mesh, TODO for parallel
      SUBROUTINE IFEM_FINDHIDDEN_ELM_FS_DIR(lMBg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

!     Find fluid node that are inside the bbxom simply, not also inside any element of the fluid mesh 

      INTEGER(KIND=IKIND) :: a, e, i, Ac, aBg, find, nbrHid, count, AcL
      INTEGER(KIND=IKIND) :: FlagBg(lMBg%nNo), FlagBgElm(lMBg%nEl)
      INTEGER(KIND=IKIND) :: FlagBgDIR(lMBg%nNo), nbrDir
      REAL(KIND=RKIND) :: poly(nsd,lMBg%eNoN), xFgCur(nsd,msh(1)%nNo)
      REAL(KIND=RKIND) :: xp(nsd), minb(nsd), maxb(nsd)
      REAL(KIND=RKIND) :: diam, maxDist
      LOGICAL :: flag = .FALSE.
      

!     Create a list of flag of the FlagBg
      FlagBg = 0
      FlagBgDIR = 0

!     Create a BBox around foreground mesh 
      DO a=1, msh(2)%nNo
         Ac = msh(2)%gN(a)
         xFgCur(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
      END DO

      CALL GETMESHDIAM(msh(2))

!     Create a bounding box around of the current solid location 
      minb = HUGE(minb)
      maxb = TINY(maxb)

      DO i=1, nsd
         minb(i) = MINVAL(xFgCur(i,:)) - msh(2)%diam 
         maxb(i) = MAXVAL(xFgCur(i,:)) + msh(2)%diam
      END DO

!     Finf hidden nodes 
!     For each node of the background mesh nMesh
      DO aBg = 1, lMBg%nNo
         Ac = lMBg%gN(aBg)
         xp(:) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)


!        we check if it is inside the augmented bbox, if no, for sure is outside 
!        if it's inside, we serach an element of the foreground mesh 
!        that contains the node 

         IF((xp(1) .LE. maxb(1)) .AND. (xp(1) .GE. minb(1)) 
     2                 .AND. (xp(2) .LE. maxb(2)) 
     3                 .AND. (xp(2) .GE. minb(2))) THEN 
!           The node is inside the Bounding Box
            flag = .TRUE.
         END IF 


         IF (flag) THEN

!           Add also mesh 2 to compute hidden 
C             DO e=1, msh(2)%nEl

C                 DO a=1, msh(2)%eNoN
C                     Ac = msh(2)%IEN(a,e)
C                     poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
C                 END DO  


C                 find = IN_POLY(xp,poly)

C !               is the solid node in this fluis element 
C                 IF (find .EQ. 1) THEN 
C                     FlagBg(aBg) = 1 
C                     EXIT
C                 END IF
C             END DO

C             IF (find .NE. 1) THEN 
        
                DO e=1, msh(3)%nEl
                    DO a=1, msh(3)%eNoN
                        Ac = msh(3)%IEN(a,e)
                        poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
                    END DO  

                    find = IN_POLY(xp,poly)

    !               is the solid node in this fluis element 
                    IF (find .EQ. 1) THEN 
                        FlagBg(aBg) = 1 
                        EXIT
                    END IF
                END DO

C             END IF

         END IF
      END DO
C       write(*,*)" FlagBg = ", FlagBg

      IF( ALLOCATED(lMBg%lstHdnElm)) DEALLOCATE(lMBg%lstHdnElm)
      ALLOCATE( lMBg%lstHdnElm(lMBg%nEl) )
      lMBg%lstHdnElm = 0
      
!     For each element of the background mesh nMesh
      DO e = 1, lMBg%nEl

!        Loop on the nodes and check if they are all hidden, if yes, that's an hidden elm 
         count = 0 

         DO a=1, lMBg%eNoN
            Ac = lMBg%IEN(a,e)
            AcL = lMBg%lN(Ac)

!           Check how many nodes are in the BBox            
            IF( FlagBg(AcL) .EQ. 1) count = count + 1
         END DO

C          IF( count .EQ. lMBg%eNoN ) lMBg%lstHdnElm(e) = 1
         IF( count .GT. 0 ) lMBg%lstHdnElm(e) = 1
         IF( count .GT. 0 ) THEN 
            write(*,*)" lMBg%lstHdnElm(", e, ") = 1"
         END IF

      END DO


!     Find BC dir node
!     Loop over hidden element from fg solid mesh  
!     For each element of the background mesh nMesh
      DO e = 1, lMBg%nEl

         IF( lMBg%lstHdnElm(e) .NE. 1) CYCLE
!        Loop on the nodes and check if they are all hidden, if yes, that's an hidden elm 
         count = 0 

         DO a=1, lMBg%eNoN
            Ac = lMBg%IEN(a,e)
            AcL = lMBg%lN(Ac)

!           Check how many nodes are in the BBox            
            IF( FlagBg(AcL) .EQ. 1) count = count + 1
         END DO

         IF( (count .LT. lMBg%eNoN) .AND. (FlagBg(AcL) .EQ. 0) ) THEN 
            FlagBgDIR(AcL) = 1
            write(*,*)" lMBg%FlagBgDIR AcL ", AcL, " = 1"
         END IF

      END DO

!     Count how many nodes are hidden, if > 0 but < eNon, it's a DIR node 
      nbrDir = SUM(FlagBgDIR)

      IF( ALLOCATED( lMBg%lstDirNd )) DEALLOCATE(lMBg%lstDirNd)
      ALLOCATE(lMBg%lstDirNd(nbrDir))
      count = 0 

      DO aBg = 1, lMBg%nNo
        IF( FlagBgDIR(aBg) .EQ. 1 ) THEN 
            write(*,*)" DIR loc node mesh 1: ", aBg
            count = count + 1
            lMBg%lstDirNd(count) = lMBg%gN(aBg)
        END IF
      END DO


      RETURN
      END SUBROUTINE IFEM_FINDHIDDEN_ELM_FS_DIR

!----------------------------------------------------
!     Find hidden nodes of background mesh wrt foreground mesh, TODO for parallel
      SUBROUTINE IFEM_FINDHIDDEN_ELM_FS(lMBg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

!     Find fluid node that are inside the bbxom simply, not also inside any element of the fluid mesh 

      INTEGER(KIND=IKIND) :: a, e, i, Ac, aBg, find, nbrHid, count, AcL
      INTEGER(KIND=IKIND) :: FlagBg(lMBg%nNo), FlagBgElm(lMBg%nEl)
      INTEGER(KIND=IKIND) :: FlagBgDIR(lMBg%nNo), nbrDir
      REAL(KIND=RKIND) :: poly(nsd,lMBg%eNoN), xFgCur(nsd,msh(1)%nNo)
      REAL(KIND=RKIND) :: xp(nsd), minb(nsd), maxb(nsd)
      REAL(KIND=RKIND) :: diam, maxDist
      LOGICAL :: flag = .FALSE.
      

!     Create a list of flag of the FlagBg
      FlagBg = 0

!     Create a BBox around foreground mesh 
      DO a=1, msh(2)%nNo
         Ac = msh(2)%gN(a)
         xFgCur(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
      END DO

      CALL GETMESHDIAM(msh(2))

!     Create a bounding box around of the current solid location 
      minb = HUGE(minb)
      maxb = TINY(maxb)

      DO i=1, nsd
         minb(i) = MINVAL(xFgCur(i,:)) - msh(2)%diam 
         maxb(i) = MAXVAL(xFgCur(i,:)) + msh(2)%diam
      END DO

!     Find hidden nodes 
!     For each node of the background mesh nMesh
      DO aBg = 1, lMBg%nNo
         Ac = lMBg%gN(aBg)
         xp(:) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)


!        we check if it is inside the augmented bbox, if no, for sure is outside 
!        if it's inside, we serach an element of the foreground mesh 
!        that contains the node 

         IF((xp(1) .LE. maxb(1)) .AND. (xp(1) .GE. minb(1)) 
     2                 .AND. (xp(2) .LE. maxb(2)) 
     3                 .AND. (xp(2) .GE. minb(2))) THEN 
!           The node is inside the Bounding Box
            flag = .TRUE.
         END IF 


         IF (flag) THEN

!           Add also mesh 2 to compute hidden 
C             DO e=1, msh(2)%nEl

C                 DO a=1, msh(2)%eNoN
C                     Ac = msh(2)%IEN(a,e)
C                     poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
C                 END DO  


C                 find = IN_POLY(xp,poly)

C !               is the solid node in this fluis element 
C                 IF (find .EQ. 1) THEN 
C                     FlagBg(aBg) = 1 
C                     EXIT
C                 END IF
C             END DO

C             IF (find .NE. 1) THEN 
        
                DO e=1, msh(3)%nEl
                    DO a=1, msh(3)%eNoN
                        Ac = msh(3)%IEN(a,e)
                        poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
                    END DO  

                    find = IN_POLY(xp,poly)

    !               is the solid node in this fluis element 
                    IF (find .EQ. 1) THEN 
                        FlagBg(aBg) = 1 
                        EXIT
                    END IF
                END DO

C             END IF

         END IF
      END DO
C       write(*,*)" FlagBg = ", FlagBg

      IF( ALLOCATED(lMBg%lstHdnElm)) DEALLOCATE(lMBg%lstHdnElm)
      ALLOCATE( lMBg%lstHdnElm(lMBg%nEl) )
      lMBg%lstHdnElm = 0
      
!     For each element of the background mesh nMesh
      DO e = 1, lMBg%nEl

!        Loop on the nodes and check if they are all hidden, if yes, that's an hidden elm 
         count = 0 

         DO a=1, lMBg%eNoN
            Ac = lMBg%IEN(a,e)
            AcL = lMBg%lN(Ac)

!           Check how many nodes are in the BBox            
            IF( FlagBg(AcL) .EQ. 1) count = count + 1
         END DO

         IF( count .EQ. lMBg%eNoN ) lMBg%lstHdnElm(e) = 1 
!        To use with hiddenDR
C          IF( count .GT. 0 ) lMBg%lstHdnElm(e) = 1

C          IF( count .GT. 0 ) THEN 
         IF( count .EQ. lMBg%eNoN  ) THEN 
            write(*,*)" lMBg%lstHdnElm(", e, ") = 1"
         END IF

      END DO

      RETURN
      END SUBROUTINE IFEM_FINDHIDDEN_ELM_FS

!####################################################################

!     Find hidden nodes of background mesh wrt solid foreground mesh, TODO for parallel
      SUBROUTINE IFEM_FINDHIDDEN_SOL(lMBg, lMFg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

      INTEGER(KIND=IKIND) :: a, e, i, Ac, aBg, find, nbrHid, count
      INTEGER(KIND=IKIND) :: FlagBg(lMBg%nNo)
      REAL(KIND=RKIND) :: poly(nsd,lMBg%eNoN), xFgCur(nsd,lMFg%nNo)
      REAL(KIND=RKIND) :: xp(nsd), minb(nsd), maxb(nsd)
      REAL(KIND=RKIND) :: diam, maxDist
      LOGICAL :: flag = .FALSE.
      

!     Create a list of flag of the FlagBg
      FlagBg = 0

!     Create a BBox around foreground mesh 
      DO a=1, lMFg%nNo
         Ac = lMFg%gN(a)
         xFgCur(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
      END DO

!     Create a bounding box around of the current solid location 
      minb = HUGE(minb)
      maxb = TINY(maxb)

      CALL GETMESHDIAM(lMFg)

      DO i=1, nsd
         minb(i) = MINVAL(xFgCur(i,:)) - lMFg%diam 
         maxb(i) = MAXVAL(xFgCur(i,:)) + lMFg%diam
      END DO

!     For each node of the background mesh nMesh
      DO aBg = 1, lMBg%nNo
         Ac = lMBg%gN(aBg)
         xp(:) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)


!        we check if it is inside the augmented bbox, if no, for sure is outside 
!        if it's inside, we serach an element of the foreground mesh 
!        that contains the node 

         IF((xp(1) .LE. maxb(1)) .AND. (xp(1) .GE. minb(1)) 
     2                 .AND. (xp(2) .LE. maxb(2)) 
     3                 .AND. (xp(2) .GE. minb(2))) THEN 
!           The node is inside the Bounding Box
            flag = .TRUE.
         END IF 


         IF (flag) THEN

            DO e=1, lMFg%nEl

                DO a=1, lMFg%eNoN
                    Ac = lMFg%IEN(a,e)
                    poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
                END DO  


                find = IN_POLY(xp,poly)

!               is the solid node in this fluis element 
                IF (find .EQ. 1) THEN 
                    FlagBg(aBg) = 1 
C                     write(*,*)" Found bg node ", aBg, " in fr elm ", e
                    EXIT
                END IF
            END DO
         END IF
      END DO
C       write(*,*)" FlagBg = ", FlagBg

      nbrHid = SUM(FlagBg)

      ALLOCATE(lMBg%lstHdnNdSol(nbrHid))
      count = 0 

      DO aBg = 1, lMBg%nNo
        IF( FlagBg(aBg) .EQ. 1 ) THEN 
            count = count + 1
            lMBg%lstHdnNdSol(count) = lMBg%gN(aBg)
        END IF
      END DO

C       write(*,*)" lMBg%lstHdnNdSol = ", lMBg%lstHdnNdSol

      RETURN
      END SUBROUTINE IFEM_FINDHIDDEN_SOL
!####################################################################

      SUBROUTINE IFEM_FINDMLSW(lMBg, lMFg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

      INTEGER(KIND=IKIND) Ac, a, is, iM, mnS, idFCls, nbrSN, idFStc, nd

      REAL(KIND=RKIND), ALLOCATABLE :: Amls(:,:), Bmls(:,:), Wmls(:)
      REAL(KIND=RKIND), ALLOCATABLE :: Pmls(:,:), PPt(:,:), q(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: QMLS(:,:), Ai(:,:), qt(:,:)
      REAL(KIND=RKIND) :: xlf(nsd), xls(nsd), rnorm, wfunc, sum
      REAL(KIND=RKIND) :: Aaux(1,nsd+1)

      ALLOCATE(Amls(nsd+1,nsd+1), Pmls(nsd+1,1), PPt(nsd+1,nsd+1), 
     2         Ai(nsd+1,nsd+1))

      mnS = SIZE(lMBg%stn%ndStn,2)
      lMFg%maxNbrST = mnS
      ALLOCATE(QMLS(mnS,lMFg%nNo))
      QMLS = 0._RKIND

!     Compute mesh space discr param if not computed yet         
      CALL GETMESHDIAM(lMBg)

!     Loop over the solid nodes 
      DO a=1, lMFg%nNo
!        Global id closest bg node
         idFCls = lMFg%clsBgNd(a)  
!        Local id bg node          
         idFCls = lMBg%lN(idFCls)


 !       Nbr of node in stencil
         nbrSN = lMBg%stn%nbrNdStn(idFCls) 
 !       Update current foreground node coord
         Ac = lMFg%gN(a)
C          write(*,*)" node foreground global ", Ac
         xls = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac) 

C          write(*,*)"foreground coord:", xls
C          write(*,*)"id closest point :", idFCls
C          write(*,*)"nmb stencil node", nbrSN 

         ALLOCATE(Bmls(nsd+1,nbrSN))
         ALLOCATE(q(1,nbrSN))
         ALLOCATE(qt(nbrSN,1))
         ALLOCATE(Wmls(nbrSN)) 

         Amls = 0._RKIND
         Bmls = 0._RKIND 

         sum = 0._RKIND
 !           Loop over the stencil 
         DO is = 1, nbrSN
C                write(*,*) "inside loop stencil"
 !          Extract fluid coordinates xlf
            idFStc = lMBg%stn%ndStn(idFCls,is)
            idFStc = lMBg%gN(idFStc)
            xlf = x(:,idFStc) + lD(nsd+2:2*nsd+1,idFStc) 

C                write(*,*)"idFStc = ", idFStc
C                write(*,*)"xlf = ", xlf 

            Pmls(1,1) = 1._RKIND
            Pmls(2,1) = xlf(1)
            Pmls(3,1) = xlf(2) 

C                write(*,*)" Pmls ", Pmls
            
 !              Compute cubic spline value
            rnorm = lMBg%diam*1.3
            Wmls(is) = WHTFUNC(xls,xlf,rnorm)
 !              Check that rnorm is good 
            rnorm = DIST(xls,xlf)
            rnorm = rnorm / (lMBg%diam*1.3)
C                write(*,*)"rnorm is ", rnorm
C                write(*,*)"Wmls is ", Wmls(is) 
 

            PPt = Wmls(is)*MATMUL(Pmls, TRANSPOSE(Pmls))
C                write(*,*)"PPt = ",PPt 

            Amls = Amls + PPt
C                write(*,*)"Matrix A inside is ",Amls 

            Bmls(:,is) = Wmls(is)*Pmls(:,1) 

            sum = sum + Wmls(is)
         END DO 

C             write(*,*)"Sum W = ", sum
         nd = nsd+1
C             write(*,*)"Matrix A is ", Amls
C             write(*,*) ""//""
         Ai = MAT_INV(Amls,nd) 

         Pmls(1,1) = 1._RKIND
         Pmls(2,1) = xls(1)
         Pmls(3,1) = xls(2) 

         Aaux = MATMUL(TRANSPOSE(Pmls),Ai)
         q = MATMUL(Aaux,Bmls)
         qt = TRANSPOSE(q) 

         QMLS(1:nbrSN,a) = qt(:,1)
         
         DEALLOCATE(Bmls)
         DEALLOCATE(q)
         DEALLOCATE(qt)
         DEALLOCATE(Wmls) 

      END DO

      IF( ALLOCATED(lMFg%QMLS) ) DEALLOCATE(lMFg%QMLS) 

      ALLOCATE(lMFg%QMLS(mnS,lMFg%nNo))
      lMFg%QMLS = QMLS

C       write(*,*)" lMFg%QMLS ", lMFg%QMLS

      DEALLOCATE(Amls, Pmls, PPt, Ai, QMLS)

      RETURN
      END SUBROUTINE IFEM_FINDMLSW

!####################################################################

      SUBROUTINE IFEM_EXCHANGE(lA, lY, lYo, lD) 
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(tDof, tnNo), lY(tDof, tnNo)
      REAL(KIND=RKIND), INTENT(INOUT) :: lYo(tDof, tnNo)
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

C       write(*,*)" Calling IFEM_VELPRE_BGtoFG "
      IF ( intrp .EQ. ifemIntrp_MLS ) THEN 
         CALL IFEM_VELPRE_BGtoFG_mls(msh(1), msh(2), lYo, msh(2)%YgFG)
      ELSE 
         CALL IFEM_VELPRE_BGtoFG_fem(msh(1), msh(2),lYo,lD,msh(2)%YgFG)
      END IF 

!     To also do the solid 
C       IF( nMsh .EQ. 3)  THEN 
C          CALL IFEM_VELPRE_BGtoFG(msh(1), msh(3), lYo, msh(3)%YgFG)
C       END IF
C       write(*,*)"  msh(2)%YgFG = ",  msh(2)%YgFG

C       write(*,*)" Calling SETBCDIR_BG "
      !CALL SETBCDIR_BG(lA, lY)

      RETURN
      END SUBROUTINE IFEM_EXCHANGE

!####################################################################
C       SUBROUTINE IFEM_EXCHANGE_WITHDIRBC(lA, lY, lYo) 
C       USE COMMOD
C       USE ALLFUN
C       IMPLICIT NONE
C       REAL(KIND=RKIND), INTENT(INOUT) :: lA(tDof, tnNo), lY(tDof, tnNo)
C       REAL(KIND=RKIND), INTENT(INOUT) :: lYo(tDof, tnNo)

C C       write(*,*)" Calling IFEM_VELPRE_BGtoFG "
C       CALL IFEM_VELPRE_BGtoFG(msh(1), msh(2), lYo, msh(2)%YgFG)
C       IF( nMsh .EQ. 3)  THEN 
C          CALL IFEM_VELPRE_BGtoFG(msh(1), msh(3), lYo, msh(3)%YgFG)
C       END IF
C C       write(*,*)"  msh(2)%YgFG = ",  msh(2)%YgFG

C C       write(*,*)" Calling SETBCDIR_BG "
C       CALL SETBCDIR_BG(lA, lY)

C       RETURN
C       END SUBROUTINE IFEM_EXCHANGE_WITHDIRBC

!####################################################################
      SUBROUTINE IFEM_EXCHANGE_BG(lA, lY, lYo, lD) 
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(tDof, tnNo), lY(tDof, tnNo)
      REAL(KIND=RKIND), INTENT(INOUT) :: lYo(tDof, tnNo)
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

C       write(*,*)" Calling IFEM_EXCHANGE_BG "

      INTEGER(KIND=IKIND) :: nbrHid, i, idHdGl, idHdLc
      REAL(KIND=RKIND) :: Vg(tDof,msh(1)%nNo), tet 

      tet = 1._RKIND

C       nbrHid = SIZE(msh(1)%lstDirNd)
      nbrHid = SIZE(msh(1)%lstHdnNd)

      write(*,*)" Calling IFEM_VEL_FGtoBG "
      IF( intrp .EQ. ifemIntrp_MLS) THEN 
         CALL IFEM_UNK_FGtoBG_mls(msh(1), msh(2), lY, Vg)
      ELSE
         CALL IFEM_UNK_FGtoBG_fem(msh(1), msh(2), lY, lD, Vg)
      END IF
C       write(*,*)" msh(1)%YgBG = ", msh(1)%YgBG
      ! put the bc into lY
      DO i = 1, nbrHid
C          idHdGl = msh(1)%lstDirNd(i)
         idHdGl = msh(1)%lstHdnNd(i)
         idHdLc = msh(1)%lN(idHdGl)

C          write(*,*)" local id hidden ", idHdLc
C          write(*,*)" global id hidden ", idHdGl

C          lY(:, idHdGl) = Vg(:,idHdLc) ! with pressure 
C          lY(1:nsd, idHdGl) = Vg(1:nsd,idHdLc) ! just vel
         lY(1:nsd, idHdGl) = (1.0_RKIND - tet) * lY(1:nsd, idHdGl) 
     2                                + tet * Vg(1:nsd,idHdLc) ! relaxed 
      END DO
C       write(*,*)" Yn after ", lY

      IF( intrp .EQ. ifemIntrp_MLS) THEN 
         CALL IFEM_UNK_FGtoBG_mls(msh(1), msh(2), lA, Vg)
      ELSE
         CALL IFEM_UNK_FGtoBG_fem(msh(1), msh(2), lA, lD, Vg)
      END IF
C       write(*,*)" msh(1)%YgBG = ", msh(1)%YgBG
      !  put the bc into lA
      DO i = 1, nbrHid
C          idHdGl = msh(1)%lstDirNd(i)
         idHdGl = msh(1)%lstHdnNd(i)
         idHdLc = msh(1)%lN(idHdGl)

C          write(*,*)" local id hidden ", idHdLc
C          write(*,*)" global id hidden ", idHdGl

C          lA(:, idHdGl) = Vg(:,idHdLc) 
         lA(1:nsd, idHdGl) = (1.0_RKIND - tet) * lA(1:nsd, idHdGl) 
     2                                + tet * Vg(1:nsd,idHdLc) ! relaxed 
      END DO



      RETURN
      END SUBROUTINE IFEM_EXCHANGE_BG
!####################################################################
!####################################################################

      SUBROUTINE SETBCDIR_BG(lA, lY) 
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(tDof, tnNo), lY(tDof, tnNo)

      INTEGER(KIND=IKIND) :: nbrHid, i, idHdGl, idHdLc

!---- From fluid 
      nbrHid = SIZE(msh(1)%lstHdnNd)

C       write(*,*)" Yn before ", lY
C       write(*,*)" Calling IFEM_VEL_FGtoBG "
      CALL IFEM_VEL_FGtoBG(msh(1), msh(2), lY, msh(1)%YgBG)
C       write(*,*)" msh(1)%YgBG = ", msh(1)%YgBG
      ! put the bc into lY
      DO i = 1, nbrHid
         idHdGl = msh(1)%lstHdnNd(i)
         idHdLc = msh(1)%lN(idHdGl)

C          write(*,*)" local id hidden ", idHdLc
C          write(*,*)" global id hidden ", idHdGl

         lY(1:nsd, idHdGl) = msh(1)%YgBG(:,idHdLc)
      END DO
C       write(*,*)" Yn after ", lY

      CALL IFEM_VEL_FGtoBG(msh(1), msh(2), lA, msh(1)%YgBG)
C       write(*,*)" msh(1)%YgBG = ", msh(1)%YgBG
      !  put the bc into lA
      DO i = 1, nbrHid
         idHdGl = msh(1)%lstHdnNd(i)
         idHdLc = msh(1)%lN(idHdGl)

C          write(*,*)" local id hidden ", idHdLc
C          write(*,*)" global id hidden ", idHdGl

         lA(1:nsd, idHdGl) = msh(1)%YgBG(:,idHdLc)
      END DO

      IF( nMsh .EQ. 3) THEN 
!---- From Solid 
          nbrHid = SIZE(msh(1)%lstHdnNdSol)

C       write(*,*)" Yn before ", lY
C       write(*,*)" Calling IFEM_VEL_FGtoBG "
          CALL IFEM_VEL_FGtoBG(msh(1), msh(3), lY, msh(1)%YgBG)
C       write(*,*)" msh(1)%YgBG = ", msh(1)%YgBG
          ! put the bc into lY
          DO i = 1, nbrHid
             idHdGl = msh(1)%lstHdnNdSol(i)
             idHdLc = msh(1)%lN(idHdGl)

C          write(*,*)" local id hidden ", idHdLc
C          write(*,*)" global id hidden ", idHdGl

             lY(1:nsd, idHdGl) = msh(1)%YgBG(:,idHdLc)
          END DO
C       write(*,*)" Yn after ", lY

          CALL IFEM_VEL_FGtoBG(msh(1), msh(3), lA, msh(1)%YgBG)
C       write(*,*)" msh(1)%YgBG = ", msh(1)%YgBG
          !  put the bc into lA
          DO i = 1, nbrHid
             idHdGl = msh(1)%lstHdnNdSol(i)
             idHdLc = msh(1)%lN(idHdGl)

C          write(*,*)" local id hidden ", idHdLc
C          write(*,*)" global id hidden ", idHdGl

             lA(1:nsd, idHdGl) = msh(1)%YgBG(:,idHdLc)
          END DO   
      END IF   

      RETURN
      END SUBROUTINE SETBCDIR_BG

!####################################################################

      SUBROUTINE IFEM_SETBCDIR_FG(lA, lY) 
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(tDof, tnNo), lY(tDof, tnNo)

      INTEGER(KIND=IKIND) :: i, idGl, idLc
      REAL(KIND=RKIND) :: YgBG(tDof, msh(2)%nNo)

      YgBG = 0._RKIND

      CALL IFEM_UNK_BGtoFG(msh(1), msh(2), lY, YgBG)
      ! put the bc into lY
      DO i = 1, msh(2)%nNo
         idGl = msh(2)%gN(i)
         idLc = i

         lY(:, idGl) = YgBG(:,idLc)
      END DO

      CALL IFEM_UNK_BGtoFG(msh(1), msh(2), lA, YgBG)
      ! put the bc into lY
      DO i = 1, msh(2)%nNo
         idGl = msh(2)%gN(i)
         idLc = i

         lA(:, idGl) = YgBG(:,idLc)
      END DO

      RETURN
      END SUBROUTINE IFEM_SETBCDIR_FG

!####################################################################
!     GET VELOCITY FROM FOREGROUND MESH TO BACKGROUND using MLS
      SUBROUTINE IFEM_VEL_FGtoBG(lMBg, lMFg, Yg, Vg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(OUT) :: Vg(nsd,lMBg%nNo)

      INTEGER(KIND=IKIND) a, is, idFgGl, idBgClsGl, idBgClsLc, nbrSN, 
     2                    idFStcLoc, idStcLoc

      Vg = 0._RKIND

!     We define an array with lenght the nbr of nodes  
!     in the backgrounf mesh, but after we will use only the nodes  
!     hidden to impose the BC Dir 

!     Loop over foreground  nodes 
      DO a=1, lMFg%nNo

!        Global Id foreground node
         idFgGl = lMFg%gN(a)         

!        Global id closest fluid node
         idBgClsGl = lMFg%clsBgNd(a) 
!        Local Id, needed because the stencil is locally defined 
         idBgClsLc = lMBg%lN(idBgClsGl)

!        Nbr of node in stencil
         nbrSN = lMBg%stn%nbrNdStn(idBgClsLc) 
!        Loop over the stencil 
         DO is = 1, nbrSN
C           write(*,*) "inside loop stencil"
            idStcLoc = lMBg%stn%ndStn(idBgClsLc,is)

            Vg(:,idStcLoc) = Vg(:,idStcLoc) + 
     2                              lMFg%QMLS(is,a) * Yg(1:nsd,idFgGl)

         END DO
      END DO

      RETURN
      END SUBROUTINE IFEM_VEL_FGtoBG 
!------------------------------------------------
!     GET VELOCITY FROM FOREGROUND MESH TO BACKGROUND using MLS
      SUBROUTINE IFEM_UNK_FGtoBG_mls(lMBg, lMFg, Yg, Vg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(OUT) :: Vg(tDof,lMBg%nNo)

      INTEGER(KIND=IKIND) a, is, idFgGl, idBgClsGl, idBgClsLc, nbrSN, 
     2                    idFStcLoc, idStcLoc

      Vg = 0._RKIND

!     We define an array with lenght the nbr of nodes  
!     in the backgrounf mesh, but after we will use only the nodes  
!     hidden to impose the BC Dir 

!     Loop over foreground  nodes 
      DO a=1, lMFg%nNo

!        Global Id foreground node
         idFgGl = lMFg%gN(a)         

!        Global id closest fluid node
         idBgClsGl = lMFg%clsBgNd(a) 
!        Local Id, needed because the stencil is locally defined 
         idBgClsLc = lMBg%lN(idBgClsGl)

!        Nbr of node in stencil
         nbrSN = lMBg%stn%nbrNdStn(idBgClsLc) 
!        Loop over the stencil 
         DO is = 1, nbrSN
C           write(*,*) "inside loop stencil"
            idStcLoc = lMBg%stn%ndStn(idBgClsLc,is)

            Vg(:,idStcLoc) = Vg(:,idStcLoc) + 
     2                              lMFg%QMLS(is,a) * Yg(:,idFgGl)

         END DO
      END DO

      RETURN
      END SUBROUTINE IFEM_UNK_FGtoBG_mls 

!------------------------------------------------
!     GET VELOCITY FROM FOREGROUND MESH TO BACKGROUND using MLS
      SUBROUTINE IFEM_UNK_FGtoBG_fem(lMBg, lMFg, Yg, lD, Vg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(OUT) :: Vg(tDof,lMBg%nNo)

      INTEGER(KIND=IKIND) a, i, nbrHid, Ac, eFg, idHdGl, idHdLc
      REAL(KIND=RKIND) :: wFem(lMBg%eNoN), xl(nsd,lMBg%eNoN), xiCur(nsd)
      REAL(KIND=RKIND) :: xiRef(nsd)

      Vg = 0._RKIND

C       nbrHid = SIZE(lMBg%lstDirNd)
      nbrHid = SIZE(lMBg%lstHdnNd)

!     Loop over hidden background nodes  
      DO i = 1, nbrHid
C          idHdGl = lMBg%lstDirNd(i)
         idHdGl = lMBg%lstHdnNd(i)
         idHdLc = lMBg%lN(idHdGl)

!        Update coordinates FG node
         xiCur(:) = x(:,idHdGl) + lD(nsd+2:2*nsd+1,idHdGl)

         eFg = lMBg%lstNdToELm(i)

         DO a=1, lMFg%eNoN
            Ac = lMFg%IEN(a,eFg)   
            xl(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
         END DO  

!        Get coord point in the parent element 
         CALL GET_PntCurToRef(lMFg, xl, xiCur, xiRef)         

!        Compute fem weights 
         CALL GETN(nsd, lMFg%eType, lMFg%eNoN, xiRef, wFem)

!        Compute unknowns at FG point via fem proj
         DO a=1, lMFg%eNoN
            Ac = lMFg%IEN(a,eFg)
            
            Vg(:,idHdLc) = Vg(:,idHdLc) + wFem(a) * Yg(1:nsd+1,Ac)
         END DO  
      END DO

      RETURN
      END SUBROUTINE IFEM_UNK_FGtoBG_fem

!####################################################################
!     GET VELOCITY AND PRSSURE FROM BACKGROUND MESH 
!     Interpolate data at IFEM nodes from background mesh using MLS
      SUBROUTINE IFEM_VELPRE_BGtoFG_mls(lMBg, lMFg, Yg, Vg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(OUT) :: Vg(nsd+1,lMFg%nNo)


      INTEGER(KIND=IKIND) a, is, idFgGl, idFgLc, idBgClsGl, idBgClsLc,  
     2                    nbrSN, idStcLoc, idStcGl


      Vg = 0._RKIND

!     Loop over foreground  nodes 
      DO idFgLc = 1, lMFg%nNo

!        Global Id foreground node
         idFgGl = lMFg%gN(idFgLc)         

!        Global id closest fluid node
         idBgClsGl = lMFg%clsBgNd(idFgLc) 
!        Local Id, needed because the stencil is locally defined 
         idBgClsLc = lMBg%lN(idBgClsGl)

!        Nbr of node in stencil
         nbrSN = lMBg%stn%nbrNdStn(idBgClsLc) 
!        Loop over the stencil 
         DO is = 1, nbrSN
C           write(*,*) "inside loop stencil"
            idStcLoc = lMBg%stn%ndStn(idBgClsLc,is)
            idStcGl = lMBg%gN(idStcLoc)

            Vg(:,idFgLc) = Vg(:,idFgLc) + 
     2                       lMFg%QMLS(is,idFgLc) * Yg(1:nsd+1,idStcGl)

         END DO
      END DO

      RETURN
      END SUBROUTINE IFEM_VELPRE_BGtoFG_mls
!--------------------------------------------------------------------
!     GET VELOCITY AND PRSSURE FROM BACKGROUND MESH 
!     Interpolate data at IFEM nodes from background mesh using MLS
      SUBROUTINE IFEM_VELPRE_BGtoFG_fem(lMBg, lMFg, Yg, lD, Vg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(OUT) :: Vg(nsd+1,lMFg%nNo)

      INTEGER(KIND=IKIND) a, is, idFgGl, idFgLc, idBgClsGl, idBgClsLc,  
     2                    nbrSN, idStcLoc, idStcGl, Ac, eBg
      REAL(KIND=RKIND) :: wFem(lMFg%eNoN), xl(nsd,lMFg%eNoN), xiCur(nsd)
      REAL(KIND=RKIND) :: xiRef(nsd)


      Vg = 0._RKIND

!     Loop over foreground nodes 
      DO idFgLc = 1, lMFg%nNo

!        Global Id foreground node
         idFgGl = lMFg%gN(idFgLc)         

!        Update coordinates FG node
         xiCur(:) = x(:,idFgGl) + lD(nsd+2:2*nsd+1,idFgGl)

!        Get coordinates of the BG elm containing FG node 
         eBg = lMFg%lstNdToELm(idFgLc)
         DO a=1, lMBg%eNoN
            Ac = lMBg%IEN(a,eBg)
            xl(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
         END DO  

!        Get coord point in the parent element 
         CALL GET_PntCurToRef(lMBg, xl, xiCur, xiRef)         

!        Compute fem weights 
         CALL GETN(nsd, lMBg%eType, lMBg%eNoN, xiRef, wFem)

!        Compute unknowns at FG point via fem proj
         DO a=1, lMBg%eNoN
            Ac = lMBg%IEN(a,eBg)
            
            Vg(:,idFgLc) = Vg(:,idFgLc) + wFem(a) * Yg(1:nsd+1,Ac)
         END DO  
      END DO

      RETURN
      END SUBROUTINE IFEM_VELPRE_BGtoFG_fem
!####################################################################

!     GET VELOCITY AND PRSSURE FROM BACKGROUND MESH 
!     Interpolate data at IFEM nodes from background mesh using MLS
      SUBROUTINE IFEM_UNK_BGtoFG(lMBg, lMFg, Yg, Vg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(OUT) :: Vg(tDof,lMFg%nNo)


      INTEGER(KIND=IKIND) a, is, idFgGl, idFgLc, idBgClsGl, idBgClsLc,  
     2                    nbrSN, idStcLoc, idStcGl


      Vg = 0._RKIND

!     Loop over foreground  nodes 
      DO idFgLc = 1, lMFg%nNo

!        Global Id foreground node
         idFgGl = lMFg%gN(idFgLc)         

!        Global id closest fluid node
         idBgClsGl = lMFg%clsBgNd(idFgLc) 
!        Local Id, needed because the stencil is locally defined 
         idBgClsLc = lMBg%lN(idBgClsGl)

!        Nbr of node in stencil
         nbrSN = lMBg%stn%nbrNdStn(idBgClsLc) 
!        Loop over the stencil 
         DO is = 1, nbrSN
C           write(*,*) "inside loop stencil"
            idStcLoc = lMBg%stn%ndStn(idBgClsLc,is)
            idStcGl = lMBg%gN(idStcLoc)

            Vg(:,idFgLc) = Vg(:,idFgLc) + 
     2                       lMFg%QMLS(is,idFgLc) * Yg(:,idStcGl)

         END DO
      END DO

      RETURN
      END SUBROUTINE IFEM_UNK_BGtoFG

!####################################################################
!####################################################################

      SUBROUTINE GET_PntCurToRef(lM, xl, xiCur, xiRef)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: xl(nsd,lM%eNoN)
      REAL(KIND=RKIND), INTENT(IN) :: xiCur(nsd)
      REAL(KIND=RKIND), INTENT(OUT) :: xiRef(nsd)

      INTEGER(KIND=IKIND) :: g
      REAL(KIND=RKIND) :: Def(nsd,nsd), invDef(nsd,nsd), rhsD(nsd) 

      IF (nsd .EQ. 2) THEN
         Def(1,1) = xl(1,1) - xl(1,3) 
         Def(1,2) = xl(1,2) - xl(1,3) 
         Def(2,1) = xl(2,1) - xl(2,3) 
         Def(2,2) = xl(2,2) - xl(2,3) 
         
      ELSE 
         write(*,*)"****** 3D implementation still TODO *****" 
!           Need to check how the basis function are ordered
!           and rewrite Def and rhsD            
      END IF 

!     Find the coordinates in the current element 
      rhsD(:) = xiCur(:) - xl(:,3) 
      invDef = MAT_INV(Def,nsd)
      xiRef = MATMUL(invDef,rhsD)
C       write(*,*)" xiCur = " , xiCur, "xiRef coord are ", xiRef

      RETURN
      END SUBROUTINE GET_PntCurToRef

!--------------------------------------------------------------------
!####################################################################
!     This routine reads IFEM options
      SUBROUTINE IFEM_READOPTS(list)
      USE COMMOD
      USE LISTMOD
      IMPLICIT NONE
      TYPE(listType), INTENT(INOUT) :: list

      CHARACTER(LEN=stdL) :: ctmp
      TYPE(listType), POINTER :: lPtr

      lPtr => list%get(ctmp, "IFEM interpolation method")
      CALL TO_LOWER(ctmp)
            write(*,*)" ctmp  ", ctmp

      SELECT CASE (ctmp)
      CASE ("fem")
         intrp = ifemIntrp_FEM
         std = " IFEM interpolation method: "//CLR("Finite Element",3)
      CASE ("mls")
         intrp = ifemIntrp_MLS
         std = " IFEM interpolation method: "
     2                     //CLR("Moving Leas-squares",3)
      CASE DEFAULT
         err = " Invalid IFEM interpolation chosen"
      END SELECT

      RETURN
      END SUBROUTINE IFEM_READOPTS
!####################################################################
!####################################################################

      SUBROUTINE PLOT_LOC_GLOB
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

!     Plot info mesh 1 
      write(*,*)" ------- PRINT msh 1 : "
      CALL PRINT_MESH(msh(1))
      write(*,*)"-------"

!     Plot info mesh 2 
      write(*,*)" ------- PRINT msh 2 : "
      CALL PRINT_MESH(msh(2))
      write(*,*)"-------"

!     Plot info mesh 3 
      write(*,*)" ------- PRINT msh 3 : "
      CALL PRINT_MESH(msh(3))
      write(*,*)"-------"

      RETURN
      END SUBROUTINE PLOT_LOC_GLOB
!####################################################################

      SUBROUTINE PRINT_MESH(lM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM

      write(*,*) " gnEl =  ", lM%gnEl
      write(*,*) " gnNo =  ", lM%gnNo
      write(*,*) " nEl =  ", lM%nEl
      write(*,*) " nFa =  ", lM%nFa
      write(*,*) " nNo =  ", lM%nNo

      write(*,*)""

      write(*,*)" gN size  = ", SIZE(lM%gN)
      write(*,*)" gN(:) = ", lM%gN(:)
      write(*,*)""
      write(*,*)" gpN size = ", SIZE(lM%gpN)
      write(*,*)" gpN(:) = ", lM%gpN(:) 
      write(*,*)""
      write(*,*)" lN size = ", SIZE(lM%lN)
      write(*,*)" lN(:) = ", lM%lN(:)
      write(*,*)""

      write(*,*)" Coord elm 1 gIEN "
      write(*,*) x(:,lM%gIEN(1,1))
      write(*,*) x(:,lM%gIEN(2,1))
      write(*,*) x(:,lM%gIEN(3,1))
      write(*,*)" gIEN(:,1) = ", lM%gIEN(:,1)
      write(*,*)""

      write(*,*)" Coord elm 1 IEN "
      write(*,*) x(:,lM%IEN(1,1))
      write(*,*) x(:,lM%IEN(2,1))
      write(*,*) x(:,lM%IEN(3,1))
      write(*,*)" IEN(:,1) = ", lM%IEN(:,1)
      write(*,*)""

      RETURN
      END SUBROUTINE PRINT_MESH

!####################################################################
