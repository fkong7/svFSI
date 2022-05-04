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
         CALL IFEM_FINDMLSW(msh(1), msh(iM), Dg)
      END DO

      CALL IFEM_FINDHIDDEN(msh(1), msh(2), Dg)

!     Allocation of arrays to store unknown in Fg and Bg meshes 
      ALLOCATE(msh(1)%YgBG(nsd,msh(1)%nNo))
      ALLOCATE(msh(2)%YgFG(nsd+1,msh(2)%nNo))

      RETURN
      END SUBROUTINE IFEM_INIT

!####################################################################
!####################################################################

      SUBROUTINE CONSTRUCT_FLUID_MM(lM, Ag, Yg, iM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(IN) :: iM

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

      IF( iM .GE. 2) ALLOCATE(ylFg(nsd+1,eNoN)) 

!     Loop over all elements of mesh
      DO e=1, lM%nEl
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
     2               fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, ylFg, lR)
               END IF

             ELSE IF (nsd .EQ. 2) THEN
               CALL FLUID2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl,
     3            bfl, lR, lK)
               IF( iM .GE. 2) THEN 
                  CALL IFEM_2DRES(fs(1)%eNoN, fs(2)%eNoN, w,
     2               fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, ylFg, lR)
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
     2   Nqx, yl, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Nw(eNoNw), Nq(eNoNq),
     2   Nwx(2,eNoNw), Nqx(2,eNoNq), yl(nsd+1,eNoNw) 
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

         u(1)    = u(1) + Nw(a)*yl(1,a)
         u(2)    = u(2) + Nw(a)*yl(2,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)

      END DO

!     Pressure 
      p  = 0._RKIND
      DO a=1, eNoNq
         p  = p + Nq(a)*yl(3,a)
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
!     Find closest lM node for each node of the foreground mesh lMFg, TODO for parallel
      SUBROUTINE IFEM_FINDCLOSEST(lMBg, lMFg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)


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


      ALLOCATE(lMFg%clsBgElm(lMFg%nNo))
      ALLOCATE(lMFg%clsBgNd(lMFg%nNo))

      lMFg%clsBgNd = nodeBg
      lMFg%clsBgElm = elmBg

!     ---------------------------
!     Printing 
      DO a=1, lMFg%nNo
         write(*,*)"closest foreground mesh id ", a,  
     2    " is node id ", lMFg%clsBgNd(a), " in elem ", 
     3    lMFg%clsBgElm(a)
      END DO

       write(*,*) "fg nodes are into bg elem: ", lMFg%clsBgElm 
      write(*,*) "closest local id is : ", lMFg%clsBgNd
!     ---------------------------

      DEALLOCATE(xFgCur)
      DEALLOCATE(elmBg)
      DEALLOCATE(nodeBg)
      DEALLOCATE(poly)

      RETURN
      END SUBROUTINE IFEM_FINDCLOSEST

!####################################################################

!     Find hidden nodes of background mesh wrt foreground mesh, TODO for parallel
      SUBROUTINE IFEM_FINDHIDDEN(lMBg, lMFg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

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
                    write(*,*)" Found bg node ", aBg, " in fr elm ", e
                    EXIT
                END IF
            END DO
         END IF
      END DO
C       write(*,*)" FlagBg = ", FlagBg

      nbrHid = SUM(FlagBg)

      ALLOCATE(lMBg%lstHdnNd(nbrHid))
      count = 0 

      DO aBg = 1, lMBg%nNo
        IF( FlagBg(aBg) .EQ. 1 ) THEN 
            count = count + 1
            lMBg%lstHdnNd(count) = lMBg%gN(aBg)
        END IF
      END DO

C       write(*,*)" lMBg%lstHdnNd = ", lMBg%lstHdnNd

      RETURN
      END SUBROUTINE IFEM_FINDHIDDEN

!####################################################################

      SUBROUTINE IFEM_FINDMLSW(lMBg, lMFg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lMBg 
      TYPE(mshType), INTENT(INOUT) :: lMFg 
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

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
            rnorm = lMBg%diam*1.2
            Wmls(is) = WHTFUNC(xls,xlf,rnorm)
 !              Check that rnorm is good 
            rnorm = DIST(xls,xlf)
            rnorm = rnorm / (lMBg%diam*1.15)
C                write(*,*)"rnorm is ", rnorm
C                write(*,*)"Wmls is ", Wmls(is) 
 

            PPt = Wmls(is)*MATMUL(Pmls, TRANSPOSE(Pmls))
C                write(*,*)"PPt = ",PPt 

            Amls = Amls + PPt
C                write(*,*)"Matrix A inside is ",Amls 

            Bmls(:,is) = Wmls(is)*Pmls(:,1) 

            sum = sum + Wmls(is)
         END DO 

            write(*,*)"Sum W = ", sum
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

      IF(.NOT.ALLOCATED(lMFg%QMLS)) ALLOCATE(lMFg%QMLS(mnS,lMFg%nNo))
      lMFg%QMLS = QMLS

C       write(*,*)" lMFg%QMLS ", lMFg%QMLS

      DEALLOCATE(Amls, Pmls, PPt, Ai, QMLS)

      RETURN
      END SUBROUTINE IFEM_FINDMLSW

!####################################################################

      SUBROUTINE IFEM_EXCHANGE(lA, lY) 
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(tDof, tnNo), lY(tDof, tnNo)


C       write(*,*)" Calling IFEM_VELPRE_BGtoFG "
      CALL IFEM_VELPRE_BGtoFG(msh(1), msh(2), lY, msh(2)%YgFG)
C       write(*,*)"  msh(2)%YgFG = ",  msh(2)%YgFG


C       write(*,*)" Calling SETBCDIR_BG "
      CALL SETBCDIR_BG(lA, lY)


      RETURN
      END SUBROUTINE IFEM_EXCHANGE

!####################################################################

      SUBROUTINE SETBCDIR_BG(lA, lY) 
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(tDof, tnNo), lY(tDof, tnNo)

      INTEGER(KIND=IKIND) :: nbrHid, i, idHdGl, idHdLc


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

      RETURN
      END SUBROUTINE SETBCDIR_BG

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


!####################################################################
!     GET VELOCITY AND PRSSURE FROM BACKGROUND MESH 
!     Interpolate data at IFEM nodes from background mesh using MLS
      SUBROUTINE IFEM_VELPRE_BGtoFG(lMBg, lMFg, Yg, Vg)
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
      END SUBROUTINE IFEM_VELPRE_BGtoFG
