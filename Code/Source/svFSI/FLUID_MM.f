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

      SUBROUTINE CONSTRUCT_FLUID_MM(lM, Ag, Yg, iM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(IN) :: iM

      LOGICAL :: vmsStab
      INTEGER(KIND=IKIND) a, e, g, l, Ac, eNoN, cPhys
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)
      TYPE(fsType) :: fs(2)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   bfl(:,:), lR(:,:), lK(:,:,:)
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

             ELSE IF (nsd .EQ. 2) THEN
               CALL FLUID2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl,
     3            bfl, lR, lK)
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

      CALL DESTROY(fs(1))
      CALL DESTROY(fs(2))


C !     Additional part for Multimesh strategy with IFEM methodology
C !     If mesh 1: need to clean the velocity row, we will impose strongly  
C !     the velocity after the system resolution 
C       IF( iM .EQ. 1) THEN 
C         CALL APPLY_DIR
C       END IF

C !     If mesh 2, use the given vel and pressure to compute FSI res forse 
C !     and add it to the residue 
C       IF( iM .EQ. 2) THEN 
C         CALL APPLY_NEU
C       END IF


      RETURN
      END SUBROUTINE CONSTRUCT_FLUID_MM
!####################################################################
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
         write(*,*)" idFCls = ", idFCls
!        Local id bg node          
         idFCls = lMBg%lN(idFCls)
         write(*,*)" idFCls = ", idFCls


 !       Nbr of node in stencil
         nbrSN = lMBg%stn%nbrNdStn(idFCls) 
 !       Update current foreground node coord
         Ac = lMFg%gN(a)
         write(*,*)" node foreground global ", Ac
         xls = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac) 

         write(*,*)"foreground coord:", xls
         write(*,*)"id closest point :", idFCls
         write(*,*)"nmb stencil node", nbrSN 

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

      write(*,*)" lMFg%QMLS ", lMFg%QMLS

      DEALLOCATE(Amls, Pmls, PPt, Ai, QMLS)

      RETURN
      END SUBROUTINE IFEM_FINDMLSW





!####################################################################
GET VELOCITY FROM FOREGROUND MESH TO BACKGROUND  

      SUBROUTINE IFEM_RASSEMBLY()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) a, is, iM, idFCls, nbrSN, idFStc
      INTEGER(KIND=IKIND) idFEl, al, Ac

!     TODO for multiple fluid mesh
      iM = 1


!     Definiamo un vettore vel lungo quando il numero di nodiu 
!     in the backgrounf mesh, but we keep the node only if they are 
!     hidden 

!        Loop over the solid nodes 
         DO a=1, ifem%tnNo
C           write(*,*)"inside loop solid node IFEM_CONSTRUCT"
!           Global id closest fluid node
            idFCls = ifem%clsFNd(a) 
!           Nbr of node in stencil
            nbrSN = msh(iM)%stn%nbrNdStn(idFCls) 
!           Loop over the stencil 
            DO is = 1, nbrSN
C              write(*,*) "inside loop stencil"
!              Extract fluid coordinates xlf
               idFStc = msh(iM)%stn%ndStn(idFCls,is)

               R(1:nsd,idFStc) = R(1:nsd,idFStc) + 
     2                              ifem%QMLS(is,a)*ifem%Rsolid(:,a)

            END DO
         END DO



!     OLD Final residual assembly 
C       DO a=1, tnNo
C          R(1:nsd,a) = R(1:nsd,a) + ifem%Rfluid(:,a)
C       END DO

      RETURN
      END SUBROUTINE IFEM_RASSEMBLY 







!####################################################################
GET VELOCITY AND PRSSURE FROM BACKGROUND MESH 

!     Interpolate data at IFEM nodes from background mesh using fem
!     nodal function 
      SUBROUTINE IFEM_FINDSOLVEL_MLS(m, Ug, Dg, Ub)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m
      REAL(KIND=RKIND), INTENT(IN) :: Ug(m,tnNo), Dg(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(OUT) :: Ub(nsd,ifem%tnNo) !new solid vel

      INTEGER(KIND=IKIND) :: a, iM, idFCls, nbrSN, idFStc, is

      Ub = 0._RKIND

      !     TODO for multiple fluid mesh
      iM = 1
!     Loop over the solid nodes 
      DO a=1, ifem%tnNo
C          write(*,*)"inside loop solid node IFEM_CONSTRUCT"
!        Global id closest fluid node
         idFCls = ifem%clsFNd(a) 
!        Nbr of node in stencil
         nbrSN = msh(iM)%stn%nbrNdStn(idFCls) 
!        Loop over the stencil 
         DO is = 1, nbrSN
C             write(*,*) "inside loop stencil"
!           Extract fluid coordinates xlf
            idFStc = msh(iM)%stn%ndStn(idFCls,is)

            Ub(:,a) = Ub(:,a) + ifem%QMLS(is,a)*Ug(1:nsd,idFStc)

         END DO
      END DO

      RETURN
      END SUBROUTINE IFEM_FINDSOLVEL_MLS
