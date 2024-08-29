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
!     Collection of subroutines for the Resistive Immersed Surface 
!     method. 
!
!--------------------------------------------------------------------
!####################################################################
!     This subroutine computes the mean pressure and flux on the 
!     ris surface 
      SUBROUTINE RIS_MEANQ
      USE TYPEMOD
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) :: iEq, nPrj, m, s, e, i, iM, iFa, iProj
      REAL(KIND=RKIND) :: tmp
      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:)

!     Let's first conpute the mean pressure at the boudary on the two meashes

      iEq = 1

      IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
      ALLOCATE (tmpV(maxnsd,tnNo))

      nPrj = RIS%nbrRIS
  
      m = 1
      s = eq(iEq)%s + nsd 
      e = s + m - 1

      tmpV(1:m,:) = Yn(s:e,:)
      
!     Future loop on all the nbgProj nPrj
      DO iProj = 1, nPrj
        RIS%meanP(iProj, :) = 0._RKIND
        RIS%meanFl(iProj) = 0._RKIND
        DO i = 1, 2 ! We always have two meshes 
           iM = RIS%lst(i,1,iProj)
           iFa = RIS%lst(i,2,iProj)

!          ERROR. HERRE FOR ALE, recompure the area with the new displacement 
!          ^FK: But that shouldn't be a problem for pressure difference
!          check since the area (change) above and below a RIS is the
!          same?           
           tmp = msh(iM)%fa(iFa)%area
           RIS%meanP(iProj, i) = Integ(msh(iM)%fa(iFa),tmpV,1)/tmp
        END DO
      END DO

!     For the velocity 
      m = nsd
      s = eq(iEq)%s
      e = s + m - 1

      DO iProj = 1, nPrj
        IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
        ALLOCATE (tmpV(maxnsd,tnNo))
        tmpV(1:m,:) = Yn(s:e,:)
        iM = RIS%lst(1,1,iProj)
        iFa = RIS%lst(1,2,iProj)

        RIS%meanFl(iProj) = Integ(msh(iM)%fa(iFa),tmpV,1,m)
        std = "For RIS projection "//iProj
        std = "    The average pressure is: "//RIS%meanP(iProj,
     2      1)//", "//RIS%meanP(iProj, 2)
        std = "    The average flow is: "//RIS%meanFl(iProj)
      END DO


      IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)

      RETURN
      END SUBROUTINE RIS_MEANQ

!####################################################################
!     Weak treatment of RIS resistance boundary conditions
      SUBROUTINE RIS_RESBC(Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) :: iFa, iM, nPrj, i, cPhys, iProj
      TYPE(bcType) :: lBc


      nPrj = RIS%nbrRIS ! for the moment 


      DO iProj = 1, nPrj
         IF(.NOT.RIS%clsFlg(iProj)) CYCLE
!        Weak Dirichlet BC for fluid/FSI equations
         lBc%weakDir = .TRUE.
         lBc%tauB = RIS%Res(iProj)
         lBc%bType = IBSET(lBc%bType,bType_Dir)
         lBc%bType = IBSET(lBc%bType,bType_std)
         lBc%bType = IBSET(lBc%bType,bType_flat)
         ALLOCATE(lBc%eDrn(nsd))
         lBc%eDrn = 0
         DO i = 1, 2 ! We always have two meshes for 1 projection
            iM = RIS%lst(i,1,iProj)
            iFa = RIS%lst(i,2,iProj)

            ALLOCATE(lBc%gx(msh(iM)%fa(iFa)%nNo))
            lBc%gx = 1._RKIND

            cPhys = eq(cEq)%dmn(cDmn)%phys
            IF (cPhys .EQ. phys_fluid) THEN
!               TO build the correct bc 
!               CALL SETBCRIS(lBc, msh(iM), msh(iM)%fa(iFa), Yg, Dg)
               CALL SETBCDIRWL(lBc, msh(iM), msh(iM)%fa(iFa), Yg, Dg)
            END IF

            DEALLOCATE(lBc%gx)

         END DO
         DEALLOCATE(lBc%eDrn)
      END DO
      

      RETURN
      END SUBROUTINE RIS_RESBC
!####################################################################
!--------------------------------------------------------------------
      SUBROUTINE SETBCRIS(lBc, lM, lFa, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(bcType), INTENT(IN) :: lBc
      TYPE(mshType), INTENT(IN) :: lM
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      LOGICAL :: flag, eDir(maxnsd)
      INTEGER(KIND=IKIND) :: a, e, i, g, Ac, Ec, ss, ee, lDof, nNo, nEl,
     2   eNoN, eNoNb, cPhys
      REAL(KIND=RKIND) :: w, Jac, xp(nsd), xi(nsd), xi0(nsd), nV(nsd),
     2   ub(nsd), tauB(2), Ks(nsd,nsd)
      TYPE(fsType) :: fs(2)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: Nw(:), Nq(:), Nwx(:,:), Nqx(:,:),
     2   Nwxi(:,:), xl(:,:), xbl(:,:), xwl(:,:), xql(:,:), yl(:,:),
     3   ubl(:,:), ubg(:,:), tmpA(:,:), tmpY(:,:), lR(:,:), lK(:,:,:)

      nNo   = lFa%nNo
      nEl   = lFa%nEl
      eNoNb = lFa%eNoN
      eNoN  = lM%eNoN
      tauB  = lBc%tauB

      IF (lM%nFs .EQ. 1) THEN
         flag = .TRUE.
      ELSE
         flag = .FALSE.
      END IF

!     Compute the Dirichlet value to be applied weakly on the face
      ss    = eq(cEq)%s
      ee    = eq(cEq)%e
      IF (eq(cEq)%dof .EQ. nsd+1) ee = ee - 1
      eDir  = .FALSE.
      lDof  = 0
      DO i=1, nsd
         IF (lBc%eDrn(i) .NE. 0) THEN
            eDir(i) = .TRUE.
            lDof = lDof + 1
         END IF
      END DO
      IF (lDof .EQ. 0) lDof = ee - ss + 1

      ALLOCATE(tmpA(lDof,nNo), tmpY(lDof,nNo))
      CALL SETBCDIRL(lBc, lFa, tmpA, tmpY, lDof)
      IF (BTEST(lBc%bType,bType_impD)) tmpY(:,:) = tmpA(:,:)

!     Transfer Dirichlet value to global numbering (lFa%nNo -> tnNo).
!     Take effective direction into account if set.
      ALLOCATE(ubg(nsd,tnNo))
      ubg = 0._RKIND
      IF (ANY(eDir)) THEN
         DO a=1, nNo
            Ac = lFa%gN(a)
            DO i=1, nsd
               lDof = 0
               IF (eDir(i)) THEN
                  lDof = lDof + 1
                  ubg(i,Ac) = tmpY(lDof,a)
               END IF
            END DO
         END DO
      ELSE
         DO a=1, nNo
            Ac = lFa%gN(a)
            ubg(:,Ac) = tmpY(:,a)
         END DO
      END IF
      DEALLOCATE(tmpA, tmpY)

      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), yl(tDof,eNoN), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN))

      ALLOCATE(xbl(nsd,eNoNb), ubl(nsd,eNoNb))

!     Loop over all the elements of the face and construct residue and
!     tangent matrices
      DO e=1, nEl
         Ec    = lFa%gE(e)
         cDmn  = DOMAIN(lM, cEq, Ec)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_fluid) err = "Weakly applied Dirichlet "//
     2      "BC is allowed for fluid phys only"

!        Initialize local residue and stiffness
         lR = 0._RKIND
         lK = 0._RKIND

!        Create local copies of fluid velocity and position vector
         DO a=1, eNoN
            Ac = lM%IEN(a,Ec)
            ptr(a)  = Ac
            xl(:,a) = x(:,Ac)
            IF (mvMsh) xl(:,a) = xl(:,a) + Dg(nsd+2:2*nsd+1,Ac)
            yl(:,a) = Yg(:,Ac)
         END DO

!        Set function spaces for velocity and pressure on mesh
         CALL GETTHOODFS(fs, lM, flag, 1)

         ALLOCATE(xwl(nsd,fs(1)%eNoN), Nw(fs(1)%eNoN),
     2      Nwx(nsd,fs(1)%eNoN), Nwxi(nsd,fs(1)%eNoN))
         ALLOCATE(xql(nsd,fs(2)%eNoN), Nq(fs(2)%eNoN),
     2      Nqx(nsd,fs(2)%eNoN))
         xwl(:,:) = xl(:,:)
         xql(:,:) = xl(:,1:fs(2)%eNoN)

!        Create local copies of the wall/solid/interface quantites
         DO a=1, eNoNb
            Ac = lFa%IEN(a,e)
            xbl(:,a) = x(:,Ac)
            ubl(:,a) = ubg(:,Ac)
            IF (mvMsh) xbl(:,a) = xbl(:,a) + Dg(nsd+2:2*nsd+1,Ac)
         END DO

!        Initialize parameteric coordinate for Newton's iterations.
!        Newton method is used to compute derivatives on the face using
!        mesh-based shape functions as an inverse problem.
         xi0 = 0._RKIND
         DO g=1, fs(1)%nG
            xi0 = xi0 + fs(1)%xi(:,g)
         END DO
         xi0 = xi0 / REAL(fs(1)%nG, KIND=RKIND)

!        Gauss integration 1
         DO g=1, lFa%nG
            CALL GNNB(lFa, e, g, nsd-1, eNoNb, lFa%Nx(:,:,g), nV)
            Jac = SQRT(NORM(nV))
            nV  = nV/Jac
            w   = lFa%w(g) * Jac

            xp = 0._RKIND
            ub = 0._RKIND
            DO a=1, eNoNb
               xp = xp + xbl(:,a)*lFa%N(a,g)
               ub = ub + ubl(:,a)*lFa%N(a,g)
            END DO

!           Compute Nw and Nwxi of the mesh at the face integration
!           point using Newton method. Then calculate Nwx.
            xi = xi0
            CALL GETNNX(fs(1)%eType, fs(1)%eNoN, xwl, fs(1)%xib,
     2         fs(1)%Nb, xp, xi, Nw, Nwxi)
            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF)
     2         CALL GNN(fs(1)%eNoN, nsd, Nwxi, xwl, Nwx, Jac, Ks)

!           Compute Nq of the mesh at the face integration point using
!           Newton method.
            xi = xi0
            CALL GETNNX(fs(2)%eType, fs(2)%eNoN, xql, fs(2)%xib,
     2         fs(2)%Nb, xp, xi, Nq, Nqx)

            IF (nsd .EQ. 3) THEN
                CALL BWFLUID3D(fs(1)%eNoN, fs(2)%eNoN, w, Nw, Nq, Nwx,
     2             yl, ub, nV, tauB, lR, lK)
            ELSE
                CALL BWFLUID2D(fs(1)%eNoN, fs(2)%eNoN, w, Nw, Nq, Nwx,
     2             yl, ub, nV, tauB, lR, lK)
            END IF
         END DO

         DEALLOCATE(xwl, xql, Nw, Nq, Nwx, Nqx, Nwxi)

!     Now doing the assembly part
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
         END IF
#endif
      END DO

      DEALLOCATE(ptr, xl, yl, lR, lK, xbl, ubl, ubg)

      CALL DESTROY(fs(1))
      CALL DESTROY(fs(2))

      RETURN
      END SUBROUTINE SETBCRIS
!####################################################################
!     This subroutine updates the resistance and activation flag for the 
!     closed and open configurations of the RIS surfaces 
      SUBROUTINE RIS_UPDATER
      USE TYPEMOD
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND) :: iProj

      DO iProj = 1, RIS%nbrRIS
         !IF((RIS%nbrIter(iProj).LE.1) .AND. (RIS%status(iProj))) CYCLE

!        The valve is closed check if it should open
         IF (RIS%clsFlg(iProj)) THEN 
!          OPENING CONDITION: Check condition on the pressure difference        
            IF( RIS%meanP(iProj, 1) .GT. RIS%meanP(iProj, 2)  ) THEN 
               RIS%clsFlg(iProj) = .FALSE.
               std="RIS Proj "//iProj//": Going from close to open."
               RIS%nbrIter(iProj) = 0 
               ! I needed to update the state variables when the valve 
               ! goes from close to open to prevent the valve goes back
               ! to close at the next iteration. This is needed only for
               ! close to open and cannot be used for open to close.
               Ao = An
               Yo = Yn
               IF (dFlag) Do = Dn
               cplBC%xo = cplBC%xn
               RIS%restoreP(iProj) = .TRUE.
               RIS%pbc(iProj)%g = 0.5 * RIS%meanP(iProj, 1) + 
     2              0.5 * RIS%meanP(iProj, 2)
            END iF
         ELSE 
!        The valve is open, check if it should close. 
!           CLOSING CONDITION: Check existence of a backflow
            IF( RIS%meanFl(iProj) .LT. 0.) THEN 
               RIS%clsFlg(iProj) = .TRUE.
               std="RIS Proj "//iProj//": Going from open to close "
               RIS%nbrIter(iProj) = 0 
            END IF
         END IF
      END DO

      RETURN
      END SUBROUTINE RIS_UPDATER

!####################################################################
!     This subroutine will check the valve status if it is admissible 
!     or not, if not admissible we recompute the iteration until it will be 
      SUBROUTINE RIS_STATUS
      USE TYPEMOD
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND) :: iProj

      DO iProj = 1, RIS%nbrRIS
         RIS%nbrIter(iProj) = RIS%nbrIter(iProj) + 1
         RIS%status(iProj) = .TRUE.
!        If the valve is closed, chech the pressure difference, 
!        if the pressure difference is negative the valve should be open
!        -> the status is then not admissible
         IF (RIS%clsFlg(iProj)) THEN 
            IF( RIS%meanP(iProj,1) .GT. RIS%meanP(iProj,2)  ) THEN 
               std= "RIS Proj "//iProj//": **** Not admissible, 
     2              it should be open **** "
               RIS%status(iProj) = .FALSE.
            END iF
         ELSE 

!        If the valve is open, chech the flow, 
!        if the flow is negative the valve should be closed
!        -> the status is then not admissible
            IF( RIS%meanFl(iProj) .LT. 0.) THEN 
               std = "RIS Proj "//iProj//": **** Not admissible, 
     2              it should be closed **** "
               RIS%status(iProj) = .FALSE.
            END IF
         END IF
      END DO

      RETURN
      END SUBROUTINE RIS_STATUS



!####################################################################
!     This subroutine assembels the element stiffness matrix into the
!     global stiffness matrix (Val sparse matrix formatted as a vector)
      SUBROUTINE DOASSEM_RIS (d, eqN, lK, lR)
      USE TYPEMOD
      USE COMMOD, ONLY: dof, rowPtr, colPtr, R, Val, nMsh, nsd,
     2  grisMapList, cm, RIS
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: d, eqN(d)
      REAL(KIND=RKIND), INTENT(IN) :: lK(dof*dof,d,d), lR(dof,d)

      INTEGER(KIND=IKIND) a, b, ptr, rowN, colN, left, right, mapIdx(2), 
     2                    jM, mapIdxC(2), iProj
      INTEGER(KIND=IKIND) :: rowNadj=0
     
      !IF(.NOT.cm%mas()) RETURN
      DO iProj=1, RIS%nbrRIS
         IF(RIS%clsFlg(iProj)) CYCLE
         DO a=1, d

            rowN = eqN(a)
            IF (rowN .EQ. 0) CYCLE
            mapIdx = FINDLOC(grisMapList(iProj)%map, rowN)
            
            IF(mapIdx(1).EQ.0) CYCLE 

            DO jM=1, 2
               IF(jM .EQ. mapIdx(1)) CYCLE
               IF(grisMapList(iProj)%map(jM, mapIdx(2)) .EQ. 0) CYCLE
               
               rowNadj = grisMapList(iProj)%map(jM, mapIdx(2))

            END DO
C             R(1:nsd,rowNadj) = R(1:nsd,rowNadj) + lR(1:nsd,a)
            IF (rowNadj .EQ. 0) CYCLE
            R(:,rowNadj) = R(:,rowNadj) + lR(:,a)

            DO b=1, d
               colN = eqN(b)

!              If colN is also a ris node, we have to connect the cooresponding 
!              rowN node with the corresponding colN node. 
               mapIdxC = FINDLOC(grisMapList(iProj)%map, colN)

               IF(mapIdxC(1).NE.0) THEN 
                  DO jM=1, 2
                     IF(jM .EQ. mapIdxC(1)) CYCLE
                     IF(grisMapList(iProj)%map(jM,mapIdxC(2)).EQ.0) THEN
                        CYCLE
                     END IF
                     
                     colN = grisMapList(iProj)%map(jM, mapIdxC(2))
                  END DO
!               ELSE
!                   CYCLE
               END IF

               IF (colN .EQ. 0) CYCLE
               left  = rowPtr(rowNadj)
               right = rowPtr(rowNadj+1)
               ptr   = (right + left)/2

               DO WHILE (colN.NE.colPtr(ptr))
                  IF (colN .GT. colPtr(ptr)) THEN
                     left  = ptr
                  ELSE
                     right = ptr
                  END IF
                  ptr = (right + left)/2
               END DO
             
               Val(:,ptr) = Val(:,ptr) + lK(:,a,b)

C                Val(1:nsd*2+2,ptr) = Val(1:nsd*2+2,ptr) + lK(1:nsd*2+2,a,b)

C                Val(1:nsd,ptr) = Val(1:nsd,ptr) + lK(1:nsd,a,b)
C                Val(nsd+2:2*nsd+2,ptr) = Val(nsd+2:2*nsd+2,ptr) +
C      2                                            lK(nsd+2:2*nsd+2,a,b)
            END DO
         END DO
      END DO
      RETURN
      END SUBROUTINE DOASSEM_RIS
!####################################################################
!     This subroutine assembels the element stiffness matrix into the
!     global stiffness matrix (Val sparse matrix formatted as a vector)
      SUBROUTINE DOASSEM_VELRIS (d, eqN, lK, lR)
      USE TYPEMOD
      USE COMMOD, ONLY: dof, rowPtr, colPtr, R, Val, nMsh, nsd,
     2      grisMapList, RIS
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: d, eqN(d)
      REAL(KIND=RKIND), INTENT(IN) :: lK(dof*dof,d,d), lR(dof,d)

      INTEGER(KIND=IKIND) a, b, ptr, rowN, colN, left, right, mapIdx(2), 
     2                    jM, rowNadj, mapIdxC(2), iProj
      DO iProj=1, RIS%nbrRIS
         IF(RIS%clsFlg(iProj)) CYCLE
         DO a=1, d

            rowN = eqN(a)
            IF (rowN .EQ. 0) CYCLE
            mapIdx = FINDLOC(grisMapList(iProj)%map, rowN)
            
            IF(mapIdx(1).EQ.0) CYCLE 

            DO jM=1, 2
               
               IF(jM .EQ. mapIdx(1)) CYCLE
               IF(grisMapList(iProj)%map(jM, mapIdx(2)) .EQ. 0) CYCLE
               
               rowNadj = grisMapList(iProj)%map(jM, mapIdx(2))
            END DO

            R(1:nsd,rowNadj) = R(1:nsd,rowNadj) + lR(1:nsd,a)
C             R(:,rowNadj) = R(:,rowNadj) + lR(:,a)

            DO b=1, d
               colN = eqN(b)

!              If colN is also a ris node, we have to connect the cooresponding 
!              rowN node woth the corresponding colN node
               mapIdxC = FINDLOC(grisMapList(iProj)%map, colN)

               IF(mapIdxC(1).NE.0) THEN 
                  DO jM=1, 2
                     IF(jM .EQ. mapIdxC(1)) CYCLE
                     IF(grisMapList(iProj)%map(jM,mapIdxC(2)).EQ.0) THEN
                        CYCLE
                     END IF
                     
                     colN = grisMapList(iProj)%map(jM, mapIdxC(2))
                  END DO
               END IF

               IF (colN .EQ. 0) CYCLE
               left  = rowPtr(rowNadj)
               right = rowPtr(rowNadj+1)
               ptr   = (right + left)/2
               DO WHILE (colN .NE. colPtr(ptr))
                  IF (colN .GT. colPtr(ptr)) THEN
                     left  = ptr
                  ELSE
                     right = ptr
                  END IF
                  ptr = (right + left)/2
               END DO
               
C                Val(:,ptr) = Val(:,ptr) + lK(:,a,b)

               Val(1:nsd*2+2,ptr) = Val(1:nsd*2+2,ptr) + 
     2              lK(1:nsd*2+2,a,b)

C                Val(1:nsd,ptr) = Val(1:nsd,ptr) + lK(1:nsd,a,b)
C                Val(nsd+2:2*nsd+2,ptr) = Val(nsd+2:2*nsd+2,ptr) +
C         2                                         lK(nsd+2:2*nsd+2,a,b)
            END DO
         END DO
      END DO
      RETURN
      END SUBROUTINE DOASSEM_VELRIS
!####################################################################
!     This subroutine assembels the element stiffness matrix into the
!     global stiffness matrix (Val sparse matrix formatted as a vector)
      SUBROUTINE CLEAN_R_RIS
      USE TYPEMOD
      USE COMMOD
      IMPLICIT NONE

      INTEGER(KIND=IKIND) iM, j, Ac, nStk, iProj

      DO iProj=1, RIS%nbrRIS
         nStk = SIZE(grisMapList(iProj)%map,2)

         DO iM=1, 2
            DO j=1, nStk
               Ac = grisMapList(iProj)%map(iM,j)
               IF(Ac .EQ. 0) CYCLE
               R(1:nsd,Ac) = 0._RKIND
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE CLEAN_R_RIS
!####################################################################

      SUBROUTINE SETBCDIR_RIS(lA, lY, lD)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(tDof, tnNo), lY(tDof, tnNo),
     2   lD(tDof, tnNo)

      INTEGER(KIND=IKIND) iM, j, Ac, nStk, iProj

      DO iProj=1, RIS%nbrRIS
         nStk = SIZE(grisMapList(iProj)%map,2)

         DO iM=1, 2 
            DO j=1, nStk
               Ac = grisMapList(iProj)%map(iM,j)
               IF(Ac .EQ. 0) CYCLE
               lA(1:nsd,Ac) = 0._RKIND
               lY(1:nsd,Ac) = 0._RKIND
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE SETBCDIR_RIS

!####################################################################
!####################################################################
!     Weak treatment of RIS0D boundary conditions
      SUBROUTINE RIS0D_BC(Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) iFa, iBc, iM
      TYPE(bcType) :: lBc

      DO iBc=1, eq(cEq)%nBc
         iFa = eq(cEq)%bc(iBc)%iFa
         iM  = eq(cEq)%bc(iBc)%iM

         IF (.NOT.BTEST(eq(cEq)%bc(iBc)%bType,bType_Ris0D)) CYCLE
         
         IF ( eq(cEq)%bc(iBc)%clsFlgRis .EQ. 1 ) THEN

!           Weak Dirichlet BC for fluid/FSI equations
            lBc%weakDir = .TRUE.
            lBc%tauB    = eq(cEq)%bc(iBc)%res
            lBc%bType   = IBSET(lBc%bType,bType_Dir)
            lBc%bType   = IBSET(lBc%bType,bType_std)
            lBc%bType   = IBSET(lBc%bType,bType_flat)

            ALLOCATE(lBc%eDrn(nsd))
            lBc%eDrn = 0

!           Apply bc Dir 
            ALLOCATE(lBc%gx(msh(iM)%fa(iFa)%nNo))
            lBc%gx = 1._RKIND

            CALL SETBCDIRWL(lBc, msh(iM), msh(iM)%fa(iFa), Yg, Dg)

            DEALLOCATE(lBc%gx)
            DEALLOCATE(lBc%eDrn)
         ELSE 

!           Apply Neu bc 
            CALL SETBCNEUL(eq(cEq)%bc(iBc), msh(iM)%fa(iFa), Yg, Dg)

         END IF
      
      END DO
      
      RETURN
      END SUBROUTINE RIS0D_BC

!####################################################################
!####################################################################
!     This subroutine computes the mean pressure and flux on the 
!     ris surface 
      SUBROUTINE RIS0D_STATUS(Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) iFa, iBc, iM, m, s, e, nbrIter
      TYPE(bcType) :: lBc
      TYPE(faceType):: lFa
      REAL(KIND=RKIND) :: meanP = 0._RKIND
      REAL(KIND=RKIND) :: meanFl = 0._RKIND
      REAL(KIND=RKIND) :: tmp, tmp_new
      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:), sA(:)

      DO iBc=1, eq(cEq)%nBc
         iFa = eq(cEq)%bc(iBc)%iFa
         iM  = eq(cEq)%bc(iBc)%iM

         IF (.NOT.BTEST(eq(cEq)%bc(iBc)%bType,bType_Ris0D)) CYCLE
         
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         ALLOCATE (tmpV(maxnsd,tnNo))

!        Compute mean Q and pressure difference 
         m = 1
         s = eq(cEq)%s + nsd 
         e = s + m - 1

         tmpV(1:m,:) = Yn(s:e,:)

C          IF( msh(iM)%fa(iFa)%nEl .EQ. 0) THEN 
C             write(*,*), " outside RIS0D_STATUS proc ", cm%id() 
C             RETURN
C          END IF

         tmp = msh(iM)%fa(iFa)%area
         !meanP = Integ(msh(iM)%fa(iFa),tmpV,1)/tmp
         ALLOCATE(sA(tnNo))
         sA   = 1._RKIND
         lFa = msh(iM)%fa(iFa)
!        such update may be not correct
         tmp_new = Integ(lFa, sA)
         print*, "area og is ", tmp
         print*, "area up is ", tmp_new
         meanP = Integ(msh(iM)%fa(iFa),tmpV,1,m)/tmp_new

!        For the velocity 
         m = nsd
         s = eq(cEq)%s
         e = s + m - 1

         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         ALLOCATE (tmpV(maxnsd,tnNo))
         tmpV(1:m,:) = Yn(s:e,:)

         meanFl = Integ(msh(iM)%fa(iFa),tmpV,1,m)

         write(*,*)" The average pressure is: ", meanP
         write(*,*)" The pressure from 0D is: ", eq(cEq)%bc(iBc)%g 
         write(*,*)" The average flow is: ", meanFl

         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)

         RisnbrIter = RisnbrIter + 1
         IF( RisnbrIter .LE. 25 .AND. cTS .GT. 1) THEN 
            IF (.NOT.cm%seq()) CALL cm%bcast(RisnbrIter)  
            RETURN
         END IF

!--- 
!--- Update RES
!        Update the resistance - determine the configuration 
!     The valve is closed check if it should open
         IF (eq(cEq)%bc(iBc)%clsFlgRis .EQ. 1) THEN 
 !       OPENING CONDITION: Check condition on the pressure difference        
            IF( eq(cEq)%bc(iBc)%g .LT. meanP  ) THEN 
               eq(cEq)%bc(iBc)%clsFlgRis  = 0
               IF (.NOT.cm%seq()) THEN 
                  CALL cm%bcast(eq(cEq)%bc(iBc)%clsFlgRis) 
               END IF
               std = "!!! -- Going from close to open "
               RisnbrIter = 0 

!            CALL restore equal velocity at the interface - mean vel at the node 
            END iF
         ELSE 
!     The valve is open, check if it should close. 
!        CLOSING CONDITION: Check existence of a backflow
            IF( meanFl .LT. 0.) THEN 
               eq(cEq)%bc(iBc)%clsFlgRis = 1
               IF (.NOT.cm%seq()) THEN 
                  CALL cm%bcast(eq(cEq)%bc(iBc)%clsFlgRis) 
               END IF
               std = "!!! -- Going from open to close "
               RisnbrIter = 0 
            END IF
         END IF

!--- 

!--- Check for the status
      ! status = .TRUE.
!     If the valve is closed, chech the pressure difference, 
!     if the pressure difference is negative the valve should be open
!     -> the status is then not admissible
         IF (eq(cEq)%bc(iBc)%clsFlgRis .EQ. 1) THEN 
            IF( eq(cEq)%bc(iBc)%g .LT. meanP  ) THEN 
               std = "** Not admissible, should be open **"
               ! status = .FALSE.
            END IF
         ELSE 
!     If the valve is open, chech the flow, 
!     if the flow is negative the valve should be closed
!     -> the status is then not admissible
            IF( meanFl .LT. 0.) THEN 
               std = "** Not admissible, it should be closed **"
               ! status = .FALSE.
            END IF
         END IF
!---- 
!---- 
      END DO

      IF (.NOT.cm%seq()) CALL cm%bcast(RisnbrIter)  

      RETURN
      END SUBROUTINE RIS0D_STATUS
