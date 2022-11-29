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
!-----------------------------------------------------------------------
!
!     This routine applies penalty-based contact model for thick 
!     structures with a fixed wall 
!
!-----------------------------------------------------------------------

      SUBROUTINE LCONTACTFORCES(Yg,Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) iM, iFa

!     Formulation for linear contact of a thick structure and a fixed 
!     wall, yL = defined form file, for the moment we consider 0 

!     The integrals are performed only on the structure boundary. 
!     We consider a relaxed penalization approach, whcih means that a gap 
!     is included in the formulation.


!     Loop on faces and call the BASSEMLCONT

      DO iM=1, nMsh
         DO iFa=1, msh(iM)%nFa
            CALL BASSEMLCONT(msh(iM)%fa(iFa), Yg, Dg)             
         END DO
      END DO

      RETURN
      END SUBROUTINE LCONTACTFORCES
!#######################################################################
!     Construct Neumann BCs
      SUBROUTINE BASSEMLCONT(lFa, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, Ec, iM, cPhys, eNoN, b
      REAL(KIND=RKIND) w, h, nV(nsd), y(tDof), Jac, hc(nsd), gap, 
     2   gamma, yw, afu, ys

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: N(:), yl(:,:), lR(:,:),
     2   lK(:,:,:), xl(:,:), dl(:,:)

      gap = 0.02_RKIND
      yw = 0._RKIND
C       yw = 0._RKIND
      hc = 0._RKIND
      gamma = 1.e5 !1.e4


      iM   = lFa%iM
      eNoN = lFa%eNoN
      DO e=1, lFa%nEl
         Ec    = lFa%gE(e)
         cDmn  = DOMAIN(msh(iM), cEq, Ec)
         cPhys = eq(cEq)%dmn(cDmn)%phys
!        This should be the coeff in from of the directional derivative 
!        containing the displacement 
         afu   = eq(cEq)%af*eq(cEq)%beta*dt*dt

C          IF (cPhys .NE. phys_struct) RETURN 
         IF ((cPhys .EQ. phys_struct) .OR. 
     2             (cPhys .EQ. phys_ustruct)) THEN 
            CONTINUE 
         ELSE 
            RETURN 
         END IF

C          IF (cPhys .EQ. phys_fluid) RETURN 

         ALLOCATE(ptr(eNoN), N(eNoN), yl(tDof,eNoN), xl(nsd,eNoN), 
     2      lR(dof,eNoN), lK(dof*dof,eNoN,eNoN), dl(tDof,eNoN))
         lK = 0._RKIND
         lR = 0._RKIND
         DO a=1, eNoN
            Ac      = lFa%IEN(a,e)
            ptr(a)  = Ac
            xl(:,a) = x(:,Ac)
            yl(:,a) = Yg(:,Ac)
            dl(:,a) = Dg(:,Ac)
         END DO

C          IF( e .EQ. 1) THEN 
C             write(*,*)" dof = ", dof
C             write(*,*)" tDof = ", tDof
C             write(*,*)" xl = ", xl(2,1)
C             write(*,*)" dl = ", dl(2,1)
C             write(*,*)" xl current = ", xl(2,1) + dl(2,1)
C             write(*,*)" yl = ", yl(2,1)
C             write(*,*)""//""
C          END IF

!        Updating the shape functions, if neccessary
         IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(iM), lFa, e)

         DO g=1, lFa%nG
            CALL GNNB(lFa, e, g, nsd-1, eNoN, lFa%Nx(:,:,g), nV)
            Jac = SQRT(NORM(nV))
            nV  = nV/Jac
            w   = lFa%w(g)*Jac
            N   = lFa%N(:,g)

            y = 0._RKIND
            ys = 0._RKIND
            DO a=1, eNoN
               y  = y + N(a)*yl(:,a)
               ys = ys + N(a)*(xl(2,a) + dl(2,a))
            END DO

            !CALL BFLUID(eNoN, w, N, y, h, nV, lR, lK)


C             wl  = w*eq(cEq)%af*eq(cEq)%gam*dt

C             udn = 0._RKIND
C             IF (mvMsh) THEN
C                DO i=1, nsd
C                   j    = i + nsd + 1
C                   u(i) = y(i) - y(j)
C                   udn  = udn + u(i)*nV(i)
C                END DO
C             ELSE
C                DO i=1, nsd
C                   u(i) = y(i)
C                   udn  = udn + u(i)*nV(i)
C                END DO
C             END IF

C             udn = 0.5_RKIND*eq(cEq)%dmn(cDmn)%prop(backflow_stab)*
C            2   eq(cEq)%dmn(cDmn)%prop(fluid_density)*(udn - ABS(udn))

!           Computation of the gap function


      !     Here the loop is started for constructing left and right hand side
C             IF (nsd .EQ. 2) THEN
C                DO a=1, eNoN
C C                   lR(1,a) = lR(1,a) - w*N(a)*hc(1)

C                   hc  = xl(2,a) + 
C                   lR(2,a) = lR(2,a) - w*N(a)*hc(2)
C                   DO b=1, eNoN
C                      T1        = wl*N(a)*N(b)*udn
C C                      lK(1,a,b) = lK(1,a,b) - T1
C C                      lK(5,a,b) = lK(5,a,b) - T1
C                   END DO
C                END DO
C             ELSE
               DO a=1, eNoN

C                   hc(2) = xl(2,a) + dl(2,a) - ( yw + gap) 
                  hc(2) = ys - ( yw + gap) 

C                   write(*,*)"xl(2,a) =  ", xl(2,a)
C                   write(*,*)"dl(2,a) =  ", dl(2,a)
C                   write(*,*)" ys is ", ys
C                   write(*,*)" gap is ", hc(2)

                  IF( hc(2) .LT. 0._RKIND ) THEN 

C                      write(*,*)" **** CONTACT ACTIVE **** "
C                      write(*,*)"xl(2,a) =  ", xl(2,a)
C                      write(*,*)"dl(2,a) =  ", dl(2,a)
C                      write(*,*)" ys is ", ys
C                      write(*,*)" gap is ", hc(2)

                     lR(1,a) = lR(1,a) + gamma*w*N(a)*hc(1)
                     lR(2,a) = lR(2,a) + gamma*w*N(a)*hc(2)
                     lR(3,a) = lR(3,a) + gamma*w*N(a)*hc(3)
                     DO b=1, eNoN
C                         lK(1,a,b)  = lK(1,a,b)  - gamma*w*afu*N(a)*N(b)
                        lK(5,a,b)  = lK(5,a,b) + gamma*w*afu*N(a)*N(b)
C                         lK(9,a,b)  = lK(9,a,b)  - gamma*w*afu*N(a)*N(b)
                     END DO
                  END IF
               END DO
C             END IF



         END DO

!        Now doing the assembly part
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
         END IF
#endif
         DEALLOCATE(ptr, N, yl, lR, lK, xl, dl)
      END DO

      RETURN
      END SUBROUTINE BASSEMLCONT
!#######################################################################




C       LOGICAL :: flag
C       INTEGER(KIND=IKIND) :: i, j, k, l, m, iM, jM, e, a, Ac, b, Bc, g,
C      2   eNoN, insd, nnb, maxNnb
C       REAL(KIND=RKIND) :: kl, hl, w, Jac, al, c, d, pk, nV1(nsd),
C      2   nV2(nsd), x1(nsd), x2(nsd), x12(nsd), xmin(nsd), xmax(nsd)

C       INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), bBox(:,:)
C       REAL(KIND=RKIND), ALLOCATABLE :: sA(:), sF(:,:), N(:), Nx(:,:),
C      2   gCov(:,:), gCnv(:,:), xl(:,:), lR(:,:)

C       IF (eq(cEq)%phys .NE. phys_shell) RETURN

C       i  = eq(cEq)%s
C       j  = i + 1
C       k  = j + 1
C       kl = cntctM%k
C       hl = cntctM%h

C !     Compute normal vectors at each node in the current configuration
C       ALLOCATE(sF(nsd,tnNo), sA(tnNo))
C       sF = 0._RKIND
C       sA = 0._RKIND
C       DO iM=1, nMsh
C          IF (.NOT.msh(iM)%lShl) CYCLE
C          eNoN = msh(iM)%eNoN
C          insd = nsd - 1
C          ALLOCATE(Nx(insd,eNoN), N(eNoN), xl(nsd,eNoN), gCov(nsd,insd),
C      2      gCnv(nsd,insd))
C          Nx = msh(iM)%Nx(:,:,1)
C          DO e=1, msh(iM)%nEl
C             DO a=1, eNoN
C                Ac = msh(iM)%IEN(a,e)
C                xl(1,a) = x(1,Ac) + Dg(i,Ac)
C                xl(2,a) = x(2,Ac) + Dg(j,Ac)
C                xl(3,a) = x(3,Ac) + Dg(k,Ac)
C             END DO
C             CALL GNNS(eNoN, Nx, xl, nV1, gCov, gCnv)
C             Jac = SQRT(NORM(nV1(:)))
C             nV1(:) = nV1(:)/Jac

C             DO g=1, msh(iM)%nG
C                w = msh(iM)%w(g) * Jac
C                N(:) = msh(iM)%N(:,g)
C                DO a=1, eNoN
C                   Ac = msh(iM)%IEN(a,e)
C                   sA(Ac) = sA(Ac) + w*N(a)
C                   sF(:,Ac) = sF(:,Ac) + w*N(a)*nV1(:)
C                END DO
C             END DO
C          END DO
C          DEALLOCATE(N, Nx, xl, gCov, gCnv)
C       END DO

C       CALL COMMU(sF)
C       CALL COMMU(sA)

C       DO Ac=1, tnNo
C          IF (.NOT.ISZERO(sA(Ac))) sF(:,Ac) = sF(:,Ac)/sA(Ac)
C          Jac = SQRT(SUM(sF(:,Ac)**2))
C          IF (.NOT.ISZERO(Jac)) sF(:,Ac) = sF(:,Ac) / Jac
C       END DO

C !     Create a bounding box around possible region of contact and bin
C !     the box with neighboring nodes
C       maxNnb = 15
C  101  maxNnb = maxNnb + 5
C       IF (ALLOCATED(bBox)) DEALLOCATE(bBox)
C       ALLOCATE(bBox(maxNnb,tnNo))
C       bBox = 0
C       DO iM=1, nMsh
C          IF (.NOT.msh(iM)%lShl) CYCLE
C          DO a=1, msh(iM)%nNo
C             Ac = msh(iM)%gN(a)
C             x1(1) = x(1,Ac) + Dg(i,Ac)
C             x1(2) = x(2,Ac) + Dg(j,Ac)
C             x1(3) = x(3,Ac) + Dg(k,Ac)

C !           Box limits for each node
C             xmin(:) = x1(:) - cntctM%c
C             xmax(:) = x1(:) + cntctM%c

C !           Load the box with neighboring nodes lying within it
C             DO jM=1, nMsh
C                IF (iM.EQ.jM .OR. .NOT.msh(jM)%lShl) CYCLE
C                DO b=1, msh(jM)%nNo
C                   Bc = msh(jM)%gN(b)
C                   x2(1) = x(1,Bc) + Dg(i,Bc)
C                   x2(2) = x(2,Bc) + Dg(j,Bc)
C                   x2(3) = x(3,Bc) + Dg(k,Bc)
C                   IF (x2(1).GE.xmin(1) .AND. x2(1).LE.xmax(1) .AND.
C      2                x2(2).GE.xmin(2) .AND. x2(2).LE.xmax(2) .AND.
C      3                x2(3).GE.xmin(3) .AND. x2(3).LE.xmax(3)) THEN
C                      DO l=1, maxNnb
C                         IF (bBox(l,Ac) .EQ. 0) THEN
C                            bBox(l,Ac) = Bc
C                            EXIT
C                         END IF
C                         IF (Bc .GT. bBox(l,Ac)) CYCLE
C                         IF (Bc .EQ. bBox(l,Ac)) EXIT
C                         IF (bBox(maxNnb,Ac) .NE. 0) GOTO 101
C                         DO m=maxNnb, l+1, -1
C                            bBox(m,Ac) = bBox(m-1,Ac)
C                         END DO
C                         bBox(l,Ac) = Bc
C                         EXIT
C                      END DO
C                      IF (l .GT. maxNnb) GOTO 101
C                   END IF
C                END DO ! b
C             END DO ! jM
C          END DO ! a
C       END DO ! iM

C !     Check if any node is strictly involved in contact and compute
C !     corresponding penalty forces assembled to the residue
C       ALLOCATE(lR(dof,tnNo), incNd(tnNo))
C       lR    = 0._RKIND
C       incNd = 0
C       DO Ac=1, tnNo
C          IF (bBox(1,Ac) .EQ. 0) CYCLE
C          x1(1)  = x(1,Ac) + Dg(i,Ac)
C          x1(2)  = x(2,Ac) + Dg(j,Ac)
C          x1(3)  = x(3,Ac) + Dg(k,Ac)
C          nV1(:) = sF(:,Ac)
C          nNb    = 0
C          DO a=1, maxNnb
C             Bc = bBox(a,Ac)
C             IF (Bc .EQ. 0) CYCLE
C             x2(1)  = x(1,Bc) + Dg(i,Bc)
C             x2(2)  = x(2,Bc) + Dg(j,Bc)
C             x2(3)  = x(3,Bc) + Dg(k,Bc)
C             nV2(:) = sF(:,Bc)

C             x12 = x1(:) - x2(:)
C             c   = SQRT(NORM(x12))
C             al  = SQRT(ABS(NORM(nV1, nV2)))

C             IF (c.LE.cntctM%c .AND. al.GE.cntctM%al) THEN
C                d = NORM(x12, nV2)
C                flag = .FALSE.
C                IF (d.GE.-cntctM%h .AND. d.LT.0._RKIND) THEN
C                   pk = 0.5_RKIND*kl/hl * (d + hl)**2._RKIND
C                   flag = .TRUE.
C                ELSE IF (d .GE. 0._RKIND) THEN
C                   pk = 0.5_RKIND*kl * (hl + d)
C                   flag = .TRUE.
C                ELSE
C                   pk = 0._RKIND
C                END IF
C                IF (flag) THEN
C                   incNd(Ac) = 1
C                   nNb = nNb + 1
C                   lR(1:nsd,Ac) = lR(1:nsd,Ac) - pk*nV1(:)
C                END IF
C             END IF
C          END DO
C          IF (nNb .NE. 0) lR(:,Ac) = lR(:,Ac) / REAL(nNb, KIND=RKIND)
C       END DO
C       DEALLOCATE(sA, sF, bBox)

C !     Return if no penalty forces are to be added
C       IF (SUM(incNd) .EQ. 0) RETURN

C !     Global assembly
C       DO Ac=1, tnNo
C          IF (incNd(Ac) .EQ. 0) CYCLE
C          R(:,Ac) = R(:,Ac) + lR(:,Ac)
C       END DO
C       DEALLOCATE(lR, incNd)

C       RETURN
C       END SUBROUTINE CONTACTFORCES
C !####################################################################
