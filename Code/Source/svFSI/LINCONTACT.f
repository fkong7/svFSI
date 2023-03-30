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

      INTEGER(KIND=IKIND) iM, iFa, cntAct

      cntAct = 0

!     Formulation for linear contact of a thick structure and a fixed 
!     wall, yL = defined form file, for the moment we consider 0 

!     The integrals are performed only on the structure boundary. 
!     We consider a relaxed penalization approach, whcih means that a gap 
!     is included in the formulation.


!     Loop on faces and call the BASSEMLCONT

      DO iM=1, nMsh
         DO iFa=1, msh(iM)%nFa
C             write(*,*)" id proc ", cm%id()," iM,iFa ", iM,iFa
            CALL BASSEMLCONT(msh(iM)%fa(iFa), Yg, Dg, cntAct)             
         END DO
      END DO

      cntAct = cm%reduce(cntAct)

      IF(cntAct.GT.0) contAct = .TRUE.
      IF(contAct.AND.cm%mas()) write(*,*)"*** Contact activated ***" 

      RETURN
      END SUBROUTINE LCONTACTFORCES
!#######################################################################
!     Construct Neumann BCs
      SUBROUTINE BASSEMLCONT(lFa, Yg, Dg, cntAct)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(OUT) :: cntAct

      INTEGER(KIND=IKIND) a, e, g, Ac, Ec, iM, cPhys, eNoN, b
      REAL(KIND=RKIND) w, h, nV(nsd), y(tDof), Jac, hc(nsd), gap, 
     2   gamma, yw, afu, ys

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: N(:), yl(:,:), lR(:,:),
     2   lK(:,:,:), xl(:,:), dl(:,:)

      gap = cntGap
      yw = 0._RKIND
      hc = 0._RKIND
      gamma = 1.e5 !1.e4

C       write(*,*)" proc: ",cm%id(),", lFa%nEl", lFa%nEl

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
                     cntAct = 1

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
!     Routine to compute vf, scalar value for each local mesh, defining the 
!     spatial distribution of K for contact in the context of NSB eqs
      SUBROUTINE CMPVfCONTACT()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) iM, cPhys, nEl, e, a, cnt, Ac, find, Acl, eNgb
      INTEGER(KIND=IKIND), ALLOCATABLE :: idTet(:)
      REAL(KIND=RKIND) sEp

      sEp=1.E-6_RKIND

      IF(.NOT.(flagLCONT.OR.flagNLCONT)) RETURN

      DO iM=1, nMsh
!        Assumption that the physics is constant for all the elements of the mesh 
!        Look at element 1 for eq 1 (it should be FSI/fluid/struct )
         cDmn  = DOMAIN(msh(iM), 1, 1)
         cPhys = eq(1)%dmn(cDmn)%phys
         nEl = msh(iM)%nEl

         IF(nEl.EQ.0) CYCLE
                
         IF(.NOT.ALLOCATED(msh(iM)%vf)) ALLOCATE(msh(iM)%vf(nEl))    

         msh(iM)%vf = 0._RKIND

         IF (cPhys.EQ.phys_struct .OR. 
     2           cPhys.EQ.phys_ustruct) msh(iM)%vf = 2._RKIND
        
         IF (cPhys.NE.phys_fluid) CYCLE 

         ALLOCATE(idTet(msh(iM)%eNoN))

!--      Layer 1 - vf = 1.
         DO e=1, msh(iM)%nEl 

            cnt = 0

            DO a=1, msh(iM)%eNoN
               Acl = msh(iM)%IEN(a,e)
               Ac = ltg(Acl)

               find = FINDLOC(listPrjID, Ac,dim=1)
               IF(find.GT.0) THEN 
                  cnt = cnt + 1
               END IF

            END DO

            IF(cnt.GE.1) msh(iM)%vf(e) = 1._RKIND

         END DO

!--      Layer 2 - vf = 0.75
         DO e=1, msh(iM)%nEl 

            IF( msh(iM)%vf(e).LT.1._RKIND-sEp) CYCLE 

            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               idTet(a) = Ac
            END DO

            DO eNgb =1, msh(iM)%nEl 

               IF(msh(iM)%vf(eNgb).GE.1._RKIND-sEp) CYCLE

               cnt = 0
               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%IEN(a,eNgb)

                  find = FINDLOC(idTet, Ac,dim=1)
                  IF(find.GT.0) THEN 
                     cnt = cnt + 1
                  END IF
               END DO

               IF(cnt.GT.0) msh(iM)%vf(eNgb)=0.75_RKIND 

            END DO

         END DO

!--      Layer 3 - vf = 0.5, looking from neig vf=0, with neigh vf=0.75
         DO e=1, msh(iM)%nEl 

!           Select only elem with vf=0.
            IF( msh(iM)%vf(e).GT.0.75_RKIND-sEp ) CYCLE 

            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               idTet(a) = Ac
            END DO

            DO eNgb =1, msh(iM)%nEl 

!              Search neigh with vf=0.5
               IF(msh(iM)%vf(eNgb).GT.0.75_RKIND+sEp) CYCLE
               IF(msh(iM)%vf(eNgb).LT.0.75_RKIND-sEp) CYCLE

               cnt = 0
               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%IEN(a,eNgb)

                  find = FINDLOC(idTet, Ac,dim=1)
                  IF(find.GT.0) THEN 
                     cnt = cnt + 1
                  END IF
               END DO

               IF(cnt.GT.0) msh(iM)%vf(e)=0.5_RKIND 

            END DO

         END DO

!--      Layer 4 - vf = 0.25, looking from neig vf=0, with neigh vf=0.5
         DO e=1, msh(iM)%nEl 

!           Select only elem with vf=0.
            IF( msh(iM)%vf(e).GT.0.5_RKIND-sEp ) CYCLE 

            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               idTet(a) = Ac
            END DO

            DO eNgb =1, msh(iM)%nEl 

!              Search neigh with vf=0.5
               IF(msh(iM)%vf(eNgb).GT.0.5_RKIND+sEp) CYCLE
               IF(msh(iM)%vf(eNgb).LT.0.5_RKIND-sEp) CYCLE

               cnt = 0
               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%IEN(a,eNgb)

                  find = FINDLOC(idTet, Ac,dim=1)
                  IF(find.GT.0) THEN 
                     cnt = cnt + 1
                  END IF
               END DO

               IF(cnt.GT.0) msh(iM)%vf(e)=0.25_RKIND 

            END DO

         END DO

         DEALLOCATE(idTet)

      END DO

      RETURN
      END SUBROUTINE CMPVfCONTACT

!####################################################################
