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
!     Main routine that contains the calls to all major routines and
!     general structure of the code.
!
!--------------------------------------------------------------------

      PROGRAM MAIN
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      LOGICAL l1, l2, l3
      INTEGER(KIND=IKIND) i, iM, iBc, ierr, iEqOld, stopTS, j, 
     2      iProj
      REAL(KIND=RKIND) timeP(3)

      INTEGER(KIND=IKIND), ALLOCATABLE :: incL(:)
      REAL(KIND=RKIND), ALLOCATABLE :: Ag(:,:), Yg(:,:), Dg(:,:), res(:)

      REAL(KIND=RKIND) v1(3), v2(3), v3(3), v4(3), p(3), V(3,2)
      REAL(KIND=RKIND)  v41(3), vp1(3), N(3), dotV4, dotP, dt_ori,
     2      dt_cum
      INTEGER(KIND=IKIND) same 

      IF (IKIND.NE.LSIP .OR. RKIND.NE.LSRP) THEN
         STOP "Incompatible datatype precision between solver and FSILS"
      END IF

      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.


      savedOnce = .FALSE.
      CALL MPI_INIT(i)
      CALL cm%new(MPI_COMM_WORLD)

!     Initiating the exception tracing
      CALL EXCEPTIONS

      resetSim  = .FALSE.
      rmsh%cntr = 0

!     Reading the user-defined parameters from foo.inp
 101  CALL READFILES
!     Doing the partitioning and distributing the data to the all
!     Processors
      CALL DISTRIBUTE

!     Initializing the solution vectors and constructing LHS matrix
!     format
      CALL INITIALIZE(timeP)
      stopTS = nTS

      dbg = 'Allocating intermediate variables'
      ALLOCATE(Ag(tDof,tnNo), Yg(tDof,tnNo), Dg(tDof,tnNo),
     2   res(nFacesLS), incL(nFacesLS))

!--------------------------------------------------------------------
!     Outer loop for marching in time. When entring this loop, all old
!     variables are completely set and satisfy BCs.
      IF (cTS .LE. nITS) dt = dt/10._RKIND
      dt_ori = dt
      dt_cum = 0._RKIND
      DO
!     Adjusting the time step size once initialization stage is over
         IF (cTS .EQ. nITS) THEN
            dt = dt*10._RKIND
            std = " New time step size: "//dt
            dt_ori = dt
         END IF
!     Incrementing time step, hence cTS will be associated with new
!     variables, i.e. An, Yn, and Dn
         cTS    = cTS + 1
         time   = time + dt
! --- RIS GOTO 1 should be here, we do not update the time but we redo all the rest          
11       cEq    = 1
         eq%itr = 0
         eq%ok  = .FALSE.

!     Compute mesh properties to check if remeshing is required
         IF (mvMsh .AND. rmsh%isReqd) THEN
            CALL CALCMESHPROPS(nMsh, msh)
            IF (resetSim) EXIT
         END IF

!     Predictor step
         CALL PICP

!     Apply Dirichlet BCs strongly
         CALL SETBCDIR(An, Yn, Dn)
!     Inner loop for iteration
         DO
            iEqOld = cEq
            IF (cplBC%coupled .AND. cEq.EQ.1) THEN
               CALL SETBCCPL
               CALL SETBCDIR(An, Yn, Dn)
            END IF

!        Initiator step (quantities at n+am, n+af)
            CALL PICI(Ag, Yg, Dg)
            IF (ALLOCATED(Rd)) THEN
               Rd = 0._RKIND
               Kd = 0._RKIND
            END IF

            dbg = 'Allocating the RHS and LHS'
            CALL LSALLOC(eq(cEq))

!        Compute body forces. If phys is shells or CMM (init), apply
!        contribution from body forces (pressure) to residue
            CALL SETBF(Dg)

            dbg = "Assembling equation <"//eq(cEq)%sym//">"
            DO iM=1, nMsh
               CALL GLOBALEQASSEM(msh(iM), Ag, Yg, Dg)
               dbg = "Mesh "//iM//" is assembled"
            END DO

!        Treatment of boundary conditions on faces
!        Apply Neumman or Traction boundary conditions
            CALL SETBCNEU(Yg, Dg)

!        Apply CMM BC conditions
            IF (.NOT.cmmInit) CALL SETBCCMM(Ag, Dg)

!        Apply weakly applied Dirichlet BCs
            CALL SETBCDIRW(Yg, Dg)
            IF (risFlag)  CALL RIS_RESBC(Yg, Dg)
            IF (ris0DFlag)  CALL RIS0D_BC(Yg, Dg)

!        Apply contact model and add its contribution to residue
            IF (iCntct) CALL CONSTRUCT_CONTACTPNLTY(Dg)

!        Synchronize R across processes. Note: that it is important
!        to synchronize residue, R before treating immersed bodies as
!        ib%R is already communicated across processes
            IF (.NOT.eq(cEq)%assmTLS) CALL COMMU(R)

!        Update residue in displacement equation for USTRUCT phys.
!        Note that this step is done only first iteration. Residue
!        will be 0 for subsequent iterations
            IF (sstEq) CALL USTRUCTR(Yg)

            IF ((eq(cEq)%phys .EQ. phys_stokes) .OR.
     2          (eq(cEq)%phys .EQ. phys_fluid)  .OR.
     3          (eq(cEq)%phys .EQ. phys_ustruct).OR.
     4          (eq(cEq)%phys .EQ. phys_fsi)) THEN
               CALL THOOD_ValRC()
            END IF

!        Update LHS for clamped BC
            CALL SETBC_CLMPD()

!        IB treatment: for explicit coupling, simply construct residue.
            IF (ibFlag) THEN
               IF (ib%cpld .EQ. ibCpld_I) THEN
                  CALL IB_IMPLICIT(Ag, Yg, Dg)
               END IF
               CALL IB_CONSTRUCT()
            END IF

            incL = 0
            IF (eq(cEq)%phys .EQ. phys_mesh) incL(nFacesLS) = 1
            IF (cmmInit) incL(nFacesLS) = 1
            DO iBc=1, eq(cEq)%nBc
               i = eq(cEq)%bc(iBc)%lsPtr
               IF (i .NE. 0) THEN
!              Resistance term for coupled Neu face tangent contribution
                  res(i)  = eq(cEq)%gam*dt*eq(cEq)%bc(iBc)%r
                  incL(i) = 1
               END IF
            END DO

            dbg = "Solving equation <"//eq(cEq)%sym//">"
            !2, incL, res
            CALL LSSOLVE(eq(cEq), incL, res)

!        Solution is obtained. Corrector and check for convergence
            CALL PICC

!        Checking for exceptions
            CALL EXCEPTIONS

!        Writing out the time passed, residual, and etc.
            IF (ALL(eq%ok)) EXIT
            CALL OUTRESULT(timeP, 2, iEqOld)

         END DO
!     End of inner loop

!     IB treatment: interpolate flow data on IB mesh from background
!     fluid mesh for explicit coupling, update old solution for implicit
!     coupling
         IF (ibFlag) THEN
            CALL IB_INTERPYU(Yn, Dn)
            IF (ib%cpld .EQ. ibCpld_I) THEN
               ib%Auo = ib%Aun
               ib%Ubo = ib%Ubn
            END IF
         END IF
! ---- Here we probably have to update the ris resistance value
! ---- If the state has to change, we recompute this time step GOTO 1
! ---- Control where if the time and the new has changed! 
         IF ( risFlag ) THEN 
            RIS%restoreP = .FALSE.
            CALL RIS_MEANQ
            CALL RIS_STATUS
            std = " Iteration: "//cTS
            dt_cum = dt_cum + dt
            DO iProj=1, RIS%nbrRIS
                std = "Status for RIS projection: "//iProj
                std = "  RIS iteration: "//RIS%nbrIter(iProj)
                std = "  Is the valve close? "//RIS%clsFlg(iProj)
                std = "  The status is "//RIS%status(iProj)
                IF (.NOT.RIS%status(iProj)) THEN
                    ! If pressure difference is too big, repeat this
                    ! time step with a much smaller dt
                    IF ((RIS%meanP(iProj,1)-RIS%meanP(iProj,2))
     2                  .GT.ABS(RIS%meanP(iProj,1)*0.01).AND.
     3                      dt_ori/dt.LE.8._RKIND) THEN
                        dt = dt/2._RKIND
                        write(*,*) "Reducing time step to", dt
                        dt_cum = 0._RKIND
                        An = Ao
                        Yn = Yo
                        IF (dFlag) Dn = Do
                        cplBC%xn = cplBC%xo
                        GOTO 11
                    END IF
                END IF
            END DO
            IF (cTS .GE. 0) THEN
                IF((.NOT.ALL(RIS%status)))THEN 
                    IF (ANY(RIS%nbrIter.LE.2).AND.(cTS.GE.2)) THEN
                        std = "Valve status just changed. Do not update"
                    ELSE
                        CALL RIS_UPDATER
                        dt_cum = 0._RKIND
                        GOTO 11
                    END IF
                END IF
            END IF
            IF (dt_cum .LT. dt_ori) THEN
                write(*,*) "Current time is: ", dt_cum, "total time is",
     2              dt_ori, ". Solve for another time step with", dt
                Ao = An
                Yo = Yn
                IF (dFlag) Do = Dn
                cplBC%xo = cplBC%xn
                CALL TXT(.FALSE.)
                time = time + dt
                write(*,*) "Set time to be", time
                GOTO 11
            ELSE
                write(*,*) "Finished iterations for RIS"
                time = time + dt - dt_cum
                dt = dt_ori
                dt_cum = 0._RKIND
                write(*,*) "Restoring time step to", dt
                write(*,*) "Restoring time to", time+dt
            END IF
         END IF

!     Saving the TXT files containing average and fluxes
         CALL TXT(.FALSE.)

         IF (rmsh%isReqd) THEN
            l1 = MOD(cTS,rmsh%cpVar) .EQ. 0
            IF (l1) THEN
               rmsh%rTS = cTS-1
               rmsh%time = time-dt
               rmsh%iNorm(:) = eq(:)%iNorm
               rmsh%A0(:,:) = Ao(:,:)
               rmsh%Y0(:,:) = Yo(:,:)
               rmsh%D0(:,:) = Do(:,:)
            END IF
         END IF

         IF (cm%mas()) THEN
            INQUIRE(FILE=stopTrigName, EXIST=l1)
            IF (l1) THEN
               OPEN(664,FILE=TRIM(stopTrigName))
               READ(664,'(I10)',ADVANCE='NO',IOSTAT=ierr) stopTS
               CLOSE(664)
               IF (ierr .EQ. -1) stopTS = cTS
            ELSE
               stopTS = nTS
            END IF
         END IF
         CALL cm%bcast(stopTS)
         l1 = cTS .GE. stopTS
         l2 = MOD(cTS,stFileIncr) .EQ. 0
!     Saving the result to restart bin file
         IF (l1 .OR. l2) CALL WRITERESTART(timeP)

!     Writing results into the disk with VTU format
         IF (saveVTK) THEN
            l2 = MOD(cTS,saveIncr) .EQ. 0
            l3 = cTS .GE. saveATS
            IF (l2 .AND. l3) THEN
               CALL OUTRESULT(timeP, 3, iEqOld)
               CALL WRITEVTUS(An, Yn, Dn, .FALSE.)
               IF (ibFlag) CALL IB_WRITEVTUS(ib%Yb, ib%Ubo)
            ELSE
               CALL OUTRESULT(timeP, 2, iEqOld)
            END IF
         ELSE
            CALL OUTRESULT(timeP, 2, iEqOld)
         END IF
         IF (pstEq) CALL OUTDNORM()

         IF (ibFlag) CALL IB_OUTCPUT()


         IF ((cEq.EQ.1) .AND. ris0DFlag ) THEN 
            CALL RIS0D_STATUS

!          Redo the fluid iteration without updating the time 
           IF( RisnbrIter .LE. 2) GOTO 11

         END IF

!        Part related to unfitted RIS
!        If the valve is active, look at the pressure difference 
         IF(urisFlag) THEN 

            IF( urisActFlag ) THEN 
               cntURIS = cntURIS + 1

               IF(.NOT.((cTS.GE.10).AND.(cntURIS .LT. 10))) THEN 
                  CALL URIS_MEANP  
               END IF

               !CALL URIS_MEANV 
               IF( cntURIS .EQ. 0) THEN 
                  IF( cntURIS .LT. 10) GOTO 11
               END IF

            ELSE 
               cntURIS = cntURIS + 1
               
               IF(.NOT.((cTS.GE.10).AND.(cntURIS .LT. 10))) THEN 
                  CALL URIS_MEANV  
               END IF

               !CALL URIS_MEANV 
               IF( cntURIS .EQ. 0) THEN 
                  IF( cntURIS .LT. 10) GOTO 11
               END IF

            END IF

C             IF( time .EQ. 0.5) urisActFlag = .FALSE.

            IF (mvMsh) THEN 
               CALL URIS_UpdateDisp !(Do,Dn)
            END IF
            CALL URIS_WRITEVTUS(uris%Yd)
         END IF
!---     end RIS/URIS stuff 

!     Solution is stored here before replacing it at next time step
         Ao = An
         Yo = Yn
         IF (dFlag) Do = Dn
         cplBC%xo = cplBC%xn
!     Exiting outer loop if l1
         IF (l1) EXIT


      END DO
!     End of outer loop

      IF (resetSim) THEN
         CALL REMESHRESTART(timeP)
         DEALLOCATE(Ag, Yg, Dg, incL, res)
         IF (ALLOCATED(tls)) THEN
            DEALLOCATE(tls%ltg, tls%W, tls%R)
            DEALLOCATE(tls)
         END IF
         GOTO 101
      END IF

      IF (l1 .AND. saveAve) CALL CALCAVE

      DEALLOCATE(Ag, Yg, Dg, incL, res)
      CALL FINALIZE()
      CALL MPI_FINALIZE(ierr)

      END PROGRAM MAIN
!####################################################################
      SUBROUTINE STOPSIM()

      CALL FINALIZE
      STOP "MPI is forced to stop by a fatal error"

      END SUBROUTINE STOPSIM
!####################################################################
