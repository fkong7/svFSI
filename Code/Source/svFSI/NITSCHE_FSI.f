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
!     This is for constructing FSI equations on fluid and solid
!     domains.
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!
!     TODO Description
!
!--------------------------------------------------------------------
      SUBROUTINE CONSTRUCT_NITSCHE_FSI(lM, Ag, Yg, Dg, lD, lY)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo), lD(tDof, tnNo),  lY(tDof,tnNo)

      LOGICAL :: vmsStab, flag, computeBoth, vmsCIPStab, CIPStab
      INTEGER(KIND=IKIND) a, e, g, l, Ac, eNoN, cPhys, iFn, nFn, N
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd), FinvT(nsd,nsd)
      TYPE(fsType) :: fs(2), fsFace, fsFaceCur, fsFaceOpp


      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:), ptrOpp(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:), 
     2   dl(:,:), bfl(:,:), fN(:,:), pS0l(:,:), pSl(:), ya_l(:),
     3   lR(:,:), lK(:,:,:), lKd(:,:,:), xfacCur(:,:), xfacOpp(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xwl(:,:), xql(:,:), Nwx(:,:),
     2   Nwxx(:,:), Nqx(:,:), xp(:), xiCur(:,:), ylo(:,:)

      INTEGER(KIND=IKIND) maxSubTri, maxQuadPnt, bng, end, cnt, is, As,
     2   find, FlagToDel, i, j, maxLevel, eNoNF, eTypeF, nGF, idSolFac

      INTEGER(KIND=IKIND), ALLOCATABLE :: FlagSubElm(:), FlagLevel(:), 
     2                     FlagQP(:)
      REAL(KIND=RKIND), ALLOCATABLE :: lstSubElm(:,:,:), 
     2                              lstQdPntRef(:,:), lstQdPntCur(:,:) 
      REAL(KIND=RKIND) :: poly(nsd,msh(2)%eNoN),
     2                    xlToDel(nsd,msh(2)%eNoN), x4(nsd), x5(nsd), 
     3                    x6(nsd), JacSub, Nw(msh(1)%eNoN), 
     4                    Nq(msh(1)%eNoN)

      INTEGER(KIND=IKIND) :: f, b, faceID(msh(1)%eNoN, msh(1)%eNoN-1), 
     2                       eOpp, iFa, iFNit, isf, maxNbrSElm, nbrSElm, 
     3                       nGFTot, dofToJump(msh(1)%eNoN), findInt
      REAL(KIND=RKIND) :: Nfac(nsd), NCur(nsd), nrm, n1,n2,n3, dNU(nsd), 
     2                    dNUOpp(nsd), dNVb, dNVa, ghsP, visc, diamFace, 
     3                    gammag, measMetric, wt, diam, nitscheP, px(2),
     4                    pxOpp(2), CipP, linfVel, velF(nsd), nvelF, 
     5                    rho, ReF
      REAL(KIND=RKIND), ALLOCATABLE :: covBas(:,:), xFacOp(:,:), 
     2                   xFace(:,:), xlOpp(:,:), NwxOp(:,:), lROpp(:,:), 
     3                   lKCurOpp(:,:,:), lKOpp(:,:,:), lKOppCur(:,:,:),
     4                   ylOp(:,:), covMetric(:,:), FSubElm(:,:,:),
     5                   qpFlu(:,:), qpSol(:,:), xFRef(:,:)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptrSOl(:)

!     Parent shape function sub element 
      REAL(KIND=RKIND), ALLOCATABLE :: Nwf(:,:), Nws(:,:), lRS(:,:), 
     2        lKS(:,:,:), lKSF(:,:,:), lKFS(:,:,:), yls(:,:)      



      INTERFACE
         SUBROUTINE GET_SElmQuadPoint(nbrSElm, eTypeF, eNoNF, nG, nGF, 
     2                         crdFace, quadPnt, SElm, qPRef )
         USE COMMOD
         IMPLICIT NONE
         INTEGER(KIND=IKIND), INTENT(IN) :: nbrSElm, eTypeF,nG,nGF,eNoNF
         REAL(KIND=RKIND), INTENT(IN) :: crdFace(nsd, eNoNF)
         REAL(KIND=RKIND), INTENT(OUT) :: quadPnt(nsd, nGF)
         REAL(KIND=RKIND),INTENT(OUT),OPTIONAL :: SElm(nsd-1,
     2                                                    eNoNF,nbrSElm)
         REAL(KIND=RKIND),INTENT(OUT),OPTIONAL :: qPRef(nsd-1, nGF)
         END SUBROUTINE GET_SElmQuadPoint
      END INTERFACE


      eNoN = lM%eNoN
      nFn  = lM%nFn
      IF (nFn .EQ. 0) nFn = 1

      IF (lM%nFs .EQ. 1) THEN
         vmsStab = .TRUE.
         vmsCIPStab = .TRUE.
      ELSE
         vmsStab = .FALSE.
         vmsCIPStab = .FALSE.
      END IF

      CIPStab = .FALSE.

!     l = 3, if nsd==2 ; else 6;
      l = nsymd

      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), bfl(nsd,eNoN), fN(nsd,nFn), pS0l(nsymd,eNoN),
     3   pSl(nsymd), ya_l(eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN),
     4   lKd(dof*nsd,eNoN,eNoN), xp(nsd), ylo(tDof,eNoN))


!     Loop over all elements of mesh
      DO e=1, lM%nEl
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
            ylo(:,a) = lY(:,Ac)
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

!        To decomment to activate CIP
         IF (cPhys.EQ.phys_fluid) THEN 
            vmsCIPStab = .FALSE.
            CIPStab = .TRUE.
            CipP = 1.E-3
         END IF

!        Definition of different integration in this element based on physics 
         IF((cPhys.EQ.phys_fluid).AND.(intFElmFlag(e) .EQ. 2)) GOTO 20

!        Jump also for normal fluid to integration in subelem, nothing should change         
!         IF((cPhys.EQ.phys_fluid).AND.(intFElmFlag(e) .EQ. 1)) GOTO 20

!        Hidden fluid element
         IF((cPhys.EQ.phys_fluid).AND.(intFElmFlag(e) .EQ. 0)) GOTO 100

!-----------------------------------------------------------------------
!        Normal fluid or structure element 

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
                  CALL FLUID3D_M(vmsCIPStab, fs(1)%eNoN, fs(2)%eNoN, w,
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
                  !CALL FLUID2D_M_NoVMS(vmsCIPStab, fs(1)%eNoN, 
                  CALL FLUID2D_M(vmsCIPStab, fs(1)%eNoN, 
     2    fs(2)%eNoN, w, ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, 
     3          Nwxx, al, yl, bfl, lR, lK)

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
                  CALL FLUID3D_C(vmsCIPStab, fs(1)%eNoN, fs(2)%eNoN, w,
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
                  CALL FLUID2D_C(vmsCIPStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               CASE (phys_ustruct)
                  CALL USTRUCT2D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl,
     3               dl, bfl, lR, lK, lKd)
               END SELECT
            END IF
         END DO ! g: loop

         GOTO 100


!-----------------------------------------------------------------------
!        Fluid element cut by the Nitsche's boundary.
!        Need to perform special integration with finite cell 
!        Let consider for the moment fs(1)%nG = fs(2)%nG 
20       CONTINUE

C          write(*,*)" For element ", e, " cut-fem "
C          write(*,*)" Element coord: "
C          write(*,*) xl(:,1)
C          write(*,*) xl(:,2)
C          write(*,*) xl(:,3)


!-------------- No Bulk term integration 

C !        Gauss integration 1
C          DO g=1, fs(1)%nG
C             IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
C                CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
C      2            ksix)
C                IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
C             END IF

C             IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
C                CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
C      2            ksix)
C                IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

C                CALL GNNxx(l, fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g),
C      2            fs(1)%Nxx(:,:,g), xwl, Nwx, Nwxx)
C             END IF
C             w = fs(1)%w(g) * Jac

C             IF (nsd .EQ. 3) THEN
C                SELECT CASE (cPhys)
C                CASE (phys_fluid)
C                   CALL FLUID3D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
C      2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
C      3               al, yl, bfl, lR, lK)

C                CASE (phys_lElas)
C                   CALL LELAS3D(fs(1)%eNoN, w, fs(1)%N(:,g), Nwx, al, dl,
C      2               bfl, pS0l, pSl, lR, lK)

C                CASE (phys_struct)
C                   CALL STRUCT3D(fs(1)%eNoN, nFn, w, fs(1)%N(:,g), Nwx,
C      2               al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK)

C                CASE (phys_ustruct)
C                   CALL USTRUCT3D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, nFn,
C      2               w, Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, al, yl,
C      3               dl, bfl, fN, ya_l, lR, lK, lKd)

C                END SELECT

C             ELSE IF (nsd .EQ. 2) THEN
C                SELECT CASE (cPhys)
C                CASE (phys_fluid)
C                   CALL FLUID2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
C      2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
C      3               al, yl, bfl, lR, lK)

C                CASE (phys_lElas)
C                   CALL LELAS2D(fs(1)%eNoN, w, fs(1)%N(:,g), Nwx, al, dl,
C      2               bfl, pS0l, pSl, lR, lK)

C                CASE (phys_struct)
C                   CALL STRUCT2D(fs(1)%eNoN, nFn, w, fs(1)%N(:,g), Nwx,
C      2               al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK)

C                CASE (phys_ustruct)
C                   CALL USTRUCT2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, nFn,
C      2               w, Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, al, yl,
C      3               dl, bfl, fN, ya_l, lR, lK, lKd)

C                END SELECT
C             END IF
C          END DO ! g: loop

C !        Set function spaces for velocity and pressure.
C          CALL GETTHOODFS(fs, lM, vmsStab, 2)

C !        Gauss integration 2
C          DO g=1, fs(2)%nG
C             IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
C                CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
C      2            ksix)
C                IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
C             END IF

C             IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
C                CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
C      2            ksix)
C                IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
C             END IF
C             w = fs(2)%w(g) * Jac

C             IF (nsd .EQ. 3) THEN
C                SELECT CASE (cPhys)
C                CASE (phys_fluid)
C                   CALL FLUID3D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
C      2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
C      3               al, yl, bfl, lR, lK)

C                CASE (phys_ustruct)
C                   CALL USTRUCT3D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
C      2               Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl,
C      3               dl, bfl, lR, lK, lKd)
C                END SELECT

C             ELSE IF (nsd .EQ. 2) THEN
C                SELECT CASE (cPhys)
C                CASE (phys_fluid)
C                   CALL FLUID2D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
C      2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
C      3               al, yl, bfl, lR, lK)

C                CASE (phys_ustruct)
C                   CALL USTRUCT2D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
C      2               Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl,
C      3               dl, bfl, lR, lK, lKd)
C                END SELECT
C             END IF
C          END DO ! g: loop


!-------------- Cut Bulk term integration  
C          write(*,*)" Beginning integration cut cell "

         nbrSElm = 4**(nbrFntCllSD-1)
         maxQuadPnt = nbrSElm*fs(1)%nG

!        Allocate lstSubTri: list of sub element coordinates in the reference element 
         ALLOCATE ( lstSubElm(nsd,fs(1)%eNoN,nbrSElm) )
         ALLOCATE ( FlagSubElm(nbrSElm) )
         ALLOCATE ( FlagQP(maxQuadPnt) )
         ALLOCATE ( FlagLevel(nbrSElm) )
!        Allocate lstQdPnt: list of quad point in the reference element   
         ALLOCATE ( lstQdPntRef(nsd,maxQuadPnt) )
         ALLOCATE ( lstQdPntCur(nsd,maxQuadPnt) )


         lstSubElm  = -1._RKIND
         lstQdPntRef   = -1._RKIND
         lstQdPntCur   = -1._RKIND
         FlagSubElm = -1
         FlagQP = -1
         FlagLevel  = -1
         maxLevel = nbrFntCllSD


         CALL GET_BulkSElmQP(nbrSElm, lM%eType, eNoN, fs(1)%nG, 
     2              maxQuadPnt, xl, lstQdPntRef, lstQdPntCur, lstSubElm)

         CALL GET_FlagBulkSElm(lD, lM, fs(1)%nG, maxQuadPnt, nbrSElm,  
     2                                          lstQdPntCur, FlagSubElm)

         CALL GET_FlagBulkQP(lD, lM, fs(1)%nG, maxQuadPnt, nbrSElm,  
     2                                          lstQdPntCur, FlagQP)

C          write(*,*)" element e ", e 
C          write(*,*)" xl = ", xl 
C          write(*,*)" Jac Khat -> K = " , Jac

         a = 0
!        Loop over sub-element 
         DO i = 1, nbrSElm
!           If the sub-elem is hidden, just cycle 
            IF ( FlagSubElm(i) .EQ. 0) a = a + 1
         END DO
         IF( a .EQ. 0) THEN 
            write(*,*)" WARNING!! we don't integrate in element ", e
         END IF


!        Loop over sub-element 
         DO i = 1, nbrSElm

!           If the sub-elem is hidden, just cycle 
            IF ( FlagSubElm(i) .EQ. 0) CYCLE

            bng = fs(1)%nG*(i-1) 

            Nw = 0._RKIND
            Nq = 0._RKIND

!           Update test function and compute J K^hat -> K^tilde 
!           Jac from ref element to subElm in the reference element
!           check this 
C             CALL COMP_SubElmJAC(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), 
C      2              lstSubElm(:,:,i), JacSub) 

!           Loop over quad point 1
            DO g=1, fs(1)%nG

!              Nx are constant, we can use the previous function               
!              If this is not the case anymore, this needs to be modified 
!              fs(2)%lShpF is true here     

!              If quad point is inside any solid element, skip 
               IF( FlagQP(bng+g) .EQ. 0) CYCLE          

!              Jac from ref element to subElm in the reference element
               CALL COMP_SubElmJAC(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), 
     2              lstSubElm(:,:,i), JacSub)   

!              Evaluate test function in the new quadrature points
!              Quad point in the sub-elm of the ref element                   
               CALL GETGNN(nsd, lM%eType, fs(1)%eNoN,  
     2             lstQdPntRef(:,bng+g), Nw, Nwx)
               CALL GETGNN(nsd, lM%eType, fs(2)%eNoN,  
     2             lstQdPntRef(:,bng+g), Nq, Nqx)

               CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, 
     2                 Jac, ksix)
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx,
     2                 Jac, ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

C                write(*,*)" Nwx= " , Nwx
C                write(*,*)" lstSubElm = ", lstSubElm(:,:,i)
C                write(*,*)" lstQdPntRef = ", lstQdPntRef(:,bng+g)
C                write(*,*)" lstQdPntCur = ", lstQdPntCur(:,bng+g)
C                write(*,*)""
C                write(*,*)""



!              Compute new weight (w * J k_ref -> K_cut * J k_ref -> K_sub )
               w = fs(1)%w(g) * Jac * JacSub        

               IF (nsd .EQ. 3) THEN
                  CALL FLUID3D_M(vmsCIPStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, Nw, Nq, Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)
               ELSE IF (nsd .EQ. 2) THEN
                  !CALL FLUID2D_M_NoVMS(vmsCIPStab, fs(1)%eNoN, 
                  CALL FLUID2D_M(vmsCIPStab, fs(1)%eNoN, 
     2               fs(2)%eNoN, w, ksix, Nw, Nq, Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)
               END IF

               IF (nsd .EQ. 3) THEN
                  CALL FLUID3D_C(vmsCIPStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, Nw, Nq, Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               ELSE IF (nsd .EQ. 2) THEN
                  CALL FLUID2D_C(vmsCIPStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, Nw, Nq, Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               END IF


            END DO

         END DO

         
         DEALLOCATE( lstSubElm, FlagSubElm, FlagLevel, lstQdPntRef) 
         DEALLOCATE( lstQdPntCur, FlagQP ) 
         


!-------------- Nitsche's term segments intersection

C          write(*,*)" Fluid element e = ", e
C          write(*,*)" coord = ", xl
C          write(*,*)" Nbr of solid faces inside "
C          write(*,*) mapFElmSElm(e, 1, :)
C          write(*,*) mapFElmSElm(e, 2, :)
C          write(*,*) mapFElmSElm(e, 3, :)
!        Set function spaces for velocity and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 1)

!        Compute fluid element diameter h
         CALL COMP_ELM_DIAM( xl, eNoN, diam)

!        Define nitsche's coefficient gamma*mu/h
         nitscheP = eq(cEq)%dmn(cDmn)%prop(nitsche_param)
         visc = eq(cEq)%dmn(cDmn)%visc%mu_i

         nitscheP = nitscheP * visc / diam

         iFNit = 0
         maxNbrSElm = SIZE(mapFElmSElm, 3) 

         iFa = 1 
         eNoNF  = msh(2)%fa(iFa)%eNoN
         eTypeF = msh(2)%fa(iFa)%eType
         nGF    = msh(2)%fa(iFa)%nG

         ALLOCATE(qpFlu(nsd,nGF), FSubElm(nsd,eNoNF,1))
         ALLOCATE(qpSol(nsd-1, nGF))

         ALLOCATE(xFace(nsd,eNoNF),xFRef(nsd,eNoNF))
         ALLOCATE(Nws(eNoNF,nGF), Nwf(eNoN,nGF))

         ALLOCATE( lRS(dof,eNoNF), lKS(dof*dof,eNoNF,eNoNF),
     2   lKFS(dof*dof,eNoN,eNoNF), lKSF(dof*dof,eNoNF,eNoN), 
     3   yls(nsd,eNoNF), ptrSOl(eNoNF))

!        Loop over solid faces 
         DO iFa = 1, msh(2)%nFa

            IF( .NOT.msh(2)%fa(iFa)%isNts ) CYCLE

!           If the face is a nitsche's face continue 
            iFNit = iFNit + 1

C             write(*,*)" Nitsche face ", iFNit, " with elem ", e
C             write(*,*)" mapFElmSElm = " ,  mapFElmSElm(e, iFNit, :)


            qpFlu = 0._RKIND
            qpSol = 0._RKIND
            FSubElm = 0._RKIND
            xFace = 0._RKIND
            xFRef = 0._RKIND

            FinvT = 0._RKIND
            Jac = 0._RKIND

!           Loop over face elements 
            DO isf = 1, maxNbrSElm

               lRS  = 0._RKIND
               lKS  = 0._RKIND
               lKFS = 0._RKIND
               lKSF = 0._RKIND
               yls  = 0._RKIND

!              If the id != 0, continue, otherwise EXIT 
               IF( mapFElmSElm(e, iFNit, isf) .EQ. 0 ) EXIT  

               idSolFac = mapFElmSElm(e, iFNit, isf)


!              Get initial face coordinates  
               DO a = 1, eNoNF
                  Ac = msh(2)%fa(iFa)%IEN(a,idSolFac)
                  ptrSOl(a) = Ac
                  xFace(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
                  yls(:,a) = Yg(1:nsd,Ac)
               END DO

C                write(*,*)" Face coord : ", xFace
C                write(*,*)" ptrSOl : ", ptrSOl

!              Update fluid gradient (constant) wrt K^hat and compute J, F^-T
               g = 1

C                write(*,*)" fs(1)%Nx(:,:,g) = ", fs(1)%Nx(:,:,g)

               CALL GNN_FACE( fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xl, Nwx, 
     2            Jac, FinvT)     

C                write(*,*)" Nwx = ", Nwx        
C                write(*,*)" Jac = ", Jac        
C                write(*,*)" FinvT = ", FinvT        

!              Get solid face finite element space 
               CALL GETFSFace(msh(2), fsFace, vmsStab)

!              Map solid face to reference fluid elememnt and compute 
!              det F * | F^{-T} n_ref |, with the normal to the segment 
!              in reference configuration
               DO a = 1, eNoNF
                  CALL GET_PntCurToRef(lM, xl, xFace(:,a), xFRef(:,a))
               END DO

C                write(*,*)"xFace =  ", xFace
C                write(*,*)"xFRef =  ", xFRef

               IF ( nsd .EQ. 2 ) THEN

!                 Normal computation in the fluid reference elm  
                  n1 =   xFRef(2,2) - xFRef(2,1) 
                  n2 = - (xFRef(1,2) - xFRef(1,1))
                  nrm = SQRT(n1*n1 + n2*n2)
                  Nfac(1) = n1/nrm
                  Nfac(2) = n2/nrm

                  Nfac = - Nfac
C                 write(*,*)" Normal in the ref is : ", Nfac

!                 Normal computation in the fluid current elm 
                  n1 =   xFace(2,2) - xFace(2,1) 
                  n2 = - (xFace(1,2) - xFace(1,1))
                  nrm = SQRT(n1*n1 + n2*n2)
                  NCur(1) = n1/nrm
                  NCur(2) = n2/nrm

                  NCur = - NCur

C                   write(*,*)" xFace = " , xFace
C                   write(*,*)" Normal in the cur is : ", NCur

               ELSE IF ( nsd .EQ. 3 ) THEN
                  ALLOCATE(covBas(2,3))

!                 Normal computation in the fluid reference elm  
                  covBas(1,:) = xFRef(:,2) - xFRef(:,1) 
                  covBas(2,:) = xFRef(:,3) - xFRef(:,1) 

                  n1 = covBas(1,2)*covBas(2,3)-covBas(1,3)*covBas(2,2);
                  n2 = covBas(1,3)*covBas(2,1)-covBas(1,1)*covBas(2,3);
                  n3 = covBas(1,1)*covBas(2,2)-covBas(1,2)*covBas(2,1);

                  nrm = SQRT(n1*n1 + n2*n2 + n3*n3)
                  Nfac(1) = n1/nrm
                  Nfac(2) = n2/nrm
                  Nfac(3) = n3/nrm

                  Nfac = - Nfac

!                 Normal computation in the fluid current elm 
                  covBas(1,:) = xFace(:,2) - xFace(:,1) 
                  covBas(2,:) = xFace(:,3) - xFace(:,1) 

                  n1 = covBas(1,2)*covBas(2,3)-covBas(1,3)*covBas(2,2);
                  n2 = covBas(1,3)*covBas(2,1)-covBas(1,1)*covBas(2,3);
                  n3 = covBas(1,1)*covBas(2,2)-covBas(1,2)*covBas(2,1);

                  nrm = SQRT(n1*n1 + n2*n2 + n3*n3)
                  NCur(1) = n1/nrm
                  NCur(2) = n2/nrm
                  NCur(3) = n3/nrm

                  NCur = - NCur

                  DEALLOCATE(covBas)
               END IF

               wt = SQRT(NORM(MATMUL(FinvT,Nfac)))
C                write(*,*)" NCur = ", NCur
C                write(*,*)" Nfac = ", Nfac
C                write(*,*)" | F^-T n | = ", wt

               wt = Jac*SQRT(NORM(MATMUL(FinvT,Nfac)))

!              Get list of internal sub element and respective quad point 
!              in the reference fluid configuration 
               CALL CompInter(eTypeF,eNoNF,nGF,xFRef, qpFlu,  
     2             FSubElm(:,:,1), qpSol, findInt)        

C                write(*,*)" findInt = ", findInt

               IF( findInt .EQ. 0) CYCLE

C                write(*,*)" xFRef = ", xFRef
C                write(*,*)""
C                write(*,*) " Intersec solid ref fluid elm = "
C                write(*,*) " P1  = ", FSubElm(:,1,:)
C                write(*,*) " P2  = ", FSubElm(:,2,:)
C                write(*,*)""
C                write(*,*) "qpFlu = ", qpFlu
C                write(*,*)""
C                write(*,*) "qpSol = ", qpSol
C                write(*,*)""
C                write(*,*)""


!              Loop over sub-element 
               DO i = 1, 1

!                 Check that the sub-element is inside the fluid element 
C                   IF( FlagSubElm(i) .EQ. 0 ) CYCLE    

C                   write(*,*)" elem  ", e
C                   write(*,*)" xl  ", xl
C                   write(*,*)" xFace  ", xFace
C                   write(*,*)" subELm ", i
C                   write(*,*)" FSubElm(:,1,i) " , FSubElm(:,1,i)
C                   write(*,*)" FSubElm(:,2,i) " , FSubElm(:,2,i)


                  Nwf = 0._RKIND         
!                 Evaluate fluid nodal finite element function at the quad point 
                  DO g = 1, nGF
                     bng = nGF*(i-1)+g
                     CALL GETN(nsd,lM%eType,eNoN,qpFlu(:,bng),Nwf(:,g))

C                      write(*,*)" Quad point fluid " , qpFlu(:,bng)
C                      write(*,*)" value fluid test func ", Nwf(:,g)
                  END DO

!                 Compute metric a=F^T F and sqrt(a) for solid-sub element 
!                 Metric and measure camputations: sqrt(det (F_(S^ -> S)^T F(S^ -> S)))
                  IF ( nsd .EQ. 2 ) THEN
                     ALLOCATE(covMetric(2,1))

                     !covMetric(:,1) = 0.5_RKIND*(xFace(:,2)-xFace(:,1))

                     covMetric(:,1) = 0.5_RKIND * ( FSubElm(:,2,i)
     2                                                - FSubElm(:,1,i) )

                     measMetric = SQRT( covMetric(1,1)*covMetric(1,1) + 
     2                                covMetric(2,1)*covMetric(2,1) )

                     DEALLOCATE(covMetric)
                  ELSE IF ( nsd .EQ. 3 ) THEN
                     err = "TODO "
                  ELSE 
                     err = "Dimension not supported in Measure comp"
                  END IF

!                 Find solid quad point in the face reference element and 
!                 evaluate the solid fe test func there 
                  DO g = 1, nGF
                     bng = nGF*(i-1)+g
                     CALL GETN(nsd-1,eTypeF,eNoNF, qpSol(:,bng), 
     2                                                         Nws(:,g))
C                      write(*,*)" Quad point solid " , qpSol(:,bng)
C                      write(*,*)" value fluid test func ", Nws(:,g)
                  END DO                  


!                 We can now integrate considering as weights
C                   write(*,*)" integration between element fluid ", e 
C                   write(*,*)" and solid nodes ", ptrSOl 
C                   write(*,*)" FSubElm = ", FSubElm(:,1,i) 
C                   write(*,*)"           ", FSubElm(:,2,i) 

!                 Loop over solid quad point and assembly 
                  DO g = 1, nGF

C                      IF( FlagQP(i) .EQ. 0 ) CYCLE   

!                    det F * | F^{-T} n  | * | alpha'| * wq (wq ref face solid quad )
                     w = fsFace%w(g) * measMetric * wt

C                      write(*,*)"  fsFace%w(g) ", fsFace%w(g)
C                      write(*,*)"   w = ", w 

                     CALL NTS_ALLOC_2D( eNoN, eNoNF, w, Nwf(:,g), Nwx, 
     2                    Nws(:,g), yl, yls, lR, lK, lRS, lKS, lKFS, 
     3                                  lKSF, NCur, nitscheP)

                  END DO



               END DO

C                write(*,*)" lRS = ", lRS
C                write(*,*)" lKS = ", lKS
C                write(*,*)" lKFS = ", lKFS
C                write(*,*)" lKSF = ", lKSF

!              Call special assembly for the coupling fluid solid terms
               CALL NTS_ASSEMBLY(eNoN, eNoNF, ptr, ptrSOl, lRS, lKS, 
     2                                                     lKFS, lKFS) 


            END DO

         END DO

         DEALLOCATE(xFace, xFRef, qPSol, qpFlu, FSubElm,
     2       lRS, lKS, lKFS, lKSF, yls, Nwf, Nws, ptrSOl)

!--------------
         GOTO 100


!-------------- Nitsche's term finite cell

C          nbrFntCllSD = 7

C          write(*,*)" Fluid element e = ", e
C          write(*,*)" coord = ", xl
C          write(*,*)" Nbr of solid faces inside "
C          write(*,*) mapFElmSElm(e, 1, :)
C          write(*,*) mapFElmSElm(e, 2, :)
C          write(*,*) mapFElmSElm(e, 3, :)
!        Set function spaces for velocity and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 1)

!        Compute fluid element diameter h
         CALL COMP_ELM_DIAM( xl, eNoN, diam)

!        Define nitsche's coefficient gamma*mu/h
         nitscheP = eq(cEq)%dmn(cDmn)%prop(nitsche_param)
         visc = eq(cEq)%dmn(cDmn)%visc%mu_i

         nitscheP = nitscheP * visc / diam

         iFNit = 0
         maxNbrSElm = SIZE(mapFElmSElm, 3) 

         iFa = 1 
         eNoNF  = msh(2)%fa(iFa)%eNoN
         eTypeF = msh(2)%fa(iFa)%eType
         nGF    = msh(2)%fa(iFa)%nG

         IF ( nsd .EQ. 2) nbrSElm = 2**(nbrFntCllSD-1)
         IF ( nsd .EQ. 3) nbrSElm = 4**(nbrFntCllSD-1)
         nGFTot = nGF*nbrSElm

         ALLOCATE(qpFlu(nsd,nGFTot), FSubElm(nsd,eNoNF,nbrSElm))
         ALLOCATE(qpSol(nsd-1, nGFTot), FlagSubElm(nbrSElm))
         ALLOCATE(FlagQP(nGFTot))
         ALLOCATE(xFace(nsd,eNoNF),xFRef(nsd,eNoNF))
         ALLOCATE(Nws(eNoNF,nGF), Nwf(eNoN,nGF))

         ALLOCATE( lRS(dof,eNoNF), lKS(dof*dof,eNoNF,eNoNF),
     2   lKFS(dof*dof,eNoN,eNoNF), lKSF(dof*dof,eNoNF,eNoN), 
     3   yls(nsd,eNoNF), ptrSOl(eNoNF))

!        Loop over solid faces 
         DO iFa = 1, msh(2)%nFa

            IF( .NOT.msh(2)%fa(iFa)%isNts ) CYCLE

!           If the face is a nitsche's face continue 
            iFNit = iFNit + 1

C             write(*,*)" Nitsche face ", iFNit, " with elem ", e
C             write(*,*)" mapFElmSElm = " ,  mapFElmSElm(e, iFNit, :)


            qpFlu = 0._RKIND
            qpSol = 0._RKIND
            FSubElm = 0._RKIND
            FlagSubElm = -1
            xFace = 0._RKIND
            xFRef = 0._RKIND

            FinvT = 0._RKIND
            Jac = 0._RKIND

!           Loop over face elements 
            DO isf = 1, maxNbrSElm

               lRS  = 0._RKIND
               lKS  = 0._RKIND
               lKFS = 0._RKIND
               lKSF = 0._RKIND
               yls  = 0._RKIND

C                write(*,*)" isf = ", isf
C                write(*,*) mapFElmSElm(e, iFNit, :)

!              If the id != 0, continue, otherwise EXIT 
               IF( mapFElmSElm(e, iFNit, isf) .EQ. 0 ) EXIT  

               idSolFac = mapFElmSElm(e, iFNit, isf)

C                write(*,*) " NEED to do Nitsche here "

!              Get initial face coordinates  
               DO a = 1, eNoNF
                  Ac = msh(2)%fa(iFa)%IEN(a,idSolFac)
                  ptrSOl(a) = Ac
                  xFace(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
                  yls(:,a) = Yg(1:nsd,Ac)
               END DO

C                write(*,*)" Face coord : ", xFace
C                write(*,*)" ptrSOl : ", ptrSOl

!              Update fluid gradient (constant) wrt K^hat and compute J, F^-T
               g = 1

C                write(*,*)" fs(1)%Nx(:,:,g) = ", fs(1)%Nx(:,:,g)

               CALL GNN_FACE( fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xl, Nwx, 
     2            Jac, FinvT)              

!              Get solid face finite element space 
               CALL GETFSFace(msh(2), fsFace, vmsStab)

!              Map solid face to reference fluid elememnt and compute 
!              det F * | F^{-T} n_ref |, with the normal to the segment 
!              in reference configuration
               DO a = 1, eNoNF
                  CALL GET_PntCurToRef(lM, xl, xFace(:,a), xFRef(:,a))
               END DO

C                write(*,*)" elem ", e
C                write(*,*)" xl = ", xl
C                write(*,*)"xFace =  ", xFace
C                write(*,*)"xFRef =  ", xFRef

               IF ( nsd .EQ. 2 ) THEN

!                 Normal computation in the fluid reference elm  
                  n1 =   xFRef(2,2) - xFRef(2,1) 
                  n2 = - (xFRef(1,2) - xFRef(1,1))
                  nrm = SQRT(n1*n1 + n2*n2)
                  Nfac(1) = n1/nrm
                  Nfac(2) = n2/nrm

                  Nfac = - Nfac
C                 write(*,*)" Normal in the ref is : ", Nfac

!                 Normal computation in the fluid current elm 
                  n1 =   xFace(2,2) - xFace(2,1) 
                  n2 = - (xFace(1,2) - xFace(1,1))
                  nrm = SQRT(n1*n1 + n2*n2)
                  NCur(1) = n1/nrm
                  NCur(2) = n2/nrm

                  NCur = - NCur

C                   write(*,*)" xFace = " , xFace
C                   write(*,*)" Normal in the cur is : ", NCur

               ELSE IF ( nsd .EQ. 3 ) THEN
                  ALLOCATE(covBas(2,3))

!                 Normal computation in the fluid reference elm  
                  covBas(1,:) = xFRef(:,2) - xFRef(:,1) 
                  covBas(2,:) = xFRef(:,3) - xFRef(:,1) 

                  n1 = covBas(1,2)*covBas(2,3)-covBas(1,3)*covBas(2,2);
                  n2 = covBas(1,3)*covBas(2,1)-covBas(1,1)*covBas(2,3);
                  n3 = covBas(1,1)*covBas(2,2)-covBas(1,2)*covBas(2,1);

                  nrm = SQRT(n1*n1 + n2*n2 + n3*n3)
                  Nfac(1) = n1/nrm
                  Nfac(2) = n2/nrm
                  Nfac(3) = n3/nrm

                  Nfac = - Nfac

!                 Normal computation in the fluid current elm 
                  covBas(1,:) = xFace(:,2) - xFace(:,1) 
                  covBas(2,:) = xFace(:,3) - xFace(:,1) 

                  n1 = covBas(1,2)*covBas(2,3)-covBas(1,3)*covBas(2,2);
                  n2 = covBas(1,3)*covBas(2,1)-covBas(1,1)*covBas(2,3);
                  n3 = covBas(1,1)*covBas(2,2)-covBas(1,2)*covBas(2,1);

                  nrm = SQRT(n1*n1 + n2*n2 + n3*n3)
                  NCur(1) = n1/nrm
                  NCur(2) = n2/nrm
                  NCur(3) = n3/nrm

                  NCur = - NCur

                  DEALLOCATE(covBas)
               END IF

               wt = SQRT(NORM(MATMUL(FinvT,Nfac)))
C                write(*,*)" | F^-T n | = ", wt

               wt = Jac*SQRT(NORM(MATMUL(FinvT,Nfac)))

C               


!              Get list of internal sub element and respective quad point 
!              in the reference fluid configuration 
               CALL GET_SElmQuadPoint(nbrSElm, eTypeF, eNoNF, nGF,  
     2                         nGFTot, xFRef, qpFlu, FSubElm, qpSol )

C                write(*,*)" before GET_InternalQP "
               CALL GET_InternalSubElm(lM, nGF, nGFTot, nbrSElm, qpFlu, 
     2                                                    FlagSubElm)

               CALL GET_InternalQP(lM, nGF, nGFTot, nbrSElm, qpFlu, 
     2                                                    FlagQP)

C                write(*,*)" elem ", e
C                write(*,*)" xFRef = ", xFRef
C                write(*,*)""
C                write(*,*) "FlagSubElm =  ", FlagSubElm
C                write(*,*)""
C                write(*,*) "FSubElm = ", FSubElm
C                write(*,*)""
C                write(*,*) "qpFlu = ", qpFlu
C                write(*,*)""
C                write(*,*) "qpSol = ", qpSol



!              Loop over sub-element 
               DO i = 1, nbrSElm

!                 Check that the sub-element is inside the fluid element 
                  IF( FlagSubElm(i) .EQ. 0 ) CYCLE    

C                   write(*,*)" elem  ", e
C                   write(*,*)" xl  ", xl
C                   write(*,*)" xFace  ", xFace
C                   write(*,*)" subELm ", i
C                   write(*,*)" FSubElm(:,1,i) " , FSubElm(:,1,i)
C                   write(*,*)" FSubElm(:,2,i) " , FSubElm(:,2,i)


                  Nwf = 0._RKIND         
!                 Evaluate fluid nodal finite element function at the quad point 
                  DO g = 1, nGF
                     bng = nGF*(i-1)+g
                     CALL GETN(nsd,lM%eType,eNoN,qpFlu(:,bng),Nwf(:,g))

C                      write(*,*)" Quad point fluid " , qpFlu(:,bng)
C                      write(*,*)" value fluid test func ", Nwf(:,g)
                  END DO

!                 Compute metric a=F^T F and sqrt(a) for solid-sub element 
!                 Metric and measure camputations: sqrt(det (F_(S^ -> S)^T F(S^ -> S)))
                  IF ( nsd .EQ. 2 ) THEN
                     ALLOCATE(covMetric(2,1))

                     !covMetric(:,1) = 0.5_RKIND*(xFace(:,2)-xFace(:,1))

                     covMetric(:,1) = 0.5_RKIND * ( FSubElm(:,2,i)
     2                                                - FSubElm(:,1,i) )

                     measMetric = SQRT( covMetric(1,1)*covMetric(1,1) + 
     2                                covMetric(2,1)*covMetric(2,1) )

                     DEALLOCATE(covMetric)
                  ELSE IF ( nsd .EQ. 3 ) THEN
                     err = "TODO "
                  ELSE 
                     err = "Dimension not supported in Measure comp"
                  END IF

!                 Find solid quad point in the face reference element and 
!                 evaluate the solid fe test func there 
                  DO g = 1, nGF
                     bng = nGF*(i-1)+g
                     CALL GETN(nsd-1,eTypeF,eNoNF, qpSol(:,bng), 
     2                                                         Nws(:,g))
C                      write(*,*)" Quad point solid " , qpSol(:,bng)
C                      write(*,*)" value fluid test func ", Nws(:,g)
                  END DO                  


!                 We can now integrate considering as weights
C                   write(*,*)" integration between element fluid ", e 
C                   write(*,*)" and solid nodes ", ptrSOl 
C                   write(*,*)" FSubElm = ", FSubElm(:,1,i) 
C                   write(*,*)"           ", FSubElm(:,2,i) 

!                 Loop over solid quad point and assembly 
                  DO g = 1, nGF

                     IF( FlagQP(i) .EQ. 0 ) CYCLE   

!                    det F * | F^{-T} n  | * | alpha'| * wq (wq ref face solid quad )
                     w = fsFace%w(g) * measMetric * wt

C                      write(*,*)"  fsFace%w(g) ", fsFace%w(g)
C                      write(*,*)"   w = ", w 

                     CALL NTS_ALLOC_2D( eNoN, eNoNF, w, Nwf(:,g), Nwx, 
     2                    Nws(:,g), yl, yls, lR, lK, lRS, lKS, lKFS, 
     3                                  lKSF, NCur, nitscheP)

                  END DO



               END DO

C                write(*,*)" lRS = ", lRS
C                write(*,*)" lKS = ", lKS
C                write(*,*)" lKFS = ", lKFS
C                write(*,*)" lKSF = ", lKSF

!              Call special assembly for the coupling fluid solid terms
               CALL NTS_ASSEMBLY(eNoN, eNoNF, ptr, ptrSOl, lRS, lKS, 
     2                                                     lKFS, lKFS) 


            END DO

         END DO

         DEALLOCATE(xFace, xFRef, qPSol, qpFlu, FSubElm, FlagSubElm,
     2       FlagQP, lRS, lKS, lKFS, lKSF, yls, Nwf, Nws, ptrSOl)

!--------------


!-----------------------------------------------------------------------
100      CONTINUE





















!-------------- Ghost Penalty stabilization 
C          write(*,*)" Beginning ghost element ", e 

         ghsP = eq(cEq)%dmn(cDmn)%prop(ghostP_param)
         ghsP = 0._RKIND
         visc = eq(cEq)%dmn(cDmn)%visc%mu_i

         IF( ghsP .LT. 0._RKIND ) GOTO 101
         IF( cPhys .NE. phys_fluid ) GOTO 101

         ALLOCATE( xlOpp(nsd,fs(1)%eNoN), NwxOp(nsd,fs(1)%eNoN) )
         ALLOCATE( xFacOp(nsd,fs(1)%eNoN-1), xFace(nsd,fs(1)%eNoN-1) )
         ALLOCATE( lROpp(dof,eNoN), lKOpp(dof*dof,eNoN,eNoN) )
         ALLOCATE( lKCurOpp(dof*dof,eNoN,eNoN), ptrOpp(eNoN) )
         ALLOCATE( lKOppCur(dof*dof,eNoN,eNoN), ylOp(tDof,eNoN) )
         Nwx    = 0._RKIND
         NwxOp  = 0._RKIND

         IF ( nsd .EQ. 2 ) THEN
            faceID(1,:) = (/2, 3/)
            faceID(2,:) = (/3, 1/)
            faceID(3,:) = (/1, 2/)

C             faceIDOp(1,:) = (/3, 2/)
C             faceIDOp(2,:) = (/1, 3/)
C             faceIDOp(3,:) = (/2, 1/)
         ELSE IF ( nsd .EQ. 3 ) THEN
            faceID(1,:) = (/2, 3, 4/)
            faceID(2,:) = (/4, 3, 1/)
            faceID(3,:) = (/1, 2, 4/)
            faceID(4,:) = (/1, 3, 2/)

C             faceIDOp(1,:) = (/2, 4, 3/)
C             faceIDOp(2,:) = (/4, 1, 3/)
C             faceIDOp(3,:) = (/1, 4, 2/)
C             faceIDOp(4,:) = (/1, 2, 3/)
         ELSE 
            err = "Dimension not supported in GETNEIGH"
         END IF

         
!        xl(:,a)  = x(:,Ac) contains the coordinates of the current element 
!        Loop over faces 
         DO f = 1, eNoN

!           Initialize residue and tangents
            lROpp     = 0._RKIND
            lKOpp     = 0._RKIND
            lKCurOpp  = 0._RKIND
            lKOppCur  = 0._RKIND

!           Define the opposite element 
            eOpp = lM%neigh(e,f)

!           If the face is a boundary face, go to next face 
            IF( eOpp .EQ. -1 ) CYCLE

            computeBoth = .FALSE.

!-----      Do nothing in the hidden fluid element 
!           Compute also opposite if it is fluid and jump fluid elem  
            IF( (intFElmFlag(eOpp) .EQ. 1) ) computeBoth = .TRUE.
            IF( (intFElmFlag(e) .EQ. 1) ) CYCLE
!           
!           Do it only in intersected part with the 1 or 2         
C             IF( (intFElmFlag(eOpp) .EQ. 0) ) CYCLE
            IF( (intFElmFlag(e) .EQ. 0) ) CYCLE
            
!           Do not assemble test physical, unknown unphysical 
             IF( (intFElmFlag(eOpp) .EQ. 0) .AND. 
     2                 (intFElmFlag(e) .EQ. 2) ) CYCLE 


!-----      Do Ghost, but remove contribution unphysical-physical  
C !           Compute also opposite if it is fluid and jump fluid elem  
C             IF( (intFElmFlag(eOpp) .EQ. 1) ) computeBoth = .TRUE.
C             IF( (intFElmFlag(e) .EQ. 1) ) CYCLE
            
C !           Do not assemble test physical, unknown unphysical 
C              IF( (intFElmFlag(eOpp) .EQ. 0) .AND. 
C      2                 (intFElmFlag(e) .EQ. 2) ) CYCLE 




!           Get coordinates and velocity opposite element        
            DO a = 1, eNoN
               Ac = lM%IEN(a,eOpp)
               ptrOpp(a) = Ac
               xlOpp(:,a)  = x(:,Ac)

C                write(*,*) xlOpp(:,a)

               dl(:,a)     = Dg(:,Ac)
               ylOp(:,a)   = Yg(:,Ac)
            END DO
            xlOpp(:,:) = xlOpp(:,:) + dl(nsd+2:2*nsd+1,:)

C             write(*,*)" xl = ", xl
C             write(*,*)" xlOpp = ", xlOpp

            
!           Define dof to jump during assembly 
            IF( intFElmFlag(e) .EQ. 0) THEN 
               DO i=1, eNoN
                  Ac = lM%IEN(i,e)
                  dofToJump(i) = dofToJumpAss(Ac)
               END DO
            END IF


            CALL GETFSFace(lM, fsFaceCur, vmsStab)
            CALL GETFSFace(lM, fsFaceOpp, vmsStab)

!           No need to locate the surface reference quad point in the  
!           reference area/volume element because we are assuming P1 fem  
!           the gradient is therefore constant in the element.
!           We evaluate the grand in both currenrt and opposite 
!           just once because it will be constant over the element 
            g = 1
            CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xl, Nwx, Jac,
     2            ksix)
            CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xlOpp, NwxOp, 
     2            Jac, ksix)


C             write(*,*)" The gradient are ", Nwx
C             write(*,*)" For the opposite ", NwxOp


!           Get coordinated current and opposite face 
            DO a = 1, eNoN-1
               xFace(:,a) = xl(:,faceID(f,a))
C                xFacOp(:,a) = xl(:,faceIDOp(f,a))
            END DO

C             write(*,*)" face coordinates : "
C             write(*,*) xFace(:,1)
C             write(*,*) xFace(:,2)

!           Computing ghost penalty coefficient 
            CALL COMP_ELM_DIAM( xFace, eNoN-1, diamFace)

C             write(*,*)" diam face = ", diamFace

            gammag = ghsP * visc * diamFace

            wt   = eq(cEq)%af * eq(cEq)%gam * dt

!           Compute face Normal current element (Nopp = - Ncur)
            IF ( nsd .EQ. 2 ) THEN

               n1 =   xFace(2,2) - xFace(2,1) 
               n2 = - (xFace(1,2) - xFace(1,1))
               nrm = SQRT(n1*n1 + n2*n2)
               Nfac(1) = n1/nrm
               Nfac(2) = n2/nrm

C                write(*,*)" Normal is : ", Nfac

            ELSE IF ( nsd .EQ. 3 ) THEN
               ALLOCATE(covBas(2,3))

               covBas(1,:) = xFace(:,2) - xFace(:,1) 
               covBas(2,:) = xFace(:,3) - xFace(:,1) 

               n1 = covBas(1,2) * covBas(2,3) - covBas(1,3)*covBas(2,2);
               n2 = covBas(1,3) * covBas(2,1) - covBas(1,1)*covBas(2,3);
               n3 = covBas(1,1) * covBas(2,2) - covBas(1,2)*covBas(2,1);

               nrm = SQRT(n1*n1 + n2*n2 + n3*n3)
               Nfac(1) = n1/nrm
               Nfac(2) = n2/nrm
               Nfac(3) = n3/nrm

               DEALLOCATE(covBas)
            ELSE 
               err = "Dimension not supported in Normal comp"
            END IF

!           Metric and measure camputations: sqrt(det (F_(S^ -> S)^T F(S^ -> S)))
            IF ( nsd .EQ. 2 ) THEN
               ALLOCATE(covMetric(2,1))

               covMetric(:,1) = 0.5_RKIND*(xFace(:,2) - xFace(:,1))

               measMetric = SQRT( covMetric(1,1)*covMetric(1,1) + 
     2                            covMetric(2,1)*covMetric(2,1) )

               DEALLOCATE(covMetric)
            ELSE IF ( nsd .EQ. 3 ) THEN
               err = "TODO "
            ELSE 
               err = "Dimension not supported in Measure comp"
            END IF


!---------- Defining vectors for known velocity current and opposite  
            dNU = 0._RKIND
            dNUOpp = 0._RKIND
            DO a=1, eNoN
               dNU(1) = dNU(1) + Nwx(1,a)*yl(1,a)*Nfac(1) + 
     2                                       Nwx(2,a)*yl(1,a)*Nfac(2)

               dNU(2) = dNU(2) + Nwx(1,a)*yl(2,a)*Nfac(1) + 
     2                                       Nwx(2,a)*yl(2,a)*Nfac(2)
            END DO

            DO a=1, eNoN
               dNUOpp(1) = dNUOpp(1) + NwxOp(1,a)*ylOp(1,a)*Nfac(1)  
     2                                  + NwxOp(2,a)*ylOp(1,a)*Nfac(2)

               dNUOpp(2) = dNUOpp(2) + NwxOp(1,a)*ylOp(2,a)*Nfac(1) 
     2                                  + NwxOp(2,a)*ylOp(2,a)*Nfac(2)
            END DO

!---------- Loop over quadrature points 
            DO g=1, fsFaceCur%nG

C                write(*,*)" num quad = ", fsFaceCur%nG

               w = fsFaceCur%w(g) * measMetric

C                write(*,*)" Original weight is ", fsFaceCur%w(g) 
C                write(*,*)" New weight is ", w

!------------- row Current - col Current/Opposite assembly   
!              Local Residue current - current and current - opposite 
               DO a = 1, fs(1)%eNoN  

                  IF( dofToJump(a) .EQ. 1 ) CYCLE

                  dNVa = Nwx(1,a)*Nfac(1) + Nwx(2,a)*Nfac(2)
                  lR(1,a) = lR(1,a) + w * dNU(1) * dNVa
                  lR(2,a) = lR(2,a) + w * dNU(2) * dNVa

                  lR(1,a) = lR(1,a) - w * dNUOpp(1) * dNVa
                  lR(2,a) = lR(2,a) - w * dNUOpp(2) * dNVa

               END DO             

!              Local current - current Tangent matrices 
               DO a = 1, fs(1)%eNoN  ! row 
                  
                  IF( dofToJump(a) .EQ. 1 ) CYCLE

                  DO b = 1, fs(1)%eNoN ! col

                     dNVb = Nwx(1,b)*Nfac(1) + Nwx(2,b)*Nfac(2)
                     dNVa = Nwx(1,a)*Nfac(1) + Nwx(2,a)*Nfac(2)

                     lK(1,a,b) = lK(1,a,b) + w*wt *  dNVb * dNVa

                     lK(5,a,b) = lK(5,a,b) + w*wt *  dNVb * dNVa

                  END DO  
               END DO  

!              Local current - opposite Tangent matrices 
               DO a = 1, fs(1)%eNoN  ! row current 

                  IF( dofToJump(a) .EQ. 1 ) CYCLE

                  DO b = 1, fs(1)%eNoN ! col opposite

                     dNVb = NwxOp(1,b)*Nfac(1) + NwxOp(2,b)*Nfac(2)
                     dNVa = Nwx(1,a)*Nfac(1) + Nwx(2,a)*Nfac(2)

                     lKCurOpp(1,a,b) = lKCurOpp(1,a,b) - w*wt*dNVb*dNVa

                     lKCurOpp(5,a,b) = lKCurOpp(5,a,b) - w*wt*dNVb*dNVa

                  END DO  
               END DO  
       
!------------- row Current - col Current/Opposite assembly   
               IF( computeBoth ) THEN

                  w = fsFaceOpp%w(g) * measMetric

C                   write(*,*)" New weight Opp is ", w

!                 Local Residue opposite - opposite and opposite - current 
                  DO a = 1, fs(1)%eNoN  

                     dNVa = NwxOp(1,a)*Nfac(1) + NwxOp(2,a)*Nfac(2)
                     lROpp(1,a) = lROpp(1,a) - w * dNU(1) * dNVa
                     lROpp(2,a) = lROpp(2,a) - w * dNU(2) * dNVa

                     lROpp(1,a) = lROpp(1,a) + w * dNUOpp(1) * dNVa
                     lROpp(2,a) = lROpp(2,a) + w * dNUOpp(2) * dNVa

                  END DO             

!                 Local opposite - opposite Tangent matrices 
                  DO a = 1, fs(1)%eNoN  ! row 
                     DO b = 1, fs(1)%eNoN ! col

                        dNVb = NwxOp(1,b)*Nfac(1) + NwxOp(2,b)*Nfac(2)
                        dNVa = NwxOp(1,a)*Nfac(1) + NwxOp(2,a)*Nfac(2)

                        lKOpp(1,a,b) = lKOpp(1,a,b) + w*wt* dNVb * dNVa

                        lKOpp(5,a,b) = lKOpp(5,a,b) + w*wt* dNVb * dNVa

                     END DO  
                  END DO  

!                 Local opposite - current Tangent matrices 
                  DO a = 1, fs(1)%eNoN  ! row opposite 
                     DO b = 1, fs(1)%eNoN ! col current

                        dNVb = Nwx(1,b)*Nfac(1) + Nwx(2,b)*Nfac(2)
                        dNVa = NwxOp(1,a)*Nfac(1) + NwxOp(2,a)*Nfac(2)

                        lKOppCur(1,a,b)=lKOppCur(1,a,b) - w*wt*dNVb*dNVa

                        lKOppCur(5,a,b)=lKOppCur(5,a,b) - w*wt*dNVb*dNVa

                     END DO  
                  END DO  

               END IF

            END DO ! quad point 

            CALL DOASSEM_RC( eNoN, ptr, ptrOpp, lKCurOpp )
            IF( computeBoth ) THEN 
               CALL DOASSEM_Opp(eNoN,ptrOpp,ptr, lKOpp, lROpp, lKOppCur)
            END IF

C             write(*,*)" lK ", lK
C             write(*,*)" lR ", lR
C             write(*,*)" lKOpp ", lKOpp
C             write(*,*)" lROpp ", lROpp
C             write(*,*)" lKOppCur ", lKOppCur
C             write(*,*)" lKCurOpp ", lKCurOpp

         END DO ! faces 


         DEALLOCATE(xlOpp, NwxOp, xFacOp, xFace, lROpp, lKOpp, lKCurOpp, 
     2              lKOppCur, ptrOpp, ylOp)   

C          write(*,*)" ending ghost penalty "   
!-------------- 





































101      CONTINUE




!-------------- Continue Interior Penalty Stab 

         IF( CIPStab ) THEN 

C             write(*,*)" Beginning CIP stab elm ", e 
C             write(*,*)" lM%neigh(e,:)  = ", lM%neigh(e,:)

            ALLOCATE( xlOpp(nsd,fs(1)%eNoN), NwxOp(nsd,fs(1)%eNoN) )
            ALLOCATE( xFacOp(nsd,fs(1)%eNoN-1), xFace(nsd,fs(1)%eNoN-1))
            ALLOCATE( lKCurOpp(dof*dof,eNoN,eNoN), ptrOpp(eNoN) )
            ALLOCATE( ylOp(tDof,eNoN) )

            Nwx    = 0._RKIND
            NwxOp  = 0._RKIND
!           0 do not jump, 1 not assemble row of this dof 
            dofToJump = 0

            IF ( nsd .EQ. 2 ) THEN
               faceID(1,:) = (/2, 3/)
               faceID(2,:) = (/3, 1/)
               faceID(3,:) = (/1, 2/)

C             faceIDOp(1,:) = (/3, 2/)
C             faceIDOp(2,:) = (/1, 3/)
C             faceIDOp(3,:) = (/2, 1/)
            ELSE IF ( nsd .EQ. 3 ) THEN
               faceID(1,:) = (/2, 3, 4/)
               faceID(2,:) = (/4, 3, 1/)
               faceID(3,:) = (/1, 2, 4/)
               faceID(4,:) = (/1, 3, 2/)

C             faceIDOp(1,:) = (/2, 4, 3/)
C             faceIDOp(2,:) = (/4, 1, 3/)
C             faceIDOp(3,:) = (/1, 4, 2/)
C             faceIDOp(4,:) = (/1, 2, 3/)
            ELSE 
               err = "Dimension not supported in GETNEIGH"
            END IF

            
!           xl(:,a)  = x(:,Ac) contains the coordinates of the current element 
!           Loop over faces 
C             write(*,*)" beginning loop faces "

!           Define dof to jump during assembly 
            IF( intFElmFlag(e) .EQ. 0) THEN 
               DO i=1, eNoN
                  Ac = lM%IEN(i,e)
                  dofToJump(i) = dofToJumpAss(Ac)
               END DO
            END IF

!           Loop over faces
            DO f = 1, eNoN

!           Initialize residue and tangents
               lKCurOpp  = 0._RKIND

!           Define the opposite element 
               eOpp = lM%neigh(e,f)

!              If the face is a boundary face, go to next face 
               IF( eOpp .EQ. -1 ) CYCLE

!-----         Only compute CIP in the physical part 
!              Jump hidden          
               IF( intFElmFlag(e) .EQ. 0 ) CYCLE

!              Do not assemble test physical, unknown unphysical 
               IF( (intFElmFlag(e) .EQ. 2) .AND. 
     2                 (intFElmFlag(eOpp) .EQ. 0) ) CYCLE 


!-----      CIP in both physical and unphysical, but remove contrinution unphysical nodes  
C                IF( (intFElmFlag(eOpp) .EQ. 0) .AND. 
C      2                           (intFElmFlag(e) .EQ. 2)) CYCLE 

!-----
C                IF( (intFElmFlag(e) .EQ. 0) .AND. 
C      2                           (intFElmFlag(eOpp) .EQ. 2)) CYCLE 



C             IF( (intFElmFlag(e) .EQ. 0) .AND. 
C      2                 (intFElmFlag(eOpp) .EQ. 2) ) CYCLE 

!           Jump fluid and cutted         
C             IF( intFElmFlag(e) .EQ. 1 ) CYCLE
C             IF( intFElmFlag(e) .EQ. 2 ) CYCLE


C             write(*,*)" Beginning loop faces for opposite ", eOpp

   !           Get coordinates and velocity opposite element   
C             write(*,*)" Coordinates opposite: "
           
               DO a = 1, eNoN
                  Ac = lM%IEN(a,eOpp)
                  ptrOpp(a) = Ac
                  xlOpp(:,a)  = x(:,Ac)

C                write(*,*) xlOpp(:,a)

                  dl(:,a)     = Dg(:,Ac)
                  ylOp(:,a)   = Yg(:,Ac)
               END DO
               xlOpp(:,:) = xlOpp(:,:) + dl(nsd+2:2*nsd+1,:)

C             write(*,*)" xl = ", xl
C             write(*,*)" xlOpp = ", xlOpp


               CALL GETFSFace(lM, fsFaceCur, vmsStab)
               CALL GETFSFace(lM, fsFaceOpp, vmsStab)

   !           No need to locate the surface reference quad point in the  
   !           reference area/volume element because we are assuming P1 fem  
   !           the gradient is therefore constant in the element.
   !           We evaluate the grand in both currenrt and opposite 
   !           just once because it will be constant over the element 
               g = 1
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xl, Nwx, Jac,
     2            ksix)
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xlOpp, NwxOp, 
     2            Jac, ksix)


C             write(*,*)" The gradient are ", Nwx
C             write(*,*)" For the opposite ", NwxOp


   !           Get coordinated current and opposite face 
               DO a = 1, eNoN-1
                  xFace(:,a) = xl(:,faceID(f,a))
C                xFacOp(:,a) = xl(:,faceIDOp(f,a))
               END DO

C                write(*,*)" face coordinates : "
C                write(*,*) xFace(:,1)
C                write(*,*) xFace(:,2)

   !           Computing ghost penalty coefficient 
               CALL COMP_ELM_DIAM( xFace, eNoN-1, diamFace)

C             write(*,*)" diam face = ", diamFace

!               CipP = 1.e-2  ! eq(cEq)%dmn(cDmn)%prop(ghostP_param)
               visc = eq(cEq)%dmn(cDmn)%visc%mu_i
               rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)

               linfVel = 0._RKIND
               DO a = 1, eNoN-1

                  velF(:) = lY( 1:nsd,lM%IEN(faceID(f,a),e) )
C                   velF(:) = Yg( 1:nsd,lM%IEN(faceID(f,a),e) )
C                   write(*,*)" vel is ", velF
                  nvelF = SQRT(NORM(velF))
                  IF ( nvelF .GT. linfVel ) linfVel = nvelF
               END DO 

C                write(*,*)" linfVel = ", linfVel 

               ReF = rho * linfVel * diamFace / visc

C                write(*,*)" ReF ", ReF

               IF( ReF .LT. 1._RKIND ) THEN 
                  CipP = CipP * (diamFace**3) *  rho / visc
               ELSE 
                  CipP = CipP * (diamFace**2) / linfVel
               END IF

C                write(*,*)" CipP = " , CipP

               wt   = eq(cEq)%af * eq(cEq)%gam * dt

C    !           Compute face Normal current element (Nopp = - Ncur)
C                IF ( nsd .EQ. 2 ) THEN

C                   n1 =   xFace(2,2) - xFace(2,1) 
C                   n2 = - (xFace(1,2) - xFace(1,1))
C                   nrm = SQRT(n1*n1 + n2*n2)
C                   Nfac(1) = n1/nrm
C                   Nfac(2) = n2/nrm

C C                write(*,*)" Normal is : ", Nfac

C                ELSE IF ( nsd .EQ. 3 ) THEN
C                   ALLOCATE(covBas(2,3))

C                   covBas(1,:) = xFace(:,2) - xFace(:,1) 
C                   covBas(2,:) = xFace(:,3) - xFace(:,1) 

C                   n1 = covBas(1,2)*covBas(2,3)-covBas(1,3)*covBas(2,2);
C                   n2 = covBas(1,3)*covBas(2,1)-covBas(1,1)*covBas(2,3);
C                   n3 = covBas(1,1)*covBas(2,2)-covBas(1,2)*covBas(2,1);

C                   nrm = SQRT(n1*n1 + n2*n2 + n3*n3)
C                   Nfac(1) = n1/nrm
C                   Nfac(2) = n2/nrm
C                   Nfac(3) = n3/nrm

C                   DEALLOCATE(covBas)
C                ELSE 
C                   err = "Dimension not supported in Normal comp"
C                END IF

   !           Metric and measure camputations: sqrt(det (F_(S^ -> S)^T F(S^ -> S)))
               IF ( nsd .EQ. 2 ) THEN
                  ALLOCATE(covMetric(2,1))

                  covMetric(:,1) = 0.5_RKIND*(xFace(:,2) - xFace(:,1))

                  measMetric = SQRT( covMetric(1,1)*covMetric(1,1) + 
     2                            covMetric(2,1)*covMetric(2,1) )

                  DEALLOCATE(covMetric)
               ELSE IF ( nsd .EQ. 3 ) THEN
                  err = "TODO "
               ELSE 
                  err = "Dimension not supported in Measure comp"
               END IF








   !---------- Defining vectors for known pressure grad   

               px = 0._RKIND
               DO a=1, eNoN
                  px(1) = px(1) + Nwx(1,a)*yl(3,a)
                  px(2) = px(2) + Nwx(2,a)*yl(3,a)
               END DO

               pxOpp = 0._RKIND
               DO a=1, eNoN
                  pxOpp(1) = pxOpp(1) + NwxOp(1,a)*ylOp(3,a)
                  pxOpp(2) = pxOpp(2) + NwxOp(2,a)*ylOp(3,a)
               END DO


   !---------- Loop over quadrature points 
               DO g=1, fsFaceCur%nG

C                write(*,*)" num quad = ", fsFaceCur%nG

                  w = fsFaceCur%w(g) * measMetric

C                write(*,*)" Original weight is ", fsFaceCur%w(g) 
C                write(*,*)" New weight is ", w

   !------------- row Current - col Current/Opposite assembly   
   !              Local Residue current - current and current - opposite 
                  DO a = 1, fs(1)%eNoN  

                     IF( dofToJump(a) .EQ. 1 ) CYCLE

                     lR(3,a) = lR(3,a) + CipP * w * px(1) * Nwx(1,a)
                     lR(3,a) = lR(3,a) + CipP * w * px(2) * Nwx(2,a)

                     lR(3,a) = lR(3,a) - CipP * w * pxOpp(1) * Nwx(1,a)
                     lR(3,a) = lR(3,a) - CipP * w * pxOpp(2) * Nwx(2,a)

                  END DO             

   !              Local current - current Tangent matrices 
                  DO a = 1, fs(1)%eNoN  ! row 

                     IF( dofToJump(a) .EQ. 1 ) CYCLE

                     DO b = 1, fs(1)%eNoN ! col

                        lK(9,a,b) = 
     2                   lK(9,a,b) + CipP * w * wt * Nwx(1,b) * Nwx(1,a)

                        lK(9,a,b) = 
     2                   lK(9,a,b) + CipP * w * wt * Nwx(2,b) * Nwx(2,a)

                     END DO  
                  END DO  

   !              Local current - opposite Tangent matrices 
                  DO a = 1, fs(1)%eNoN  ! row current 
                     
                     IF( dofToJump(a) .EQ. 1 ) CYCLE

                     DO b = 1, fs(1)%eNoN ! col opposite

                        lKCurOpp(9,a,b) = 
     2                  lKCurOpp(9,a,b) - CipP*w*wt*NwxOp(1,b)*Nwx(1,a)

                        lKCurOpp(9,a,b) = 
     2                  lKCurOpp(9,a,b) - CipP*w*wt*NwxOp(2,b)*Nwx(2,a)

                     END DO  
                  END DO  
          
               END DO ! quad point 

C                write(*,*)" before assembly"
C                write(*,*)" ptr " , ptr
C                write(*,*)" ptrOpp " , ptrOpp
C                write(*,*)" lKCurOpp " , lKCurOpp
               
               CALL DOASSEM_RC( eNoN, ptr, ptrOpp, lKCurOpp )


C                IF( computeBoth ) THEN 
C                   CALL DOASSEM_Opp(eNoN,ptrOpp,ptr, lKOpp, lROpp, lKOppCur)
C                END IF

C              write(*,*)" lK ", lK
C              write(*,*)" lR ", lR

            END DO ! faces 


            DEALLOCATE(xlOpp, NwxOp, xFacOp, xFace,lKCurOpp,ptrOpp,ylOp)   





         END IF

!-------------- 
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
     2   lKd, ylo)

      CALL DESTROY(fs(1))
      CALL DESTROY(fs(2))

!     This need to be removed once we have a ghost penalty stabilization 
!      CALL CLEAN_GHOSTPNT

      RETURN
      END SUBROUTINE CONSTRUCT_NITSCHE_FSI
!####################################################################
      SUBROUTINE NTS_ALLOC_2D(eNoN, eNoNF, w, Nw, Nwx, Nws, yl, yls, lR, 
     2                    lK, lRS, lKS, lKFS, lKSF, N, coef)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, eNoNF
      
      REAL(KIND=RKIND), INTENT(IN) :: w, Nw(eNoN), Nws(eNoNF),
     2  Nwx(nsd,eNoN), yl(tDof,eNoN), yls(nsd,eNoNF), N(nsd), coef
      
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2  lK(dof*dof,eNoN,eNoN), lRS(dof,eNoNF), lKS(dof*dof,eNoNF,eNoNF), 
     3  lKFS(dof*dof,eNon,eNoNF), lKSF(dof*dof,eNoNF,eNoN)


      INTEGER(KIND=IKIND) :: a, b
      REAL(KIND=RKIND) :: u(nsd), p, dd(nsd), Rnts(nsd), ux(2,2), 
     2                    es(2,2), esN(2), vis, wt, st(2), T1, uN, rho

C       write(*,*)" Inside NTS_ALLOC_2D "
C       write(*,*)" the normal is N ", N

      u   = 0._RKIND
      ux  = 0._RKIND
      p   = 0._RKIND
      p   = 0._RKIND
      dd  = 0._RKIND
      st  = 0._RKIND
      vis = eq(1)%dmn(cDmn)%visc%mu_i
      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)

      wt  = eq(1)%af * eq(1)%gam * dt

      DO a=1, eNoN
         u(1) = u(1) + Nw(a)*yl(1,a)
         u(2) = u(2) + Nw(a)*yl(2,a)

         p = p + Nw(a)*yl(3,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(1,a)*yl(2,a)
         ux(1,2) = ux(1,2) + Nwx(2,a)*yl(1,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)
      END DO

C       write(*,*)" u = ", u

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,1) = ux(2,1) + ux(1,2)
      es(2,2) = ux(2,2) + ux(2,2)
      es(1,2) = es(2,1)

      esN(1) = vis*(es(1,1)*N(1) + es(1,2)*N(2))
      esN(2) = vis*(es(2,1)*N(1) + es(2,2)*N(2))

      uN = (u(1)*N(1) + u(2)*N(2))*0.5_RKIND

      st(1) = esN(1) - p*N(1) 
      st(2) = esN(2) - p*N(2)

C       write(*,*)" second loop solid vel " 

      DO a=1, eNoNF
         dd(1) = dd(1) + Nws(a)*yls(1,a)
         dd(2) = dd(2) + Nws(a)*yls(2,a)
      END DO

!------------- row Fluid - col Fluid  
C       write(*,*)" Assembly F - F "
!     Local Residue 
      DO a = 1, eNoN  

         Rnts(1) = coef*(u(1)-dd(1)) - st(1)
         Rnts(2) = coef*(u(2)-dd(2)) - st(2)
         
         lR(1,a) = lR(1,a) + w * Rnts(1) * Nw(a)
         lR(2,a) = lR(2,a) + w * Rnts(2) * Nw(a)  

!        Adding extra terms for Temam's trick  
C          lR(1,a) = lR(1,a) - rho * w * uN * u(1) * Nw(a)            
C          lR(2,a) = lR(2,a) - rho * w * uN * u(2) * Nw(a)             

      END DO

!     Local Tangent Matrix 
      DO a = 1, eNoN ! row 
         DO b = 1, eNoN ! col

!           Nitsche penalty terms (u1.v1 + u2.v2)
            lK(1,a,b) = lK(1,a,b) + w*wt * coef * Nw(a)*Nw(b)  
            lK(5,a,b) = lK(5,a,b) + w*wt * coef * Nw(a)*Nw(b)  

!           Stability Nitsche's terms ( sigma(u,p) n . v ) 
!           Part related to v1
            T1 = 2*vis*Nwx(1,b)*N(1) + vis*Nwx(2,b)*N(2)  
            lK(1,a,b) = lK(1,a,b) - w*wt * T1 * Nw(a) 

            T1 = vis*Nwx(1,b)*N(2) 
            lK(2,a,b) = lK(2,a,b) - w*wt * T1 * Nw(a)

            T1 = - Nw(b)*N(1) 
            lK(3,a,b) = lK(3,a,b) - w*wt * T1 * Nw(a)
!           Part related to v2
            T1 = vis*Nwx(2,b)*N(1) 
            lK(4,a,b) = lK(4,a,b) - w*wt * T1 * Nw(a) 

            T1 = 2*vis*Nwx(2,b)*N(2) + vis*Nwx(1,b)*N(1)              
            lK(5,a,b) = lK(5,a,b) - w*wt * T1 * Nw(a)

            T1 = - Nw(b)*N(2) 
            lK(6,a,b) = lK(6,a,b) - w*wt * T1 * Nw(a)


!           Adding extra terms for Temam's trick  
C             lK(1,a,b) = lK(1,a,b) - rho * w*wt * uN * Nw(b) * Nw(a)            
C             lK(5,a,b) = lK(5,a,b) - rho * w*wt * uN * Nw(b) * Nw(a)  

         END DO
      END DO
!-------------  

!------------- row Fluid - col Solid  
C       write(*,*)" Assembly F - S "
!     Local Tangent Matrix 
      DO a = 1, eNoN ! row 
         DO b = 1, eNoNF ! col

!           Nitsche penalty terms (-ddot1.v1 - ddot2.v2)
            lKFS(1,a,b) = lKFS(1,a,b) - w*wt * coef * Nw(a)*Nws(b)  
            lKFS(5,a,b) = lKFS(5,a,b) - w*wt * coef * Nw(a)*Nws(b)              


         END DO
      END DO

!-------------  

!------------- row Solid - col Solid  
C       write(*,*)" Assembly S - S "

!     Local Residue 
      DO a = 1, eNoNF  
         Rnts(1) = - coef*(u(1)-dd(1)) + st(1)
         Rnts(2) = - coef*(u(2)-dd(2)) + st(2)
         
         lRS(1,a) = lRS(1,a) + w * Rnts(1) * Nws(a)
         lRS(2,a) = lRS(2,a) + w * Rnts(2) * Nws(a)       
      END DO

!     Local Tangent Matrix 
      DO a = 1, eNoNF ! row 
         DO b = 1, eNoNF ! col

!           Nitsche penalty terms (ddot1.w1 + ddot2.w2)
            lKS(1,a,b) = lKS(1,a,b) + w*wt * coef * Nws(a)*Nws(b)  
            lKS(5,a,b) = lKS(4,a,b) + w*wt * coef * Nws(a)*Nws(b)    

            
         END DO
      END DO
!-------------  

!------------- row Solid - col Fluid  
C       write(*,*)" Assembly S - F "

!     Local Tangent Matrix 
      DO a = 1, eNoNF ! row 
         DO b = 1, eNoN ! col

!           Nitsche penalty terms (u1.w1 + u2.w2)
            lKSF(1,a,b) = lKSF(1,a,b) - w*wt * coef * Nws(a)*Nw(b)  
            lKSF(5,a,b) = lKSF(5,a,b) - w*wt * coef * Nws(a)*Nw(b) 

!           Stability Nitsche's terms ( sigma(u,p) n . w ) 
!           Part related to w1
            T1 = 2*vis*Nwx(1,b)*N(1) + vis*Nwx(2,b)*N(2)  
            lKSF(1,a,b) = lKSF(1,a,b) + w*wt * T1 * Nws(a) 

            T1 = vis*Nwx(1,b)*N(2) 
            lKSF(2,a,b) = lKSF(2,a,b) + w*wt * T1 * Nws(a)

            T1 = - Nw(b)*N(1) 
            lKSF(3,a,b) = lKSF(3,a,b) + w*wt * T1 * Nws(a)
!           Part related to w2
            T1 = vis*Nwx(2,b)*N(1) 
            lKSF(4,a,b) = lKSF(4,a,b) + w*wt * T1 * Nws(a) 

            T1 = 2*vis*Nwx(2,b)*N(2) + vis*Nwx(1,b)*N(1)              
            lKSF(5,a,b) = lKSF(5,a,b) + w*wt * T1 * Nws(a)

            T1 = - Nw(b)*N(2) 
            lKSF(6,a,b) = lKSF(6,a,b) + w*wt * T1 * Nws(a)


         END DO
      END DO
!-------------  

      RETURN
      END SUBROUTINE NTS_ALLOC_2D

!####################################################################

      SUBROUTINE GET_InternalSubElm(lM, nGF, nGFTot, nbrSElm, xFRef, 
     2    FlagSubElm)
      USE COMMOD
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: nGF, nGFTot, nbrSElm
      REAL(KIND=RKIND), INTENT(IN) :: xFRef(nsd, nGFTot)
      INTEGER(KIND=IKIND), INTENT(OUT) :: FlagSubElm(nbrSElm)


      REAL(KIND=RKIND) :: poly(nsd,lM%eNoN)
      INTEGER(KIND=IKIND) :: bng, cnt, find, g, i


!     Store only info of internal (to the reference fluid element) quadrature point
      SELECT CASE(lM%eType)
      CASE(eType_TRI3)
         poly(1,1) = 1._RKIND
         poly(2,1) = 0._RKIND

         poly(1,2) = 0._RKIND
         poly(2,2) = 1._RKIND

         poly(1,3) = 0._RKIND
         poly(2,3) = 0._RKIND

      CASE(eType_TET4)
         poly(1,1) = 1._RKIND
         poly(2,1) = 0._RKIND
         poly(3,1) = 0._RKIND

         poly(1,2) = 0._RKIND
         poly(2,2) = 1._RKIND
         poly(3,2) = 0._RKIND

         poly(1,3) = 0._RKIND
         poly(2,3) = 0._RKIND
         poly(3,3) = 1._RKIND

         poly(1,4) = 0._RKIND
         poly(2,4) = 0._RKIND
         poly(3,4) = 0._RKIND

      END SELECT


      DO i = 1, nbrSElm

!        If outside just delete the sub element
!        If inisde keep it, and go to the next one 
!        If partially inside subdivide the element 
         find = 0
         cnt = 0
         DO g = 1, nGF 
            bng = nGF*(i-1)+g 

!           Is the fluid node inside any solid element?
            CALL IN_POLY(xFRef(:,bng),poly, find)

            IF (find .EQ. 1) THEN 
               cnt = cnt + 1
            END IF

         END DO 

C          write(*,*)" nbr inside quad point = ", cnt 

!        FlagSubElm = 1: inside fluid element, 0: outside fluid elm, 
!        2: to be divided
         IF ( cnt .EQ. nGF ) THEN 
            FlagSubElm(i) = 1
C              write(*,*)" Flag sub elem ", i , "is inside "
         ELSE IF ( cnt .EQ. 0 ) THEN 
            FlagSubElm(i) = 0
C              write(*,*)" Flag sub elem ", i , "is outside "
         ELSE 
            FlagSubElm(i) = 1 !we consider outside the seg with only one qp inside 
C              write(*,*)" Flag sub elem ", i , "is partially inside "
         END IF

      END DO


      RETURN
      END SUBROUTINE GET_InternalSubElm

!----------------------------------------

      SUBROUTINE GET_InternalQP(lM, nGF, nGFTot, nbrSElm, xFRef, 
     2    FlagQP)
      USE COMMOD
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: nGF, nGFTot, nbrSElm
      REAL(KIND=RKIND), INTENT(IN) :: xFRef(nsd, nGFTot)
      INTEGER(KIND=IKIND), INTENT(OUT) :: FlagQP(nGFTot)


      REAL(KIND=RKIND) :: poly(nsd,lM%eNoN)
      INTEGER(KIND=IKIND) :: bng, cnt, find, g, i

      FlagQP = 0

!     Store only info of internal (to the reference fluid element) quadrature point
      SELECT CASE(lM%eType)
      CASE(eType_TRI3)
         poly(1,1) = 1._RKIND
         poly(2,1) = 0._RKIND

         poly(1,2) = 0._RKIND
         poly(2,2) = 1._RKIND

         poly(1,3) = 0._RKIND
         poly(2,3) = 0._RKIND

      CASE(eType_TET4)
         poly(1,1) = 1._RKIND
         poly(2,1) = 0._RKIND
         poly(3,1) = 0._RKIND

         poly(1,2) = 0._RKIND
         poly(2,2) = 1._RKIND
         poly(3,2) = 0._RKIND

         poly(1,3) = 0._RKIND
         poly(2,3) = 0._RKIND
         poly(3,3) = 1._RKIND

         poly(1,4) = 0._RKIND
         poly(2,4) = 0._RKIND
         poly(3,4) = 0._RKIND

      END SELECT


      DO i = 1, nbrSElm

!        If outside just delete the sub element
!        If inisde keep it, and go to the next one 
!        If partially inside subdivide the element 
         find = 0
         cnt = 0
         DO g = 1, nGF 
            bng = nGF*(i-1)+g 

!           Is the fluid node inside any solid element?
            CALL IN_POLY(xFRef(:,bng),poly, find)

            IF (find .EQ. 1) THEN 
               FlagQP(bng) = 1
            END IF

         END DO 

      END DO


      RETURN
      END SUBROUTINE GET_InternalQP

!####################################################################

      SUBROUTINE GET_FlagBulkSElm(lD, lM, nG, nGTot, nbrSElm, xQPCur, 
     2    FlagSubElm)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof, tnNo)
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: nG, nGTot, nbrSElm
      REAL(KIND=RKIND), INTENT(IN) :: xQPCur(nsd, nGTot)
      INTEGER(KIND=IKIND), INTENT(OUT) :: FlagSubElm(nbrSElm)


      REAL(KIND=RKIND) :: poly(nsd, msh(2)%eNoN)
      INTEGER(KIND=IKIND) :: bng, cnt, find, g, i, a, As, is


      DO i = 1, nbrSElm

C          write(*,*)" subelem ", i

         cnt = 0
         DO g = 1, nG 

            !write(*,*)" quad point ", g
            bng = nG*(i-1) + g 

            DO a = 1, msh(2)%nEl

               DO is = 1, msh(2)%eNoN
                  As = msh(2)%IEN(is,a)
                  poly(:,is) = x(:,As) + lD(nsd+2:2*nsd+1,As)
               END DO

C                write(*,*)" Localize quad point ", xQPCur(:,bng)
!                 Is the fluid node inside any solid element?
               CALL IN_POLY(xQPCur(:,bng),poly, find)

               !write(*,*)" find it in elm sol ", a , " ? ", find

               IF (find .EQ. 1) THEN 
                  cnt = cnt + 1
C                   write(*,*)" quad point ", g, " inside elm solid ", a
                  EXIT
               END IF

            END DO 
         END DO 

!        FlagSubElm = 1: fluid element, 0: under solid 
         IF ( cnt .GE. 3 ) THEN 
            FlagSubElm(i) = 0
C             write(*,*)" cnt = ", cnt
C           write(*,*)" Flag sub elem ", i , "is inside "
         ELSE 
            FlagSubElm(i) = 1
C             write(*,*)" cnt = ", cnt
         END IF

      END DO


      RETURN
      END SUBROUTINE GET_FlagBulkSElm

!-----------------------------------------------------------------------         
      SUBROUTINE GET_FlagBulkQP(lD, lM, nG, nGTot, nbrSElm, xQPCur, 
     2    FlagQP)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof, tnNo)
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: nG, nGTot, nbrSElm
      REAL(KIND=RKIND), INTENT(IN) :: xQPCur(nsd, nGTot)
      INTEGER(KIND=IKIND), INTENT(OUT) :: FlagQP(nGTot)


      REAL(KIND=RKIND) :: poly(nsd, msh(2)%eNoN)
      INTEGER(KIND=IKIND) :: bng, cnt, find, g, i, a, As, is


      FlagQP = 1

      DO g = 1, nGTot 

         DO a = 1, msh(2)%nEl

            DO is = 1, msh(2)%eNoN
               As = msh(2)%IEN(is,a)
               poly(:,is) = x(:,As) + lD(nsd+2:2*nsd+1,As)
            END DO

C           write(*,*)" Localize quad point ", xQPCur(:,bng)
!           Is the fluid node inside any solid element?
            CALL IN_POLY(xQPCur(:,g),poly, find)

            !write(*,*)" find it in elm sol ", a , " ? ", find

            IF (find .EQ. 1) FlagQP(g) = 0
               
         END DO 
      END DO 

      RETURN
      END SUBROUTINE GET_FlagBulkQP


!####################################################################
!####################################################################

      SUBROUTINE CLEAN_GHOSTPNT
      USE TYPEMOD
      USE COMMOD
      IMPLICIT NONE

      INTEGER(KIND=IKIND) a, Ac, ptr, rowN, colN, left, right

      DO a = 1, msh(1)%nNo
         Ac = msh(1)%gN(a)
         IF( ghostFNd(a) .EQ. 0) THEN
            R(:,Ac) = 0._RKIND

            left  = rowPtr(Ac)
            right = rowPtr(Ac+1)

            Val(:,left:right) = 0._RKIND

            ptr   = (right + left)/2
            DO WHILE (Ac .NE. colPtr(ptr))
               IF (Ac .GT. colPtr(ptr)) THEN
                  left  = ptr
               ELSE
                  right = ptr
               END IF
               ptr = (right + left)/2
            END DO
            Val(:,ptr) = 1._RKIND
         END IF
      END DO

      RETURN
      END SUBROUTINE CLEAN_GHOSTPNT

!####################################################################
!####################################################################
!--------------------------------------------------------------------
!
!     Giving the coordinates of the reference quad point in the 
!     current element with vertices xl
!
!--------------------------------------------------------------------
      SUBROUTINE GET_CURQUAD(lM, fs, xl, xiCur)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      TYPE(fsType), INTENT(IN) :: fs
      REAL(KIND=RKIND), INTENT(IN) :: xl(nsd,lM%eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: xiCur(nsd, fs%nG)

      INTEGER(KIND=IKIND) :: g
      REAL(KIND=RKIND) :: Def(nsd,nsd), rhsD(nsd) 

!     Loop over Gauss quad point
      DO g=1, fs%nG

         IF (nsd .EQ. 2) THEN
            Def(1,1) = xl(1,1) - xl(1,3) 
            Def(1,2) = xl(1,2) - xl(1,3) 
            Def(2,1) = xl(2,1) - xl(2,3) 
            Def(2,2) = xl(2,2) - xl(2,3) 
            rhsD(1) = lM%xi(1,g)
            rhsD(2) = lM%xi(2,g)
         ELSE 
            write(*,*)"****** 3D implementation still TODO *****" 
!           Need to check how the basis function are ordered
!           and rewrite Def and rhsD            
         END IF 

!        Find the coordinates in the current element 
         xiCur(:,g) = MATMUL(Def,rhsD)
         xiCur(:,g) = xiCur(:,g) + xl(:,3) 
!         write(*,*)" xiCur coord are ", xiCur(:,g) 

      END DO

      RETURN
      END SUBROUTINE GET_CURQUAD

!--------------------------------------------------------------------
      SUBROUTINE GET_PntRefToCur(lM, fs, xl, xi, xiCur)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      TYPE(fsType), INTENT(IN) :: fs
      REAL(KIND=RKIND), INTENT(IN) :: xl(nsd,lM%eNoN)
      REAL(KIND=RKIND), INTENT(IN) :: xi(nsd, fs%nG)
      REAL(KIND=RKIND), INTENT(OUT) :: xiCur(nsd, fs%nG)

      INTEGER(KIND=IKIND) :: g
      REAL(KIND=RKIND) :: Def(nsd,nsd), rhsD(nsd) 

C       write(*,*)" Inside GET_PntRefToCur "

!     Loop over Gauss quad point
      DO g=1, fs%nG

         IF (nsd .EQ. 2) THEN
            Def(1,1) = xl(1,1) - xl(1,3) 
            Def(1,2) = xl(1,2) - xl(1,3) 
            Def(2,1) = xl(2,1) - xl(2,3) 
            Def(2,2) = xl(2,2) - xl(2,3) 
            rhsD(1) = xi(1,g)
            rhsD(2) = xi(2,g)
         ELSE 
            write(*,*)"****** 3D implementation still TODO *****" 
!           Need to check how the basis function are ordered
!           and rewrite Def and rhsD            
         END IF 

!        Find the coordinates in the current element 
         xiCur(:,g) = MATMUL(Def,rhsD)
         xiCur(:,g) = xiCur(:,g) + xl(:,3) 
!         write(*,*)" xiCur coord are ", xiCur(:,g) 

      END DO

      RETURN
      END SUBROUTINE GET_PntRefToCur
!--------------------------------------------------------------------
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


!####################################################################
!####################################################################
!--------------------------------------------------------------------
!
!     Compute Jac from subElement (in the reference element) 
!     and reference element 
!
!--------------------------------------------------------------------
      PURE SUBROUTINE COMP_SubElmJAC(eNoN, insd, Nxi, x, Jac)
      USE COMMOD, ONLY: nsd
      USE UTILMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, insd
      REAL(KIND=RKIND), INTENT(IN) :: Nxi(insd,eNoN), x(nsd,eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: Jac 

      INTEGER(KIND=IKIND) a
      REAL(KIND=RKIND) xXi(nsd,insd), xiX(insd,nsd)

      xXi = 0._RKIND
      Jac = 0._RKIND
      IF (insd .EQ. 1) THEN
         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
         END DO

         Jac = SQRT(NORM(xXi)) + 1.E+3_RKIND*eps

      ELSE IF (insd .EQ. 2) THEN
         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
         END DO

         Jac = xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1)

      ELSE IF (insd .EQ. 3) THEN
         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + x(:,a)*Nxi(3,a)
         END DO

         Jac = xXi(1,1)*xXi(2,2)*xXi(3,3)
     2       + xXi(1,2)*xXi(2,3)*xXi(3,1)
     3       + xXi(1,3)*xXi(2,1)*xXi(3,2)
     4       - xXi(1,1)*xXi(2,3)*xXi(3,2)
     5       - xXi(1,2)*xXi(2,1)*xXi(3,3)
     6       - xXi(1,3)*xXi(2,2)*xXi(3,1)

      END IF

      RETURN
      END SUBROUTINE COMP_SubElmJAC

!####################################################################
!####################################################################
!--------------------------------------------------------------------
!
!     TODO
!
!--------------------------------------------------------------------

      SUBROUTINE NTS_INIT(lA, lY, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(tDof, tnNo), lY(tDof, tnNo),
     2   lD(tDof, tnNo)

      INTEGER(KIND=IKIND) :: e, a, Ac, nbrSN, is, As, i, j, count, eOpp
      REAL(KIND=RKIND) :: xp(nsd), poly(nsd,msh(1)%eNoN), xs(nsd)
      REAL(KIND=RKIND) :: minXs(nsd), maxXs(nsd)
      REAL(KIND=RKIND) :: minb(nsd), maxb(nsd)
      INTEGER(KIND=IKIND) :: find
      INTEGER(KIND=IKIND) :: maxNbrSNd = 0
      INTEGER(KIND=IKIND), ALLOCATABLE :: FNdFlag(:), cnt(:), 
     2                                                      FElmSNd(:,:)
      INTEGER(KIND=IKIND) :: f, faceID(msh(1)%eNoN, msh(1)%eNoN-1)

C       write(*,*)" Calling NTS_INIT "

!     We will make the assumption that the first mesh is the fluid and the 
!     second the solid, this can be generalized 

!     This is now done in READFILES if Nitsche flag is active 
! !     Define fluid mesh stencil 
!       CALL GETNSTENCIL(msh(1))
!       CALL GETNSTENCIL(msh(2))!
!       write(*,*)" Stencil done "!
! !     Define neighbour structure necessary for ghost penalty stabilization 
!       CALL GETNEIGH(msh(1))

      ALLOCATE( ghostFNd(msh(1)%nNo), FElmSNd(msh(1)%nEl,20), 
     2          FNdFlag(msh(1)%nNo),  mapSNdFElm(msh(2)%nNo), 
     3          cnt(msh(1)%nEl),      intFElmFlag(msh(1)%nEl) )

      intFElmFlag = 0
      FElmSNd = 0
      cnt = 0
      mapSNdFElm = 0
      ghostFNd = 0
      FNdFlag = 0

!     Compute meshes diameter 
!     CALL COMP_DIAM(msh(1),msh(1)%mDiam)
      CALL COMP_DIAM(msh(2),msh(2)%mDiam)
C       write(*,*)" Diam solid computed "

!     Create a bounding box around of the current solid location 
      minb = HUGE(minb)
      maxb = TINY(maxb)

      minXs = HUGE(minXs) 
      maxXs = TINY(maxXs)

!     Find max/min solid coord in each direction
      DO e=1, msh(2)%nEl
         DO a=1, msh(2)%eNoN
            Ac = msh(2)%IEN(a,e)
            xp(:) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)

            DO i=1, nsd 
               minXs(i) = MIN(xp(i),minXs(i))
               maxXs(i) = MAX(xp(i),maxXs(i))
            END DO
         END DO 
      END DO

      DO i=1, nsd
         minb(i) = minXs(i) - msh(2)%mDiam 
         maxb(i) = maxXs(i) + msh(2)%mDiam 
      END DO

C       write(*,*)"S BBOX NTS_INIT x-dir is: ", minb(1), " and ", maxb(1)
C       write(*,*)"S BBOX NTS_INITy-dir is: ", minb(2), " and ", maxb(2)

!     Check that the solid is inside the fluid deformed mesh 
C       DO iM=1, nMsh
C          DO a=1, msh(iM)%nNo
C             Ac = msh(iM)%gN(a)
C             xp(:) = x(:,Ac) + lD(:,Ac)
C             DO i=1, nsd
C                IF (minb(i) .GT. xp(i)) minb(i) = xp(i)
C                IF (maxb(i) .LT. xp(i)) maxb(i) = xp(i)
C             END DO
C          END DO
C       END DO

!     loop over the fluid element to search if any of the solid node is inside 
      DO a = 1, msh(1)%nNo

         find = 0

         Ac = msh(1)%gN(a)

         xp(:) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)

!        Is the fluid node in the solid BBox? 
         IF((xp(1) .LE. maxb(1)) .AND. (xp(1) .GE. minb(1)) 
     2                   .AND. (xp(2) .LE. maxb(2)) 
     3                   .AND. (xp(2) .GE. minb(2))) THEN

            DO e = 1, msh(2)%nEl 

               DO is = 1, msh(2)%eNoN
                  As = msh(2)%IEN(is,e)
                  poly(:,is) = x(:,As) + lD(nsd+2:2*nsd+1,As)
               END DO

!              Is the fluid node inside any solid element?
               CALL IN_POLY(xp,poly, find)

               IF (find .EQ. 1) EXIT
            END DO 
         END IF

!        Fluid element is not intersected by the solid mesh 
         IF (find .EQ. 0) THEN 

!           Set the fluid node as outside node 
            FNdFlag(Ac) = 1

!           Loop on stencil and set to 1 all the fluid node in ghostFNd
            nbrSN = msh(1)%stn%nbrNdStn(Ac)

            DO is = 1, nbrSN
               ghostFNd(msh(1)%stn%ndStn(Ac,is)) = 1
            END DO
         END IF 

      END DO

!---- Filling intFElmFlag: 1 fluid, 2 intersected from Nitsche Bnd, 0 hidden element 
      DO e = 1, msh(1)%nEl 

         count = 0

         DO a = 1, msh(1)%eNoN
            Ac = msh(1)%IEN(a,e)
            IF ( FNdFlag(Ac) .EQ. 1 ) count = count + 1
         END DO

         IF( count .EQ. msh(1)%eNoN) THEN 
!           The fluid element is outside the solid region 
            intFElmFlag(e) = 1
         ELSE IF ( count .EQ. 0) THEN 
!           The fluid element is fully inside the solid region 
            intFElmFlag(e) = 0 
         ELSE 
!           The fluid element is intersected from the Nitsche boundary  
            intFElmFlag(e) = 2
         END IF

      END DO

!---- Filling mapSNdFElm, find the fluid element that contains the corresponding solid node
      DO a = 1, msh(2)%nNo

         find = 0

         As = msh(2)%gN(a)
         xp(:) = x(:,As) + lD(nsd+2:2*nsd+1,As)

!        Loop over fluid element 
         DO e = 1, msh(1)%nEl 

!           Is the fluid node in the solid BBox? 
            IF((xp(1) .LE. maxb(1)) .AND. (xp(1) .GE. minb(1)) 
     2                   .AND. (xp(2) .LE. maxb(2)) 
     3                   .AND. (xp(2) .GE. minb(2))) THEN

               DO is = 1, msh(1)%eNoN
                  Ac = msh(1)%IEN(is,e)
                  poly(:,is) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
               END DO

   !           Is the fluid node inside any solid element?
               CALL IN_POLY(xp,poly, find)

               IF (find .EQ. 1) THEN 
                  mapSNdFElm(a) = e

!                 Check that we have labeled it correctly
                  IF( intFElmFlag(e) .EQ. 1) THEN 
                     intFElmFlag(e) = 2
                     j = e-1 
!                     write(*,*)" Modifying Fluid Elm ", j, " Flag to 2 "
                  END IF

                  EXIT
               END IF
            END IF
         END DO 
      END DO

!---- Filling FElmSNd, loop over mapSNdFElm and add the analogous contribution 
      DO a = 1, msh(2)%nNo

C          Ac = msh(2)%gN(a)

         e = mapSNdFElm(a)
         cnt(e) = cnt(e) + 1

         FElmSNd( e, cnt(e)) = a

      END DO

C       write(*,*)" cnt = ", cnt

      maxNbrSNd = MAXVAL(cnt)

!---- Filling mapFElmSNd from FElmSNd 
      ALLOCATE(mapFElmSNd(msh(1)%nEl ,maxNbrSNd))
      mapFElmSNd = 0
      mapFElmSNd(:,:) = FElmSNd(:,1:maxNbrSNd)


!---- Filling dofToJumpAss
      ALLOCATE( dofToJumpAss(msh(1)%nNo) ) 
      dofToJumpAss = 0

      IF ( nsd .EQ. 2 ) THEN
         faceID(1,:) = (/2, 3/)
         faceID(2,:) = (/3, 1/)
         faceID(3,:) = (/1, 2/)

      ELSE IF ( nsd .EQ. 3 ) THEN
         faceID(1,:) = (/2, 3, 4/)
         faceID(2,:) = (/4, 3, 1/)
         faceID(3,:) = (/1, 2, 4/)
         faceID(4,:) = (/1, 3, 2/)

      ELSE 
         err = "Dimension not supported in GETNEIGH"
      END IF

      DO e=1, msh(1)%nEl 
!        Define dof to jump dering assembly 
         DO f = 1, msh(1)%eNoN
            eOpp = msh(1)%neigh(e,f)
            IF( eOpp .EQ. -1 ) CYCLE

            IF( (intFElmFlag(e) .EQ. 2) .AND. 
     2                 (intFElmFlag(eOpp) .EQ. 0) ) THEN

C                write(*,*)"intFElmFlag elem  ", e , " is 2"

               DO i=1, msh(1)%eNoN-1
                  Ac = msh(1)%IEN(faceID(f,i),e)
C                   write(*,*)" Ac = ", Ac
                  dofToJumpAss(Ac) = 1
               END DO
            END IF 
         END DO
      END DO

C       DO e = 1, msh(1)%nNo
C C          write(*,*)" dofToJumpAss(",e,") = ", dofToJumpAss(e)
C          IF( dofToJumpAss(e) .EQ. 1) write(*,*)" jumping dof ", e
C       END DO




!------------------------------
!     Plots useful for debug 
!     Plot ghostFNd and FNdFlag
C       write(*,*)" Beginning printing stuff "
C       DO i = 1, msh(1)%nNo
C          j = i-1
C          IF( FNdFlag(i) .EQ. 0 ) THEN 
C             IF( ghostFNd(i) .EQ. 0 ) THEN 
C                write(*,*)" ghostFNd and FNdFlag 0 for fluid node ", j
C             ELSE 
C             write(*,*)" Only FNdFlag 0 for fluid node ", j
C             END IF
C          END IF
C       END DO
C !     Plot intFElmFlag
C       DO e = 1, msh(1)%nEl 
C          j = e-1
C          IF ( intFElmFlag(e) .EQ. 0) THEN 
C             write(*,*)" Fluid element ", j , " is behind solid"  
C          ELSE IF ( intFElmFlag(e) .EQ. 2) THEN 
C             write(*,*)" Fluid element ", j , " intersected Nit Bnd "  
C          END IF
C       END DO
C !     Plot mapSNdFElm 
C       DO a = 1, msh(2)%nNo
C          j = a-1
C          i = mapSNdFElm(a) - 1
C          write(*,*)" Solid ", j, " is fluid elm ", i
C       END DO
C !     Plot mapFElmSNd, loop over mapSNdFElm and add the analogous contribution 
C       write(*,*)" maxNbrSNd = ", maxNbrSNd
C       mapFElmSNd = mapFElmSNd - 1
C       DO a = 1, msh(1)%nEl
C          IF ( intFElmFlag(a) .NE. 1) THEN 
C             j = a -1 
C             write(*,*)" Fluid Elm ", j, ": internal Solid nodes "
C             write(*,*) " " , mapFElmSNd(a,1:cnt(a))
C          END IF
C       END DO
C       mapFElmSNd = mapFElmSNd + 1


      DEALLOCATE(FElmSNd, FNdFlag)

      RETURN
      END SUBROUTINE NTS_INIT
!####################################################################
!--------------------------------------------------------------------
!
!     TODO
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTR_mapFElmSElm(lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

      INTEGER(KIND=IKIND) :: iFa, eNoNF, idFaceElm, a, Ac, nbrNdF, 
     2                       cntElm, iFNit, e, eTypeF, nG, nGF, idFlElm, 
     3                       nbrSElm, maxNbrSElm, maxNbr, i, findInF, 
     4                       idF, iq, find, eOpp, f, Ac1, Ac2, inside, 
     5                       int, onFace, cntInt, onVertex
      INTEGER(KIND=IKIND) :: faceID(msh(1)%eNoN,msh(1)%eNoN-1)
      INTEGER(KIND=IKIND), ALLOCATABLE :: FElmSElm(:,:,:), cnt(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xFace(:,:), ds(:,:), QPnt(:,:)
      REAL(KIND=RKIND) :: xq(nsd), poly(nsd,msh(1)%eNoN), crd1(nsd), 
     2             crdFace(nsd,msh(1)%eNoN-1), segm(nsd,msh(1)%eNoN-1), 
     3             Pint(nsd), crd2(nsd), nrm, tan(nsd), tanN(nsd)
      LOGICAL :: isNeigh

      INTERFACE
         SUBROUTINE GET_SElmQuadPoint(nbrSElm, eTypeF, eNoNF, nG, nGF, 
     2                         crdFace, quadPnt, SElm, qPRef )
         USE COMMOD
         IMPLICIT NONE
         INTEGER(KIND=IKIND), INTENT(IN) :: nbrSElm, eTypeF,nG,nGF,eNoNF
         REAL(KIND=RKIND), INTENT(IN) :: crdFace(nsd, eNoNF)
         REAL(KIND=RKIND), INTENT(OUT) :: quadPnt(nsd, nGF)
         REAL(KIND=RKIND),INTENT(OUT),OPTIONAL :: SElm(nsd-1,
     2                                                    eNoNF,nbrSElm)
         REAL(KIND=RKIND),INTENT(OUT),OPTIONAL :: qPRef(nsd-1, nGF)
         END SUBROUTINE GET_SElmQuadPoint
      END INTERFACE



C       write(*,*) " Inside CONSTR_mapFElmSElm "

!---- Filling FElmSElm, loop over mapSNdFElm and add the analogous contribution 
      ALLOCATE(FElmSElm(msh(1)%nEl, msh(2)%nNtsFa, 10), 
     2                                  cnt(msh(1)%nEl, msh(2)%nNtsFa) )

      nbrNdF = msh(1)%nNo
      cnt = 0
      iFNit = 0
      FElmSElm = 0

      faceID(1,:) = (/2, 3/)
      faceID(2,:) = (/3, 1/)
      faceID(3,:) = (/1, 2/)

!---  Initial way of fulling FElmSElm using quad point, failing for 
!     intersection purposes 

C !     Loop over solid faces 
C       DO iFa = 1, msh(2)%nFa 
C C       DO iFa = 1, 1

C !        We cycle if this is not a Nitsche-face 
C          IF( .NOT.msh(2)%fa(iFa)%isNts ) CYCLE

C !        Update counter Nitsche's face          
C          iFNit = iFNit + 1

C C          write(*,*)" Lookinfg at Nitsche solid face ", iFa

C          eNoNF  = msh(2)%fa(iFa)%eNoN
C          eTypeF = msh(2)%fa(iFa)%eType
C          nG     = msh(2)%fa(iFa)%nG

C !        Define number of face sub-element/tot quad point in the reference face element 
C          IF ( nsd .EQ. 2) nbrSElm = 2**(nbrFntCllSD-1)
C          IF ( nsd .EQ. 3) nbrSElm = 4**(nbrFntCllSD-1)
C          nGF = nG*nbrSElm

C          ALLOCATE(QPnt(nsd,nGF))

C C          write(*,*)"Solid face ",iFa,"nbr of elem = ",msh(2)%fa(iFa)%nEl

C !        Loop over face element 
C          DO idFaceElm = 1, msh(2)%fa(iFa)%nEl

C C             write(*,*)" local idFaceElm = ", idFaceElm

C !           Are all element's nodes inside the same fluid elem? 
C             cntElm = 0
C             Ac1 = msh(2)%fa(iFa)%IEN(1,idFaceElm) - nbrNdF

C             DO a = 1, eNoNF
C                Ac = msh(2)%fa(iFa)%IEN(a,idFaceElm)
C C                Ac = msh(2)%fa(iFa)%gN(Ac) ! check this for parallel 
C                Ac = Ac - nbrNdF
C C                write(*,*)" Ac = ", Ac, " inside fluid ", mapSNdFElm(Ac) 

C                IF ( mapSNdFElm(Ac) .EQ. mapSNdFElm(Ac1) ) THEN 
C                   cntElm = cntElm + 1
C                END IF

C             END DO

C !---------- Check if the solid element nodes are in the same fluid elm 
C             IF (cntElm .EQ. eNoNF) THEN 
               
C !-----         The solid face is in only one fluid elm
C C                 write(*,*)" The solid face is in only one fluid elm "

C !              Set fluid element in which the solid face is                
C                idFlElm = mapSNdFElm(Ac)
C C                write(*,*)" idFlElm = ", idFlElm

C !              Update counter for nmr of solids in fluid elm                
C                cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) + 1 
C                e = cnt(idFlElm,iFNit) 

C !              Store local (wrt the face) id face into idFlElm map  
C !              Check if we have already add this solid for this fluid   
C                DO a = 1, e           
C                   IF(FElmSElm(idFlElm,iFNit, a) .EQ. idFaceElm) THEN
C                      cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) - 1  
C                      EXIT
C                   END IF
C                   IF(a .EQ. e) FElmSElm(idFlElm,iFNit, e) = idFaceElm
C                END DO    

C             ELSE 
C !-----         The solid extremes are in different flu elem
C C                write(*,*)" The solid extremes are in different flu elem"

C !-----         Store id fluid elem for each solid extreme
C                DO i = 1, eNoNF
C                   Ac = msh(2)%fa(iFa)%IEN(i,idFaceElm)
C C                 Ac = msh(2)%fa(iFa)%gN(Ac) ! check this for parallel 
C                   Ac = Ac - nbrNdF

C !                 Set fluid element in which the solid face is                
C                   idFlElm = mapSNdFElm(Ac)

C !                 Update counter for nmr of solids in fluid elm                
C                   cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) + 1 
C                   e = cnt(idFlElm,iFNit) 

C !                 Store local (wrt the face) id face into idFlElm map  
C !                 Check if we have already add this solid for this fluid   
C                   DO a = 1, e           
C                      IF(FElmSElm(idFlElm,iFNit, a) .EQ. idFaceElm) THEN
C                         cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) - 1  
C                         EXIT
C                      END IF
C                      IF(a .EQ. e) FElmSElm(idFlElm,iFNit, e) = idFaceElm
C                   END DO 
C                END DO


C !-----         check the quad point 
C                ALLOCATE(xFace(nsd,eNoNF), ds(tDof,eNoNF))

C                DO a = 1, eNoNF
C                   Ac = msh(2)%fa(iFa)%IEN(a,idFaceElm)
C                   ds(:,a) = lD(:,Ac)
C                   xFace(:,a) = x(:,Ac) + ds(nsd+2:2*nsd+1,a)
C                END DO

C C                write(*,*)" xFace solid is ", xFace

C                CALL GET_SElmQuadPoint(nbrSElm, eTypeF, eNoNF, nG, nGF, 
C      2                                            xFace, QPnt)

C !--            Given the quadrature points in the current solid face element 
C !              we have to find the fluid element that contains the quad point                

C                DO iq = 1, nGF 
C                   xq = QPnt(:,iq)

C C                   write(*,*)" quad point: ", iq
C C                   write(*,*) xq

C                   find = 0

C                   DO a = 1, eNoNF

C                      IF( find .EQ. 1 ) EXIT 

C                      Ac = msh(2)%fa(iFa)%IEN(a,idFaceElm)
C C                    Ac = msh(2)%fa(iFa)%gN(Ac) ! check this for parallel 
C                      Ac = Ac - nbrNdF

C                      idF = mapSNdFElm(Ac) 

C C                      write(*,*)" looking into fluid elm ", idF

C                      find = 0

C                      CALL SrcSTENCIL(msh(1),lD, xq, idF, idFlElm, find)
                   
C !--                  If we have found the quadrature point in the fluid mesh 
C !                    we update the map                      
C                      IF( find .EQ. 1 ) THEN 
C !                       Update counter for nmr of solids in fluid elm                
C                         cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) + 1 
C                         e = cnt(idFlElm,iFNit) 

C !                       Store local (wrt the face) id face into idFlElm map  
C !                       Check if we have already add this solid for this fluid   
C                         DO i = 1, e           
C                            IF(FElmSElm(idFlElm,iFNit, i) 
C      2                                         .EQ.idFaceElm) THEN

C                               cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit)-1  
C                               EXIT
C                            END IF
C                            IF(i .EQ. e) THEN 
C                               FElmSElm(idFlElm,iFNit, e) = idFaceElm
C                            END IF
C                         END DO

C                      END IF

C                   END DO 

C !--               If at the end of the loop we didn't find it, loop 
C !                 over the whole fluid mesh 
C                   IF( find .EQ. 0 ) THEN 
C                      write(*,*)"!!!! @@@ NEED TO LOOP OVER MESH "
C                      CALL SrcMESH(msh(1),lD, xq, idFlElm, find)

C                      IF( find .EQ. 1 ) THEN 
C !                       Update counter for nmr of solids in fluid elm                
C                         cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) + 1 
C                         e = cnt(idFlElm,iFNit) 

C !                       Store local (wrt the face) id face into idFlElm map  
C !                       Check if we have already add this solid for this fluid   
C                         DO i = 1, e           
C                            IF(FElmSElm(idFlElm,iFNit, i) 
C      2                                         .EQ.idFaceElm) THEN

C                               cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit)-1  
C                               EXIT
C                            END IF
C                            IF(i .EQ. e) THEN 
C                               FElmSElm(idFlElm,iFNit, e) = idFaceElm
C                            END IF
C                         END DO

C                      ELSE
C                         write(*,*)" quad solid node not in fluid mesh "
C                         !err = " quad solid node not in fluid mesh "
C                      END IF
C                   END IF


C                END DO

C                DEALLOCATE(xFace, ds)


C             END IF 

C          END DO 

C          DEALLOCATE(QPnt)

C C          write(*,*)""
C C          write(*,*)""
C C          write(*,*)""

C       END DO 


!---  New way using intersection 


!     Loop over solid faces 
      DO iFa = 1, msh(2)%nFa 
C       DO iFa = 1, 1

!        We cycle if this is not a Nitsche-face 
         IF( .NOT.msh(2)%fa(iFa)%isNts ) CYCLE

!        Update counter Nitsche's face          
         iFNit = iFNit + 1

C          write(*,*)" Lookinfg at Nitsche solid face ", iFa

         eNoNF  = msh(2)%fa(iFa)%eNoN
         eTypeF = msh(2)%fa(iFa)%eType
         nG     = msh(2)%fa(iFa)%nG

!        Loop over face element 
         DO idFaceElm = 1, msh(2)%fa(iFa)%nEl

C             write(*,*)""
C             write(*,*)""
C             write(*,*)" idFaceElm solid = ", idFaceElm

            Ac1 = mapSNdFElm(msh(2)%fa(iFa)%IEN(1,idFaceElm) - nbrNdF)
            Ac2 = mapSNdFElm(msh(2)%fa(iFa)%IEN(2,idFaceElm) - nbrNdF)

1              CONTINUE

            write(*,*)" solid extreme are in flui elm ", Ac1, Ac2

            CALL IS_NEIGH( msh(1), Ac1, Ac2, isNeigh)

            IF (Ac1 .EQ. Ac2) THEN                
!              The solid face is in only one fluid elm

!              Set fluid element in which the solid face is                
               idFlElm = Ac1
C                write(*,*)" adding CASE 1 idFlElm = ", idFlElm

!              Update counter for nmr of solids in fluid elm                
               cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) + 1 
               e = cnt(idFlElm,iFNit) 

!              Store local (wrt the face) id face into idFlElm map  
!              Check if we have already add this solid for this fluid   
               DO a = 1, e           
                  IF(FElmSElm(idFlElm,iFNit, a) .EQ. idFaceElm) THEN
                     cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) - 1  
                     EXIT
                  END IF
                  IF(a .EQ. e) FElmSElm(idFlElm,iFNit, e) = idFaceElm
               END DO  

            ELSE IF (isNeigh) THEN 

C                write(*,*)" adding CASE 2 idFlElm = ", Ac1
C                write(*,*)" adding CASE 2 idFlElm = ", Ac2

!              Set fluid element in which the solid face is                
               idFlElm = Ac1
!              Update counter for nmr of solids in fluid elm                
               cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) + 1 
               e = cnt(idFlElm,iFNit) 
!              Store local (wrt the face) id face into idFlElm map  
!              Check if we have already add this solid for this fluid   
               DO a = 1, e           
                  IF(FElmSElm(idFlElm,iFNit, a) .EQ. idFaceElm) THEN
                     cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) - 1  
                     EXIT
                  END IF
                  IF(a .EQ. e) FElmSElm(idFlElm,iFNit, e) = idFaceElm
               END DO  


!              Set fluid element in which the solid face is                
               idFlElm = Ac2
!              Update counter for nmr of solids in fluid elm                
               cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) + 1 
               e = cnt(idFlElm,iFNit) 
!              Store local (wrt the face) id face into idFlElm map  
!              Check if we have already add this solid for this fluid   
               DO a = 1, e           
                  IF(FElmSElm(idFlElm,iFNit, a) .EQ. idFaceElm) THEN
                     cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) - 1  
                     EXIT
                  END IF
                  IF(a .EQ. e) FElmSElm(idFlElm,iFNit, e) = idFaceElm
               END DO 

            ELSE
 
C                write(*,*)" adding CASE 3 aa idFlElm = ", Ac2
C !              Set fluid element in which the solid face is                
C                idFlElm = Ac2
C !              Update counter for nmr of solids in fluid elm                
C                cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) + 1 
C                e = cnt(idFlElm,iFNit) 
C !              Store local (wrt the face) id face into idFlElm map  
C !              Check if we have already add this solid for this fluid   
C                DO a = 1, e           
C                   IF(FElmSElm(idFlElm,iFNit, a) .EQ. idFaceElm) THEN
C                      cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) - 1  
C                      EXIT
C                   END IF
C                   IF(a .EQ. e) FElmSElm(idFlElm,iFNit, e) = idFaceElm
C                END DO 

               inside = 0
               cntInt = 0

               Ac = msh(2)%fa(iFa)%IEN(1,idFaceElm)
               crd1(:) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
C                write(*,*)" crd1 ", crd1

               Ac = msh(2)%fa(iFa)%IEN(2,idFaceElm)
               crd2(:) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
C                write(*,*)" crd2 ", crd2
            
               DO WHILE ( inside .EQ. 0)
                  write(*,*)" inside while idFaceElm solid ", idFaceElm
!                 Is Ac2 inside? yes, we are done and we insert it. 
!                 No we insert Ac1 and we search for an interscetion on the edges != from crd1

!                 Define crd first point and poly
C                   write(*,*)" Ac1 = ", Ac1
                  DO a = 1, msh(1)%eNoN
                     Ac = msh(1)%IEN(a,Ac1)
C                      write(*,*)" Ac ", Ac
                     poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
                  END DO
C                   write(*,*)" poly = ", poly 

                  DO i = 1, msh(1)%eNoN
                     crdFace(:,1) = poly(:,faceID(i,1))
                     crdFace(:,2) = poly(:,faceID(i,2))

                     segm(:,1) = crd1
                     segm(:,2) = crd2

C                      write(*,*)" searching inter crdFace ", i ," = ",
C      2                                              crdFace
C                      write(*,*)" searching inter segm ", segm

                     CALL SegSegInt(crdFace, segm, int, Pint)
C                      write(*,*)" int ", int, " and Pint = ", Pint
                     
                     IF( int .EQ. 0 ) CYCLE

                     cntInt = 1

C                    write(*,*)" adding CASE 3 aa idFlElm = ", Ac1
!                    Set fluid element in which the solid face is                
                     idFlElm = Ac1
!                    Update counter for nmr of solids in fluid elm                
                     cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) + 1 
                     e = cnt(idFlElm,iFNit) 
!                    Store local (wrt the face) id face into idFlElm map  
!                    Check if we have already add this solid for this fluid   
                     DO a = 1, e           
                        IF(FElmSElm(idFlElm,iFNit, a).EQ.idFaceElm) THEN
                           cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) - 1  
                           EXIT
                        END IF
                        IF(a .EQ. e) FElmSElm(idFlElm,iFNit,e)=idFaceElm
                     END DO 


!                    Once we have the intersection we update crd1 and Ac1 with the neigh to 
!                    the face in which we have found the intersection
                     
                     crd1 = Pint
C                      write(*,*)" new intersection is ", crd1
                     CALL IN_POLY_FV(Pint, poly, int, onFace, onVertex)
C                      write(*,*)" int ", int
C                      write(*,*)" on face ", onFace
                     idFlElm = Ac1
                     Ac1 = msh(1)%neigh(idFlElm,onFace)
C                      write(*,*)" and with neigh ", Ac1

!                    Set fluid element in which the solid face is                
                     idFlElm = Ac1
!                    Update counter for nmr of solids in fluid elm                
                     cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) + 1 
                     e = cnt(idFlElm,iFNit) 
!                    Store local (wrt the face) id face into idFlElm map  
!                    Check if we have already add this solid for this fluid   
                     DO a = 1, e           
                        IF(FElmSElm(idFlElm,iFNit, a).EQ.idFaceElm) THEN
                           cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) - 1  
                           EXIT
                        END IF
                        IF(a .EQ. e) THEN 
                           FElmSElm(idFlElm,iFNit, e) = idFaceElm
                           write(*,*)" adding CASE 3 b idFlElm = ", Ac1
                        END IF
                     END DO 



                     EXIT
                  END DO

                  IF( cntInt .EQ. 0) THEN 
                     write(*,*)" we are starting from a wrong element "

!                    Find out where is the node in the element 
                     CALL IN_POLY_FV(crd1, poly, int, onFace, onVertex)

                     tan = crd2 - crd1 
                     nrm = SQRT(tan(1)**2 + tan(2)**2)
                     tanN(1) = tan(1)/nrm
                     tanN(2) = tan(2)/nrm
                     Pint = crd1 + tanN * msh(1)%mDiam * 1.E-1_RKIND

!                    If on vertex: move an eps from the point and seach 
!                    in the stencil 
                     IF( onVertex .GT. 0) THEN 
                        Ac = msh(1)%IEN(onVertex,Ac1)
                        CALL SrcSTENCIL_NODE(msh(1), lD, Pint, Ac, 
     2                                                       Ac1, int)
                        write(*,*)" node new starting elm ", Ac1
                        GOTO 1
                     END IF 


!                    If on an edge, look the neigh wrt that face 
                     IF( onFace .GT. 0) THEN 
                        Ac1 = msh(1)%neigh(Ac1,onFace)

                        write(*,*)" face new starting elm ", Ac1

                        GOTO 1
                     END IF                      

                  END IF


                  CALL IN_POLY(crd2, poly, inside)
C                   write(*,*)" is crd 2 inside ", inside 
                  
               END DO

C                write(*,*)" adding CASE 3 c idFlElm = ", Ac2

!              Add neigh for Ac2  
!              Set fluid element in which the solid face is                
               idFlElm = Ac2
!              Update counter for nmr of solids in fluid elm                
               cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) + 1 
               e = cnt(idFlElm,iFNit) 
!              Store local (wrt the face) id face into idFlElm map  
!              Check if we have already add this solid for this fluid   
               DO a = 1, e           
                  IF(FElmSElm(idFlElm,iFNit, a) .EQ. idFaceElm) THEN
                     cnt(idFlElm,iFNit) = cnt(idFlElm,iFNit) - 1  
                     EXIT
                  END IF
                  IF(a .EQ. e) FElmSElm(idFlElm,iFNit, e) = idFaceElm
               END DO              

            END IF 
         END DO
      END DO



!---      

      maxNbrSElm = MAXVAL(cnt(:,1))
      DO i=1, iFNit
         maxNbr = MAXVAL(cnt(:,i))
         IF(maxNbr .GT. maxNbrSElm) maxNbrSElm = maxNbr
      END DO
C       write(*,*)" maxNbrSElm =  ", maxNbrSElm


!     Copy FElmSElm into  mapFElmSElm
!---- Filling mapFElmSElm  
      
      IF (ALLOCATED(mapFElmSElm)) DEALLOCATE(mapFElmSElm)
      ALLOCATE( mapFElmSElm(msh(1)%nEl, msh(2)%nNtsFa, maxNbrSElm) )
      mapFElmSElm = 0
      mapFElmSElm(:,:,:) = FElmSElm(:,:,1:maxNbrSElm)

!-------------------------
!     Update intFElmFlag, if we find out that there are no solid quad point
!     into the fluid element  
      DO e = 1, msh(1)%nEl 

         IF( intFElmFlag(e) .EQ. 2 ) THEN 
            intFElmFlag(e) = 1

            DO a = 1, msh(2)%nNtsFa  
               IF( mapFElmSElm(e,a,1) .NE. 0 ) intFElmFlag(e) = 2
            END DO

C             IF( intFElmFlag(e) .NE. 2 ) THEN 
C                write(*,*)"  intFElmFlag(e) changed to ",  intFElmFlag(e)
C             END IF

         END IF
      END DO

      DO e = 1, msh(1)%nEl 
         IF( intFElmFlag(e) .EQ. 1 ) THEN 

            DO f = 1, msh(1)%eNoN
               eOpp = msh(1)%neigh(e,f)
               
               IF( eOpp .EQ. -1 ) CYCLE
               
               IF(intFElmFlag(eOpp) .EQ. 0) THEN 
                  intFElmFlag(e) = 0
                  EXIT
               END IF
            END DO
         END IF
      END DO


!-------------------------
!     Useful printing part 
C       write(*,*)" PRINTING PART mapFElmSElm "
      
      iFNit =0

      DO e = 1, msh(2)%nFa
         IF( .NOT.msh(2)%fa(e)%isNts ) CYCLE

         iFNit = iFNit + 1

         DO a = 1, msh(1)%nEl
            IF( mapFElmSElm(a,iFNit,1) .EQ. 0 ) CYCLE
            write(*,*)" solid face nische ", e, " element fluid: ", a 

            DO i=1, maxNbrSElm

               IF( mapFElmSElm(a,iFNit,i) .GT. msh(2)%fa(e)%nEl) THEN 
                  write(*,*)" error here!!!!! %%%%$$$$%@^ "
               END IF
               write(*,*) mapFElmSElm(a,iFNit,i) 
            END DO
         END DO
      END DO


C       write(*,*)" END PRINTING PART mapFElmSElm "


      DEALLOCATE(FElmSElm, cnt)


      RETURN
      END SUBROUTINE CONSTR_mapFElmSElm

!####################################################################
!     Search a point xq inside the element stencil of node idNodeStc
      SUBROUTINE SrcSTENCIL(lM, lD, xq, idElmStc, findIn, find)
      USE COMMOD
      IMPLICIT NONE
      
      TYPE(mshType),  INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(IN) :: xq(nsd)
      INTEGER(KIND=IKIND), INTENT(IN) :: idElmStc 
     
      INTEGER(KIND=IKIND), INTENT(OUT) :: findIn, find
      
      INTEGER(KIND=IKIND) :: nbrElmStc, is, e, a, Ac, idNodeStc, in
      REAL(KIND=RKIND), ALLOCATABLE :: poly(:,:)

      ALLOCATE(poly(nsd,lM%eNoN))

      find = 0
      findIn = 0

      DO in = 1, lM%eNoN

         IF(find .EQ. 1) EXIT

         idNodeStc = lM%IEN(in,idElmStc) 
C          write(*,*)" Node for stencil is ", idNodeStc

         nbrElmStc = lM%stn%nbrNdStn(idNodeStc)-1

C       write(*,*)" beginning loop "
         DO is = 1, nbrElmStc

            e = lM%stn%elmStn(idNodeStc,is)

            IF( e .EQ. 0 ) EXIT
            
C          write(*,*)" inside loop stencil elem, fliud el id = ", e

            DO a = 1, lM%eNoN
               Ac = lM%IEN(a,e)
C             poly(:,a) = x(:,Ac) + lD(:,Ac)
C                write(*,*)" Ac = ", Ac
               poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
            END DO 

C             write(*,*)" poly ", poly
C          write(*,*)" xq ", xq

            CALL IN_POLY(xq,poly, find)

            IF(find .EQ. 1) THEN 

C             write(*,*)" inside loop stencil elem, fliud el id = ", e
               findIn = e

C                write(*,*)" find it in fluid element ", findIn
               EXIT
            END IF
         END DO
            

      END DO 

C       write(*,*)" "

      DEALLOCATE(poly)

      RETURN

      END SUBROUTINE SrcSTENCIL
    
!--------
!     Search a point xq inside the element stencil of node idNodeStc
      SUBROUTINE SrcSTENCIL_NODE(lM, lD, xq, idNodeStc, findIn, find)
      USE COMMOD
      IMPLICIT NONE
      
      TYPE(mshType),  INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(IN) :: xq(nsd)
      INTEGER(KIND=IKIND), INTENT(IN) :: idNodeStc 
     
      INTEGER(KIND=IKIND), INTENT(OUT) :: findIn, find
      
      INTEGER(KIND=IKIND) :: nbrElmStc, is, e, a, Ac, in
      REAL(KIND=RKIND), ALLOCATABLE :: poly(:,:)

      ALLOCATE(poly(nsd,lM%eNoN))

      find = 0
      findIn = 0

      nbrElmStc = lM%stn%nbrNdStn(idNodeStc)-1
      write(*,*)" Stencil is ", lM%stn%elmStn(idNodeStc,:)

C     write(*,*)" beginning loop "
      DO is = 1, nbrElmStc

         e = lM%stn%elmStn(idNodeStc,is)

         IF( e .EQ. 0 ) EXIT
         
         DO a = 1, lM%eNoN
            Ac = lM%IEN(a,e)
            poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
         END DO 


         CALL IN_POLY(xq,poly, find)

         IF(find .EQ. 1) THEN 

            findIn = e

            EXIT
         END IF
      END DO
            



C       write(*,*)" "

      DEALLOCATE(poly)

      RETURN

      END SUBROUTINE SrcSTENCIL_NODE
!--------------------------------------------------------------------      
!     Search a point xq inside the element stencil of node idNodeStc
      SUBROUTINE SrcMESH(lM, lD, xq, findIn, find)
      USE COMMOD
      IMPLICIT NONE
      
      TYPE(mshType),  INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(IN) :: xq(nsd)
     
      INTEGER(KIND=IKIND), INTENT(OUT) :: findIn, find
      
      REAL(KIND=RKIND), ALLOCATABLE :: poly(:,:)
      INTEGER(KIND=IKIND) :: e, a, Ac

      ALLOCATE(poly(nsd,lM%eNoN))

      find = 0
      findIn = 0

      DO e = 1, lM%nEl

         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            poly(:,a) = x(:,Ac) + lD(nsd+2:2*nsd+1,Ac)
         END DO 

         CALL IN_POLY(xq,poly, find)

         IF(find .EQ. 1) THEN 
            findIn = e
            EXIT
         END IF
      END DO

      DEALLOCATE(poly)
      
      RETURN

      END SUBROUTINE SrcMESH
!####################################################################
!####################################################################
!--------------------------------------------------------------------
!
!     TODO
!
!--------------------------------------------------------------------
      SUBROUTINE NTS_UPDATE(lA, lY, lD)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(tDof, tnNo), lY(tDof, tnNo),
     2   lD(tDof, tnNo)

!     TODO


      RETURN
      END SUBROUTINE NTS_UPDATE
!####################################################################
!####################################################################
!--------------------------------------------------------------------
!
!     Creating the data structure and assembling LHS sparse matrix
!     in case of unfitted fsi with Nitsche type mortaring 
!
!--------------------------------------------------------------------
      SUBROUTINE NTS_LHSA(nnz)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(OUT) :: nnz

      LOGICAL :: flag, isInsie
      INTEGER(KIND=IKIND) a, b, e, i, j, rowN, colN, iM, iFa, masN,
     2   mnnzeic, jM, tmpr, tmpc, aopp, eOpp

      INTEGER(KIND=IKIND), ALLOCATABLE :: uInd(:,:)

      INTEGER(KIND=IKIND) :: Ac
      REAL(KIND=RKIND) :: xp(nsd)
      REAL(KIND=RKIND) :: minXs(nsd), maxXs(nsd)
      REAL(KIND=RKIND) :: minb(nsd), maxb(nsd)

      ALLOCATE(idMap(tnNo))

!     Compute B-Box around solid with fluid mesh size
      CALL COMP_DIAM(msh(1),msh(1)%mDiam)
C       write(*,*)" fluid mesh diam ", msh(1)%mDiam

!     Create a bounding box around of the current solid location 
      minb = HUGE(minb)
      maxb = TINY(maxb)

      minXs = HUGE(minXs) 
      maxXs = TINY(maxXs)

!     Find max/min solid coord in each direction
      DO e=1, msh(2)%nEl
         DO a=1, msh(2)%eNoN
            Ac = msh(2)%IEN(a,e)
            xp(:) = x(:,Ac) !+ lD(nsd+2:2*nsd+1,Ac)

            DO i=1, nsd 
               minXs(i) = MIN(xp(i),minXs(i))
               maxXs(i) = MAX(xp(i),maxXs(i))
            END DO
         END DO 
      END DO

      DO i=1, nsd
         minb(i) = minXs(i) - 5._RKIND*msh(1)%mDiam 
         maxb(i) = maxXs(i) + 5._RKIND*msh(1)%mDiam 
      END DO

      write(*,*)"The solid BBOX x-dir is: ", minb(1), " and ", maxb(1)
      write(*,*)"The solid BBOX y-dir is: ", minb(2), " and ", maxb(2)


!     Select fluid-fluid, solid-solid and fluid-solid node to connect            
      DO a=1, tnNo
         idMap(a) = a
      END DO

      mnnzeic = 20*MAXVAL(msh%eNoN)

!     First fill uInd array depending on mesh connectivity as is.
!     TODO also add the stencil nodes, we need it for ghost penalty stab

      ALLOCATE (uInd(mnnzeic,tnNo))
      uInd = 0
      DO iM=1, nMsh
!        Treat shell with triangular elements separately
         IF (shlEq .AND. msh(iM)%eType.EQ.eType_TRI3) CYCLE
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               rowN = msh(iM)%IEN(a,e)
               DO b=1, msh(iM)%eNoN
                  colN = msh(iM)%IEN(b,e)
                  CALL ADDCOL(rowN, colN)
               END DO
            END DO
         END DO
      END DO

!     Adding opposite Neigh nodes if a fluid element is inside the bigger 
!     solid bounding-box)
      write(*,*)" Adding opposite Neigh nodes "
      DO iM=1, nMsh
!        Select fluid mesh          
         IF( msh(iM)%nNtsFa .EQ. 0) THEN 
            DO e=1, msh(iM)%nEl

!              Check if the element is inside the BBox
               isInsie = .FALSE.

               DO a=1, msh(iM)%eNoN
                  rowN = msh(iM)%IEN(a,e)

!                 Check that the fluid node is in the BBox, if yes, add connection 
                  xp(:) = x(:,rowN) !+ lD(:,Ac)

!                 Is the fluid node in the solid BBox? 
                  IF((xp(1) .LE. maxb(1)) .AND. (xp(1) .GE. minb(1)) 
     2                          .AND. (xp(2) .LE. maxb(2)) 
     3                          .AND. (xp(2) .GE. minb(2))) THEN

                     isInsie = .TRUE.
                  END IF
               END DO

C                write(*,*)" Using CIP "
C                write(*,*)" isInsie = .TRUE. "
               isInsie = .TRUE.

               IF( isInsie ) THEN
C                   write(*,*)" Element ", e , " inside bbox"
                  DO a=1, msh(iM)%eNoN
                     rowN = msh(iM)%IEN(a,e)
                     eOpp = msh(iM)%neigh(e,a)

                     IF( eOpp .EQ. -1 ) CYCLE 

                     DO aopp = 1, msh(iM)%eNoN
                        colN = msh(iM)%IEN(aopp,eOpp)
!                        write(*,*)" Adding connection ", rowN, 
!     2                           " with ", colN
                        CALL ADDCOL(rowN, colN)
                        CALL ADDCOL(colN, rowN)
                     END DO
                  END DO
               END IF
            END DO
         END IF
      END DO


!     Adding Nitsche coupling connections (only Nitsche faces with fluid 
!     nodes inside a bigger solid bounding-box)
      write(*,*)" Adding Nitsche coupling connections "
      DO iM=1, nMsh
         IF( msh(iM)%nNtsFa .EQ. 0) CYCLE         

         DO iFa = 1, msh(iM)%nFa 
            IF( .NOT.msh(iM)%fa(iFa)%isNts ) CYCLE

!           Only for Nitsche's type solid faces            
            DO a=1, msh(iM)%fa(iFa)%nNo
               rowN = msh(iM)%fa(iFa)%gN(a)

               DO jM = 1, nMsh
                  IF( msh(jM)%nNtsFa .NE. 0) CYCLE  ! solid mesh 

                  DO b = 1, msh(jM)%nNo
                     colN = msh(jM)%gN(b)

!                    useful in case we need to print                      
!                     tmpr = rowN-msh(1)%nNo-1
!                     tmpc = colN-1

!                    Check that the fluid node is in the BBox, if yes, add connection 
                     xp(:) = x(:,colN) !+ lD(:,Ac)

!                    Is the fluid element in the solid BBox? 
                     IF((xp(1) .LE. maxb(1)) .AND. (xp(1) .GE. minb(1)) 
     2                          .AND. (xp(2) .LE. maxb(2)) 
     3                          .AND. (xp(2) .GE. minb(2))) THEN

!                        write(*,*)" Adding connection ", rowN, 
!     2                                                  " with ", colN
                        CALL ADDCOL(rowN, colN)
                        CALL ADDCOL(colN, rowN)
                     END IF

                  END DO
               END DO
            END DO
         END DO
      END DO

!     Treat shells with triangular elements here
      DO iM=1, nMsh
         IF (.NOT.shlEq .OR. .NOT.msh(iM)%lShl) CYCLE
         IF (msh(iM)%eType .EQ. eType_NRB) CYCLE
         DO e=1, msh(iM)%nEl
            DO a=1, 2*msh(iM)%eNoN
               IF (a .LE. msh(iM)%eNoN) THEN
                  rowN = msh(iM)%IEN(a,e)
               ELSE
                  rowN = msh(iM)%eIEN(a-msh(iM)%eNoN,e)
               END IF
               IF (rowN .EQ. 0) CYCLE
               DO b=1, 2*msh(iM)%eNoN
                  IF (b .LE. msh(iM)%eNoN) THEN
                     colN = msh(iM)%IEN(b,e)
                  ELSE
                     colN = msh(iM)%eIEN(b-msh(iM)%eNoN,e)
                  END IF
                  IF (colN .EQ. 0) CYCLE
                  CALL ADDCOL(rowN, colN)
               END DO
            END DO
         END DO
      END DO

!     Now reset idMap for undeforming Neumann BC faces. Then insert
!     master node as a column entry in each row for all the slave nodes.
!     This step is performed even for ghost master nodes where the idMap
!     points to the ghost master node.
      flag = .FALSE.
      DO i=1, nEq
         DO j=1, eq(i)%nBc
            iM  = eq(i)%bc(j)%iM
            iFa = eq(i)%bc(j)%iFa
            IF (BTEST(eq(i)%bc(j)%bType, bType_undefNeu)) THEN
               masN = eq(i)%bc(j)%masN
               IF (masN .EQ. 0) CYCLE
               DO a=1, msh(iM)%fa(iFa)%nNo
                  rowN = msh(iM)%fa(iFa)%gN(a)
                  IF (rowN .EQ. masN) CYCLE
                  idMap(rowN) = masN
!                 Insert master to the row if not already present
                  CALL ADDCOL(rowN, masN)
               END DO
               flag = .TRUE.
            END IF
         END DO
      END DO

!     Change uInd if idMap has been changed
      IF (flag) THEN
         DO a=1, tnNo
            rowN = idMap(a)
!           If the mapping is not changed, examine the mapping of the
!           column entries of uInd. Don't do anything if the mapping is
!           unchanged. If the mapping is changed, then add to the column\
!           indices of uInd if the entry is not already present.
            IF (rowN .EQ. a) THEN
               i = 0
               DO
                  i = i + 1
                  b = uInd(i,rowN)
!                 Terminate if end of column entries are reached
                  IF (b .EQ. 0) EXIT
                  colN = idMap(b)
!                 Ignore if the column entry mapping is not changed.
!                 This entry is already present and will be used to
!                 assemble mass (D) matrix
                  IF (b .EQ. colN) CYCLE
!                 As the column entry is now mapped to a new node,
!                 search all column entries and insert the new node if
!                 it is not present. This step is performed to assemble
!                 the Divergence (C) matrix
                  CALL ADDCOL(rowN, colN)
                  IF (i .EQ. mnnzeic) EXIT
               END DO
            ELSE
!           If the row mapping is changed, insert the mapped/unmapped
!           column entries of the old row into the new row if not
!           already present.
               i = 0
               DO
                  i = i + 1
                  b = uInd(i,a)
!                 Terminate if end of column entries are reached
                  IF (b .EQ. 0) EXIT
!                 Add unmapped column to assemble gradient matrix
                  CALL ADDCOL(rowN, b)
!                 If column is mapped, add the mapped column to assemble
!                 stiffness matrix
                  colN = idMap(b)
                  IF (b .NE. colN) CALL ADDCOL(rowN, colN)
                  IF (i .EQ. mnnzeic) EXIT
               END DO
            END IF
         END DO
      END IF

!--------------------------------------------------------------------
!     Finding number of non-zeros in colPtr vector
      nnz = 0
      DO rowN=1, tnNo
         IF (uInd(1,rowN) .EQ. 0) THEN
            err = "Node "//rowN//" is isolated"
         END IF
         DO i = 1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               nnz = nnz + 1
            END IF
         END DO
      END DO

!--------------------------------------------------------------------
!     Now constructing compact form of rowPtr and colPtr
      ALLOCATE (colPtr(nnz), rowPtr(tnNo+1))
      j  = 1
      rowPtr(1) = 1
      DO rowN=1, tnNo
         DO i=1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               colPtr(j) = uInd(i,rowN)
               j = j + 1
            END IF
         END DO
         rowPtr(rowN+1) = j
      END DO
      DEALLOCATE (uInd)

      RETURN
      CONTAINS
!--------------------------------------------------------------------
         SUBROUTINE ADDCOL(row, col)
         IMPLICIT NONE
         INTEGER(KIND=IKIND), INTENT(IN) :: row, col

         INTEGER(KIND=IKIND) i, j

         i = 0
         DO
            i = i + 1
            IF (i .EQ. mnnzeic) CALL RESIZ()

!           If current entry is zero, then  fill it up
            IF (uInd(i,row) .EQ. 0) THEN
               uInd(i,row) = col
               EXIT
            END IF

!           If current entry is still smaller, keep going
            IF (col .GT. uInd(i,row)) CYCLE

!           If column entry already exists, exit
            IF (col .EQ. uInd(i,row)) EXIT

!           If we are this point, then then the current entry is bigger.
!           Shift all the entries from here to the end of the list. If
!           list is full, we request a larger list, otherwise we shift
!           and add the item at the current entry position.
            IF (uInd(mnnzeic,row) .NE. 0) CALL RESIZ()
            DO j=mnnzeic, i+1, -1
               uInd(j,row) = uInd(j-1,row)
            END DO
            uInd(i,row) = col
            EXIT
         END DO

         RETURN
         END SUBROUTINE ADDCOL
!--------------------------------------------------------------------
         SUBROUTINE RESIZ()
         IMPLICIT NONE

         INTEGER(KIND=IKIND) n
         INTEGER(KIND=IKIND), ALLOCATABLE :: tmp(:,:)

         n = mnnzeic
         ALLOCATE(tmp(n,tnNo))
         tmp(:,:) = uInd(:,:)
         DEALLOCATE(uInd)
         mnnzeic = n + MAX(5,n/5)
         ALLOCATE(uInd(mnnzeic,tnNo))
         uInd(:,:)   = 0
         uInd(1:n,:) = tmp(:,:)
         DEALLOCATE(tmp)

         RETURN
         END SUBROUTINE RESIZ
!--------------------------------------------------------------------
      END SUBROUTINE NTS_LHSA
!####################################################################
!####################################################################
!     This routine computes the diam of a mesh with TRI type elements in 2D  
      PURE SUBROUTINE COMP_DIAM(lM, maxDist)
      USE COMMOD, ONLY : mshType, x
      USE ALLFUN
      USE UTILMOD
      IMPLICIT NONE 

      TYPE(mshType), INTENT(IN) :: lM 
      REAL(KIND=RKIND), INTENT(OUT) :: maxDist
      
      REAL(KIND=RKIND) :: diam
      INTEGER(KIND=IKIND) :: e, Ac, Acn

      maxDist = TINY(maxDist)
!     compute solid mesh diamter (in the current configuration)
      DO e=1, lM%nEl

         Ac = lM%IEN(1,e)
         Acn = lM%IEN(2,e)
         diam = DIST(x(:,Ac), x(:,Acn))
         IF(diam .GT. maxDist) maxDist = diam

         Ac = lM%IEN(2,e)
         Acn = lM%IEN(3,e)
         diam = DIST(x(:,Ac),x(:,Acn))
         IF(diam .GT. maxDist) maxDist = diam

         Ac = lM%IEN(1,e)
         Acn = lM%IEN(3,e)
         diam = DIST(x(:,Ac),x(:,Acn))
         IF(diam .GT. maxDist) maxDist = diam

      END DO

      RETURN
      END SUBROUTINE COMP_DIAM
!####################################################################
!     This routine computes the diam of a mesh with TRI type elements in 2D  
      PURE SUBROUTINE COMP_ELM_DIAM(poly,eNoN,maxDist)
      USE COMMOD, ONLY : nsd, mshType
      USE ALLFUN
      USE UTILMOD
      IMPLICIT NONE 

      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: poly(nsd,eNoN) 
      REAL(KIND=RKIND), INTENT(OUT) :: maxDist
      
      REAL(KIND=RKIND) :: diam 
      INTEGER(KIND=IKIND) :: e, Ac, Acn

      maxDist = TINY(maxDist)

      IF((nsd .EQ. 2) .AND. (eNoN .EQ. 3)) THEN

         diam = DIST(poly(:,1),poly(:,2))
         IF(diam .GT. maxDist) maxDist = diam

         diam = DIST(poly(:,2),poly(:,3))
         IF(diam .GT. maxDist) maxDist = diam

         diam = DIST(poly(:,1),poly(:,3))
         IF(diam .GT. maxDist) maxDist = diam

      ELSE IF((nsd .EQ. 2) .AND. (eNoN .EQ. 2)) THEN

         diam = DIST(poly(:,1),poly(:,2))
         IF(diam .GT. maxDist) maxDist = diam

      ELSE IF ((nsd .EQ. 3) .AND. (eNoN .EQ. 4)) THEN

         diam = DIST(poly(:,1),poly(:,2))
         IF(diam .GT. maxDist) maxDist = diam

         diam = DIST(poly(:,1),poly(:,3))
         IF(diam .GT. maxDist) maxDist = diam

         diam = DIST(poly(:,1),poly(:,4))
         IF(diam .GT. maxDist) maxDist = diam

         diam = DIST(poly(:,2),poly(:,3))
         IF(diam .GT. maxDist) maxDist = diam

         diam = DIST(poly(:,2),poly(:,4))
         IF(diam .GT. maxDist) maxDist = diam

         diam = DIST(poly(:,3),poly(:,4))
         IF(diam .GT. maxDist) maxDist = diam

      ELSE IF ((nsd .EQ. 3) .AND. (eNoN .EQ. 3)) THEN

         diam = DIST(poly(:,1),poly(:,2))
         IF(diam .GT. maxDist) maxDist = diam

         diam = DIST(poly(:,2),poly(:,3))
         IF(diam .GT. maxDist) maxDist = diam

         diam = DIST(poly(:,1),poly(:,3))
         IF(diam .GT. maxDist) maxDist = diam

      ELSE 
!        *** rule for elementwise distance not defined ***
         maxDist = 0._RKIND
      END IF

      RETURN
      END SUBROUTINE COMP_ELM_DIAM
!####################################################################
!--------------------------------------------------------------------
      SUBROUTINE GETFSFace(lM, fs, lStab)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(fsType), INTENT(OUT) :: fs
      TYPE(mshType), INTENT(IN) :: lM
      LOGICAL, INTENT(IN) :: lStab

      INTEGER i, g

      CALL DESTROY(fs)
      
      IF (lStab) THEN
         fs%nG    = lM%fa(1)%fs(1)%nG
         fs%eType = lM%fa(1)%fs(1)%eType
         fs%lShpF = lM%fa(1)%fs(1)%lShpF
         fs%eNoN  = lM%fa(1)%fs(1)%eNoN
         CALL ALLOCFS(fs, nsd-1)
         fs%w   = lM%fa(1)%fs(1)%w
         fs%xi  = lM%fa(1)%fs(1)%xi
         fs%N   = lM%fa(1)%fs(1)%N
         fs%Nx  = lM%fa(1)%fs(1)%Nx
         fs%xib = lM%fa(1)%fs(1)%xib
         fs%Nb  = lM%fa(1)%fs(1)%Nb
         IF (ALLOCATED(fs%Nxx)) THEN
            fs%Nxx = lM%fa(1)%fs(1)%Nxx
         END IF
      ELSE 
         err=" Non linear FEM not accepted for Nitshe-unfitted "
      END IF

      RETURN
      END SUBROUTINE GETFSFace

!####################################################################
!     

      SUBROUTINE GET_BulkSElmQP(nbrSElm, eType, eNoN, nG, nGTot, crdElm, 
     2                                quadPntRef, quadPntCur, lstSbElm)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nbrSElm, eType, eNoN, nG, nGTot 
      REAL(KIND=RKIND), INTENT(IN) :: crdElm(nsd, eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: quadPntRef(nsd, nGTot)
      REAL(KIND=RKIND), INTENT(OUT) :: quadPntCur(nsd, nGTot)
      REAL(KIND=RKIND),INTENT(OUT) :: lstSbElm(nsd,eNon,nbrSElm)

      INTEGER(KIND=IKIND) :: FlagToDel, i, j, bng, end
      INTEGER(KIND=IKIND), ALLOCATABLE :: FlagLevel(:)                                     
      REAL(KIND=RKIND), ALLOCATABLE :: x6(:), x4(:), x5(:), xlToDel(:,:)      

C       write(*,*)" Inside GET_SElmQuadPoint "

!     Define number of face sun-element in the reference face element 
C       IF ( nsd .EQ. 2) nbrSElm = 2**(nbrFntCllSD-1)
C       IF ( nsd .EQ. 3) nbrSElm = 4**(nbrFntCllSD-1)

C       write(*,*)" nbr max sub-element = ", nbrSElm

      ALLOCATE( FlagLevel(nbrSElm), 
     2           xlToDel(nsd, eNoN) )
      FlagLevel = 0
      lstSbElm = -2
      
      IF ( nsd .EQ. 2) THEN 

         lstSbElm(1,1,1) = 1._RKIND
         lstSbElm(2,1,1) = 0._RKIND
         
         lstSbElm(1,2,1) = 0._RKIND
         lstSbElm(2,2,1) = 1._RKIND
         
         lstSbElm(1,3,1) = 0._RKIND
         lstSbElm(2,3,1) = 0._RKIND

         FlagLevel(1) = 1

         i = 1

         ALLOCATE( x4(nsd), x5(nsd), x6(nsd) )         

         DO WHILE (i .LE. nbrSElm)

            IF ( FlagLevel(i) .EQ. nbrFntCllSD ) THEN 
!              We have reached the maxminum number of subdivision
               i = i + 1
C                write(*,*)" sub elem ", i, " end sundivision"

C                write(*,*)" Flag level = ", FlagLevel
C                write(*,*)" Sub element list : "
C                write(*,*) lstSbElm

               CYCLE
            END IF

!           Store into of initial element to subdivide 
            xlToDel = lstSbElm(:,:,i)  
C             write(*,*)" xlToDel = ", xlToDel
            FlagToDel = FlagLevel(i)


            IF (lstSbElm(1,1,i+1) .LT. -1.5_RKIND) THEN 
               j = i
!              Nothing to copy, just divide
               GOTO 13
            END IF

!           Need to copy on the right the previously cretaed sub-elements 

            IF (lstSbElm(1,1,i+1) .GT. -1.5_RKIND) THEN 
               j = nbrSElm - 4
               DO WHILE (j .GE. i+1)
C                   write(*,*)" j = ", j 
C                   write(*,*)" move j to  = ", j + 3
                  lstSbElm(:,:,j+3) = lstSbElm(:,:,j)
                  FlagLevel(j+3) = FlagLevel(j)
                  j = j - 1
               END DO 
               
C                write(*,*)" Flag level = ", FlagLevel
C                write(*,*)" Sub element list : "
C                write(*,*) lstSbElm

            END IF


13          CONTINUE
!           We have 4 new subTriangles
            x4(1) = (xlToDel(1,1) + xlToDel(1,2))*0.5_RKIND
            x4(2) = (xlToDel(2,1) + xlToDel(2,2))*0.5_RKIND

            x5(1) = (xlToDel(1,2) + xlToDel(1,3))*0.5_RKIND
            x5(2) = (xlToDel(2,2) + xlToDel(2,3))*0.5_RKIND

            x6(1) = (xlToDel(1,3) + xlToDel(1,1))*0.5_RKIND
            x6(2) = (xlToDel(2,3) + xlToDel(2,1))*0.5_RKIND


            lstSbElm(:,1,i) = xlToDel(:,1)
            lstSbElm(:,2,i) = x4
            lstSbElm(:,3,i) = x6

            lstSbElm(:,1,i+1) = x4
            lstSbElm(:,2,i+1) = xlToDel(:,2)
            lstSbElm(:,3,i+1) = x5

            lstSbElm(:,1,i+2) = x6
            lstSbElm(:,2,i+2) = x5
            lstSbElm(:,3,i+2) = xlToDel(:,3)

            lstSbElm(:,1,i+3) = x6
            lstSbElm(:,2,i+3) = x4
            lstSbElm(:,3,i+3) = x5

!           Update just created level by adding 1 to previous parent level 
            FlagLevel(i:i+3) = FlagToDel + 1

C             write(*,*)" Flag level = ", FlagLevel
C             write(*,*)" Sub element list : "
C             write(*,*) lstSbElm

         END DO

         DEALLOCATE(x6, x4, x5)

      ELSE IF( nsd .EQ. 3) THEN 


!        To do sub-element dub-division 

      END IF


!     Given lstSbElm, we can now define the quadrature point in the 
!     reference divided face element 

      DO i = 1, nbrSElm

         bng = nG*(i-1)+1
         end = nG*i
         CALL GETGIP_SELM(nsd, eType, eNoN, nG, lstSbElm(:,:,i), 
     2                                           quadPntRef(:,bng:end))

      END DO

C       write(*,*)" !!!!!! Quad point ref face element = ", QpFacRef

      DO i = 1, nbrSElm

         bng = nG*(i-1)+1
         end = nG*i
         CALL GETP_SELMB(nsd, eType, eNoN, nG, crdElm,  
     2                    quadPntRef(:,bng:end), quadPntCur(:,bng:end))

      END DO

      RETURN
      END SUBROUTINE GET_BulkSElmQP

!####################################################################
!     Returns intersection point between two segments
      SUBROUTINE SegSegInt(Seg1, Seg2, int, Pint)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(OUT) :: int
      REAL(KIND=RKIND), INTENT(IN)  :: Seg1(nsd,2), Seg2(nsd,2)
      REAL(KIND=RKIND), INTENT(OUT) :: Pint(nsd)


      REAL(KIND=RKIND) a1,a2,b1,b2,c1,c2, det, dis, ep

      ep = 1.E-8_RKIND

      a1 = Seg1(2,2) - Seg1(2,1)
      b1 = Seg1(1,1) - Seg1(1,2)
      c1 = Seg1(1,1)*a1 + Seg1(2,1)*b1

      a2 = Seg2(2,2) - Seg2(2,1)
      b2 = Seg2(1,1) - Seg2(1,2)
      c2 = Seg2(1,1)*a2 + Seg2(2,1)*b2
  
      det = a1*b2 - a2*b1

      Pint = 0._RKIND
      IF( det .EQ. 0) THEN 
         int = 0
      ELSE 

         int = 1 
C          write(*,*)" we have found an intersection ", int

         Pint(1) = (b2*c1 - b1*c2)/det
         Pint(2) = (a1*c2 - a2*c1)/det

!        Check that is inside both segments 
         IF(Pint(1) .LT. MIN(Seg1(1,1),Seg1(1,2))-ep ) int = 0
         IF(Pint(1) .LT. MIN(Seg2(1,1),Seg2(1,2))-ep ) int = 0
C          write(*,*)" case 1 ", int

         IF(Pint(1) .GT. MAX(Seg1(1,1),Seg1(1,2))+ep ) int = 0
         IF(Pint(1) .GT. MAX(Seg2(1,1),Seg2(1,2))+ep ) int = 0
C          write(*,*)" case 2 ", int

         IF(Pint(2) .LT. MIN(Seg1(2,1),Seg1(2,2))-ep ) int = 0
         IF(Pint(2) .LT. MIN(Seg2(2,1),Seg2(2,2))-ep ) int = 0
C          write(*,*)" case 3 ", int

         IF(Pint(2) .GT. MAX(Seg1(2,1),Seg1(2,2))+ep ) int = 0
         IF(Pint(2) .GT. MAX(Seg2(2,1),Seg2(2,2))+ep) int = 0
C          write(*,*)" case 4 ", int

!        Check that the point is not an extreme of Seg2
         dis = DIST(Pint, Seg2(:,1))
         IF( dis .LT. ep ) int = 0
C          write(*,*)" case 5 ", int

         dis = DIST(Pint, Seg2(:,2))
         IF( dis .LT. ep ) int = 0
C          write(*,*)" case 6 ", int

      END IF 

      RETURN
      END SUBROUTINE SegSegInt
!--------------------------------------------------------------------
!     Returns Gauss integration points in current coordinates from 
!     reference face (dimFace) to current face (dimFace+1)
!     nGF = nGF total  

      SUBROUTINE CompInter( eTypeF,eNoNF, nGF, xFaceRef, quadPnt, SElm, 
     2                         qPRef, findInt)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eTypeF, nGF, eNoNF
      REAL(KIND=RKIND), INTENT(IN) :: xFaceRef(nsd, eNoNF)
      REAL(KIND=RKIND), INTENT(OUT) :: quadPnt(nsd, nGF)
      REAL(KIND=RKIND), INTENT(OUT) :: SElm(nsd,eNoNF)
      REAL(KIND=RKIND), INTENT(OUT) :: qPRef(nsd-1, nGF)
      INTEGER(KIND=IKIND), INTENT(OUT) :: findInt


      REAL(KIND=RKIND) :: crdFace(nsd, eNoNF), Pint(nsd), Pint2(nsd),
     2                    SElmRef(nsd-1,eNoNF), poly(nsd,eNoNF+1), dRef,
     3                    PintRef, dis, disP2
                       
      INTEGER(KIND=IKIND) :: i, int, find, faceID(eNoNF+1,eNoNF), 
     2                       onFace, P1in, P2in, cnt, onVertex

      findInt = 0
      SElmRef = 0._RKIND

C       write(*,*)" inside CompInter"

!     Store only info of internal (to the reference fluid element) quadrature point
      SELECT CASE(eTypeF)
      CASE(eType_LIN1)
         poly(1,1) = 1._RKIND
         poly(2,1) = 0._RKIND

         poly(1,2) = 0._RKIND
         poly(2,2) = 1._RKIND

         poly(1,3) = 0._RKIND
         poly(2,3) = 0._RKIND

         faceID(1,:) = (/2, 3/)
         faceID(2,:) = (/3, 1/)
         faceID(3,:) = (/1, 2/)


      CASE(eType_TRI3)
         poly(1,1) = 1._RKIND
         poly(2,1) = 0._RKIND
         poly(3,1) = 0._RKIND

         poly(1,2) = 0._RKIND
         poly(2,2) = 1._RKIND
         poly(3,2) = 0._RKIND

         poly(1,3) = 0._RKIND
         poly(2,3) = 0._RKIND
         poly(3,3) = 1._RKIND

         poly(1,4) = 0._RKIND
         poly(2,4) = 0._RKIND
         poly(3,4) = 0._RKIND

         faceID(1,:) = (/2, 3, 4/)
         faceID(2,:) = (/4, 3, 1/)
         faceID(3,:) = (/1, 2, 4/)
         faceID(4,:) = (/1, 3, 2/)

      END SELECT

      CALL IN_POLY(xFaceRef(:,1), poly, P1in)
      CALL IN_POLY(xFaceRef(:,2), poly, P2in)

      IF( (P1in .EQ. 1) .AND. (P2in .EQ. 1) ) THEN 

         findInt = 1
         SElm(:,1) = xFaceRef(:,1)
         SElm(:,2) = xFaceRef(:,2)

         SElmRef(:,1) = -1._RKIND
         SElmRef(:,2) = 1._RKIND
      ELSE IF( (P1in .EQ. 1) .OR. (P2in .EQ. 1) ) THEN 

!        done for the moment in 2d
!        Loop om face 
         DO i = 1, eNoNF+1

            crdFace(:,1) = poly(:,faceID(i,1))
            crdFace(:,2) = poly(:,faceID(i,2))

            CALL SegSegInt(crdFace, xFaceRef, int, Pint)

C             write(*,*)" Face ref elm", crdFace
C             write(*,*)" solid seg in ref elm  ", xFaceRef
C             write(*,*)" int ", int
C             write(*,*)" Pint ", Pint

            IF( int .EQ. 0 ) CYCLE

C             write(*,*)" Is xFaceRef(:,1) = ", xFaceRef(:,1), " inside?"
            CALL IN_POLY_FV(xFaceRef(:,1), poly, find, onFace, onVertex)
C             write(*,*) find

!           We jump the face that contains a solid vertex 
C            IF( i .EQ. onFace) THEN 
C               write(*,*)" node on face ", i 
C               CYCLE
C            END IF

!           Compute Pint in the solid face reference element
            dRef = 2._RKIND*( SQRT(NORM(Pint - xFaceRef(:,1))) / 
     2                    SQRT(NORM(xFaceRef(:,2) - xFaceRef(:,1))) ) 

            PintRef = -1._RKIND + dRef

            IF( find .EQ. 1 ) THEN 
               SElm(:,1) = xFaceRef(:,1)
               SElm(:,2) = Pint

               SElmRef(:,1) = -1._RKIND
               SElmRef(:,2) = PintRef
            ELSE 
               SElm(:,1) = Pint 
               SElm(:,2) = xFaceRef(:,2)

               SElmRef(:,1) = PintRef 
               SElmRef(:,2) = 1._RKIND
            END IF
            
            findInt = 1
            EXIT

         END DO

      ELSE

         cnt = 0
         DO i = 1, eNoNF+1

            crdFace(:,1) = poly(:,faceID(i,1))
            crdFace(:,2) = poly(:,faceID(i,2))

            CALL SegSegInt(crdFace, xFaceRef, int, Pint)  

            IF( int .EQ. 0 ) CYCLE  

            cnt = cnt + 1 

            SElm(:,cnt) = Pint

         END DO

!        Check that SElm is not a point and that the order is correct

         dis = DIST(SElm(:,1),SElm(:,2))

         IF( dis .GT. 1.E-8_RKIND ) THEN 
            findInt = 1

            dis   = DIST(SElm(:,1),xFaceRef(:,1))
            disP2 = DIST(SElm(:,1),xFaceRef(:,2))

            IF( dis .GT. disP2 ) THEN 
               Pint2     = SElm(:,2)
               SElm(:,2) = SElm(:,1)
               SElm(:,1) = Pint2
            END IF

!           Compute SElmRef in the solid face reference element
            dRef = 2._RKIND*( SQRT(NORM(SElm(:,1) - xFaceRef(:,1))) / 
     2                      SQRT(NORM(xFaceRef(:,2) - xFaceRef(:,1))) ) 

            PintRef = -1._RKIND + dRef    
            SElmRef(:,1) = PintRef       

            dRef = 2._RKIND*( SQRT(NORM(SElm(:,2) - xFaceRef(:,1))) / 
     2                      SQRT(NORM(xFaceRef(:,2) - xFaceRef(:,1))) ) 

            PintRef = -1._RKIND + dRef    
            SElmRef(:,2) = PintRef    

         END IF

      END IF
      
!     Given SElmRef, we can now define the quadrature point in the 
!     reference divided face element 

      CALL GETGIP_SELM(nsd-1, eTypeF, eNoNF, nGF, SElmRef, qPRef)

      CALL GETP_SELM(nsd-1, eTypeF, eNoNF, nGF, xFaceRef, qPRef,quadPnt)


      RETURN
      END SUBROUTINE CompInter


!####################################################################
!     Returns Gauss integration points in current coordinates from 
!     reference face (dimFace) to current face (dimFace+1)
!     nGF = nGF total  

      SUBROUTINE GET_SElmQuadPoint(nbrSElm, eTypeF, eNoNF, nG, nGF, 
     2                         crdFace, quadPnt, SElm, qPRef )
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nbrSElm, eTypeF, nG, nGF, eNoNF
      REAL(KIND=RKIND), INTENT(IN) :: crdFace(nsd, eNoNF)
      REAL(KIND=RKIND), INTENT(OUT) :: quadPnt(nsd, nGF)
      REAL(KIND=RKIND),INTENT(OUT),OPTIONAL :: SElm(nsd,eNoNF,nbrSElm)
      REAL(KIND=RKIND),INTENT(OUT),OPTIONAL :: qPRef(nsd-1, nGF)


      INTEGER(KIND=IKIND) :: FlagToDel, i, j, bng, end
      INTEGER(KIND=IKIND), ALLOCATABLE :: FlagLevel(:)                                     
      REAL(KIND=RKIND), ALLOCATABLE :: lstSbElm(:,:,:), x3(:), 
     2                                                 xlToDel(:,:)
      REAL(KIND=RKIND) :: QpFacRef(nsd-1, nGF)
      

C       write(*,*)" Inside GET_SElmQuadPoint "

!     Define number of face sun-element in the reference face element 
C       IF ( nsd .EQ. 2) nbrSElm = 2**(nbrFntCllSD-1)
C       IF ( nsd .EQ. 3) nbrSElm = 4**(nbrFntCllSD-1)

C       write(*,*)" nbr max sub-element = ", nbrSElm

      ALLOCATE( lstSbElm(nsd-1, eNoNF, nbrSElm), FlagLevel(nbrSElm), 
     2           xlToDel(nsd-1, eNoNF) )
      FlagLevel = 0
      lstSbElm = -2
      
      IF ( nsd .EQ. 2) THEN 
         IF ( eTypeF .NE. eType_LIN1) err=" eTypeF not supported "

         ALLOCATE(x3(nsd-1))

         lstSbElm(1,1,1) = -1._RKIND
         lstSbElm(1,2,1) = 1._RKIND

         FlagLevel(1) = 1

         i = 1
         
C          write(*,*)" Flag level = ", FlagLevel
C          write(*,*)" Sub element list : "
C          write(*,*) lstSbElm
         

         DO WHILE (i .LE. nbrSElm)

            IF ( FlagLevel(i) .EQ. nbrFntCllSD ) THEN 
!              We have reached the maxminum number of subdivision
               i = i + 1
C                write(*,*)" sub elem ", i, " end sundivision"

C                write(*,*)" Flag level = ", FlagLevel
C                write(*,*)" Sub element list : "
C                write(*,*) lstSbElm

               CYCLE
            END IF

!           Store into of initial element to subdivide 
            xlToDel = lstSbElm(:,:,i)  
C             write(*,*)" xlToDel = ", xlToDel
            FlagToDel = FlagLevel(i)


            IF (lstSbElm(1,1,i+1) .LT. -1.5_RKIND) THEN 
               j = i
!              Nothing to copy, just divide
               GOTO 130
            END IF

!           Need to copy on the right the previously cretaed sub-elements 

            IF (lstSbElm(1,1,i+1) .GT. -1.5_RKIND) THEN 
               j = nbrSElm - 1
               DO WHILE (j .GE. i+1)
                  lstSbElm(:,:,j+1) = lstSbElm(:,:,j)
                  FlagLevel(j+1) = FlagLevel(j)
                  j = j - 1
               END DO 
               
C                write(*,*)" Flag level = ", FlagLevel
C                write(*,*)" Sub element list : "
C                write(*,*) lstSbElm

            END IF


130         CONTINUE
!           We have 2 new subFaces 
            x3(1) = (xlToDel(1,1) + xlToDel(1,2))*0.5_RKIND

            lstSbElm(:,1,i) = xlToDel(:,1)
            lstSbElm(:,2,i) = x3

            lstSbElm(:,1,i+1) = x3 
            lstSbElm(:,2,i+1) = xlToDel(:,2)

!           Update just created level by adding 1 to previous parent level 
            FlagLevel(i:i+1) = FlagToDel + 1

C             write(*,*)" Flag level = ", FlagLevel
C             write(*,*)" Sub element list : "
C             write(*,*) lstSbElm

C             write(*,*)""
C             write(*,*)""
C             write(*,*)""
C             write(*,*)""

         END DO

C          write(*,*)" For the solid reference element we have " 
C          write(*,*)" Flag level = ", FlagLevel
C          write(*,*)" Sub element list : "
C          write(*,*) lstSbElm

         DEALLOCATE(x3)

      ELSE IF( nsd .EQ. 3) THEN 
         IF ( eTypeF .NE. eType_TRI3) err=" eTypeF not supported "

         lstSbElm(1,1,1) = 1._RKIND
         lstSbElm(2,1,1) = 0._RKIND
         
         lstSbElm(1,2,1) = 0._RKIND
         lstSbElm(2,2,1) = 1._RKIND
         
         lstSbElm(1,3,1) = 0._RKIND
         lstSbElm(2,3,1) = 0._RKIND

         FlagLevel(1) = 1

!        To do sub-element dub-division 

      END IF


!     Given lstSbElm, we can now define the quadrature point in the 
!     reference divided face element 

      DO i = 1, nbrSElm

         bng = nG*(i-1)+1
         end = nG*i
         CALL GETGIP_SELM(nsd-1, eTypeF, eNoNF, nG, lstSbElm(:,:,i), 
     2                                              QpFacRef(:,bng:end))

      END DO

C        write(*,*)" lstSbElm = ", lstSbElm
C        write(*,*)" !!!!!! Quad point ref face element = ", QpFacRef

      DO i = 1, nbrSElm

         bng = nG*(i-1)+1
         end = nG*i
         CALL GETP_SELM(nsd-1, eTypeF, eNoNF, nG, crdFace,  
     2                          QpFacRef(:,bng:end), quadPnt(:,bng:end))

      END DO

C       write(*,*)" !!!!!! Quad point in current elemem ", quadPnt

C       write(*,*)" End GET_SElmQuadPoint"

      IF (PRESENT(SElm)) THEN 

         DO i = 1, nbrSElm
            CALL GETP_SELM(nsd-1, eTypeF, eNoNF, nG, crdFace,  
     2                          lstSbElm(:,:, i), SElm(:,:,i))
         END DO
      END IF

C       write(*,*)" crdFace = ", crdFace
C       write(*,*)" SElm = ", SElm
 
      IF (PRESENT(qPRef)) qPRef = QpFacRef

      DEALLOCATE ( lstSbElm, FlagLevel, xlToDel )

      RETURN
      END SUBROUTINE GET_SElmQuadPoint

!####################################################################
!     Returns Gauss integration points in current coordinates from 
!     reference face (dimFace) to current face (dimFace)
      SUBROUTINE GETGIP_SELM(insd, eType, eNoNF, nG, crdFace, xi)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: insd, eType, nG, eNoNF
      REAL(KIND=RKIND), INTENT(IN) :: crdFace(insd,eNoNF)
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
         CALL GETNFace(insd, eType, eNoNF, xiRef(:,in), N(:,in))
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
      END SUBROUTINE GETGIP_SELM
!--------------------------------------------------------------------
!     Returns shape functions at given point (reference coords)
      SUBROUTINE GETNFace(insd, eType, eNoN, xi, N)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: insd, eType, eNoN
      REAL(KIND=RKIND), INTENT(IN)  :: xi(insd)
      REAL(KIND=RKIND), INTENT(OUT) :: N(eNoN)

      IF (eType .EQ. eType_NRB) RETURN

C       write(*,*)" Inside GETNFace "

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

      END SELECT

      RETURN
      END SUBROUTINE GETNFace
!--------------------------------------------------------------------
!     Returns quad points in current coordinates (current face, dimFace) 
!     from points (nbr of point = nbr of quad points) in the reference 
!     face (dimFace-1) 
      SUBROUTINE GETP_SELM(insd, eType, eNoNF, nG, crdFace, 
     2                                                        xin, xout)
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
         CALL GETNFace(insd, eType, eNoNF, xin(:,in), N(:,in))
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
      END SUBROUTINE GETP_SELM

!--------------------------------------------------------------------
!     Returns quad points in current coordinates (current subElm, dim) 
!     from points (nbr of point = nbr of quad points) in the reference 
!     elm (dim) 
      SUBROUTINE GETP_SELMB(insd, eType, eNoN, nG, crdFace, xin, xout)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: insd, eType, nG, eNoN
      REAL(KIND=RKIND), INTENT(IN) :: crdFace(insd, eNoN), 
     2              xin(insd, nG)
      REAL(KIND=RKIND), INTENT(OUT) :: xout(insd, nG)

      REAL(KIND=RKIND) :: in, jq, N(eNoN,nG)

C       write(*,*)" face sub elemem = ", crdFace

      DO in = 1, eNoN
!        Get nodal shape function in the reference face element 
         CALL GETNFace(insd, eType, eNoN, xin(:,in), N(:,in))
      END DO

!     N_ij is the weight that we will use to define the quad point 
!     in the new subElement in the reference bulk element with 
!     coordinates = crdFace
      xout = 0._RKIND

!     Loop over quad point          
      DO jq = 1, nG
!        Loop over element nodes 
         DO in = 1, eNoN
            xout(:, jq) = xout(:, jq) + crdFace(:,in)*N(in,jq)
         END DO
      END DO  

C       write(*,*)" quad point sub elemem ", xout

      RETURN
      END SUBROUTINE GETP_SELMB
!####################################################################
!     This subroutine assembels the element stiffness matrix into the
!     global stiffness matrix (Val sparse matrix formatted as a vector)
      SUBROUTINE DOASSEM_RC (d, eqNRow, eqNCol, lK )
      USE TYPEMOD
      USE COMMOD, ONLY: dof, rowPtr, colPtr, R, Val, nsd
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: d, eqNRow(d), eqNCol(d)
      REAL(KIND=RKIND), INTENT(IN) :: lK(dof*dof,d,d)

      INTEGER(KIND=IKIND) a, b, ptr, rowN, colN, left, right

      DO a=1, d
         rowN = eqNRow(a)
         IF (rowN .EQ. 0) CYCLE

         DO b=1, d
            colN = eqNCol(b)
            IF (colN .EQ. 0) CYCLE
            left  = rowPtr(rowN)
            right = rowPtr(rowN+1)
            ptr   = (right + left)/2
            DO WHILE (colN .NE. colPtr(ptr))
               IF (colN .GT. colPtr(ptr)) THEN
                  left  = ptr
               ELSE
                  right = ptr
               END IF
               ptr = (right + left)/2
            END DO
            Val(:,ptr) = Val(:,ptr) + lK(:,a,b)
         END DO
      END DO

      RETURN
      END SUBROUTINE DOASSEM_RC
!####################################################################
!     This subroutine assembels the element stiffness matrix into the
!     global stiffness matrix (Val sparse matrix formatted as a vector)
      SUBROUTINE DOASSEM_Opp (d, eqNRow, eqNCol, lK, lR, lKrc)
      USE TYPEMOD
      USE COMMOD, ONLY: dof, rowPtr, colPtr, R, Val
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: d, eqNRow(d), eqNCol(d)
      REAL(KIND=RKIND), INTENT(IN) :: lK(dof*dof,d,d), lR(dof,d)
      REAL(KIND=RKIND), INTENT(IN) :: lKrc(dof*dof,d,d)

      INTEGER(KIND=IKIND) a, b, ptr, rowN, colN, left, right

      DO a=1, d
         rowN = eqNRow(a)
         IF (rowN .EQ. 0) CYCLE
         R(:,rowN) = R(:,rowN) + lR(:,a)

!        Loop for lK         
         DO b=1, d
            colN = eqNRow(b)
            IF (colN .EQ. 0) CYCLE
            left  = rowPtr(rowN)
            right = rowPtr(rowN+1)
            ptr   = (right + left)/2
            DO WHILE (colN .NE. colPtr(ptr))
               IF (colN .GT. colPtr(ptr)) THEN
                  left  = ptr
               ELSE
                  right = ptr
               END IF
               ptr = (right + left)/2
            END DO
            Val(:,ptr) = Val(:,ptr) + lK(:,a,b)
         END DO

!        Loop for lKrc
         DO b=1, d
            colN = eqNCol(b)
            IF (colN .EQ. 0) CYCLE
            left  = rowPtr(rowN)
            right = rowPtr(rowN+1)
            ptr   = (right + left)/2
            DO WHILE (colN .NE. colPtr(ptr))
               IF (colN .GT. colPtr(ptr)) THEN
                  left  = ptr
               ELSE
                  right = ptr
               END IF
               ptr = (right + left)/2
            END DO
            Val(:,ptr) = Val(:,ptr) + lKrc(:,a,b)
         END DO

      END DO

      RETURN
      END SUBROUTINE DOASSEM_Opp
!####################################################################
!     This subroutine assembels the element tangent matrices and red
!     for Nitsche into the global stiffness matrix
      SUBROUTINE NTS_ASSEMBLY (eNoN, eNoNF, ptrFlu, ptrSOl, lRS, lKS, 
     2                                                      lKFS, lKSF)
      USE TYPEMOD
      USE COMMOD, ONLY: dof, rowPtr, colPtr, R, Val
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, eNoNF 
      INTEGER(KIND=IKIND), INTENT(IN) :: ptrFlu(eNoN), ptrSOl(eNoNF)
      REAL(KIND=RKIND), INTENT(IN) :: lRS(dof,eNoNF)
      REAL(KIND=RKIND), INTENT(IN) :: lKS(dof*dof,eNoNF,eNoNF)
      REAL(KIND=RKIND), INTENT(IN) :: lKFS(dof*dof,eNoN,eNoNF)
      REAL(KIND=RKIND), INTENT(IN) :: lKSF(dof*dof,eNoNF,eNoN)

      INTEGER(KIND=IKIND) a, b, ptr, rowN, colN, left, right


C       write(*,*)" Starting assembly"

!     Adding Solid-Solid contribution    
      DO a=1, eNoNF
         rowN = ptrSOl(a)
         IF (rowN .EQ. 0) CYCLE
         R(:,rowN) = R(:,rowN) + lRS(:,a)

         DO b=1, eNoNF
            colN = ptrSOl(b)
            IF (colN .EQ. 0) CYCLE
            left  = rowPtr(rowN)
            right = rowPtr(rowN+1)
            ptr   = (right + left)/2
            DO WHILE (colN .NE. colPtr(ptr))
               IF (colN .GT. colPtr(ptr)) THEN
                  left  = ptr
               ELSE
                  right = ptr
               END IF
               ptr = (right + left)/2
            END DO

            Val(:,ptr) = Val(:,ptr) + lKS(:,a,b)
         END DO
      END DO

C       write(*,*)" end Solid-Solid "

!     Adding Fluid-Solid contribution  lKFS
      DO a=1, eNoN
         rowN = ptrFlu(a)
         IF (rowN .EQ. 0) CYCLE

         DO b=1, eNoNF
            colN = ptrSOl(b)
            IF (colN .EQ. 0) CYCLE
            left  = rowPtr(rowN)
            right = rowPtr(rowN+1)
            ptr   = (right + left)/2
            DO WHILE (colN .NE. colPtr(ptr))
               IF (colN .GT. colPtr(ptr)) THEN
                  left  = ptr
               ELSE
                  right = ptr
               END IF
               ptr = (right + left)/2
            END DO

            Val(:,ptr) = Val(:,ptr) + lKFS(:,a,b)


         END DO
      END DO
C       write(*,*)" end Fluid-Solid "

!     Adding Solid-Fluid contribution  lKSF
      DO a=1, eNoNF
         rowN = ptrSOl(a)
         IF (rowN .EQ. 0) CYCLE

         DO b=1, eNoN
            colN = ptrFlu(b) 
            IF (colN .EQ. 0) CYCLE
            left  = rowPtr(rowN)
            right = rowPtr(rowN+1)
            ptr   = (right + left)/2
            DO WHILE (colN .NE. colPtr(ptr))
               IF (colN .GT. colPtr(ptr)) THEN
                  left  = ptr
               ELSE
                  right = ptr
               END IF
               ptr = (right + left)/2
            END DO

            Val(:,ptr) = Val(:,ptr) + lKSF(:,a,b)


         END DO
      END DO

C       write(*,*)" end Solid-Fluid  "

      RETURN
      END SUBROUTINE NTS_ASSEMBLY