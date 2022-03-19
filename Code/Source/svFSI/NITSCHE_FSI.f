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
!     TODO
!
!--------------------------------------------------------------------
      SUBROUTINE CONSTRUCT_NITSCHE_FSI(lM, Ag, Yg, Dg, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo), lD(tDof, tnNo)

      LOGICAL :: vmsStab, flag
      INTEGER(KIND=IKIND) a, e, g, l, Ac, eNoN, cPhys, iFn, nFn, N
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)
      TYPE(fsType) :: fs(2)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   dl(:,:), bfl(:,:), fN(:,:), pS0l(:,:), pSl(:), ya_l(:),
     3   lR(:,:), lK(:,:,:), lKd(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xwl(:,:), xql(:,:), Nwx(:,:),
     2   Nwxx(:,:), Nqx(:,:), xp(:), xiCur(:,:)

      INTEGER(KIND=IKIND) maxSubTri, maxQuadPnt, bgn, end, cnt, is, As,
     2   find, FlagToDel, i, j, maxLevel
      INTEGER(KIND=IKIND), ALLOCATABLE :: FlagSubTri(:), FlagLevel(:)
      REAL(KIND=RKIND), ALLOCATABLE :: lstSubTri(:,:,:), lstQdPnt(:,:)
      REAL(KIND=RKIND) :: poly(nsd,msh(2)%eNoN),
     2                    xlToDel(nsd,msh(2)%eNoN), x4(nsd), x5(nsd), 
     3                    x6(nsd), JacSub

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
     4   lKd(dof*nsd,eNoN,eNoN), xp(nsd))

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

!        Definition of the integration to use in this element based on physics 
         IF((cPhys.EQ.phys_fluid).AND.(intFElmFlag(e) .EQ. 2)) GOTO 20
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

         ALLOCATE( xiCur(nsd,fs(1)%nG) )

         maxSubTri = 16 
         maxQuadPnt = 16*fs(1)%nG
!        Allocate lstSubTri: list of sub element coordinates in the reference element 
         ALLOCATE ( lstSubTri(nsd,fs(1)%eNoN,maxSubTri) )
         ALLOCATE ( FlagSubTri(maxSubTri) )
         ALLOCATE ( FlagLevel(maxSubTri) )
!        Allocate lstQdPnt: list of quad point in the reference element   
         ALLOCATE ( lstQdPnt(nsd,maxQuadPnt) )

         lstSubTri  = -1._RKIND
         lstQdPnt   = -1._RKIND
         FlagSubTri = -1
         FlagLevel  = -1
         maxLevel = 3

!        Initialize lists with reference element 
         lstSubTri(:,1,1) = (/ 1._RKIND , 0._RKIND /)
         lstSubTri(:,2,1) = (/ 0._RKIND , 1._RKIND /)
         lstSubTri(:,3,1) = (/ 0._RKIND , 0._RKIND /)
         !FlagSubTri(1) = -1
         FlagLevel(1) = 1
         lstQdPnt(:,1:fs(1)%nG) = lM%xi(:,:)

!        Loop over subTriList   
         i = 1
         DO WHILE (i .LE. maxSubTri)
            
            bgn = i*3-(fs(1)%nG-1)
            end = i*3

            IF (lstSubTri(1,1,i) .LT. -0.5_RKIND) THEN 
C                write(*,*)" We stop here for fluid elem ", e
               EXIT
            END IF

!           Plots 
C             write(*,*)"FlagSubTri: ", FlagSubTri(i)
C             write(*,*)"FlagLevel: ", FlagLevel(i)
C             write(*,*)"lstSubTri: "
C             write(*,*) lstSubTri(:,1,i)
C             write(*,*) lstSubTri(:,2,i)
C             write(*,*) lstSubTri(:,3,i)
C             write(*,*)"lstQdPnt: "
C             write(*,*) lstQdPnt(:,bgn)
C             write(*,*) lstQdPnt(:,bgn+1)
C             write(*,*) lstQdPnt(:,bgn+2)
C             write(*,*)""
C             write(*,*)""

            CALL GET_PntRefToCur(lM, fs(1), xl, lstQdPnt(:,bgn:end), 
     2                                                          xiCur)

C             write(*,*)" GET_PntRefToCur done, xiCur = "
C             write(*,*) xiCur(:,1)
C             write(*,*) xiCur(:,2)
C             write(*,*) xiCur(:,3)
            
            cnt = 0
            DO g = 1, fs(1)%nG 
               DO a = 1, msh(2)%nEl
                  DO is = 1, msh(2)%eNoN
                     As = msh(2)%IEN(is,a)
                     poly(:,is) = x(:,As) + lD(nsd+2:2*nsd+1,As)
                  END DO

!                 Is the fluid node inside any solid element?
                  find = IN_POLY(xiCur(:,g),poly)

                  IF (find .EQ. 1) THEN 
                     cnt = cnt + 1
                     EXIT
                  END IF

               END DO 
            END DO 

C             write(*,*)" seach on solid msh done, cnt = ", cnt

!           FlagSubTri = 1: fluid element, 0: solid element, 2: to be divided
            IF ( cnt .EQ. 0 ) THEN 
               FlagSubTri(i) = 1
C                write(*,*)" Flag sub elem ", i , "as fluid "
            ELSE IF ( cnt .EQ. fs(1)%nG ) THEN 
               FlagSubTri(i) = 0
C                write(*,*)" Flag sub elem ", i , "as solid "
            ELSE 
               FlagSubTri(i) = 2
C                write(*,*)" Flag sub elem ", i , "as to divide "
            END IF

!           Check if we have reached the maximum level of subdivision
            IF ( FlagLevel(i) .EQ. maxLevel ) THEN 
               IF ( cnt .EQ. 1 ) FlagSubTri(i) = 1
               IF ( cnt .EQ. 2 ) FlagSubTri(i) = 1
C                IF ( cnt .EQ. 0 ) write(*,*)"!! ATTENTION, we do not 
C      2                                                     int here !! "
!              This should be better, but we have to check that we integrate somewhere                
!              or we will have conditioning issue                 
C              IF ( cnt .EQ. 2 ) FlagSubTri(i) = 0
            END IF

            IF ( FlagSubTri(i) .EQ. 2 ) THEN

C                write(*,*)" Dividing subElm " , i , " of element ", e

!              Destroy lstSubTri(i), move all the next to left if the next is negative           
!              and at the end of the list add the new sub tri 
               xlToDel = lstSubTri(:,:,i)  
               FlagToDel = FlagLevel(i)

               IF (lstSubTri(1,1,i+1) .LT. -0.5_RKIND) THEN 
                  j = i
C                   write(*,*)" Nothing to copy, just divide"
                  GOTO 30
               END IF

               DO j = i+1, maxSubTri
                  IF (lstSubTri(1,1,j) .LT. -0.5_RKIND) THEN 
                     EXIT 
                  END IF

C                   write(*,*)" Need to copy subtri ", j
                  
                  bgn = j*3-(fs(1)%nG-1)
                  end = j*3

                  lstSubTri(:,:,j-1) = lstSubTri(:,:,j)
                  FlagSubTri(j-1) = FlagSubTri(j)
                  FlagLevel(j-1) = FlagLevel(j)
                  lstQdPnt(:,bgn-3:end-3) = lstQdPnt(:,bgn:end)

               END DO
!              We add to j-1 the new subTri and quad point 
               j = j-1 

30             CONTINUE
!              We have 4 new subTriangles
               x4(1) = (xlToDel(1,1) + xlToDel(1,2))*0.5_RKIND
               x4(2) = (xlToDel(2,1) + xlToDel(2,2))*0.5_RKIND

               x5(1) = (xlToDel(1,2) + xlToDel(1,3))*0.5_RKIND
               x5(2) = (xlToDel(2,2) + xlToDel(2,3))*0.5_RKIND

               x6(1) = (xlToDel(1,3) + xlToDel(1,1))*0.5_RKIND
               x6(2) = (xlToDel(2,3) + xlToDel(2,1))*0.5_RKIND

               lstSubTri(:,1,j) = xlToDel(:,1)
               lstSubTri(:,2,j) = x4
               lstSubTri(:,3,j) = x6
C                write(*,*)" SubTri 1 "
C                write(*,*) lstSubTri(:,1,j)
C                write(*,*) lstSubTri(:,2,j)
C                write(*,*) lstSubTri(:,3,j)

               lstSubTri(:,1,j+1) = x4
               lstSubTri(:,2,j+1) = xlToDel(:,2)
               lstSubTri(:,3,j+1) = x5
C                write(*,*)" SubTri 2 "
C                write(*,*) lstSubTri(:,1,j+1)
C                write(*,*) lstSubTri(:,2,j+1)
C                write(*,*) lstSubTri(:,3,j+1)

               lstSubTri(:,1,j+2) = x6
               lstSubTri(:,2,j+2) = x5
               lstSubTri(:,3,j+2) = xlToDel(:,3)
C                write(*,*)" SubTri 3 "
C                write(*,*) lstSubTri(:,1,j+2)
C                write(*,*) lstSubTri(:,2,j+2)
C                write(*,*) lstSubTri(:,3,j+2)

               lstSubTri(:,1,j+3) = x6
               lstSubTri(:,2,j+3) = x4
               lstSubTri(:,3,j+3) = x5
C                write(*,*)" SubTri 4 "
C                write(*,*) lstSubTri(:,1,j+3)
C                write(*,*) lstSubTri(:,2,j+3)
C                write(*,*) lstSubTri(:,3,j+3)

!              and 4*3 new subQuad points, we need the sub-element coord in the ref element 
               bgn = j*3-(fs(1)%nG-1)
               end = j*3
               CALL GET_CURQUAD(lM, fs(1), lstSubTri(:,:,j), 
     2                                              lstQdPnt(:,bgn:end))

               bgn = (j+1)*3-(fs(1)%nG-1)
               end = (j+1)*3
               CALL GET_CURQUAD(lM, fs(1), lstSubTri(:,:,j+1), 
     2                                              lstQdPnt(:,bgn:end))

               bgn = (j+2)*3-(fs(1)%nG-1)
               end = (j+2)*3
               CALL GET_CURQUAD(lM, fs(1), lstSubTri(:,:,j+2), 
     2                                              lstQdPnt(:,bgn:end))

               bgn = (j+3)*3-(fs(1)%nG-1)
               end = (j+3)*3
               CALL GET_CURQUAD(lM, fs(1), lstSubTri(:,:,j+3), 
     2                                              lstQdPnt(:,bgn:end))


C                i = i - 1 

               FlagSubTri(j:j+3) = -1
               FlagLevel(j:j+3) = FlagToDel + 1

            ELSE
               i = i + 1
            END IF
         
         END DO   

         maxSubTri = i - 1 
C          write(*,*)" element fluid ", e, " nbr sub elm ", maxSubTri
C          write(*,*)" FlagSubTri: ", FlagSubTri(1:maxSubTri)
C          write(*,*)" FlagLevel: ", FlagLevel(1:maxSubTri)


!        Loop over subElm
         DO i = 1, maxSubTri
!           If it is a full solid, cycle
            IF ( FlagSubTri(i) .EQ. 0) THEN 
C                write(*,*)" sub element ", i, " hidden "
               CYCLE
            END IF

            bgn = j*3-(fs(1)%nG-1)-1

!           Loop over quad point 1
            DO g=1, fs(1)%nG

!              Nx are constant, we can use the previous function               
!              If this is not the case anymore, this needs to be modified 
!              fs(2)%lShpF is true here                
               IF (g.EQ.1) THEN

!                 Evaluate test function in the new quadrature points
!                 Quad point in the sub-elm of the ref element                   
                  CALL GETGNN(nsd, lM%eType, fs(1)%eNoN,  
     2                 lstQdPnt(:,bgn+g), fs(1)%N(:,g), fs(1)%Nx(:,:,g))
                  CALL GETGNN(nsd, lM%eType, fs(2)%eNoN,  
     2                 lstQdPnt(:,bgn+g), fs(2)%N(:,g), fs(2)%Nx(:,:,g))

                  CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, 
     2                 Jac, ksix)
                  CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx,
     2                 Jac, ksix)
                  IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
               END IF

!              Jac from ref element to subElm in the reference element
               CALL COMP_SubElmJAC(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), 
     2              lstSubTri(:,:,i), JacSub)   

!              Compute new weight (w * J k_ref -> K_cut * J k_ref -> K_sub )
               w = fs(1)%w(g) * Jac  *   JacSub        

               IF (nsd .EQ. 3) THEN
                  CALL FLUID3D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)
               ELSE IF (nsd .EQ. 2) THEN
                  CALL FLUID2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)
               END IF
            END DO

!           Set function spaces for velocity and pressure.
            CALL GETTHOODFS(fs, lM, vmsStab, 2)

!           Loop over quad point 2
            DO g=1, fs(2)%nG
               IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
                  CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, 
     2            Jac, ksix)
                  IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
               END IF

               IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
                  CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, 
     2            Jac, ksix)
                  IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
               END IF
               w = fs(2)%w(g) * Jac

               IF (nsd .EQ. 3) THEN
                  CALL FLUID3D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               ELSE IF (nsd .EQ. 2) THEN
                  CALL FLUID2D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w,
     2               ksix, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx,
     3               al, yl, bfl, lR, lK)

               END IF
            END DO

         END DO

         DEALLOCATE( xiCur, lstSubTri, FlagSubTri, FlagLevel, lstQdPnt)

!-----------------------------------------------------------------------
100      CONTINUE

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

!     This need to be removed once we have a ghost penalty stabilization 
      CALL CLEAN_GHOSTPNT

      RETURN
      END SUBROUTINE CONSTRUCT_NITSCHE_FSI
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

      INTEGER(KIND=IKIND) :: e, a, Ac, nbrSN, is, As, i, j, count
      REAL(KIND=RKIND) :: xp(nsd), poly(nsd,msh(1)%eNoN), xs(nsd)
      REAL(KIND=RKIND) :: minXs(nsd), maxXs(nsd)
      REAL(KIND=RKIND) :: minb(nsd), maxb(nsd)
      INTEGER(KIND=IKIND) :: find
      INTEGER(KIND=IKIND) :: maxNbrSNd = 0
      INTEGER(KIND=IKIND), ALLOCATABLE :: FNdFlag(:), cnt(:), 
     2                                                      FElmSNd(:,:)

      write(*,*)" Calling NTS_INIT "

!     We will make the assumption that the first mesh is the fluid and the 
!     second the solid, this can be generalized 


!     Define fluid mesh stencil 
      CALL GETNSTENCIL(msh(1))
      CALL GETNSTENCIL(msh(2))

      write(*,*)" Stencil done "

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
!     msh(1)%mDiam = COMP_DIAM(msh(1))
      msh(2)%mDiam = COMP_DIAM(msh(2))
      write(*,*)" Diam solid computed "

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

      write(*,*)"S BBOX NTS_INIT x-dir is: ", minb(1), " and ", maxb(1)
      write(*,*)"S BBOX NTS_INITy-dir is: ", minb(2), " and ", maxb(2)

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
               find = IN_POLY(xp,poly)

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

!     Filling intFElmFlag: 1 normal, 2 intersected from Nitsche Bnd, 0 hidden element 
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

!     Filling mapSNdFElm, find the fluid element that contains the corresponding solid node
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
               find = IN_POLY(xp,poly)

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

!     Filling FElmSNd, loop over mapSNdFElm and add the analogous contribution 
      DO a = 1, msh(2)%nNo

C          Ac = msh(2)%gN(a)

         e = mapSNdFElm(a)
         cnt(e) = cnt(e) + 1

         FElmSNd( e, cnt(e)) = a

      END DO

C       write(*,*)" cnt = ", cnt

      maxNbrSNd = MAXVAL(cnt)

!     Filling mapFElmSNd from FElmSNd 
      ALLOCATE(mapFElmSNd(msh(1)%nEl ,maxNbrSNd))
      mapFElmSNd = 0
      mapFElmSNd(:,:) = FElmSNd(:,1:maxNbrSNd)

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

      LOGICAL flag
      INTEGER(KIND=IKIND) a, b, e, i, j, rowN, colN, iM, iFa, masN,
     2   mnnzeic, jM, tmpr, tmpc

      INTEGER(KIND=IKIND), ALLOCATABLE :: uInd(:,:)

      INTEGER(KIND=IKIND) :: Ac
      REAL(KIND=RKIND) :: xp(nsd)
      REAL(KIND=RKIND) :: minXs(nsd), maxXs(nsd)
      REAL(KIND=RKIND) :: minb(nsd), maxb(nsd)

      ALLOCATE(idMap(tnNo))



!     Compute B-Box around solid with fluid mesh size
      msh(1)%mDiam = COMP_DIAM(msh(1))
      write(*,*)" fluid mesh diam ", msh(1)%mDiam

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
         minb(i) = minXs(i) - 10._RKIND*msh(1)%mDiam 
         maxb(i) = maxXs(i) + 10._RKIND*msh(1)%mDiam 
      END DO

      write(*,*)"The solid BBOX x-dir is: ", minb(1), " and ", maxb(1)
      write(*,*)"The solid BBOX y-dir is: ", minb(2), " and ", maxb(2)


!     Select fluid-fluid, solid-solid and fluid-solid node to connect            
      write(*,*)" tnNo = ", tnNo
      DO a=1, tnNo
         idMap(a) = a
      END DO

      mnnzeic = 10*MAXVAL(msh%eNoN)

!     First fill uInd array depending on mesh connectivity as is.
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

!     Adding Nitsche coupling connections (only Nitsche faces with all fluid nodes)
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
