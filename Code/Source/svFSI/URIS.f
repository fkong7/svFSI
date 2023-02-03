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
!     Collection of subroutines for the Unfitted version of the 
!     Resistive Immersed Surface method. 
!
!--------------------------------------------------------------------
!####################################################################
!     This subroutine computes the mean pressure and flux on the 
!     ris surface 
      SUBROUTINE URIS_MEANP
      USE TYPEMOD
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) :: iEq, nPrj, m, s, e, i, iM, iFa
      REAL(KIND=RKIND) :: tmp
      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:)

!     Let's conpute the mean pressure in the two regions of the fluid mesh 
!     For the moment let's define a flag IdSubDmn(size the number of elements)
      INTEGER(KIND=IKIND) :: IdSubDmn(msh(1)%nEl),  a, Ac, cntM, cntP
      REAL(KIND=RKIND) :: volU, volD, Deps, zSurf, meanPU, meanPD, 
     2                    sDST(1,tnNo), sUPS(1,tnNo), sImm(1,tnNo)


      iM = 1

      Deps = 0.13_RKIND
      zSurf = 1.25_RKIND

      IdSubDmn = 0

!     For each element, select the nodes, check the z-comp
!     If z < zSurf - Deps -> Flag is -1      
!     If z > zSurf - Deps -> Flag is +1      
!     Loop over all elements of mesh
      DO e=1, msh(iM)%nEl
         cntM = 0
         cntP = 0
         DO a=1, msh(iM)%eNoN
            Ac = msh(iM)%IEN(a,e)
            
            IF( x(3,Ac) .LT. zSurf - Deps ) THEN 
               cntM = cntM + 1
            ELSE IF ( x(3,Ac) .GT. zSurf + Deps ) THEN 
               cntP = cntP + 1
            END IF
         END DO

         IF( cntM .EQ. 3 ) IdSubDmn(e) = -1
         IF( cntP .EQ. 3 ) IdSubDmn(e) = +1
      END DO


!     Let's compute left side 
      sUPS = 0._RKIND
      DO e=1, msh(iM)%nEl
         DO a=1, msh(iM)%eNoN
            Ac = msh(iM)%IEN(a,e)
            
            IF( x(3,Ac) .LT. zSurf - Deps ) THEN 
               sUPS(1,Ac) = 1._RKIND
            ELSE IF ( x(3,Ac) .GT. zSurf + Deps ) THEN 
               sUPS(1,Ac) = 0._RKIND
            END IF
         END DO
      END DO
      volU = Integ(iM,sUPS)
      write(*,*)" volume upstream ", volU

!     Let's compute right side 
      sDST = 0._RKIND
      DO e=1, msh(iM)%nEl
         DO a=1, msh(iM)%eNoN
            Ac = msh(iM)%IEN(a,e)
            
            IF( x(3,Ac) .LT. zSurf - Deps ) THEN 
               sDST(1,Ac) = 0._RKIND
            ELSE IF ( x(3,Ac) .GT. zSurf + Deps ) THEN 
               sDST(1,Ac) = 1._RKIND
            END IF
         END DO
      END DO
      volD = Integ(iM,sDST)
      write(*,*)" volume downstream ", volD


!     Let's compute immersed side 
C       sImm = 1._RKIND
C       DO e=1, msh(iM)%nEl
C          DO a=1, msh(iM)%eNoN
C             Ac = msh(iM)%IEN(a,e)
            
C             IF( x(3,Ac) .LT. zSurf - Deps ) THEN 
C                sImm(1,Ac) = 0._RKIND
C             ELSE IF ( x(3,Ac) .GT. zSurf + Deps ) THEN 
C                sImm(1,Ac) = 0._RKIND
C             END IF
C          END DO
C       END DO
C       volD = Integ(iM,sImm)
C       write(*,*)" volume inside ", volD



!     Now we can compute the pressure mean on each subdomain
      iEq = 1

      IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
      ALLOCATE (tmpV(maxnsd,tnNo))

      meanPU = 0._RKIND
      meanPD = 0._RKIND

      m = 1
      s = eq(iEq)%s + nsd 
      e = s + m - 1

      tmpV(1:m,:) = Yn(s:e,:)*sUPS
      
      meanPU = Integ(iM,tmpV)/volU

      tmpV(1:m,:) = Yn(s:e,:)*sDST
      meanPD = Integ(iM,tmpV)/volD

      write(*,*)" mean P upstream ", meanPU
      write(*,*)" mean P downstream ", meanPD

!     If the uris was active, check the 

      IF( (meanPU .GT. meanPD) .AND. (cntURIS .GT. 10 ) ) THEN 
         urisActFlag = .FALSE.
         write(*,*)" Change urisActFlag to FALSE, deactivate URIS "
         cntURIS = 0
      END IF

      IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)

      RETURN
      END SUBROUTINE URIS_MEANP
!####################################################################
!--------------------------------------------------------------------
!     This subroutine computes the mean velocity in the fluid elements near the 
!     immersed surface  
      SUBROUTINE URIS_MEANV
      USE TYPEMOD
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) :: iEq, nPrj, m, s, e, i, iM, iFa
      REAL(KIND=RKIND) :: tmp
      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:)

!     Let's conpute the mean pressure in the two regions of the fluid mesh 
!     For the moment let's define a flag IdSubDmn(size the number of elements)
      INTEGER(KIND=IKIND) :: IdSubDmn(msh(1)%nEl),  a, Ac, cntM, cntP
      REAL(KIND=RKIND) :: volI, Deps, zSurf, meanV(nsd), 
     2                    sImm(1,tnNo)

      

      iM = 1

      Deps = 0.13_RKIND
      zSurf = 1.25_RKIND

      IdSubDmn = 0

!     For each element, select the nodes, check the z-comp
!     If z < zSurf - Deps -> Flag is -1      
!     If z > zSurf - Deps -> Flag is +1      

!     Let's compute immersed side 
      sImm = 1._RKIND
      DO e=1, msh(iM)%nEl
         DO a=1, msh(iM)%eNoN
            Ac = msh(iM)%IEN(a,e)
            
            IF( x(3,Ac) .LT. zSurf - Deps ) THEN 
               sImm(1,Ac) = 0._RKIND
            ELSE IF ( x(3,Ac) .GT. zSurf + Deps ) THEN 
               sImm(1,Ac) = 0._RKIND
            END IF
         END DO
      END DO
      volI = Integ(iM,sImm)
      write(*,*)" volume inside ", volI


!     Now we can compute the vel mean for component 3
      iEq = 1

      IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
      ALLOCATE (tmpV(maxnsd,tnNo))

      meanV = 0._RKIND

      m = 1
      s = eq(iEq)%s + (nsd-1) 
      e = s + m - 1

      tmpV(1:m,:) = Yn(s:e,:)*sImm
      
      meanV(3) = Integ(iM,tmpV)/volI




C !     Now we can compute the vel mean on each subdomain
C       iEq = 1

C       IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
C       ALLOCATE (tmpV(maxnsd,tnNo))

C       meanV = 0._RKIND

C       m = nsd
C       s = eq(iEq)%s 
C       e = s + m - 1

C       tmpV(1:m,:) = Yn(s:e,:)*sImm
      
C       meanV = Integ(iM,tmpV,1,nsd)/volI

      write(*,*)" mean Vel ", meanV

!     If the uris was active, check the 
      IF( (meanV(3) .LT. 0._RKIND) .AND. (cntURIS .GT. 5 )) THEN 
         urisActFlag = .TRUE.
         write(*,*)" Change urisActFlag to TRUE, activate URIS "
         cntURIS = 0
      END IF

      IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)

      RETURN
      END SUBROUTINE URIS_MEANV
!####################################################################
!--------------------------------------------------------------------
!     This subroutine compute the disp of the immersed surface with 
!     fem projection  
      SUBROUTINE URIS_UpdateDisp
      USE TYPEMOD
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

C       REAL(KIND=RKIND), INTENT(IN) :: lDo(nsd,tnNo)
C       REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      INTEGER(KIND=IKIND) :: flag, iM, jM, iEln, a, nd, Ac
      REAL(KIND=RKIND) :: xp(nsd), xi(nsd), d(nsd)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), N(:), Nxi(:,:)
      LOGICAL :: ultra, fl

      IF(.NOT.ALLOCATED(uris%Yd)) THEN 
         ALLOCATE(uris%Yd(nsd,uris%msh(1)%nNo))
         uris%Yd(:,nd) = 0._RKIND
      END IF

!     For each point in the immersed surface we need to localize it 
!     = find the fluid element that contains the node
      iM = 1

      flag = 1
      ultra = .FALSE.
      xi = 0.5_RKIND
      
C       DO nd = 1, uris%msh(iM)%nNo
      nd = 0
      DO WHILE(nd .LE. uris%msh(iM)%nNo)
         
!        Check if we were able to find the tetra, if not run the search with 
!        an approximation 
         IF(flag .EQ. 0) THEN 
            IF( ultra ) THEN 
               write(*,*)"*** ERROR, tet not found for para nd ",nd-1
               nd = nd + 1
               ultra = .FALSE.
            ELSE 
               ultra = .TRUE.
               GOTO 123
            END IF
         ELSE
            nd = nd + 1
            ultra = .FALSE.
         END IF

C          write(*,*)" looking node ", nd,": ", xp
         xp = uris%x(:,nd) !+ uris%Yd(:,nd)

123      CONTINUE
         flag = 1

         DO jM=1, nMsh
            ALLOCATE(xl(nsd,msh(jM)%eNoN), N(msh(jM)%eNoN), 
     2                                          Nxi(nsd,msh(jM)%eNoN))

C             write(*,*)" here 1 "
            DO iEln = 1, msh(jM)%nEl
C                flag = 0
C                write(*,*)" here 2 "
               DO a=1, msh(jM)%eNoN
                  Ac = msh(jM)%IEN(a,iEln)
C                   xl(:,a)  = ( x(:,Ac) + disp)
C                   write(*,*)" here 3 " 
C                   IF(mvMsh) THEN 
C                      xl(:,a) = x(:,Ac) + Dn(nsd+2:2*nsd+1,Ac) ! problem here 
C                   ELSE 
                     xl(:,a) = x(:,Ac) 
C                   END IF
               END DO 

C                write(*,*)" here 4 "
C !               Check if it is inside 
C                write(*,*)" xp ", xp 
C                write(*,*)" here 5 "
C                write(*,*)" xl ", xl(:,1)
C                write(*,*)" xl ", xl(:,2)
C                write(*,*)" xl ", xl(:,3)
C                write(*,*)" xl ", xl(:,4)
C                write(*,*)" ultra ", ultra 
               CALL insideTet(msh(jM)%eNoN,xp,xl,flag, ultra)

               IF( flag .EQ. 1) THEN 
C                   write(*,*)" nd parav ",nd-1," inside ",iEln,": ",xl 
C                   write(*,*)" nd ",nd," inside ",iEln,": ",xl 
C                   write(*,*)" nd parav ",nd-1 ," inside ",iEln-1
C                   write(*,*)" with ultra ", ultra 

!                 Get displacement  
!                 Localize p inside the parent element  
                  CALL GETXI(msh(jM)%eType,msh(jM)%eNoN, xl, xp, xi,fl)  
                  IF( .NOT.fl) write(*,*)" GETXI not converging "

!                 evaluate N at xi 
                  CALL GETGNN(nsd,msh(jM)%eType,msh(jM)%eNoN,xi,N,Nxi)

!                 use this to compute disp al node xp 
                  d = 0._RKIND
                  DO a=1, msh(jM)%eNoN
                     Ac = msh(jM)%IEN(a,iEln)
C                      d(1) = d(1) - N(a)*Dn(nsd+2,Ac) 
C                      d(2) = d(2) - N(a)*Dn(nsd+3,Ac) 
C                      d(3) = d(3) - N(a)*Dn(nsd+4,Ac) 

C                      d(1) = d(1) - N(a)*Dn(1,Ac) 
C                      d(2) = d(2) - N(a)*Dn(2,Ac) 
C                      d(3) = d(3) - N(a)*Dn(3,Ac) 
!                    We have to use Do because Dn contains the result 
!                    coming from the solid 
                     d(1) = d(1) - N(a)*Dn(nsd+2,Ac) 
                     d(2) = d(2) - N(a)*Dn(nsd+3,Ac) 
                     d(3) = d(3) - N(a)*Dn(nsd+4,Ac) 
                  END DO

!                 update uris disp                                                   
                  uris%Yd(:,nd) = d
C                   write(*,*)" disp is ", d

                  DEALLOCATE(xl, Nxi, N)
                  GOTO 120
               END IF
            
            END DO
            DEALLOCATE(xl, N, Nxi)
         END DO

120      CONTINUE
      END DO


!     Compute the projection from the disp of the element 

!     Use this function to GETXI(eType, eNoN, xl, xp, xi, flag)
      
      
      RETURN
      END SUBROUTINE URIS_UpdateDisp
!####################################################################
!--------------------------------------------------------------------
!     This subroutine check if a node is inside a tetrahedron  
      SUBROUTINE insideTet(eNoN,xp,xl,flag, ext)
      USE TYPEMOD
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=RKIND), INTENT(IN) :: xp(nsd), xl(nsd,eNoN)
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      LOGICAL, INTENT(IN) :: ext
      INTEGER(KIND=IKIND), INTENT(OUT) :: flag

      REAL(KIND=RKIND) :: minb(nsd), maxb(nsd)
      INTEGER(KIND=IKIND) :: i

!     Create a bounding box around of the current solid location 
      minb = HUGE(minb)
      maxb = TINY(maxb)

      flag = 0

      DO i=1, nsd
         minb(i) = MINVAL(xl(i,:) - 0.01)  
         maxb(i) = MAXVAL(xl(i,:) + 0.01)
      END DO

!     Is the node inside the BBox? 
      IF((xp(1) .LE. maxb(1)) 
     2             .AND. (xp(1) .GE. minb(1)) 
     3             .AND. (xp(2) .LE. maxb(2)) 
     4             .AND. (xp(2) .GE. minb(2))
     5             .AND. (xp(3) .LE. maxb(3)) 
     6             .AND. (xp(3) .GE. minb(3)) ) THEN 
!        The node is inside the Bounding Box
         flag = IN_POLY(xp, xl, ext)
!         CALL IN_POLY2(xp,xl, flag)
      END IF 

      RETURN
      END SUBROUTINE insideTet
!####################################################################
!--------------------------------------------------------------------
!     Read the URIS mesh separately 
      SUBROUTINE URIS_READMSH(list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(listType), INTENT(INOUT) :: list

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: i, j, iM, iFa, a, b, Ac, e
      REAL(KIND=RKIND) :: fibN(nsd), rtmp
      CHARACTER(LEN=stdL) :: ctmp, fExt
      TYPE(listType), POINTER :: lPtr, lPM
      TYPE(fileType) :: fTmp

      REAL(KIND=RKIND), ALLOCATABLE :: tmpX(:,:), gX(:,:)

      uris%nMsh  = list%srch("Add URIS mesh",ll=1)
      std = " Number of immersed surfaces for uris: "//uris%nMsh
      ALLOCATE (uris%msh(uris%nMsh), gX(0,0))

      uris%tnNo = 0
      DO iM=1, uris%nMsh
         lPM => list%get(uris%msh(iM)%name, "Add URIS mesh", iM)
         lPtr => lPM%get(uris%msh(iM)%lShl, "Set mesh as shell")
         IF (.NOT.uris%msh(iM)%lShl) err = "Only shells are allowed"

         std  = "Reading URIS mesh <"//CLR(TRIM(uris%msh(iM)%name))//">"
         CALL READSV(lPM, uris%msh(iM))
         IF (uris%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READCCNE(lPM, uris%msh(iM))
         END IF
         IF (uris%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READNRB(lPM, uris%msh(iM))
         END IF
         IF (uris%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READGAMBIT(lPM, uris%msh(iM))
         END IF
         IF (uris%msh(iM)%eType .EQ. eType_NA) THEN
            err = " Failed to identify format of the uris mesh"
         END IF

         std = " Number of uris nodes: "//uris%msh(iM)%gnNo
         std = " Number of uris elements: "//uris%msh(iM)%gnEl

!     Making sure face names are unique
         DO iFa=1, uris%msh(iM)%nFa
            uris%msh(iM)%fa(iFa)%iM = iM
            ctmp = uris%msh(iM)%fa(iFa)%name
            DO i=1, iM
               DO j=1, uris%msh(i)%nFa
                  IF (ctmp .EQ. uris%msh(i)%fa(j)%name .AND.
     2               (i.NE.iM .OR. j.NE.iFa)) THEN
                     err = " Repeating face names is not allowed"
                  END IF
               END DO
            END DO
         END DO

!     To scale the mesh, while attaching x to gX
         uris%msh(iM)%scF = 1._RKIND
         lPtr => lPM%get(uris%msh(iM)%scF,"Mesh scale factor",
     2                                                     lb=0._RKIND)
         a = uris%tnNo + uris%msh(iM)%gnNo
         IF (iM .GT. 1) THEN
            ALLOCATE(tmpX(nsd,uris%tnNo))
            tmpX = gX
            DEALLOCATE(gX)
            ALLOCATE(gX(nsd,a))
            gX(:,1:uris%tnNo) = tmpX
            DEALLOCATE(tmpX)
         ELSE
            DEALLOCATE(gX)
            ALLOCATE(gX(nsd,a))
         END IF
         gX(:,uris%tnNo+1:a) = uris%msh(iM)%x * uris%msh(iM)%scF
         uris%tnNo           = a
         DEALLOCATE(uris%msh(iM)%x)

C          lPtr => lPM%get(uris%msh(iM)%dx,"Mesh global edge size",1)
      END DO
      ALLOCATE(uris%x(nsd,uris%tnNo))
      ALLOCATE(uris%xCu(nsd,uris%tnNo))
C       ALLOCATE(uris%xCuo(nsd,uris%tnNo))
      uris%x = gX
      uris%xCu = gX
C       uris%xCuo = gX
      DEALLOCATE(gX)

!     Setting msh%gN, msh%lN parameter
      b = 0
      DO iM=1, uris%nMsh
         uris%msh(iM)%nNo = uris%msh(iM)%gnNo
         ALLOCATE(uris%msh(iM)%gN(uris%msh(iM)%nNo), 
     2                             uris%msh(iM)%lN(uris%tnNo))
         uris%msh(iM)%gN = 0
         uris%msh(iM)%lN = 0
         DO a=1, uris%msh(iM)%nNo
            b = b + 1
            uris%msh(iM)%gN(a) = b
            uris%msh(iM)%lN(b) = a
         END DO
      END DO
      IF (b .NE. uris%tnNo) err =
     2   " Mismatch in uris%tnNo. Correction needed"

!     Remap msh%gIEN array
      DO iM=1, uris%nMsh
         uris%msh(iM)%nEl = uris%msh(iM)%gnEl
         ALLOCATE(uris%msh(iM)%IEN(uris%msh(iM)%eNoN,uris%msh(iM)%nEl))
         DO e=1, uris%msh(iM)%nEl
            DO a=1, uris%msh(iM)%eNoN
               Ac = uris%msh(iM)%gIEN(a,e)
               Ac = uris%msh(iM)%gN(Ac)
               uris%msh(iM)%IEN(a,e) = Ac
            END DO
         END DO
         DEALLOCATE(uris%msh(iM)%gIEN)
      END DO

!     Re-arranging fa structure - %gN, %lN, %IEN
      b = 0
      DO iM=1, uris%nMsh
         DO iFa=1, uris%msh(iM)%nFa
            ALLOCATE(uris%msh(iM)%fa(iFa)%lN(uris%tnNo))
            uris%msh(iM)%fa(iFa)%lN = 0
            DO a=1, uris%msh(iM)%fa(iFa)%nNo
               Ac = uris%msh(iM)%fa(iFa)%gN(a)
               Ac = uris%msh(iM)%gN(Ac)
               uris%msh(iM)%fa(iFa)%gN(a)  = Ac
               uris%msh(iM)%fa(iFa)%lN(Ac) = a
            END DO
            DO e=1, uris%msh(iM)%fa(iFa)%nEl
               DO a=1, uris%msh(iM)%fa(iFa)%eNoN
                  Ac = uris%msh(iM)%fa(iFa)%IEN(a,e)
                  Ac = uris%msh(iM)%gN(Ac)
                  uris%msh(iM)%fa(iFa)%IEN(a,e) = Ac
               END DO
            END DO
         END DO
      END DO

!     Setting dmnId parameter here, if there is at least one mesh that
!     has defined eId.
C       DO iM=1, nMsh
C          lPM => list%get(uris%msh(iM)%name,"Add uris",iM)

C          lPtr => lPM%get(i,"Domain (uris)",ll=0,
C      2                                     ul=BIT_SIZE(uris%dmnId)-1)
C          IF (ASSOCIATED(lPtr)) CALL SETDMNID(uris%msh(iM), i)

C          lPtr => lPM%get(fTmp,"Domain (uris) file path")
C          IF (ASSOCIATED(lPtr)) THEN
C             i = LEN(TRIM(fTmp%fname))
C             fExt = fTmp%fname(i-2:i)
C             IF (TRIM(fExt).EQ."vtp" .OR. TRIM(fExt).EQ."vtu") THEN
C                CALL SETDMNIDVTK(uris%msh(iM), fTmp%fname, "DOMAIN_ID")
C             ELSE
C                CALL SETDMNIDFF(uris%msh(iM), fTmp%open())
C             END IF
C          END IF

C          IF (.NOT.ALLOCATED(uris%msh(iM)%eId)) err =
C      2      " Immersed bodies require domain ID parameter to be set"
C       END DO
C       ALLOCATE(uris%dmnId(uris%tnNo))
C       uris%dmnId = 0
C       DO iM=1, uris%nMsh
C          DO e=1, uris%msh(iM)%nEl
C             DO a=1, uris%msh(iM)%eNoN
C                Ac = uris%msh(iM)%IEN(a,e)
C                uris%dmnId(Ac) = IOR(uris%dmnId(Ac),uris%msh(iM)%eId(e))
C             END DO
C          END DO
C       END DO

!     Read fiber orientation
C       flag = .FALSE.
C       DO iM=1, uris%nMsh
C          lPM => list%get(uris%msh(iM)%name,"Add uris",iM)
C          j = lPM%srch("Fiber direction file path")
C          IF (j .EQ. 0) j = lPM%srch("Fiber direction")
C          IF (j .NE. 0) THEN
C             flag = .TRUE.
C             EXIT
C          END IF
C       END DO

C       IF (flag) THEN
C          DO iM=1, uris%nMsh
C             lPM => list%get(uris%msh(iM)%name,"Add uris",iM)

C             uris%msh(iM)%nFn = lPM%srch("Fiber direction file path")
C             j = uris%msh(iM)%nFn
C             IF (uris%msh(iM)%nFn .NE. 0) THEN
C                ALLOCATE(uris%msh(iM)%fN(j*nsd,uris%msh(iM)%nEl))
C                uris%msh(iM)%fN = 0._RKIND
C                DO i=1, uris%msh(iM)%nFn
C                   lPtr => lPM%get(cTmp, "Fiber direction file path", i)
C                   CALL READFIBNFF(uris%msh(iM), cTmp, "Furis_DIR", i)
C                END DO
C             ELSE
C                uris%msh(iM)%nFn = lPM%srch("Fiber direction")
C                j = uris%msh(iM)%nFn
C                IF (uris%msh(iM)%nFn .NE. 0) THEN
C                   ALLOCATE(uris%msh(iM)%fN(j*nsd,uris%msh(iM)%nEl))
C                   uris%msh(iM)%fN = 0._RKIND
C                   DO i=1, uris%msh(iM)%nFn
C                      lPtr => lPM%get(fibN, "Fiber direction", i)
C                      rtmp = SQRT(NORM(fibN))
C                      IF (.NOT.ISZERO(rtmp)) fibN(:) = fibN(:)/rtmp
C                      DO e=1, uris%msh(iM)%nEl
C                         uris%msh(iM)%fN((i-1)*nsd+1:i*nsd,e) = 
C      2                                                      fibN(1:nsd)
C                      END DO
C                   END DO
C                END IF
C             END IF
C          END DO
C       ELSE
C          uris%msh(:)%nFn = 0
C       END IF

      IF (uris%nMsh .GT. 1) THEN
         std = " Total number of uris nodes: "//uris%tnNo
         std = " Total number of uris elements: "//SUM(uris%msh%nEl)
      END IF

      std = CLR(" uris mesh data imported successfully",3)

      RETURN
      END SUBROUTINE uris_READMSH
!####################################################################
!     Write IB solution to a vtu file
      SUBROUTINE URIS_WRITEVTUS(lU)
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lU(nsd,uris%tnNo)

      TYPE(dataType) :: d(uris%nMsh)
      TYPE(vtkXMLType) :: vtu

      INTEGER(KIND=IKIND) :: iStat, iEq, iOut, iM, a, e, Ac, Ec, nNo,
     2   nEl, s, l, ie, is, nSh, oGrp, outDof, nOut, cOut
      CHARACTER(LEN=stdL) :: fName

      INTEGER(KIND=IKIND), ALLOCATABLE :: outS(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: tmpI(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: outNames(:)

      IF (cm%slv()) THEN
         uris%savedOnce = .TRUE.
         RETURN
      END IF

!     we plot coord + displ
      nOut   = 2
      outDof = 2*nsd
C       DO iEq=1, nEq
C          DO iOut=1, eq(iEq)%nOutIB
C             IF (.NOT.eq(iEq)%outIB(iOut)%wtn(1)) CYCLE
C             nOut   = nOut + 1
C             outDof = outDof + eq(iEq)%outIB(iOut)%l
C          END DO
C       END DO

      ALLOCATE(outNames(nOut), outS(nOut+1))

!     Prepare all solultions in to dataType d
      nNo = 0
      nEl = 0
      DO iM=1, uris%nMsh
         cOut           = 1
         outS(cOut)     = 1
         outS(cOut+1)   = nsd + 1
         outNames(cOut) = ""

         IF (uris%msh(iM)%eType .EQ. eType_NRB) err =
     2      " Outputs for NURBS data is under development"

         d(iM)%nNo     = uris%msh(iM)%nNo
         d(iM)%nEl     = uris%msh(iM)%nEl
         d(iM)%eNoN    = uris%msh(iM)%eNoN
         d(iM)%vtkType = uris%msh(iM)%vtkType

         ALLOCATE(d(iM)%x(outDof,d(iM)%nNo),
     2      d(iM)%IEN(d(iM)%eNoN,d(iM)%nEl))
         DO a=1, uris%msh(iM)%nNo
            Ac = uris%msh(iM)%gN(a)
            d(iM)%x(1:nsd,a) = uris%x(:,Ac)
         END DO

         DO e=1, uris%msh(iM)%nEl
            d(iM)%IEN(:,e) = uris%msh(iM)%IEN(:,e)
         END DO

C          DO iEq=1, nEq
C             DO iOut=1, 2 !eq(iEq)%nOutIB
C                IF (.NOT.eq(iEq)%outIB(iOut)%wtn(1)) CYCLE
         l  = nsd !eq(iEq)%outIB(iOut)%l
         s  = 1 !eq(iEq)%s + eq(iEq)%outIB(iOut)%o
         e  = s + l - 1

         cOut = cOut + 1
         is   = outS(cOut)
         ie   = is + l - 1
         outS(cOut+1)   = ie + 1
C          outNames(1) = "URIS_"//TRIM(eq(iEq)%outIB(iOut)%name)
         outNames(1) = "coordinates"
         outNames(2) = "URIS_displacement"

C                oGrp = eq(iEq)%outIB(iOut)%grp
C                SELECT CASE (oGrp)
C                   CASE (outGrp_NA)
C                   err = "Undefined output grp in VTK"
C                CASE (outGrp_Y)
C                   DO a=1, uris%msh(iM)%nNo
C                      Ac = uris%msh(iM)%gN(a)
C                      d(iM)%x(is:ie,a) = lY(s:e,Ac)
C                   END DO
C                CASE (outGrp_D)
         DO a=1, uris%msh(iM)%nNo
            Ac = uris%msh(iM)%gN(a)
            d(iM)%x(is:ie,a) = lU(s:e,Ac)
         END DO
C                CASE DEFAULT
C                   err = "Undefined output "//
C      2               TRIM(eq(iEq)%outIB(iOut)%name)
C C                END SELECT
C             END DO
C          END DO

         ALLOCATE(d(iM)%xe(d(iM)%nEl,1))
         IF (.NOT.savedOnce) THEN
            IF (ALLOCATED(uris%dmnID)) THEN
               d(iM)%xe(:,1) = REAL(uris%msh(iM)%eId(:), KIND=RKIND)
            ELSE
               d(iM)%xe(:,1) = 1._RKIND
            END IF
         END  IF
         nNo = nNo +  d(iM)%nNo
         nEl = nEl +  d(iM)%nEl
      END DO

      ALLOCATE(tmpV(maxnsd,nNo))

C       write(*,*)" d%x ", d(1)%x
C       write(*,*)" d%x ", SIZE(d(1)%x,1)
C       write(*,*)" d%x ", SIZE(d(1)%x,2)


!     Writing to vtu file (master only)
      IF (cTS .GE. 1000) THEN
         fName = STR(cTS)
      ELSE
         WRITE(fName,'(I3.3)') cTS
      END IF

      fName = TRIM(saveName)//"_uris_"//TRIM(ADJUSTL(fName))//".vtu"
      dbg = "Writing VTU"

      CALL vtkInitWriter(vtu, TRIM(fName), iStat)
      IF (iStat .LT. 0) err = "VTU file write error (init)"

!     Writing the position data
      iOut = 1
      s    = outS(iOut)
      e    = outS(iOut+1)-1
      nSh  = 0
      tmpV = 0._RKIND
      DO iM=1, uris%nMsh
         DO a=1, d(iM)%nNo
            tmpV(1:nsd,a+nSh) = d(iM)%x(s:e,a)
         END DO
         nSh = nSh + d(iM)%nNo
      END DO
      CALL putVTK_pointCoords(vtu, tmpV(1:nsd,:), iStat)
      IF (iStat .LT. 0) err = "VTU file write error (coords)"

!     Writing the connectivity data
      nSh = -1
      DO iM=1, uris%nMsh
         ALLOCATE(tmpI(d(iM)%eNoN,d(iM)%nEl))
         DO e=1, d(iM)%nEl
            tmpI(:,e) = d(iM)%IEN(:,e) + nSh
         END DO
         CALL putVTK_elemIEN(vtu, tmpI, d(iM)%vtkType, iStat)
         IF (iStat .LT. 0) err = "VTU file write error (ien)"
         DEALLOCATE(tmpI)
         nSh = nSh + d(iM)%nNo
      END DO

C       write(*,*)" outS ",  outS

!     Writing all solutions
      DO iOut=1, nOut
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         s = outS(iOut)
         e = outS(iOut+1) - 1
         l = e - s + 1

         ALLOCATE(tmpV(l, nNo))
         nSh = 0
         DO iM=1, uris%nMsh
            DO a=1, d(iM)%nNo
               tmpV(:,a+nSh) = d(iM)%x(s:e,a)
            END DO
            nSh = nSh + d(iM)%nNo
         END DO
         CALL putVTK_pointData(vtu, outNames(iOut), tmpV, iStat)
         IF (iStat .LT. 0) err = "VTU file write error (point data)"
      END DO

!     Write element-based variables
C       IF (.NOT.savedOnce .OR. mvMsh) THEN
C          uris%savedOnce = .TRUE.
C          ALLOCATE(tmpI(1,nEl))
C !     Write the domain ID
C          IF (ALLOCATED(uris%dmnID)) THEN
C             Ec = 0
C             DO iM=1, uris%nMsh
C                DO e=1, d(iM)%nEl
C                   Ec = Ec + 1
C                   tmpI(1,Ec) = INT(d(iM)%xe(e,1), KIND=IKIND)
C                END DO
C             END DO
C             CALL putVTK_elemData(vtu, 'Domain_ID', tmpI, iStat)
C             IF (iStat .LT. 0) err = "VTU file write error (dom id)"
C          END IF

C !     Write the mesh ID
C          IF (uris%nMsh .GT. 1) THEN
C             Ec = 0
C             DO iM=1, uris%nMsh
C                DO e=1, d(iM)%nEl
C                   Ec = Ec + 1
C                   tmpI(1,Ec) = iM
C                END DO
C             END DO
C             CALL putVTK_elemData(vtu, 'Mesh_ID', tmpI, iStat)
C             IF (iStat .LT. 0) err = "VTU file write error (mesh id)"
C          END IF
C          DEALLOCATE(tmpI)
C       END IF

      DO iM=1, uris%nMsh
         CALL DESTROY(d(iM))
      END DO

      CALL vtkWriteToFile(vtu, iStat)
      IF (iStat .LT. 0) err = "VTU file write error"

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE URIS_WRITEVTUS
!####################################################################
!--------------------------------------------------------------------

