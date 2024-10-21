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
      SUBROUTINE URIS_MEANP(iUris)
      USE TYPEMOD
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iUris
      INTEGER(KIND=IKIND) :: iEq, nPrj, m, s, e, i, iM, iFa
      REAL(KIND=RKIND) :: tmp
      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:)

!     Let's conpute the mean pressure in the two regions of the fluid mesh 
!     For the moment let's define a flag IdSubDmn(size the number of elements)
      INTEGER(KIND=IKIND) :: a, Ac, cntM, cntP
      REAL(KIND=RKIND) :: volU, volD, Deps, meanPU, meanPD, 
     2                    sDST(1,tnNo), sUPS(1,tnNo), sImm(1,tnNo)


!     Now we can compute the pressure mean on each subdomain
      iEq = 1

      IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
      ALLOCATE (tmpV(maxnsd,tnNo))
      ! FK: What if there are multiple meshes??
      ! TO-DO: We need to have a sdf array for each mesh
      Deps = uris(iUris)%sdf_deps
      volU = 0._RKIND
      volD = 0._RKIND
!     Let's compute left side 
      sUPS = 0._RKIND
      WHERE(uris(iUris)%sdf.GE.0.AND.uris(iUris)%sdf.LE.Deps) 
     2     sUPS(1,:) = 1._RKIND
      DO iM=1, nMsh
         volU = volU + Integ(iM,sUPS)
      END DO

!     Let's compute right side 
      sDST = 0._RKIND
      WHERE(uris(iUris)%sdf.LT.0.AND.uris(iUris)%sdf.GE.-Deps) 
     2    sDST(1,:) = 1._RKIND
      DO iM=1, nMsh
         volD = volD + Integ(iM,sDST)
      END DO
      write(*,*)" volume upstream ", volU
      write(*,*)" volume downstream ", volD

      meanPU = 0._RKIND
      meanPD = 0._RKIND

      m = 1
      s = eq(iEq)%s + nsd 
      e = s + m - 1

      tmpV(1:m,:) = Yn(s:e,:)*sUPS
      DO iM = 1, nMsh
          meanPU = meanPU + Integ(iM,tmpV)
      END DO
      meanPU = meanPU/volU

      tmpV(1:m,:) = Yn(s:e,:)*sDST
      DO iM = 1, nMsh
          meanPD = meanPD + Integ(iM,tmpV)
      END DO
      meanPD = meanPD/volD

      write(*,*)" mean P upstream ", meanPU
      write(*,*)" mean P downstream ", meanPD

!     If the uris has passed the closing state
      IF (uris(iUris)%cnt.GT.SIZE(uris(iUris)%DxClose,1)) THEN
          IF( (meanPD .GT. meanPU)) THEN
              uris(iUris)%cnt = 1
              uris(iUris)%clsFlg = .FALSE.
              urisActFlag = .TRUE.
              write(*,*) "Set urisOpenFlag to TRUE for uris ", iUris
          END IF
      END IF
      IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)

      RETURN
      END SUBROUTINE URIS_MEANP
!####################################################################
!--------------------------------------------------------------------
!     This subroutine computes the mean velocity in the fluid elements near the 
!     immersed surface  
      SUBROUTINE URIS_MEANV(iUris)
      USE TYPEMOD
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND), INTENT(IN) :: iUris
      INTEGER(KIND=IKIND) :: iEq, nPrj, m, s, e, i, iM, iFa
      REAL(KIND=RKIND) :: tmp
      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:)

!     Let's conpute the mean pressure in the two regions of the fluid mesh 
!     For the moment let's define a flag IdSubDmn(size the number of elements)
      INTEGER(KIND=IKIND) :: a, Ac, cntM, cntP
      REAL(KIND=RKIND) :: volI, Deps, zSurf, meanV, 
     2                    sImm(1,tnNo), tmpVNrm(1,tnNo)

!     Let's compute the neighboring region below the valve normal. When
!     the valve is open, this region should roughly be valve oriface.
      iEq = 1

      IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
      ALLOCATE (tmpV(maxnsd,tnNo))
     
      Deps = uris(iUris)%sdf_deps
      sImm = 0._RKIND
      volI = 0._RKIND
      WHERE(uris(iUris)%sdf.LE.-Deps) sImm(1,:) = 1._RKIND
      DO iM=1, nMsh
        volI = volI + Integ(iM,sImm)
      END DO
      write(*,*)" volume inside ", volI

      m = nsd
      s = eq(iEq)%s 
      e = s + m - 1

      DO i=1, nsd
        tmpV(i,:) = Yn(s+i-1,:) * sImm(1,:)
      END DO
      DO i = 1, tnNo
        !tmpVNrm(1,i) = NORM(tmpV(1:m,i),uris(iUris)%nrm)
        tmpVNrm(1,i) = tmpV(1,i)*uris(iUris)%nrm(1) + 
     2          tmpV(2,i)*uris(iUris)%nrm(2) + 
     3          tmpV(3,i)*uris(iUris)%nrm(3)
      END DO
     
      meanV = 0._RKIND
      DO iM=1, nMsh 
        meanV = meanV + Integ(iM,tmpVNrm)/volI
        write(*,*) "!!!!WHY???", iM, Integ(iM,tmpVNrm)
      END DO
      write(*,*)" mean Vel ", meanV

!     If the uris has passed the open state
      IF (uris(iUris)%cnt.GT.SIZE(uris(iUris)%DxOpen,1)) THEN
          IF (meanV.LE.0._RKIND) THEN
              uris(iUris)%cnt = 1
              urisActFlag = .TRUE.
              uris(iUris)%clsFlg = .TRUE.
              write(*,*) "Set urisCloseFlag to TRUE"
          END IF
      END IF
      write(*,*) "urisCloseFlag is ", uris(iUris)%clsFlg

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

      INTEGER(KIND=IKIND) :: flag, iM, jM, iEln, a, nd, Ac, iUris
      REAL(KIND=RKIND) :: xp(nsd), xi(nsd), d(nsd)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), N(:), Nxi(:,:)
      LOGICAL :: ultra, fl

!     For each point in the immersed surface we need to localize it 
!     = find the fluid element that contains the node
!     Since the fluid element could be on another processor, we need to
!     gather the displacement values at the end      
!     FK: it's probably better to save the element ids so that we don't
!     have to run the search every time step, only during open or close      
      DO iUris=1, nUris
         DO iM=1, uris(iUris)%nFa
           flag = 1
           ultra = .FALSE.
           xi = 0.5_RKIND
           
           nd = 0
           DO WHILE(nd .LE. uris(iUris)%tnNo)
!             Check if we were able to find the tetra.
!             FK: if not, the tetra is on another processor 
              IF(flag .EQ. 0) THEN 
                 IF( ultra ) THEN 
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

              xp = uris(iUris)%x(:,nd) !+ uris(iUris)%Yd(:,nd)
123           CONTINUE
              flag = 1
              
              DO jM=1, nMsh
                 ALLOCATE(xl(nsd,msh(jM)%eNoN), N(msh(jM)%eNoN), 
     2                                           Nxi(nsd,msh(jM)%eNoN))
                 DO iEln = 1, msh(jM)%nEl
                    DO a=1, msh(jM)%eNoN
                       Ac = msh(jM)%IEN(a,iEln)
                          xl(:,a) = x(:,Ac) 
                    END DO 
                    CALL insideTet(msh(jM)%eNoN,xp,xl,flag, ultra)

                    IF( flag .EQ. 1) THEN 
!                      Get displacement  
!                      Localize p inside the parent element  
                       CALL GETXI(msh(jM)%eType,msh(jM)%eNoN, xl, 
     2                         xp, xi,fl)  
                       IF( .NOT.fl) write(*,*)" GETXI not converging "
!                      evaluate N at xi 
                       CALL GETGNN(nsd,msh(jM)%eType,msh(jM)%eNoN,xi,
     2                      N,Nxi)
!                      use this to compute disp al node xp 
                       d = 0._RKIND
                       DO a=1, msh(jM)%eNoN
                          Ac = msh(jM)%IEN(a,iEln)
!                         We have to use Do because Dn contains the result 
!                         coming from the solid 
                          d(1) = d(1) - N(a)*Dn(nsd+2,Ac) 
                          d(2) = d(2) - N(a)*Dn(nsd+3,Ac) 
                          d(3) = d(3) - N(a)*Dn(nsd+4,Ac) 
                       END DO
!                      update uris disp                                                   
                       uris(iUris)%Yd(:,nd) = d
                       DEALLOCATE(xl, Nxi, N)
                       GOTO 120
                    END IF
                 END DO
                 DEALLOCATE(xl, N, Nxi)
              END DO

120           CONTINUE
           END DO
         END DO
      END DO
      
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
    
      ! FK: Hard coded BBox?? This is going to cause problem if scale
      ! changes
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
         !FK: What if the node is inside the BBox but not poly??
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
      INTEGER(KIND=IKIND) :: fid, dispNtClose,dispNtOpen,dispNn,
     2  t, ioStatus, iUris
      REAL(KIND=RKIND) :: fibN(nsd), rtmp
      CHARACTER(LEN=stdL) :: ctmp, fExt
      TYPE(listType), POINTER :: lPtr, lPM, lPN
      TYPE(fileType) :: fTmp

      REAL(KIND=RKIND), ALLOCATABLE :: tmpX(:,:), gX(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: tmpD(:,:,:),dispClose(:,:,:), 
     2   dispOpen(:,:,:)

      nUris  = list%srch("Add URIS mesh",ll=1)
      std = " Number of immersed surfaces for uris: "//nUris
      ALLOCATE(uris(nUris))
      DO iUris=1, nUris
         lPM => list%get(uris(iUris)%name, "Add URIS mesh", iUris)
         std  = "Reading URIS mesh <"//CLR(TRIM(uris(iUris)%name))//">"
         uris(iUris)%scF = 1._RKIND
         lPtr => lPM%get(uris(iUris)%scF,"Mesh scale factor",
     2          lb=0._RKIND)

         uris(iUris)%nFa = lPM%srch("Add URIS face")
         ALLOCATE (uris(iUris)%msh(uris(iUris)%nFa), gX(0,0), 
     2      uris(iUris)%nrm(nsd))

         lPtr => lPM%get(fTmp, "Positive flow normal file")
         fid = fTmp%open()
         READ (fid,*) uris(iUris)%nrm(:)
         CLOSE (fid)

         lPtr => lPM%get(uris(iUris)%sdf_deps, "Thickness")
         uris(iUris)%tnNo = 0
         DO iM=1, uris(iUris)%nFa
            ! Set as shell
            uris(iUris)%msh(iM)%lShl = .TRUE.
            lPN => lPM%get(uris(iUris)%msh(iM)%name, "Add URIS face",
     2          iM)
            std = "Reading URIS face <"
     2              //CLR(TRIM(uris(iUris)%msh(iM)%name))//">"
            CALL READSV(lPN, uris(iUris)%msh(iM))
            IF (uris(iUris)%msh(iM)%eType .EQ. eType_NA) THEN
               CALL READCCNE(lPN, uris(iUris)%msh(iM))
            END IF
            IF (uris(iUris)%msh(iM)%eType .EQ. eType_NA) THEN
               CALL READNRB(lPN, uris(iUris)%msh(iM))
            END IF
            IF (uris(iUris)%msh(iM)%eType .EQ. eType_NA) THEN
               CALL READGAMBIT(lPN, uris(iUris)%msh(iM))
            END IF
            IF (uris(iUris)%msh(iM)%eType .EQ. eType_NA) THEN
               err = " Failed to identify format of the uris mesh"
            END IF

            std = " Number of uris nodes: "//uris(iUris)%msh(iM)%gnNo
            std = " Number of uris elements: "//uris(iUris)%msh(iM)%gnEl

!        Read valve motion: note that this motion is defined on the
!        reference configuration         
            lPtr => lPN%get(fTmp, "Open motion file path")
            fid = fTmp%open()
            READ (fid,*) dispNtOpen, dispNn
            IF (dispNn .NE. uris(iUris)%msh(iM)%gnNo) THEN
                write(*,*) dispNn, uris(iUris)%msh(iM)%gnNo
                err = "Mismatch in node numbers between URIS mesh and
     2                displacements."
            END IF
            ALLOCATE(dispOpen(dispNtOpen, nsd, dispNn))
            DO t=1, dispNtOpen
               DO a=1, dispNn
                   READ (fid,*, IOSTAT=ioStatus) dispOpen(t,:,a)
               END DO
            END DO
            CLOSE (fid)
            lPtr => lPN%get(fTmp, "Close motion file path")
            fid = fTmp%open()
            READ (fid,*) dispNtClose, dispNn
            IF (dispNn .NE. uris(iUris)%msh(iM)%gnNo) THEN
                write(*,*) dispNn, uris(iUris)%msh(iM)%gnNo
                err = "Mismatch in node numbers between URIS mesh and
     2                displacements."
            END IF
            ALLOCATE(dispClose(dispNtClose, nsd, dispNn))
            DO t=1, dispNtClose
               DO a=1, dispNn
                   READ (fid,*, IOSTAT=ioStatus) dispClose(t,:,a)
               END DO
            END DO
            CLOSE (fid)
!        To scale the mesh, while attaching x to gX
            a = uris(iUris)%tnNo + uris(iUris)%msh(iM)%gnNo
            IF (iM .GT. 1) THEN
               ALLOCATE(tmpX(nsd,uris(iUris)%tnNo))
               tmpX = gX
               DEALLOCATE(gX)
               ALLOCATE(gX(nsd,a))
               gX(:,1:uris(iUris)%tnNo) = tmpX
               DEALLOCATE(tmpX)
               ! Move data for open
               ALLOCATE(tmpD(dispNtOpen,nsd,uris(iUris)%tnNo))
               tmpD = uris(iUris)%DxOpen
               DEALLOCATE(uris(iUris)%DxOpen)
               ALLOCATE(uris(iUris)%DxOpen(dispNtOpen,nsd,a))
               uris(iUris)%DxOpen(:,:,1:uris(iUris)%tnNo)=tmpD
               DEALLOCATE(tmpD)
               ! Move data for close
               ALLOCATE(tmpD(dispNtClose,nsd,uris(iUris)%tnNo))
               tmpD = uris(iUris)%DxClose
               DEALLOCATE(uris(iUris)%DxClose)
               ALLOCATE(uris(iUris)%DxClose(dispNtClose,nsd,a))
               uris(iUris)%DxClose(:,:,1:uris(iUris)%tnNo)=tmpD
               DEALLOCATE(tmpD)
            ELSE
               DEALLOCATE(gX)
               ALLOCATE(gX(nsd,a))
               ALLOCATE(uris(iUris)%DxOpen(dispNtOpen,nsd,a))
               ALLOCATE(uris(iUris)%DxClose(dispNtClose,nsd,a))
            END IF
            gX(:,uris(iUris)%tnNo+1:a) = uris(iUris)%msh(iM)%x * 
     2          uris(iUris)%scF
            uris(iUris)%DxOpen(:,:,uris(iUris)%tnNo+1:a)=dispOpen
            uris(iUris)%DxClose(:,:,uris(iUris)%tnNo+1:a)=dispClose
            uris(iUris)%tnNo           = a
            DEALLOCATE(uris(iUris)%msh(iM)%x)
            DEALLOCATE(dispOpen, dispClose)
         END DO
         ALLOCATE(uris(iUris)%x(nsd,uris(iUris)%tnNo))
         uris(iUris)%x = gX
         ALLOCATE(uris(iUris)%Yd(nsd,uris(iUris)%tnNo))
         uris(iUris)%Yd = 0._RKIND
         DEALLOCATE(gX)

!        Setting msh%gN, msh%lN parameter
         b = 0
         DO iM=1, uris(iUris)%nFa
            uris(iUris)%msh(iM)%nNo = uris(iUris)%msh(iM)%gnNo
            ALLOCATE(uris(iUris)%msh(iM)%gN(uris(iUris)%msh(iM)%nNo), 
     2                        uris(iUris)%msh(iM)%lN(uris(iUris)%tnNo))
            uris(iUris)%msh(iM)%gN = 0
            uris(iUris)%msh(iM)%lN = 0
            DO a=1, uris(iUris)%msh(iM)%nNo
               b = b + 1
               uris(iUris)%msh(iM)%gN(a) = b
               uris(iUris)%msh(iM)%lN(b) = a
            END DO
         END DO
         IF (b .NE. uris(iUris)%tnNo) err =
     2      " Mismatch in uris(iUris)%tnNo. Correction needed"

!        Remap msh%gIEN array
         DO iM=1, uris(iUris)%nFa
            uris(iUris)%msh(iM)%nEl = uris(iUris)%msh(iM)%gnEl
            ALLOCATE(uris(iUris)%msh(iM)%IEN(uris(iUris)%msh(iM)%eNoN,
     2         uris(iUris)%msh(iM)%nEl))
            DO e=1, uris(iUris)%msh(iM)%nEl
               DO a=1, uris(iUris)%msh(iM)%eNoN
                  Ac = uris(iUris)%msh(iM)%gIEN(a,e)
                  Ac = uris(iUris)%msh(iM)%gN(Ac)
                  uris(iUris)%msh(iM)%IEN(a,e) = Ac
               END DO
            END DO
            DEALLOCATE(uris(iUris)%msh(iM)%gIEN)
         END DO

         IF (uris(iUris)%nFa .GT. 1) THEN
            std = " Total number of uris nodes: "//uris(iUris)%tnNo
            std = " Total number of uris elements: "//
     2              SUM(uris(iUris)%msh%nEl)
         END IF
      END DO
      std = CLR(" uris mesh data imported successfully",3)

      RETURN
      END SUBROUTINE uris_READMSH
!####################################################################
!     Write IB solution to a vtu file
      SUBROUTINE URIS_WRITEVTUS
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE

      TYPE(vtkXMLType) :: vtu

      INTEGER(KIND=IKIND) :: iStat, iEq, iOut, iM, a, e, Ac, Ec, nNo,
     2   nEl, s, l, ie, is, nSh, oGrp, outDof, nOut, cOut, iUris
      CHARACTER(LEN=stdL) :: fName

      INTEGER(KIND=IKIND), ALLOCATABLE :: outS(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: tmpI(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: outNames(:)
      TYPE(dataType), ALLOCATABLE:: d(:)

!     we plot coord + displ
      nOut   = 2
      outDof = nOut * nsd
      ALLOCATE(outNames(nOut), outS(nOut+1))
!     Prepare all solultions in to dataType d
      DO iUris=1, nUris
         ALLOCATE(d(uris(iUris)%nFa))
         nNo = 0
         nEl = 0
         DO iM=1, uris(iUris)%nFa
            cOut           = 1
            outS(cOut)     = 1
            outS(cOut+1)   = nsd + 1
            outNames(cOut) = ""

            IF (uris(iUris)%msh(iM)%eType .EQ. eType_NRB) err =
     2         " Outputs for NURBS data is under development"

            d(iM)%nNo     = uris(iUris)%msh(iM)%nNo
            d(iM)%nEl     = uris(iUris)%msh(iM)%nEl
            d(iM)%eNoN    = uris(iUris)%msh(iM)%eNoN
            d(iM)%vtkType = uris(iUris)%msh(iM)%vtkType
        
            !IF (ALLOCATED(d(iM)%x)) DEALLOCATE(d(iM)%x) 
            !IF (ALLOCATED(d(iM)%IEN)) DEALLOCATE(d(iM)%IEN) 
            ALLOCATE(d(iM)%x(outDof,d(iM)%nNo))
            ALLOCATE(d(iM)%IEN(d(iM)%eNoN,d(iM)%nEl))
            DO a=1, uris(iUris)%msh(iM)%nNo
               Ac = uris(iUris)%msh(iM)%gN(a)
               d(iM)%x(1:nsd,a) = uris(iUris)%x(:,Ac)
            END DO
            DO e=1, uris(iUris)%msh(iM)%nEl
               d(iM)%IEN(:,e) = uris(iUris)%msh(iM)%IEN(:,e)
            END DO

            l  = nsd !eq(iEq)%outIB(iOut)%l
            s  = 1 !eq(iEq)%s + eq(iEq)%outIB(iOut)%o
            e  = s + l - 1

            cOut = cOut + 1
            is   = outS(cOut)
            ie   = is + l - 1
            outS(cOut+1)   = ie + 1
            outNames(1) = "coordinates"
            outNames(2) = "URIS_displacement"

            DO a=1, uris(iUris)%msh(iM)%nNo
               Ac = uris(iUris)%msh(iM)%gN(a)
               d(iM)%x(is:ie,a) = uris(iUris)%Yd(s:e,Ac)
            END DO

            nNo = nNo +  d(iM)%nNo
            nEl = nEl +  d(iM)%nEl
         END DO

         ALLOCATE(tmpV(maxnsd,nNo))

!        Writing to vtu file (master only)
         IF (cTS .GE. 1000) THEN
            fName = STR(cTS)
         ELSE
            WRITE(fName,'(I3.3)') cTS
         END IF

         fName = TRIM(saveName)//"_uris_"// TRIM(uris(iUris)%name)//
     2          "_"//TRIM(ADJUSTL(fName))//".vtu"
         dbg = "Writing URIS VTU"
         CALL vtkInitWriter(vtu, TRIM(fName), iStat)
         IF (iStat .LT. 0) err = "VTU file write error (init)"

!        Writing the position data
         iOut = 1
         s    = outS(iOut)
         e    = outS(iOut+1)-1
         nSh  = 0
         tmpV = 0._RKIND
         DO iM=1, uris(iUris)%nFa
            DO a=1, d(iM)%nNo
               tmpV(1:nsd,a+nSh) = d(iM)%x(s:e,a)
            END DO
            nSh = nSh + d(iM)%nNo
         END DO
         CALL putVTK_pointCoords(vtu, tmpV(1:nsd,:), iStat)
         IF (iStat .LT. 0) err = "VTU file write error (coords)"

!        Writing the connectivity data
         DO iM=1, uris(iUris)%nFa
            ALLOCATE(tmpI(d(iM)%eNoN,d(iM)%nEl))
            DO e=1, d(iM)%nEl
               tmpI(:,e) = d(iM)%IEN(:,e) - 1
            END DO
            CALL putVTK_elemIEN(vtu, tmpI, d(iM)%vtkType, iStat)
            IF (iStat .LT. 0) err = "VTU file write error (ien)"
            DEALLOCATE(tmpI)
         END DO

!        Writing all solutions
         DO iOut=1, nOut
            IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
            s = outS(iOut)
            e = outS(iOut+1) - 1
            l = e - s + 1

            ALLOCATE(tmpV(l, nNo))
            nSh = 0
            DO iM=1, uris(iUris)%nFa
               DO a=1, d(iM)%nNo
                  tmpV(:,a+nSh) = d(iM)%x(s:e,a)
               END DO
               nSh = nSh + d(iM)%nNo
            END DO
            CALL putVTK_pointData(vtu, outNames(iOut), tmpV, iStat)
            IF (iStat .LT. 0) err = "VTU file write error (point data)"
         END DO

         DO iM=1, uris(iUris)%nFa
            CALL DESTROY(d(iM))
         END DO

         CALL vtkWriteToFile(vtu, iStat)
         IF (iStat .LT. 0) err = "VTU file write error"

         CALL flushVTK(vtu)
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         IF (ALLOCATED(d)) DEALLOCATE(d)
      END DO
      RETURN
      END SUBROUTINE URIS_WRITEVTUS
!####################################################################
!--------------------------------------------------------------------
!     Checks if a probe lies inside or outside an immersed boundary
      SUBROUTINE URIS_CALCSDF
      USE COMMOD
      USE ALLFUN
      
      IMPLICIT NONE
      REAL(KIND=RKIND) :: xp(nsd)
      REAL(KIND=RKIND) :: ALLOCATABLE 
      LOGICAL :: flag

      INTEGER(KIND=IKIND) :: i, ca, a, e, Ac, Ec, iM, jM,iUris, cnt
      REAL(KIND=RKIND) :: dS, minS, Jac, nV(nsd), xb(nsd), dotP
      REAL(KIND=RKIND), ALLOCATABLE :: lX(:,:), xXi(:,:)

      REAL(KIND=RKIND) :: minb(nsd), maxb(nsd), extra(nsd)
      ALLOCATE(xXi(nsd, nsd-1))

      DO iUris=1, nUris
        ! We need to check if the valve needs to move 
        ! This cut off threshold needs to be specified from the input
        ! file!
        IF (.NOT.uris(iUris)%clsFlg) THEN
            cnt = MIN(uris(iUris)%cnt, SIZE(uris(iUris)%DxOpen,1))
            uris(iUris)%x = uris(iUris)%DxOpen(cnt, :, :)
        ELSE
            cnt = MIN(uris(iUris)%cnt, SIZE(uris(iUris)%DxClose,1))
            uris(iUris)%x = uris(iUris)%DxClose(cnt, :, :)
        END IF 
        write(*,*) "C CNT: ", cnt, uris(iUris)%cnt, uris(iUris)%clsFlg
        IF (ALLOCATED(uris(iUris)%sdf).AND.cnt.LT.uris(iUris)%cnt) CYCLE

        ALLOCATE(lX(nsd, uris(iUris)%msh(1)%eNoN))
        IF (.NOT. ALLOCATED(uris(iUris)%sdf)) THEN
            ALLOCATE(uris(iUris)%sdf(tnNo))
        END IF
        write(*,*) "!!!RECOMPUTING SDF for", iUris
        uris(iUris)%sdf = uris(iUris)%sdf_default
    
        ! FK:
        ! Each time when the URIS moves, we need to recompute the signed
        ! distance function.
        ! DO WE NEED TO RECOMPUTE WHEN THE MESH MOVES? Ideally not

        ! Now we assume that for each uris, we only have one valve
        ! We will need to work on generalization later
        ! Find the bounding box of the valve, the BBox will be 10% larger
        ! than the actual valve.
        minb = HUGE(minb)
        maxb = TINY(maxb)
        DO i=1, nsd
           minb(i) = MINVAL(uris(iUris)%x(i,:))  
           maxb(i) = MAXVAL(uris(iUris)%x(i,:))
           extra(i) = (maxb(i) - minb(i)) * 0.1
        END DO

        ! Should this be computed on the reference or current
        ! configuration?
        DO ca=1, tnNo
          minS = HUGE(minS)
          xp = x(:, ca) + Do(nsd+2:2*nsd+1,ca)
!         Is the node inside the BBox? 
          IF (ALL(xp.GE.minb-extra).AND.ALL(xp.LE.maxb+extra)) THEN
              ! This point is inside the BBox
              ! Find the closest URIS face centroid
              DO iM=1, uris(iUris)%nFa
                  ! Here we are using original coordinates, plus 
                  ! the displacements
                  DO e=1, uris(iUris)%msh(iM)%nEl
                      xb = 0._RKIND
                      DO a=1, uris(iUris)%msh(iM)%eNoN
                          Ac = uris(iUris)%msh(iM)%IEN(a,e)
                          xb = xb + uris(iUris)%x(:,Ac) + 
     2                         uris(iUris)%Yd(:,Ac) 
                      END DO
                      xb = xb/REAL(uris(iUris)%msh(iM)%eNoN, KIND=RKIND)
                      dS = SQRT( SUM( (xp(:)-xb(:))**2._RKIND ) )
                      IF (dS .LT. minS) THEN
                          minS = dS
                          Ec = e
                          jM = iM
                      END IF
                  END DO
              END DO
              ! We also need to compute the sign (above or below
              ! the valve).
              ! Compute the element normal
              xXi = 0._RKIND
              lX = 0._RKIND
              xb = 0._RKIND
              DO a=1, uris(iUris)%msh(jM)%eNoN
                 Ac = uris(iUris)%msh(jM)%IEN(a,Ec)
                 xb = xb + uris(iUris)%x(:,Ac) + uris(iUris)%Yd(:,Ac)
                 lX(:,a) = uris(iUris)%x(:,Ac) + uris(iUris)%Yd(:,Ac)
              END DO
              xb   = xb / REAL(uris(iUris)%msh(jM)%eNoN, KIND=RKIND)
              DO a = 1, uris(iUris)%msh(jM)%eNoN
                  DO i = 1, nsd-1
                      xXi(:,i) = xXi(:,i) + 
     2                      lX(:,a)*uris(iUris)%msh(jM)%Nx(i,a,1)
                  END DO
              END DO
              nV(:) = CROSS(xXi)
              nV(:) = nV(:) / SQRT(NORM(nV))
              dotP = NORM(xp-xb, nV)
    
              uris(iUris)%sdf(ca) = dotP
 
          END IF
        END DO
        write(*,*) "any nan in sdf? ", ANY(ISNAN(uris(iUris)%sdf))
        IF (ALLOCATED(lX)) DEALLOCATE(lX)
      END DO 
      RETURN
      END SUBROUTINE URIS_CALCSDF
