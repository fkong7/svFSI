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
!     Functions for FSI with immersed structures with IFEM
!
!--------------------------------------------------------------------

!     This routine reads IFEM mesh data
      SUBROUTINE IFEM_READMSH(list)
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

      ifem%nMsh  = list%srch("Add IFEM",ll=1)
      std = " Number of immersed boundaries: "//ifem%nMsh
      ALLOCATE (ifem%msh(ifem%nMsh), gX(0,0))

      ifem%tnNo = 0
      DO iM=1, ifem%nMsh
         lPM => list%get(ifem%msh(iM)%name, "Add IFEM", iM)
         lPtr => lPM%get(ifem%msh(iM)%lShl, "Set mesh as shell")
         IF (ifem%msh(iM)%lShl) err = "Immersed shells are not allowed"

         std  = "Reading IFEM mesh <"//CLR(TRIM(ifem%msh(iM)%name))//">"
         CALL READSV(lPM, ifem%msh(iM))
         IF (ifem%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READCCNE(lPM, ifem%msh(iM))
         END IF
         IF (ifem%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READNRB(lPM, ifem%msh(iM))
         END IF
         IF (ifem%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READGAMBIT(lPM, ifem%msh(iM))
         END IF
         IF (ifem%msh(iM)%eType .EQ. eType_NA) THEN
            err = " Failed to identify format of the IFEM mesh"
         END IF

         std = " Number of IFEM nodes: "//ifem%msh(iM)%gnNo
         std = " Number of IFEM elements: "//ifem%msh(iM)%gnEl

!     Making sure face names are unique
         DO iFa=1, ifem%msh(iM)%nFa
            ifem%msh(iM)%fa(iFa)%iM = iM
            ctmp = ifem%msh(iM)%fa(iFa)%name
            DO i=1, iM
               DO j=1, ifem%msh(i)%nFa
                  IF (ctmp .EQ. ifem%msh(i)%fa(j)%name .AND.
     2               (i.NE.iM .OR. j.NE.iFa)) THEN
                     err = " Repeating face names is not allowed"
                  END IF
               END DO
            END DO
         END DO

!     To scale the mesh, while attaching x to gX
         ifem%msh(iM)%scF = 1._RKIND
         lPtr => lPM%get(ifem%msh(iM)%scF,"Mesh scale factor",
     2                                                     lb=0._RKIND)
         a = ifem%tnNo + ifem%msh(iM)%gnNo
         IF (iM .GT. 1) THEN
            ALLOCATE(tmpX(nsd,ifem%tnNo))
            tmpX = gX
            DEALLOCATE(gX)
            ALLOCATE(gX(nsd,a))
            gX(:,1:ifem%tnNo) = tmpX
            DEALLOCATE(tmpX)
         ELSE
            DEALLOCATE(gX)
            ALLOCATE(gX(nsd,a))
         END IF
         gX(:,ifem%tnNo+1:a) = ifem%msh(iM)%x * ifem%msh(iM)%scF
         ifem%tnNo           = a
         DEALLOCATE(ifem%msh(iM)%x)

C          lPtr => lPM%get(ifem%msh(iM)%dx,"Mesh global edge size",1)
      END DO
      ALLOCATE(ifem%x(nsd,ifem%tnNo))
      ALLOCATE(ifem%xCu(nsd,ifem%tnNo))
      ALLOCATE(ifem%xCuo(nsd,ifem%tnNo))
      ifem%x = gX
      ifem%xCu = gX
      ifem%xCuo = gX
      DEALLOCATE(gX)

!     Setting msh%gN, msh%lN parameter
      b = 0
      DO iM=1, ifem%nMsh
         ifem%msh(iM)%nNo = ifem%msh(iM)%gnNo
         ALLOCATE(ifem%msh(iM)%gN(ifem%msh(iM)%nNo), 
     2                             ifem%msh(iM)%lN(ifem%tnNo))
         ifem%msh(iM)%gN = 0
         ifem%msh(iM)%lN = 0
         DO a=1, ifem%msh(iM)%nNo
            b = b + 1
            ifem%msh(iM)%gN(a) = b
            ifem%msh(iM)%lN(b) = a
         END DO
      END DO
      IF (b .NE. ifem%tnNo) err =
     2   " Mismatch in ifem%tnNo. Correction needed"

!     Remap msh%gIEN array
      DO iM=1, ifem%nMsh
         ifem%msh(iM)%nEl = ifem%msh(iM)%gnEl
         ALLOCATE(ifem%msh(iM)%IEN(ifem%msh(iM)%eNoN,ifem%msh(iM)%nEl))
         DO e=1, ifem%msh(iM)%nEl
            DO a=1, ifem%msh(iM)%eNoN
               Ac = ifem%msh(iM)%gIEN(a,e)
               Ac = ifem%msh(iM)%gN(Ac)
               ifem%msh(iM)%IEN(a,e) = Ac
            END DO
         END DO
         DEALLOCATE(ifem%msh(iM)%gIEN)
      END DO

!     Re-arranging fa structure - %gN, %lN, %IEN
      b = 0
      DO iM=1, ifem%nMsh
         DO iFa=1, ifem%msh(iM)%nFa
            ALLOCATE(ifem%msh(iM)%fa(iFa)%lN(ifem%tnNo))
            ifem%msh(iM)%fa(iFa)%lN = 0
            DO a=1, ifem%msh(iM)%fa(iFa)%nNo
               Ac = ifem%msh(iM)%fa(iFa)%gN(a)
               Ac = ifem%msh(iM)%gN(Ac)
               ifem%msh(iM)%fa(iFa)%gN(a)  = Ac
               ifem%msh(iM)%fa(iFa)%lN(Ac) = a
            END DO
            DO e=1, ifem%msh(iM)%fa(iFa)%nEl
               DO a=1, ifem%msh(iM)%fa(iFa)%eNoN
                  Ac = ifem%msh(iM)%fa(iFa)%IEN(a,e)
                  Ac = ifem%msh(iM)%gN(Ac)
                  ifem%msh(iM)%fa(iFa)%IEN(a,e) = Ac
               END DO
            END DO
         END DO
      END DO

!     Setting dmnId parameter here, if there is at least one mesh that
!     has defined eId.
      DO iM=1, nMsh
         lPM => list%get(ifem%msh(iM)%name,"Add IFEM",iM)

         lPtr => lPM%get(i,"Domain (IFEM)",ll=0,
     2                                     ul=BIT_SIZE(ifem%dmnId)-1)
         IF (ASSOCIATED(lPtr)) CALL SETDMNID(ifem%msh(iM), i)

         lPtr => lPM%get(fTmp,"Domain (IFEM) file path")
         IF (ASSOCIATED(lPtr)) THEN
            i = LEN(TRIM(fTmp%fname))
            fExt = fTmp%fname(i-2:i)
            IF (TRIM(fExt).EQ."vtp" .OR. TRIM(fExt).EQ."vtu") THEN
               CALL SETDMNIDVTK(ifem%msh(iM), fTmp%fname, "DOMAIN_ID")
            ELSE
               CALL SETDMNIDFF(ifem%msh(iM), fTmp%open())
            END IF
         END IF

         IF (.NOT.ALLOCATED(ifem%msh(iM)%eId)) err =
     2      " Immersed bodies require domain ID parameter to be set"
      END DO
      ALLOCATE(ifem%dmnId(ifem%tnNo))
      ifem%dmnId = 0
      DO iM=1, ifem%nMsh
         DO e=1, ifem%msh(iM)%nEl
            DO a=1, ifem%msh(iM)%eNoN
               Ac = ifem%msh(iM)%IEN(a,e)
               ifem%dmnId(Ac) = IOR(ifem%dmnId(Ac),ifem%msh(iM)%eId(e))
            END DO
         END DO
      END DO

!     Read fiber orientation
      flag = .FALSE.
      DO iM=1, ifem%nMsh
         lPM => list%get(ifem%msh(iM)%name,"Add IFEM",iM)
         j = lPM%srch("Fiber direction file path")
         IF (j .EQ. 0) j = lPM%srch("Fiber direction")
         IF (j .NE. 0) THEN
            flag = .TRUE.
            EXIT
         END IF
      END DO

      IF (flag) THEN
         DO iM=1, ifem%nMsh
            lPM => list%get(ifem%msh(iM)%name,"Add IFEM",iM)

            ifem%msh(iM)%nFn = lPM%srch("Fiber direction file path")
            j = ifem%msh(iM)%nFn
            IF (ifem%msh(iM)%nFn .NE. 0) THEN
               ALLOCATE(ifem%msh(iM)%fN(j*nsd,ifem%msh(iM)%nEl))
               ifem%msh(iM)%fN = 0._RKIND
               DO i=1, ifem%msh(iM)%nFn
                  lPtr => lPM%get(cTmp, "Fiber direction file path", i)
                  CALL READFIBNFF(ifem%msh(iM), cTmp, "FIFEM_DIR", i)
               END DO
            ELSE
               ifem%msh(iM)%nFn = lPM%srch("Fiber direction")
               j = ifem%msh(iM)%nFn
               IF (ifem%msh(iM)%nFn .NE. 0) THEN
                  ALLOCATE(ifem%msh(iM)%fN(j*nsd,ifem%msh(iM)%nEl))
                  ifem%msh(iM)%fN = 0._RKIND
                  DO i=1, ifem%msh(iM)%nFn
                     lPtr => lPM%get(fibN, "Fiber direction", i)
                     rtmp = SQRT(NORM(fibN))
                     IF (.NOT.ISZERO(rtmp)) fibN(:) = fibN(:)/rtmp
                     DO e=1, ifem%msh(iM)%nEl
                        ifem%msh(iM)%fN((i-1)*nsd+1:i*nsd,e) = 
     2                                                      fibN(1:nsd)
                     END DO
                  END DO
               END IF
            END IF
         END DO
      ELSE
         ifem%msh(:)%nFn = 0
      END IF

      IF (ifem%nMsh .GT. 1) THEN
         std = " Total number of IFEM nodes: "//ifem%tnNo
         std = " Total number of IFEM elements: "//SUM(ifem%msh%nEl)
      END IF

      std = CLR(" IFEM mesh data imported successfully",3)

      RETURN
      END SUBROUTINE IFEM_READMSH
!####################################################################
!     This routine reads IFEM options
C       SUBROUTINE IFEM_READOPTS(list)
C       USE COMMOD
C       USE LISTMOD
C       IMPLICIT NONE
C       TYPE(listType), INTENT(INOUT) :: list

C       CHARACTER(LEN=stdL) :: ctmp
C       TYPE(listType), POINTER :: lPtr

C       lPtr => list%get(ctmp, "IFEM interpolation method")
C       CALL TO_LOWER(ctmp)
C       SELECT CASE (ctmp)
C       CASE ("direct", "nodal")
C          ifem%intrp = ibIntrp_DI
C          std = " IFEM interpolation method: "//CLR("Direct",3)
C       CASE ("l2", "gauss")
C          ifem%intrp = ibIntrp_L2
C          std = " IFEM interpolation method: "//CLR("L2 projection",3)
C       CASE DEFAULT
C          err = " Invalid IFEM interpolation chosen"
C       END SELECT

C       RETURN
C       END SUBROUTINE IFEM_READOPTS
!####################################################################
!     This routine reads IFEM domain properties and BCs in a given Eq
      SUBROUTINE IFEM_READEQ(lEq, list, eqName)
      USE COMMOD
      USE ALLFUN
      USE LISTMOD
      IMPLICIT NONE
      TYPE(eqType), INTENT(INOUT) :: lEq
      TYPE(listType), INTENT(INOUT) :: list
      CHARACTER(LEN=stdL), INTENT(IN) :: eqName

      INTEGER(KIND=IKIND), PARAMETER :: maxOutput = 5

      INTEGER(KIND=IKIND) i, iBc, propL(maxNProp), iDmn, iProp, prop,
     2   nDOP(4), outputs(maxOutput)
      REAL(KIND=RKIND) rtmp
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr, lPD, lPBC

      IF (.NOT.ALLOCATED(ifem%msh)) err = "No IFEM mesh is read yet"
      IF (.NOT.ALLOCATED(ifem%dmnId)) err = "No IFEM domain is read yet"

      lEq%nDmnIB = list%srch("Domain (IFEM)")
      IF(lEq%nDmnIB .GT. 0) write(*,*)": read Domain (IFEM)"
      IF (lEq%nDmnIB .EQ. 0) THEN
         err = "IFEM domain properties not specified"
      ELSE
         ALLOCATE(lEq%dmnIB(lEq%nDmnIB))
      END IF

!--------------------------------------------------------------------
!     Searching through each IFEM domain properties
      DO iDmn=1, lEq%nDmnIB
         lPD => list%get(lEq%dmnIB(iDmn)%Id,"Domain (IFEM)",iDmn,
     2      ll=0,ul=(BIT_SIZE(ifem%dmnId)-1))
         DO i=1, iDmn-1
            IF (lEq%dmnIB(iDmn)%Id .EQ. lEq%dmnIB(i)%Id) THEN
               err = TRIM(list%ping("Domain (IFEM)",lPD))//
     2            " Repeated IFEM domain ID"
            END IF
         END DO

!        IFEM Equation being solved: struct
!        Other potential equations include vms_struct/lElas/shell
         lPtr => lPD%get(ctmp, "Equation", 1)
         propL = prop_NA
         SELECT CASE(TRIM(ctmp))
         CASE("struct")
            lEq%dmnIB(iDmn)%phys = phys_struct
            propL(1) = solid_density
            propL(2) = solid_viscosity
            propL(3) = elasticity_modulus
            propL(4) = poisson_ratio
            propL(5) = ctau_M
            propL(6) = ctau_C
            propL(7) = f_x
            propL(8) = f_y
            IF (nsd .EQ. 3) propL(9) = f_z

            nDOP = (/4,1,0,0/)
            outPuts(1) = out_displacement
            outPuts(2) = out_integ
            outPuts(3) = out_velocity
            outPuts(4) = out_pressure

         CASE DEFAULT
            err = TRIM(lPD%ping("Equation (IFEM)",lPtr))//
     2         "IFEM must be solved using struct equation only"
         END SELECT

!        Domain properties
         DO iProp=1, maxNProp
            rtmp = 0._RKIND
            prop = propL(iProp)
            SELECT CASE (prop)
            CASE (prop_NA)
               EXIT
            CASE (solid_density)
               lPtr => lPD%get(rtmp,"Density",1,ll=0._RKIND)
            CASE (solid_viscosity)
               lPtr => lPD%get(rtmp,"Viscosity",ll=0._RKIND)
            CASE (elasticity_modulus)
               lPtr => lPD%get(rtmp,"Elasticity modulus",1,lb=0._RKIND)
            CASE (poisson_ratio)
               lPtr => lPD%get(rtmp,"Poisson ratio",1,ll=0._RKIND,
     2            ul=0.5_RKIND)
            CASE (f_x)
               lPtr => lPD%get(rtmp,"Force_X")
            CASE (f_y)
               lPtr => lPD%get(rtmp,"Force_Y")
            CASE (f_z)
               lPtr => lPD%get(rtmp,"Force_Z")
            CASE (ctau_M)
               rtmp = 1e-3_RKIND
               lPtr => lPD%get(rtmp,
     2            "Momentum stabilization coefficient")
            CASE (ctau_C)
               rtmp = 1e-3_RKIND
               lPtr => lPD%get(rtmp,
     2            "Continuity stabilization coefficient")
            CASE DEFAULT
               err = "Undefined properties (IFEM)"
            END SELECT
            lEq%dmnIB(iDmn)%prop(prop) = rtmp
         END DO

         IF (lEq%dmnIB(iDmn)%phys .EQ. phys_struct) THEN
            CALL READMATMODEL(lEq%dmnIB(iDmn), lPD)
         END IF
      END DO

!     Read IFEM outputs
      CALL IFEM_READOUTPUTS(lEq, nDOP, outPuts, list)

!     Set number of function spaces
      DO i=1, nMsh
         ifem%msh(i)%nFs = 1
      END DO

!--------------------------------------------------------------------
!     Searching for BCs on immersed bodies
      lEq%nBcIB = list%srch("Add BC (IFEM)")
      ALLOCATE(lEq%bcIB(lEq%nBcIB))
      std = " Number of imposed BC for equation <"//TRIM(eqName)//
     2   ">: "//lEq%nBcIB
      DO iBc=1, lEq%nBcIB
         lPBC => list%get(ctmp,"Add BC (IFEM)",iBc)
         CALL IFEM_FINDFACE(ctmp, lEq%bcIB(iBc)%iM, lEq%bcIB(iBc)%iFa)
         write(*,*) "IFEM_FINDFACE done"
         CALL IFEM_READBC(lEq%bcIB(iBc), lPBC)
      END DO

      write(*,*) "IFEM_FINDFACE and IFEM_READBC done"

      RETURN
      END SUBROUTINE IFEM_READEQ
!####################################################################
!     This subroutine is to read from input file on how to process the
!     output quantities: for VTK files or for surface/volume integrated
!     quantities. nDOP(1) is the total number of outputs, nDOP(2) is the
!     default number of outputs for VTK files, nDOP(3) is for boundaries
!     nDOP(4) is for volume
      SUBROUTINE IFEM_READOUTPUTS(lEq, nDOP, outputs, list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eqType), INTENT(INOUT) :: lEq
      INTEGER(KIND=IKIND), INTENT(IN) :: nDOP(4), outputs(nDOP(1))
      TYPE(listType), INTENT(INOUT) :: list

      INTEGER(KIND=IKIND) nOut, iOut, i, j
      CHARACTER(LEN=stdL) ctmp, stmp
      TYPE(listType), POINTER :: lPtr, lPO

      lEq%nOutIB = nDOP(1)
      ALLOCATE(lEq%outIB(nDOP(1)))

      DO iOut=1, nDOP(1)
         SELECT CASE (outputs(iOut))
         CASE (out_displacement)
            lEq%outIB(iOut)%grp  = outGrp_D
            lEq%outIB(iOut)%o    = 0
            lEq%outIB(iOut)%l    = nsd
            lEq%outIB(iOut)%name = "Displacement"
         CASE (out_integ)
            lEq%outIB(iOut)%grp  = outGrp_I
            lEq%outIB(iOut)%o    = 0
            lEq%outIB(iOut)%l    = 1
            IF (nsd .EQ. 2) THEN
               lEq%outIB(iOut)%name = "Area"
            ELSE
               lEq%outIB(iOut)%name = "Volume"
            END IF
         CASE (out_velocity)
            lEq%outIB(iOut)%grp  = outGrp_Y
            lEq%outIB(iOut)%o    = 0
            lEq%outIB(iOut)%l    = nsd
            lEq%outIB(iOut)%name = "Velocity"
         CASE (out_pressure)
            lEq%outIB(iOut)%grp  = outGrp_Y
            lEq%outIB(iOut)%o    = nsd
            lEq%outIB(iOut)%l    = 1
            lEq%outIB(iOut)%name = "Pressure"
         CASE DEFAULT
            err = "IFEM output undefined"
         END SELECT
      END DO

!     These are the default values, we use the first nDef/nBDef outputs
      DO j=1, 3
         lEq%outIB(1:nDOP(j+1))%wtn(j) = .TRUE.
      END DO

!     First reading the outputs for VTK files and then for boundaries
!     and last for the volume
      nOut = list%srch("Output (IFEM)")
      DO iOut=1, nOut
         lPO => list%get(ctmp,"Output (IFEM)",iOut)
         SELECT CASE(TRIM(ctmp))
         CASE("Spatial")
            j = 1
         CASE("B_INT", "Boundary_integral")
            j = 2
         CASE("V_INT", "Volume_integral")
            j = 3
         CASE("Alias")
            CYCLE
         CASE DEFAULT
            j = -1
            err = TRIM(list%ping("Output (IFEM)",lPO))//
     2         " Undefined keyword"
         END SELECT
         DO i=1, lEq%nOutIB
            lPtr => lPO%get(lEq%outIB(i)%wtn(j),lEq%outIB(i)%name)
         END DO
      END DO

!     Read any alias names for outputs
      DO iOut=1, nOut
         lPO => list%get(ctmp,"Output (IFEM)",iOut)
         SELECT CASE (TRIM(ctmp))
         CASE ("Alias")
            DO i=1, lEq%nOutIB
               lPtr => lPO%get(stmp, TRIM(lEq%output(i)%name))
               IF (ASSOCIATED(lPtr)) lEq%outIB(i)%name = TRIM(stmp)
            END DO
            EXIT
         CASE DEFAULT
            CYCLE
         END SELECT
      END DO

      RETURN
      END SUBROUTINE IFEM_READOUTPUTS
!####################################################################
!     Finding the face ID and mesh ID based on the face name
      SUBROUTINE IFEM_FINDFACE(faName, iM, iFa)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(OUT) :: iM, iFa
      CHARACTER(LEN=stdL) :: faName

      iFa = 0
      MY_LOOP : DO iM=1, ifem%nMsh
         IF (ifem%msh(iM)%lShl) THEN
!        If the BC is applied on the shell itself
            IF (ifem%msh(iM)%name .EQ. faName) THEN
               iFa = 0
               EXIT
            END IF
         END IF
         DO iFa=1, ifem%msh(iM)%nFa
            IF (ifem%msh(iM)%fa(iFa)%name .EQ. faName) EXIT MY_LOOP
         END DO
      END DO MY_LOOP

      IF (iM .GT. ifem%nMsh) err =
     2   " Unable to find face (IFEM) <"//TRIM(faName)//">"

      RETURN
      END SUBROUTINE IFEM_FINDFACE
!####################################################################
!     This routine reads BC for immersed bodies
      SUBROUTINE IFEM_READBC(lBc, list)
      USE COMMOD
      USE ALLFUN
      USE LISTMOD
      IMPLICIT NONE
      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(listType), INTENT(INOUT) :: list

      LOGICAL ltmp
      INTEGER(KIND=IKIND) a, b, i, j, Ac, iM, iFa, fid, nNo
      REAL(KIND=RKIND) rtmp
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr
      TYPE(fileType) fTmp

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)

      iM  = lBc%iM
      iFa = lBc%iFa
      nNo = ifem%msh(iM)%fa(iFa)%nNo

!     Reading the type: Dir/Neu
      lPtr => list%get(ctmp,"Type")
      SELECT CASE (ctmp)
      CASE ("Dirichlet","Dir")
         lBc%bType = IBSET(lBc%bType,bType_Dir)
      CASE ("Neumann","Neu")
         lBc%bType = IBSET(lBc%bType,bType_Neu)
      CASE DEFAULT
         err = TRIM(list%ping("Type",lPtr))//" Unexpected IFEM BC type"
      END SELECT

      ALLOCATE(lBc%eDrn(nsd))
      lBc%eDrn = 0
      lPtr => list%get(lBc%eDrn,"Effective direction")

      ctmp = "Steady"
      lPtr => list%get(ctmp,"Time dependence")
      SELECT CASE (ctmp)
      CASE ('Steady')
         lBc%bType = IBSET(lBc%bType,bType_std)
         lPtr => list%get(lBc%g,"Value",1)
      CASE ('Unsteady')
         lBc%bType = IBSET(lBc%bType,bType_ustd)
         ALLOCATE(lBc%gt)
         lBc%gt%d = 1
         lPtr => list%get(fTmp,"Temporal values file path")
         IF (ASSOCIATED(lPtr)) THEN
            ltmp = .FALSE.
            lPtr => list%get(ltmp,"Ramp function")
            lBc%gt%lrmp = ltmp
            fid = fTmp%open()
            READ(fid,*) i, j
            IF (i .LT. 2) THEN
               std = "Enter nPnts nFCoef; nPts*(t Q)"
               err = "Wrong format in: "//fTmp%fname
            END IF
            lBc%gt%n = j
            IF (lBc%gt%lrmp) lBc%gt%n = 1
            ALLOCATE(lBc%gt%qi(1), lBc%gt%qs(1))
            ALLOCATE(lBc%gt%r(1,j))
            ALLOCATE(lBc%gt%i(1,j))
            CALL FFT(fid, i, lBc%gt)
            CLOSE(fid)
         ELSE
            lPtr => list%get(fTmp,"Fourier coefficients file path",1)
            IF (.NOT.ASSOCIATED(lPtr)) err = "Undefined inputs for "//
     2         "unsteady type BC"
            lBc%gt%lrmp = .FALSE.
            ALLOCATE(lBc%gt%qi(1), lBc%gt%qs(1))
            fid = fTmp%open()
            READ (fid,*) lBc%gt%ti
            READ (fid,*) lBc%gt%T
            READ (fid,*) lBc%gt%qi(1)
            READ (fid,*) lBc%gt%qs(1)
            READ (fid,*) j
            lBc%gt%n = j
            ALLOCATE(lBc%gt%r(1,j))
            ALLOCATE(lBc%gt%i(1,j))
            DO i=1, j
               READ (fid,*) lBc%gt%r(1,i), lBc%gt%i(1,i)
            END DO
            CLOSE(fid)
         END IF
      CASE ('Coupled')
         err = " Cannot apply Coupled BCs for immersed bodies"
      CASE ('Resistance')
         err = " Cannot apply Resistance BCs for immersed bodies"
      CASE ('RCR', 'Windkessel')
         err = " Cannot apply RCR BCs for immersed bodies"
      CASE ('Spatial')
         err = " Cannot apply Spatial Neu BCs for immersed bodies"
      CASE ('General')
         lBc%bType = IBSET(lBc%bType,bType_gen)
         lPtr => list%get(ftmp,"BCT file path")
         IF (ASSOCIATED(lPtr)) THEN
            ALLOCATE(lBc%gm)
            CALL READBCT(lBc%gm, ifem%msh(iM), ifem%msh(iM)%fa(iFa),
     2         ftmp%fname)
         ELSE
            lPtr =>list%get(fTmp,
     2         "Temporal and spatial values file path",1)
            IF (.NOT.ASSOCIATED(lPtr)) err = "Input file for General"//
     2         " BC (bct_ib.vtp) is not provided"
            fid = fTmp%open()
            READ (fid,*) i, j, a
            IF (a .NE. nNo) THEN
               err = "Number of nodes does not match between "//
     2            TRIM(ifem%msh(iM)%fa(iFa)%name)//" and "//fTmp%fname
            END IF
            IF (i.LT.1 .OR. i.GT.nsd) err = "0 < dof <= "//nsd//
     2         " is violated in "//fTmp%fname

            ALLOCATE(lBc%gm)
            ALLOCATE(lBc%gm%t(j),lBc%gm%d(i,a,j),ptr(ifem%msh(iM)%nNo))
!     I am seting all the nodes to zero just in case a node is not set
            lBc%gm%d   = 0._RKIND
            lBc%gm%dof = i
            lBc%gm%nTP = j
            ptr        = 0
!     Preparing the pointer array
            DO a=1, nNo
               Ac = ifem%msh(iM)%fa(iFa)%gN(a)
               Ac = ifem%msh(iM)%lN(Ac)
               IF (Ac .EQ. 0) err = "Incorrect global node number "//
     2            "detected for BC. Mesh: "//TRIM(ifem%msh(iM)%name)//
     3            ", Face: "//TRIM(ifem%msh(iM)%fa(iFa)%name)//
     4            ", Node: "//STR(a)//" gN: "//
     5            STR(ifem%msh(iM)%fa(iFa)%gN(a))
               ptr(Ac) = a
            END DO
            DO i=1, j
               READ (fid,*) rtmp
               lBc%gm%t(i) = rtmp
               IF (i .EQ. 1) THEN
                  IF (.NOT.ISZERO(rtmp)) err = "First time step"//
     2               " should be zero in <"//TRIM(ftmp%fname)//">"
               ELSE
                  rtmp = rtmp - lBc%gm%t(i-1)
                  IF (ISZERO(rtmp) .OR. rtmp.LT.0._RKIND) err =
     2               "Non-increasing time trend is found in <"//
     3               TRIM(ftmp%fname)//">"
               END IF
            END DO

            lBc%gm%period = lBc%gm%t(j)
            DO b=1, nNo
               READ (fid,*) Ac
               IF (Ac.GT.ifem%msh(iM)%nNo .OR. Ac.LE.0) THEN
                  err = "Entry "//b//" is out of bound in "//ftmp%fname
               END IF
               a = ptr(Ac)
               IF (a .EQ. 0) THEN
                  err = "Entry "//b//" not found in "//ftmp%fname
               END IF
               DO i=1, j
                  READ (fid,*) lBc%gm%d(:,a,i)
               END DO
            END DO
            CLOSE(fid)
            DEALLOCATE(ptr)
         END IF
      CASE DEFAULT
         err = TRIM(list%ping("Time dependence",lPtr))//
     2      " Unexpected type"
      END SELECT

!     To zero-out perimeter or not. Default is .true. for Dir
      ltmp = .FALSE.
      IF (BTEST(lBc%bType,bType_Dir)) ltmp = .TRUE.
      lPtr => list%get(ltmp,"Zero out perimeter")
      lBc%bType = IBCLR(lBc%bType,bType_zp)
      IF (ltmp) lBc%bType = IBSET(lBc%bType,bType_zp)

!     Impose BC on the state variable or its integral
      ltmp = .TRUE.
      lPtr => list%get(ltmp,"Impose on state variable integral")
      lBc%bType = IBCLR(lBc%bType,bType_impD)
      IF (ltmp) lBc%bType = IBSET(lBc%bType,bType_impD)

!     Reading the spatial profile: flat/para/ud
      ctmp = "Flat"
      lPtr => list%get(ctmp,"Profile")
      SELECT CASE (ctmp)
      CASE ('Flat')
         lBc%bType = IBSET(lBc%bType,bType_flat)
      CASE ('Parabolic')
         lBc%bType = IBSET(lBc%bType,bType_para)
      CASE ('User_defined')
         lBc%bType = IBSET(lBc%bType,bType_ud)
         lPtr => list%get(fTmp,"Spatial profile file path",1)
         fid = fTmp%open()
         ALLOCATE(lBc%gx(nNo), ptr(ifem%msh(iM)%nNo))
!     I am seting all the nodes to zero just in case a node is not set
         ptr = 0
!     Preparing the pointer array
         DO a=1, nNo
            lBc%gx(a) = 0._RKIND
            Ac = ifem%msh(iM)%fa(iFa)%gN(a)
            Ac = ifem%msh(iM)%lN(Ac)
            IF (Ac .EQ. 0) err = "Incorrect global node number "//
     2         "detected for BC. Mesh: "//TRIM(ifem%msh(iM)%name)//
     3         ", Face: "//TRIM(ifem%msh(iM)%fa(iFa)%name)//
     4         ", Node: "//STR(a)//" gN: "//
     5         STR(ifem%msh(iM)%fa(iFa)%gN(a))
            ptr(Ac) = a
         END DO

         DO b=1, nNo
            READ (fid,*) Ac, rtmp
            IF (Ac.GT.ifem%msh(iM)%nNo .OR. Ac.LE.0) THEN
               err = "Entry "//b//" is out of bound in "//fTmp%fname
            END IF
            a = ptr(Ac)
            IF (a .EQ. 0) THEN
               err = "Entry <"//b//" not found in face "//fTmp%fname
            END IF
            lBc%gx(a) = rtmp
         END DO
         CLOSE(fid)
         DEALLOCATE(ptr)
      CASE DEFAULT
         err = TRIM(list%ping("Profile",lPtr))//" Unexpected profile"
      END SELECT

      lBc%weakDir = .FALSE.
      lBc%flwP = .FALSE.

      RETURN
      END SUBROUTINE IFEM_READBC
!####################################################################
!     Allocates memory for IFEM data structures
      SUBROUTINE IFEM_MEMALLOC()
      USE COMMOD
      IMPLICIT NONE

!     What is what: 
!     Yb = unknown  
!     Auo = time derivative of displacement old 
!     Ubo = displacement old
!     Rfluid = FSI force spread to fluid mesh 
!     Rsolid = whole FSI force  
! ?? TODO
      ALLOCATE(ifem%Yb(nsd,ifem%tnNo)) ! ?? why nsd+1
      ALLOCATE(ifem%Auo(nsd,ifem%tnNo))
      ALLOCATE(ifem%Auoo(nsd,ifem%tnNo))
      ALLOCATE(ifem%Ubo(nsd,ifem%tnNo))
      ALLOCATE(ifem%Rfluid(nsd,tnNo))
      ALLOCATE(ifem%Rsolid(nsd,ifem%tnNo))

C       IF (ifem%cpld .EQ. ibCpld_I) THEN
C          ALLOCATE(ifem%Aun(nsd,ifem%tnNo))
C          ALLOCATE(ifem%Auk(nsd,ifem%tnNo))
C          ALLOCATE(ifem%Ubn(nsd,ifem%tnNo))
C          ALLOCATE(ifem%Ubk(nsd,ifem%tnNo))
C          ALLOCATE(ifem%Uo(nsd,tnNo))
C          ALLOCATE(ifem%Un(nsd,tnNo))
C          ALLOCATE(ifem%Ru(nsd,tnNo))
C          ALLOCATE(ifem%Rub(nsd,ifem%tnNo))
C          ALLOCATE(ifem%Ku((nsd+1)*nsd,lhs%nnz))
C       END IF

      RETURN
      END SUBROUTINE IFEM_MEMALLOC
!####################################################################
!     This routine initializes IFEM solution and FSILS data structures
!     from the fluid displacement old Dg = Do 
      SUBROUTINE IFEM_INIT(Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, iM, iFa, iEq, iDmn, iBc, nnz, e
      REAL(KIND=RKIND) tt, lD(nsd,tnNo), s(1,ifem%tnNo)
      CHARACTER(LEN=stdL) sOut

      tt = CPUT()
      ifem%callD = 0._RKIND

      std = " ================================="
      std = " Initializing IFEM data structures.."

      ifem%cEq = 0
      DO iEq=1, nEq
         IF (eq(iEq)%phys .EQ. phys_fluid) THEN
            ifem%cEq = iEq
            EXIT ! exit the DO loop 
         END IF
      END DO
      IF (ifem%cEq .EQ. 0) err = "Fluid equation not detected "//
     2   "to treat immersed bodies"

      lD = 0._RKIND
      IF (mvMsh) THEN
         DO a=1, tnNo
            lD(:,a) = Dg(nsd+2:2*nsd+1,a)
         END DO
      END IF

!     Initialize function spaces
      DO iM=1, ifem%nMsh
         CALL INITFSMSH(ifem%msh(iM))
         DO iFa=1, ifem%msh(iM)%nFa
            CALL INITFSFACE(ifem%msh(iM), ifem%msh(iM)%fa(iFa))
         END DO
      END DO

!     Calculating the volume of each domain 
!     still TODO
      s = 1._RKIND
      DO iEq=1, nEq
         DO iDmn=1, eq(iEq)%nDmnIB
            eq(iEq)%dmnIB(iDmn)%v =
     2         Integ(eq(iEq)%dmnIB(iDmn)%Id, s, 1, 1)
            std = " Volume of domain IFEM "//STR(eq(iEq)%dmnIB(iDmn)%v)
            IF (ISZERO(eq(iEq)%dmnIB(iDmn)%v)) wrn = "Volume of "//
     2         "domain "//iDmn//" of equation "//iEq//" is zero"
         END DO
      END DO

!     Initialize IFEM face normals
      write(*,*) "::: Calling IFEM_FACEINI :::"
      DO iM=1, ifem%nMsh
         DO iFa=1, ifem%msh(iM)%nFa
            ifem%msh(iM)%fa(iFa)%iM = iM
            CALL IFEM_FACEINI(ifem%msh(iM), ifem%msh(iM)%fa(iFa))
         END DO
      END DO

!     Initialize IFEM face BC profile
      write(*,*) "::: Calling IFEM_BCINI :::"
      DO iEq=1, nEq
         DO iBc=1, eq(iEq)%nBcIB
            iFa = eq(iEq)%bcIB(iBc)%iFa
            iM  = eq(iEq)%bcIB(iBc)%iM
            CALL IFEM_BCINI(eq(iEq)%bcIB(iBc), ifem%msh(iM)%fa(iFa))
         END DO
      END DO

!     Create sparse matrix data structures 
!     ?? I don't see why do we need this
C       CALL IFEM_LHSA(nnz)
C       std = "    Non-zeros in LHS matrix (IFEM): "//nnz

!     Set IFEM Dirichlet BCs
      write(*,*) "::: Calling IFEM_SETBCDIR :::"
      CALL IFEM_SETBCDIR(ifem%Yb, ifem%Ubo)

!     Compute the nodal stencil/adjacency for each fluid node 
      DO iM=1, nMsh
         msh(iM)%iGC = 0
         write(*,*) "::: Calling GETNSTENCIL :::"
         CALL GETNSTENCIL(msh(iM))
      END DO

!     Find the closest fluid node for each solid node 
!     and search in which fluid element the solid node belongs
      DO iM=1, ifem%nMsh
         CALL IFEM_FINDCLOSEST(ifem%msh(iM), lD)
      END DO

!     DO something for paralle here??
!     Initialize IFEM communication structure
!      CALL IFEM_SETCOMMU()

!     Identify ghost cells for immersed boundaries
!      CALL IFEM_SETIGHOST()

      tt = CPUT() - tt
      WRITE(sOut,'(F6.2)') tt
      WRITE(sOut,'(A)') "    IFEM setting time: "//TRIM(sOut)//" sec"
      std = TRIM(sOut)

      std = " IFEM initialization complete.."
      std = " ================================="

      RETURN
      END SUBROUTINE IFEM_INIT
!####################################################################
!     Update iblank, tracers and communication structure
      SUBROUTINE IFEM_UPDATE(Dg)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) :: a, iM
      REAL(KIND=RKIND) tt
      REAL(KIND=RKIND), ALLOCATABLE :: lD(:,:)

      ifem%callD = 0._RKIND
      tt = CPUT()

      ALLOCATE(lD(nsd,tnNo))
      lD = 0._RKIND
      IF (mvMsh) THEN
         DO a=1, tnNo
            lD(:,a) = Dg(nsd+2:2*nsd+1,a)
         END DO
      END IF

!     TO DO FOR PARALLEL Reset IFEM communicator
!       CALL IFEM_SETCOMMU()

!     Reompute the nodal stencil/adjacency for each fluid node 
!     in case of remeshing of the fluid mesh 
!      DO iM=1, nMsh
!         msh(iM)%iGC = 0
!         write(*,*) "::: Calling GETNSTENCIL :::"
!         CALL GETNSTENCIL(msh(iM))
!      END DO

!     Find the closest fluid node for each solid node 
!     and search in which fluid element the solid node belongs
      DO iM=1, ifem%nMsh
         CALL IFEM_UPDATECLS(ifem%msh(iM), lD)
      END DO

      write(*,*)"END IFEM_UPDATECLS "

      ifem%callD(2) = CPUT() - tt
      ifem%callD(4) = ifem%callD(3)
      ifem%callD(3) = 0._RKIND

      RETURN
      END SUBROUTINE IFEM_UPDATE
!####################################################################
!     Initializing immersed boundary faces
      SUBROUTINE IFEM_FACEINI(lM, lFa)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa

      LOGICAL flag
      INTEGER(KIND=IKIND) a, b, e, g, Ac, Bc, Ec
      REAL(KIND=RKIND) area, tmp, xi0(nsd), xi(nsd), xp(nsd), nV(nsd)
      TYPE(fsType) :: fs, fsb

      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: sV(:,:), sVl(:,:), xbl(:,:),
     2   xl(:,:), N(:), Nxi(:,:)

!     Calculating face area
      ALLOCATE(N(ifem%tnNo))
      N    = 1._RKIND
      area = Integ(lFa, N)
      std  = "    Area of face <"//TRIM(lFa%name)//"> is "//STR(area)
      IF (ISZERO(area)) THEN
         IF (cm%mas()) wrn = " <"//TRIM(lFa%name)//"> area is zero"
      END IF
      lFa%area = area
      DEALLOCATE(N)

!     Compute face normals at nodes
      IF (ALLOCATED(lFa%nV)) DEALLOCATE(lFa%nV)
      ALLOCATE(lFa%nV(nsd,lFa%nNo), sV(nsd,ifem%tnNo))

      flag = .FALSE.
      IF (lM%eType.EQ.eType_TRI6  .OR. lM%eType.EQ.eType_QUD8  .OR.
     2    lM%eType.EQ.eType_QUD9  .OR. lM%eType.EQ.eType_TET10 .OR.
     3    lM%eType.EQ.eType_HEX20 .OR. lM%eType .EQ. eType_HEX27) THEN
         flag =.TRUE.
      END IF

      IF (.NOT.flag) THEN
!        For linear elements or NURBS, we simply project element normals
!        to nodes
         DO e=1, lFa%nEl
            IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(lM, lFa, e)
            DO g=1, lFa%nG
               CALL GNNIFEM(lFa, e, g, nV)
               DO a=1, lFa%eNoN
                  Ac       = lFa%IEN(a,e)
                  sV(:,Ac) = sV(:,Ac) + nV*lFa%N(a,g)*lFa%w(g)
               END DO
            END DO
         END DO

      ELSE
!        For higher order elements, use reduced order basis on mesh
!        to project element normals. Lumping method is used to project
!        to face corners. Normals at edge nodes are computed by simple
!        interpolation from reduced basis. Standard lumping using higher
!        order basis could lead to spurious errors
         CALL SETTHOODFS(fs, lM%eType)
         CALL ALLOCFS(fs, nsd)
         CALL SETTHOODFS(fsb, lFa%eType)
         CALL GETGIP(nsd, fs%eType, fs%nG, fs%w, fs%xi)
         DO g=1, fs%nG
            CALL GETGNN(nsd, fs%eType, fs%eNoN, fs%xi(:,g), fs%N(:,g),
     2         fs%Nx(:,:,g))
         END DO
         CALL GETNNBNDS(fs%eType, fs%eNoN, fs%xib, fs%Nb)

         xi0 = 0._RKIND
         DO g=1, fs%nG
            xi0 = xi0 + fs%xi(:,g)
         END DO
         xi0 = xi0 / REAL(fs%nG, KIND=RKIND)

         ALLOCATE(sVl(nsd,lFa%eNoN), xbl(nsd,lFa%eNoN), xl(nsd,fs%eNoN),
     2      N(fs%eNoN), Nxi(nsd,fs%eNoN), ptr(fs%eNoN))
         DO e=1, lFa%nEl
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               xbl(:,a) = ifem%x(:,Ac) + ifem%Ubo(:,Ac)
            END DO

            Ec  = lFa%gE(e)
            ptr = 0
            DO a=1, fs%eNoN
               Ac = lM%IEN(a,Ec)
               xl(:,a) = ifem%x(:,Ac) + ifem%Ubo(:,Ac)
               DO b=1, fsb%eNoN
                  Bc = lFa%IEN(b,e)
                  IF (Ac .EQ. Bc) THEN
                     ptr(a) = b
                     EXIT
                  END IF
               END DO
            END DO

            sVl(:,:) = 0._RKIND
            DO g=1, lFa%nG
               xp = 0._RKIND
               DO a=1, lFa%eNoN
                  xp = xp + lFa%N(a,g)*xbl(:,a)
               END DO

               xi = xi0
               CALL GETNNX(fs%eType, fs%eNoN, xl, fs%xib, fs%Nb, xp,
     2            xi, N, Nxi)

               CALL GNNIB(lFa, e, g, nV)

               DO a=1, fs%eNoN
                  b = ptr(a)
                  IF (b .EQ. 0) CYCLE
                  Ac = lM%IEN(a,Ec)
                  sVl(:,b) = sVl(:,b) + lFa%w(g)*N(a)*nV(:)
                  sV(:,Ac) = sV(:,Ac) + sVl(:,b)
               END DO
            END DO

            DO b=fsb%eNoN+1, lFa%eNoN
               xp = xbl(:,b)
               xi = xi0
               CALL GETNNX(fs%eType, fs%eNoN, xl, fs%xib, fs%Nb, xp,
     2            xi, N, Nxi)

               DO a=1, fs%eNoN
                  IF (ptr(a) .EQ. 0) CYCLE
                  sVl(:,b) = sVl(:,b) + N(a)*sVl(:,ptr(a))
               END DO

               Ac = lFa%IEN(b,e)
               sV(:,Ac) = sV(:,Ac) + sVl(:,b)
            END DO
         END DO
         DEALLOCATE(sVl, xbl, xl, N, Nxi, ptr)
         CALL DESTROY(fs)
         CALL DESTROY(fsb)
      END IF

      flag = .TRUE.
      DO a=1, lFa%nNo
         Ac  = lFa%gN(a)
         tmp = SQRT(NORM(sV(:,Ac)))
         IF (ISZERO(tmp)) THEN
            IF (flag) THEN
               wrn = " Skipping normal calculation of node "//a//
     2            " in face <"//TRIM(lFa%name)//">"
               flag = .FALSE.
            END IF
            lFa%nV(:,a) = 0._RKIND
            lFa%nV(1,a) = 1._RKIND
            CYCLE
         END IF
         lFa%nV(:,a) = sV(:,Ac)/tmp
      END DO
      DEALLOCATE(sV)

      RETURN
      END SUBROUTINE IFEM_FACEINI
!--------------------------------------------------------------------
!     Set BC spatial profile on the face
      SUBROUTINE IFEM_BCINI(lBc, lFa)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER(KIND=IKIND) iM, iFa, jFa, i, a, b, Ac, j
      REAL(KIND=RKIND) tmp, nV(nsd), center(nsd), maxN

      INTEGER(KIND=IKIND), ALLOCATABLE :: gNodes(:)
      REAL(KIND=RKIND), ALLOCATABLE :: s(:), sV(:,:), sVl(:,:)

      IF (BTEST(lBc%bType,bType_gen)) THEN
         IF (BTEST(lBc%bType,bType_Neu) .AND. lBc%gm%dof.NE.1) THEN
            err = " Only DOF=1 is accepted for Neu general BCs"
         END IF
         RETURN
      END IF

      iM  = lFa%iM
      iFa = lBc%iFa
      IF (.NOT.ALLOCATED(lBc%gx)) ALLOCATE(lBc%gx(lFa%nNo))

      ALLOCATE(s(ifem%tnNo))
      s = 0._RKIND
      IF (BTEST(lBc%bType,bType_flat)) THEN
!     Just a constant value for Flat profile
         DO a=1, lFa%nNo
            Ac    = lFa%gN(a)
            s(Ac) = 1._RKIND
         END DO
      ELSE IF (BTEST(lBc%bType,bType_para)) THEN
!     Here is the method that is used for imposing parabolic profile:
!     1- Find the coordinate of the points on the boundary 2- find unit
!     vector from center to each of points on the boundary: ew
!     3- maximize ew(i).e where e is the unit vector from current
!     point to the center 4- Use the point i as the diam here
         center = 0._RKIND
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            center(:) = center(:) + ifem%x(:,Ac)
         END DO
         center(:) = center(:)/REAL(lFa%nNo, KIND=RKIND)
         ALLOCATE(gNodes(ifem%tnNo), sV(nsd,ifem%tnNo))
!     gNodes is one if a node located on the boundary (beside iFa)
         gNodes = 0
         DO jFa=1, ifem%msh(iM)%nFa
            IF (jFa .EQ. iFa) CYCLE
            DO a=1, ifem%msh(iM)%fa(jFa)%nNo
               Ac         = ifem%msh(iM)%fa(jFa)%gN(a)
               gNodes(Ac) = 1
            END DO
         END DO
!     "j" is a counter for the number of nodes that are located on the
!     boundary of lFa and sVl contains the list of their coordinates
         j  = 0
         sV = 0._RKIND
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            IF (gNodes(Ac) .EQ. 1) THEN
               j = j + 1
               sV(:,j) = ifem%x(:,Ac)
            END IF
         END DO
         IF (cm%mas() .AND. j.EQ.0) err = " No perimeter"//
     2      " found for face "//lFa%name
!     sVl will keep the normal unit vector from center to perimeter
         ALLOCATE(sVl(nsd,j))
         DO a=1, j
            sV(:,a)  = sV(:,a) - center
            sVl(:,a) = sV(:,a)/SQRT(NORM(sV(:,a)))
         END DO
!     "s" is going to keep the ew.e value
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            nV = ifem%x(:,Ac) - center
            maxN = NORM(nV, sVl(:,1))
            i = 1
            DO b=2, j
               tmp = NORM(nV, sVl(:,b))
               IF (tmp .GT. maxN) THEN
                  maxN = tmp
                  i = b
               END IF
            END DO
            s(Ac) = 1._RKIND - NORM(nV)/NORM(sV(:,i))
         END DO
      ELSE IF (BTEST(lBc%bType,bType_ud)) THEN
         DO a=1, lFa%nNo
            Ac    = lFa%gN(a)
            s(Ac) = lBc%gx(a)
         END DO
      END IF

!     Now correcting the inlet BC for the inlet ring
      IF (BTEST(lBc%bType,bType_zp)) THEN
         DO jFa=1, ifem%msh(iM)%nFa
            IF (jFa .EQ. iFa) CYCLE
            DO a=1, ifem%msh(iM)%fa(jFa)%nNo
               Ac    = ifem%msh(iM)%fa(jFa)%gN(a)
               s(Ac) = 0._RKIND
            END DO
         END DO
      END IF

      DO a=1, lFa%nNo
         Ac        = lFa%gN(a)
         lBc%gx(a) = s(Ac)
      END DO

      RETURN
      END SUBROUTINE IFEM_BCINI
!####################################################################
!     Set iblank field
!        iblank is set only immersed solids and not set for thin shells
!        iblank(A) = 1   =>   node A is inside the immersed solid
!        iblank(A) = 0   =>   node A is outside the immersed solid
      SUBROUTINE IFEM_SETIBLANK(lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      LOGICAL :: incL
      INTEGER(KIND=IKIND) :: a, Ac, i, iM, itmp
      REAL(KIND=RKIND) :: minb(nsd), maxb(nsd), xp(nsd), dx

      LOGICAL, ALLOCATABLE :: chck(:)

!     Initialize iblank field
      iblank(:) = 0

!     Fit a bounding box around IFEM and probe only those nodes lying
!     inside the box
      dx = TINY(dx)
      DO iM=1, ifem%nMsh
         IF (dx .LT. ifem%msh(iM)%dx) THEN
            dx = ifem%msh(iM)%dx
         END IF
      END DO

      DO i=1, nsd
         minb(i) = MINVAL(ifem%x(i,:) + ifem%Ubo(i,:)) - dx
         maxb(i) = MAXVAL(ifem%x(i,:) + ifem%Ubo(i,:)) + dx
      END DO

      ALLOCATE(chck(tnNo))
      chck = .FALSE.
      DO Ac=1, tnNo
         xp(:) = x(:,Ac) + lD(:,Ac)
         itmp = 0
         DO i=1, nsd
            IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i)) THEN
               itmp = itmp + 1
            END IF
         END DO
         IF (itmp .EQ. nsd) chck(Ac) = .TRUE.
      END DO

      DO iM=1, nMsh
!        Probe each node if it is inside or outside the immersed body
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            IF (.NOT.chck(Ac)) CYCLE

            xp(:) = x(:,Ac) + lD(:,Ac)
            incL  = .FALSE.
            CALL IFEM_CHECKINOUT(xp, incL)
            IF (incL) iblank(Ac) = 1
            chck(Ac) = .FALSE.
         END DO
      END DO

      RETURN
      END SUBROUTINE IFEM_SETIBLANK
!--------------------------------------------------------------------
!     Update iblank field around the neighborhood of previous iblank. In
!     case iblank is 0 entirely on the current process, redo iblank
!     initialization to look for any fresh iblank node
      SUBROUTINE IFEM_UPDATEIBLANK(lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      LOGICAL flag
      INTEGER(KIND=IKIND) a, b, i, Ac, Bc, iswp, iM, nNo
      REAL(KIND=RKIND) xp(nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), ptr(:)

!     Just in case a fresh iblank appears in the current process
      IF (SUM(iblank) .EQ. 0) CALL IFEM_SETIBLANK(lD)

      ALLOCATE(incNd(tnNo))
      DO a=1, tnNo
         incNd(a) = iblank(a)
      END DO

!     Do multiple sweeps to get the neighborhood of iblank field. Probe
!     is performed only on these neighborhood nodes.
      iswp = 0
      DO WHILE (iswp .LE. 2)
         DO iM=1, nMsh
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               IF (incNd(Ac) .EQ. 1) THEN
                  DO i=msh(iM)%nAdj%prow(a), msh(iM)%nAdj%prow(a+1)-1
                     b = msh(iM)%nAdj%pcol(i)
                     Bc = msh(iM)%gN(b)
                     incNd(Bc) = 1
                  END DO
               END IF
            END DO
         END DO
         iswp = iswp + 1
      END DO

      a = 0
      DO Ac=1, tnNo
         IF (incNd(Ac) .EQ. 1) a = a + 1
      END DO
      nNo = a

!     Transfer all the nodes to be probed to a local pointer array
      ALLOCATE(ptr(nNo))
      a = 0
      DO Ac=1, tnNo
         IF (incNd(Ac) .EQ. 1) THEN
            a = a + 1
            ptr(a) = Ac
         END IF
      END DO

!     Probe each node if it is inside or outside the immersed body
      iblank = 0
      DO a=1, nNo
         Ac = ptr(a)
         xp = x(:,Ac) + lD(:,Ac)
         flag  = .FALSE.
         CALL IFEM_CHECKINOUT(xp, flag)
         IF (flag) iblank(Ac) = 1
      END DO
      DEALLOCATE(incNd, ptr)

      RETURN
      END SUBROUTINE IFEM_UPDATEIBLANK
!--------------------------------------------------------------------
!     Checks if a probe lies inside or outside an immersed boundary
      SUBROUTINE IFEM_CHECKINOUT(xp, flag)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: xp(nsd)
      LOGICAL, INTENT(INOUT) :: flag

      INTEGER(KIND=IKIND) :: a, e, Ac, Ec, iM, iFa, jM, jFa
      REAL(KIND=RKIND) :: dS, minS, Jac, nV(nsd), xb(nsd), dotP

!     Find the closest immersed face centroid from the probe
      minS = HUGE(minS)
      DO iM=1, ifem%nMsh
         DO iFa=1, ifem%msh(iM)%nFa
            DO e=1, ifem%msh(iM)%fa(iFa)%nEl
               xb = 0._RKIND
               DO a=1, ifem%msh(iM)%fa(iFa)%eNoN
                  Ac = ifem%msh(iM)%fa(iFa)%IEN(a,e)
                  xb = xb + ifem%x(:,Ac) + ifem%Ubo(:,Ac)
               END DO
               xb = xb / REAL(ifem%msh(iM)%fa(iFa)%eNoN, KIND=RKIND)
               dS = SQRT( SUM( (xp(:)-xb(:))**2._RKIND ) )
               IF (dS .LT. minS) THEN
                  minS = dS
                  Ec   = e
                  jM   = iM
                  jFa  = iFa
               END IF
            END DO
         END DO
      END DO

!     Compute the normal of the face element
      CALL GNNIB(ifem%msh(jM)%fa(jFa), Ec, 1, nV)
      Jac = SQRT(NORM(nV))
      nV  = nV(:)/Jac

!     Check the sign of op.n
      xb = 0._RKIND
      DO a=1, ifem%msh(jM)%fa(jFa)%eNoN
         Ac = ifem%msh(jM)%fa(jFa)%IEN(a,Ec)
         xb = xb + ifem%x(:,Ac) + ifem%Ubo(:,Ac)
      END DO
      xb   = xb / REAL(ifem%msh(jM)%fa(jFa)%eNoN, KIND=RKIND)
      dotP = NORM(xp-xb, nV)

      IF (dotP .LT. -1.E-9_RKIND) THEN
!        probe lies inside IFEM
         flag = .TRUE.
      ELSE IF (dotP .GT. 1.E-9_RKIND) THEN
!        probe lies outside IFEM
         flag = .FALSE.
      ELSE
!        if the probe is along the tangent, perform sign check with one
!        of the vertices of the closest node instead of the centroid
         Ac   = ifem%msh(jM)%fa(jFa)%IEN(1,Ec)
         xb   = ifem%x(:,Ac) + ifem%Ubo(:,Ac)
         dotP = NORM(xp-xb,nV)
         IF (ABS(dotP) .LT. 1.E-9_RKIND) THEN
            flag = .TRUE.
         ELSE
            flag = .FALSE.
         END IF
      END IF

      RETURN
      END SUBROUTINE IFEM_CHECKINOUT
!####################################################################
!     Find closest fluid node for each solid node, TODO for parallel
      SUBROUTINE IFEM_FINDCLOSEST(lM, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM ! ifem mesh 
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      INTEGER(KIND=IKIND) :: a, e, i, j, iM, Ac, Acn
      REAL(KIND=RKIND), ALLOCATABLE :: xSCur(:,:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: elmF(:), nodeF(:)
      REAL(KIND=RKIND), ALLOCATABLE :: poly(:,:)
      REAL(KIND=RKIND) :: xp(nsd), xs(nsd), minb(nsd), maxb(nsd)
      REAL(KIND=RKIND) :: diam, maxDist
      LOGICAL :: flag = .FALSE.
      LOGICAL :: find = .FALSE.      


!     Update current solid (from Lag to Eul) and fluid location (in case of ALE)
      ALLOCATE(xSCur(nsd,lM%nNo))
      ALLOCATE(elmF(lM%nNo))
      ALLOCATE(nodeF(lM%nNo))
      ALLOCATE(poly(nsd,lM%eNoN))

      xsCur = 0._RKIND
      elmF = 0
      nodeF = 0

      DO a=1, lM%nNo
         Ac = lM%gN(a)
         ! ?? check this for parallel 
C          xSCur(:,a) = ifem%x(:,Ac) + ifem%Ubo(:,Ac)
         xSCur(:,a) = ifem%x(:,a) + ifem%Ubo(:,a)
      END DO

C       write(*,*)"The solid ifem%Ubo is: ", ifem%Ubo
C       write(*,*)"The solid xSCur is: ", xSCur

      maxDist = TINY(maxDist)
!     compute solid mesh diamter (in the current configuration)
      DO e=1, lM%nEl

         Ac = lM%IEN(1,e)
         Acn = lM%IEN(2,e)
         diam = DIST(ifem%x(:,Ac),ifem%x(:,Acn))
         IF(diam .GT. maxDist) maxDist = diam

         Ac = lM%IEN(2,e)
         Acn = lM%IEN(3,e)
         diam = DIST(ifem%x(:,Ac),ifem%x(:,Acn))
         IF(diam .GT. maxDist) maxDist = diam

         Ac = lM%IEN(1,e)
         Acn = lM%IEN(3,e)
         diam = DIST(ifem%x(:,Ac),ifem%x(:,Acn))
         IF(diam .GT. maxDist) maxDist = diam

      END DO

C       write(*,*)"The solid diameter is : ", maxDist

!     Create a bounding box around of the current solid location 
      minb = HUGE(minb)
      maxb = TINY(maxb)

      DO i=1, nsd
         minb(i) = MINVAL(xSCur(i,:)) - maxDist
         maxb(i) = MAXVAL(xSCur(i,:)) + maxDist
      END DO

C       write(*,*)"The solid BBOX x-dir is: ", minb(1), " and ", maxb(1)
C       write(*,*)"The solid BBOX y-dir is: ", minb(2), " and ", maxb(2)

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
      DO iM=1, nMsh
         DO e=1, msh(iM)%nEl ! fluid elem id
C             write(*,*)"Testing fluid element: ", e
            flag = .FALSE.

            DO a=1, msh(iM)%eNoN

               Ac = msh(iM)%IEN(a,e)

               xp(:) = x(:,Ac) + lD(:,Ac)

!              is the fluid element in the solid BBox? 
               IF((xp(1) .LE. maxb(1)) .AND. (xp(1) .GE. minb(1)) 
     2                   .AND. (xp(2) .LE. maxb(2)) 
     3                   .AND. (xp(2) .GE. minb(2))) THEN 
!                 The node is inside the Bounding Box
                  flag = .TRUE.
                  EXIT
               END IF 
            END DO

            IF (flag) THEN
!              Loop over the solid node, and search if is inside the fluid element 
C                write(*,*) "Element ", e, " inside BBox, with vertices: "

               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%IEN(a,e)
                  poly(:,a) = x(:,Ac) + lD(:,Ac)
C                   write(*,*) poly(:,a)
               END DO  

               DO a=1, lM%nNo
                  xs = xSCur(:,a)

                  find = IN_POLY(xs,poly)
C                   write(*,*)"find node solid ", a, " is ", find

!                 is the solid node in this fluis element 
                  IF (find) THEN 
                     elmF(a) = e 
!                    If inside, serach for the closest point   
!                    local (element) id of the closest point
                     nodeF(a) = CLOSEST_POINT(xs,poly) 
                     nodeF(a) = msh(iM)%IEN(nodeF(a),e) 
C                      write(*,*)"Closes id is ", msh(iM)%IEN(nodeF(a),e)
C                      EXIT
                  END IF
               END DO
            END IF

         END DO
      END DO

      ALLOCATE(ifem%clsFElm(lM%nNo))
      ALLOCATE(ifem%clsFNd(lM%nNo))

      ifem%clsFNd = nodeF
      ifem%clsFElm = elmF

C       DO a=1, lM%nNo
C          e = a+1595
C          write(*,*)"closest solid id (mergefile notation)", e,  
C      2    " is fluid id ", ifem%clsFNd(a), " in fluid elem ", 
C      3    ifem%clsFElm(a)
C       END DO

C       write(*,*) "solid nodes are into flui elem: ", ifem%clsFElm 
C       write(*,*) "closest local id is : ", ifem%clsFNd

      DEALLOCATE(xSCur)
      DEALLOCATE(elmF)
      DEALLOCATE(nodeF)
      DEALLOCATE(poly)

      RETURN
      END SUBROUTINE IFEM_FINDCLOSEST
!####################################################################
!     Update closest fluid node for each solid node, searching between the 
!     fluid elem's neighbors
      SUBROUTINE IFEM_UPDATECLS(lM, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM ! ifem mesh 
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      INTEGER(KIND=IKIND) :: snd, a, js, iM, Ac, jf
      REAL(KIND=RKIND), ALLOCATABLE :: poly(:,:)

      REAL(KIND=RKIND) :: xs(nsd)
      INTEGER(KIND=IKIND) :: nodeF, idFEl, idFNd, idNFEl, nbSE
      LOGICAL :: find = .FALSE.      


!     Update current solid (from Lag to Eul) and fluid location (in case of ALE)
      ALLOCATE(poly(nsd,lM%eNoN))

      nodeF = 0

      iM = 1

!     Check if each solid node is still inside the same fluid element
      DO snd=1, ifem%tnNo
!        Find fluid elem in which is it in the def config 
         idFEl = ifem%clsFElm(snd)
         idFNd = ifem%clsFNd(snd)

!        Built poly for inside search
         DO a=1, msh(iM)%eNoN
            Ac = msh(iM)%IEN(a,idFEl)
            poly(:,a) = x(:,Ac) + lD(:,Ac)
         END DO  

!        Get solid node coord in the current config
         xs = ifem%xCu(:,snd)
         find = IN_POLY(xs,poly)
C          write(*,*) "find node solid ", snd, " is ", find

!        Is the solid node in this fluis element? yes, nothing to do
         IF (.NOT. find) THEN 
!           Search in the elm stencil of idFEl   
            nbSE = SIZE(msh(iM)%stn%elmStn,2)
C             write(*,*)"Begin search neigh" , nbSE
            
            DO js = 1, nbSE
C                write(*,*)" id neigh is ", msh(iM)%stn%elmStn(idFNd,:)
               idNFEl = msh(iM)%stn%elmStn(idFNd,js)

               IF (idNFEl .EQ. 0) EXIT

               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%IEN(a,idNFEl)
                  poly(:,a) = x(:,Ac) + lD(:,Ac)
               END DO 

               find = IN_POLY(xs,poly)

               IF (find) THEN
C                   write(*,*)"Find node in elem ", idNFEl
                  ifem%clsFElm(snd) = idNFEl

                  nodeF = CLOSEST_POINT(xs,poly) 
                  ifem%clsFNd(snd) = msh(iM)%IEN(nodeF,idNFEl) 
                  EXIT
               END IF   
            END DO

            IF (.NOT. find) THEN 
               ! loop over all the fluid elem 
               DO idNFEl = 1, msh(iM)%nEl

                  DO a=1, msh(iM)%eNoN
                     Ac = msh(iM)%IEN(a,idNFEl)
                     poly(:,a) = x(:,Ac) + lD(:,Ac)
                  END DO 

                  find = IN_POLY(xs,poly)

                  IF (find) THEN
C                      write(*,*)"Find node in elem ", idNFEl
                     ifem%clsFElm(snd) = idNFEl

                     nodeF = CLOSEST_POINT(xs,poly) 
                     ifem%clsFNd(snd) = msh(iM)%IEN(nodeF,idNFEl) 
                     EXIT
                  END IF   
               END DO
            END IF

            IF (.NOT. find) THEN 
               write(*,*)"ERROR FIND CLSEST PNT NEIGH" 
               write(*,*)"OHHH NOOOO!"
               CALL EXIT(1)
            END IF

         END IF
      
      END DO

      DEALLOCATE(poly)

      RETURN
      END SUBROUTINE IFEM_UPDATECLS
!--------------------------------------------------------------------
!     Find traces of IFEM nodes and integration points
      SUBROUTINE IFEM_FINDTRACES(lM, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

!     Find nodal traces
      CALL IFEM_FINDNDTRACES(lM, lD)

!     Find parametric coordinate for nodal traces
      CALL IFEM_FINDXINDTRC(lM, lD)

!     Find Gauss point traces
      CALL IFEM_FINDGPTRACES(lM, lD)

!     Find parametric coordinate for Gauss point traces
      CALL IFEM_FINDXIGPTRC(lM, lD)

      RETURN
      END SUBROUTINE IFEM_FINDTRACES
!--------------------------------------------------------------------
!     Find traces of IFEM nodes on the background mesh elements and those
!     elements are tagged as ghost elements.
      SUBROUTINE IFEM_FINDNDTRACES(lM, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      INTEGER(KIND=IKIND) :: a, b, e, i, j, iM, Ac, Ec, nNe, ne, ierr
      REAL(KIND=RKIND) :: xp(nsd), xi(nsd), minb(nsd), maxb(nsd), tt
      TYPE(queueType) :: probeNdQ

      LOGICAL, ALLOCATABLE :: ichk(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: nPtr(:,:), incNd(:), eList(:),
     2   masEList(:), part(:), rootEl(:), sCount(:), disps(:), ptr(:),
     3   gptr(:), tmpI(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xpL(:,:)

!     Create a list of all the mesh nodes
      ALLOCATE(xpL(nsd,lM%nNo))
      xpL = 0._RKIND
      DO a=1, lM%nNo
         Ac = lM%gN(a)
         xpL(:,a) = ifem%x(:,Ac) + ifem%Ubo(:,Ac)
      END DO

!     Create a bounding box around the intersection of immersed body
!     and background mesh
      minb = HUGE(minb)
      maxb = TINY(maxb)
      DO iM=1, nMsh
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            xp(:) = x(:,Ac) + lD(:,Ac)
            DO i=1, nsd
               IF (minb(i) .GT. xp(i)) minb(i) = xp(i)
               IF (maxb(i) .LT. xp(i)) maxb(i) = xp(i)
            END DO
         END DO
      END DO

      DO i=1, nsd
         minb(i) = MAX(MINVAL(xpL(i,:)), minb(i))
         maxb(i) = MIN(MAXVAL(xpL(i,:)), maxb(i))
      END DO
      minb(:) = minb(:) - lM%dx
      maxb(:) = maxb(:) + lM%dx

!     Loop over all possible background mesh and find traces of each IFEM
!     node onto the elements of the background mesh.
      ALLOCATE(nPtr(2,lM%nNo))
      nPtr(:,:) = 0
      DO iM=1, nMsh
!        Initialize all the mesh nodes whose traces are to be found
         ALLOCATE(incNd(lM%nNo))
         incNd(:) = 1

!        Identify the IFEM nodes which are owned by the current process.
!        Search for traces is performed only these elements
         DO a=1, lM%nNo
            xp = xpL(:,a)
            j  = 0
            DO i=1, nsd
               IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2            j = j + 1
            END DO
            IF (j .EQ. nsd) incNd(a) = 0
         END DO

!        Create a master list of elements of the background mesh
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
!           Reset eList if the there is no overlap with iblank field
            b = 0
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               b  = b + iblank(Ac)
            END DO
            IF (b .EQ. 0) eList(e) = 0
         END DO

!        Save to master list
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         ne = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               ne = ne + 1
               masElist(ne) = e
            END IF
         END DO

         IF (.NOT.cm%seq()) THEN
            ALLOCATE(part(lM%nNo))
            CALL IFEM_PARTMSH(lM, msh(iM), lD, eList, part)
            DO e=1, lM%nEl
               b = 0
               DO a=1, lM%eNoN
                  Ac = lM%IEN(a,e)
                  Ac = lM%lN(Ac)
                  IF (part(Ac) .NE. cm%tF()) b = b + 1
               END DO
               IF (b .EQ. lM%eNoN) THEN
                  DO a=1, lM%eNoN
                     Ac = lM%IEN(a,e)
                     Ac = lM%lN(Ac)
                     incNd(Ac) = 0
                  END DO
               END IF
            END DO
            DEALLOCATE(part)
         END IF
         DEALLOCATE(eList)
         IF (SUM(incNd) .EQ. 0) THEN
            DEALLOCATE(incNd, masEList)
            CYCLE
         END IF

!        At this stage, only the processes involved in search continue.
!        We have also identified the IFEM nodes local to the process and
!        prepared a master list of background mesh elements involved in
!        the trace search process.
         ALLOCATE(ichk(lM%nNo), rootEl(lM%nNo))
         ichk   = .FALSE.
         rootEl = 0

!        Identify a seed IFEM node using master element trace search and
!        the corresponding trace element is chosen as root element.
!        Identify the neighbors of the seed node and form a queue for
!        subsequent search. The root element for the seed node is
!        assigned as an initial search element for the que'd neighboring
!        nodes
         Ec = 0
         E_LOOP: DO a=1, lM%nNo
            ichk(a) = .TRUE.
            IF (incNd(a) .EQ. 0) CYCLE
            xp = xpL(:,a)
            CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec, xi)
            IF (Ec .EQ. 0) CYCLE
            ichk(a) = .FALSE.
            EXIT E_LOOP
         END DO E_LOOP

         IF (ALL(ichk(:))) THEN
            DEALLOCATE(incNd, masElist, ichk, rootEl)
            CYCLE
         END IF

         nPtr(1,a) = Ec
         nPtr(2,a) = iM
         rootEl(a) = Ec
         DO i=lM%nAdj%prow(a), lM%nAdj%prow(a+1)-1
            b = lM%nAdj%pcol(i)
            IF (incNd(b) .EQ. 1) THEN
               CALL ENQUEUE(probeNdQ, b)
               rootEl(b) = Ec
            END IF
         END DO

!        The nonlinear grid-grid search begins here.
         DO WHILE (DEQUEUE(probeNdQ, a))
            IF (ALL(ichk(:))) EXIT
            IF (ichk(a)) CYCLE
            ichk(a) = .TRUE.

            Ec = rootEl(a)
            ne = msh(iM)%eAdj%prow(Ec+1) - msh(iM)%eAdj%prow(Ec)
            IF (ALLOCATED(eList)) DEALLOCATE(eList)
            ALLOCATE(eList(ne))
            j = 0
            DO i=msh(iM)%eAdj%prow(Ec), msh(iM)%eAdj%prow(Ec+1)-1
               j = j + 1
               eList(j) = msh(iM)%eAdj%pcol(i)
            END DO

            xp = xpL(:,a)
            CALL FINDE(xp, msh(iM), x, lD, tnNo, ne, eList, Ec, xi)

!           If a trace is not found, then include the neighborhood
!           elements for search. Repeat this twice. If not found yet,
!           search using master element list. If master list search
!           also fails, continue
            IF (Ec .EQ. 0) THEN
               CALL IFEM_FPSRCH(xp, msh(iM), lD, ne, eList, 2, Ec, xi)
               IF (Ec .EQ. 0) THEN
                  CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist,
     2               Ec, xi)
                  IF (Ec .EQ. 0) CYCLE
               END IF
            END IF

!           At this point, the trace is definitely found. Save it.
            nPtr(1,a) = Ec
            nPtr(2,a) = iM

!           Use the trace as a seed search element for the neigboring
!           IFEM elements in the following trace search
            rootEl(a)  = Ec
            DO i=lM%nAdj%prow(a), lM%nAdj%prow(a+1)-1
               b = lM%nAdj%pcol(i)
               IF (incNd(b).EQ.1 .AND. .NOT.ichk(b)) THEN
                  CALL ENQUEUE(probeNdQ, b)
                  rootEl(b) = Ec
               END IF
            END DO
         END DO
         IF (ALLOCATED(eList)) DEALLOCATE(eList)
         CALL DESTROY(probeNdQ)

!        Perform a brute search on any missed elements
         DO a=1, lM%nNo
            IF (.NOT.ichk(a) .AND. incNd(a).EQ.1) THEN
               xp = xpL(:,a)
               CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec,
     2            xi)
               IF (Ec .EQ. 0) CYCLE
               nPtr(1,a) = Ec
               nPtr(2,a) = iM
            END IF
         END DO

         DEALLOCATE(incNd, masElist, ichk, rootEl)
      END DO

!     Transfer nPtr to trace data structure
      IF (ALLOCATED(lM%trc%gN))   DEALLOCATE(lM%trc%gN)
      IF (ALLOCATED(lM%trc%nptr)) DEALLOCATE(lM%trc%nptr)
      i = 0
      DO a=1, lM%nNo
         IF (nPtr(1,a) .NE. 0) i = i + 1
      END DO
      lM%trc%n = i
      ALLOCATE(lM%trc%gN(i), lM%trc%nptr(2,i))
      i = 0
      DO a=1, lM%nNo
         IF (nPtr(1,a) .NE. 0) THEN
            i = i + 1
            Ec = nPtr(1,a)
            iM = nPtr(2,a)
            lM%trc%gN(i)     = a
            lM%trc%nptr(:,i) = nPtr(:,a)
            msh(iM)%iGC(Ec)  = 1
         END IF
      END DO

!     Check if all traces are found and return if yes
      ALLOCATE(ptr(lM%nNo))
      ptr = 0
      DO a=1, lM%nNo
         ptr = nPtr(1,a)
      END DO

      tt = CPUT()
      gptr = cm%reduce(ptr, MPI_MAX)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      i = 0
      DO a=1, lM%nNo
         IF (gptr(a) .NE. 0) i = i + 1
      END DO
      IF (i .EQ. lM%nNo) THEN
         DEALLOCATE(nPtr, ptr, gptr, xpL)
         RETURN
      END IF

!     If all traces are not found, add a foolproof search on master
!     list. If the search fails here, the element is perhaps distorted
!     and hence, the simulation aborts.
!     First, create a list of all successfully found traces on master
      ALLOCATE(sCount(cm%np()), disps(cm%np()))
      sCount = 0
      disps  = 0
      i  = lM%trc%n

      tt = CPUT()
      CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
     2   ierr)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      j = SUM(sCount(:))
      sCount = 3*sCount(:)
      DO i=2, cm%np()
         disps(i) = disps(i-1) + sCount(i-1)
      END DO

      DEALLOCATE(ptr, gptr)
      ALLOCATE(ptr(3*lM%trc%n), gptr(3*j))
      DO i=1, lM%trc%n
         ptr(3*i-2) = lM%trc%gN(i)
         ptr(3*i-1) = lM%trc%nptr(1,i)
         ptr(3*i)   = lM%trc%nptr(2,i)
      END DO

      i  = 3*lM%trc%n

      tt = CPUT()
      CALL MPI_GATHERV(ptr, i, mpint, gptr, sCount, disps, mpint,
     2   master, cm%com(), ierr)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      DEALLOCATE(ptr, disps, sCount)

      IF (cm%mas()) THEN
         ALLOCATE(tmpI(2,lM%nNo), ptr(lM%nNo))
         tmpI = 0
         ptr  = 0
         DO i=1, j
            a  = gptr(3*i-2)
            Ec = gptr(3*i-1)
            iM = gptr(3*i)
            tmpI(1,a) = Ec
            tmpI(2,a) = iM
         END DO

         DO a=1, lM%nNo
            Ec = tmpI(1,a)
            iM = tmpI(2,a)
            IF (Ec.EQ.0 .OR. iM.EQ.0) ptr(a) = 1
         END DO

         j = SUM(ptr)
         ALLOCATE(incNd(j))
         i = 0
         DO a=1, lM%nNo
            IF (ptr(a) .EQ. 1) THEN
               i = i + 1
               incNd(i) = a
            END IF
         END DO
         DEALLOCATE(tmpI, ptr)
      END IF
      DEALLOCATE(gptr)

!     Share the element list to all processes
      tt = CPUT()
      CALL cm%bcast(j)
      IF (cm%slv()) ALLOCATE(incNd(j))
      CALL cm%bcast(incNd)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

!     Loop over all the background mesh and find traces of each left out
!     integration point onto the elements of the background mesh
      ALLOCATE(ichk(j))
      ichk = .FALSE.
      DO iM=1, nMsh
!        Create a master list of elements of the background mesh based
!        on IFEM bounding box
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
         END DO
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         i = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               i = i + 1
               masElist(i) = e
            END IF
         END DO

!        Now perform search for each integration point of an element
!        whose trace was not determined earlier
         DO i=1, j
            IF (ichk(i)) CYCLE
            a = incNd(i)
            xp = xpL(:,a)
            CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec, xi)
            IF (Ec .NE. 0) THEN
               nPtr(1,a) = Ec
               nPtr(2,a) = iM
               ichk(i) = .TRUE.
            END IF
         END DO
         DEALLOCATE(eList, masElist)
      END DO
      DEALLOCATE(incNd, ichk, xpL)

!     Transfer ePtr to trace data structure
      IF (ALLOCATED(lM%trc%gN))   DEALLOCATE(lM%trc%gN)
      IF (ALLOCATED(lM%trc%nptr)) DEALLOCATE(lM%trc%nptr)
      i = 0
      DO a=1, lM%nNo
         IF (nPtr(1,a) .NE. 0) i = i + 1
      END DO
      lM%trc%n = i
      ALLOCATE(lM%trc%gN(i), lM%trc%nptr(2,i))
      i = 0
      DO a=1, lM%nNo
         IF (nPtr(1,a) .NE. 0) THEN
            i = i + 1
            Ec = nPtr(1,a)
            iM = nPtr(2,a)
            lM%trc%gN(i)     = a
            lM%trc%nptr(:,i) = nPtr(:,a)
            msh(iM)%iGC(Ec) = 1
         END IF
      END DO

!     Abort simulation if all traces are still not found
!     Check if all traces are found and return if yes
      ALLOCATE(ptr(lM%nNo))
      ptr = 0
      DO a=1, lM%nNo
         ptr = nPtr(1,a)
      END DO

      tt = CPUT()
      gptr = cm%reduce(ptr, MPI_MAX)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      i = 0
      DO a=1, lM%nNo
         IF (gptr(a) .NE. 0) i = i + 1
      END DO
      DEALLOCATE(nPtr, ptr, gptr)
      IF (i .EQ. lM%nNo) RETURN

      CALL DEBUGIBNDTRCS(lM)

      END SUBROUTINE IFEM_FINDNDTRACES
!--------------------------------------------------------------------
!     Find traces of IFEM integration points on the background mesh
!     elements and those elements are tagged as ghost elements.
      SUBROUTINE IFEM_FINDGPTRACES(lM, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      INTEGER(KIND=IKIND) :: a, b, e, g, i, j, iM, Ac, Ec, nNe, ne, ierr
      REAL(KIND=RKIND) :: xp(nsd), xi(nsd), minb(nsd), maxb(nsd), tt
      TYPE(queueType) :: probeElQ

      LOGICAL, ALLOCATABLE :: ichk(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: ePtr(:,:,:), incEl(:),
     2   eList(:), masEList(:), part(:), rootEl(:), sCount(:), disps(:),
     3   ptr(:), gptr(:), tmpI(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), xpL(:,:)

!     Create a list of all the mesh nodes
      ALLOCATE(xpL(nsd,lM%nNo))
      xpL = 0._RKIND
      DO a=1, lM%nNo
         Ac = lM%gN(a)
         xpL(:,a) = ifem%x(:,Ac) + ifem%Ubo(:,Ac)
      END DO

!     Create a bounding box around the intersection of immersed body
!     and background mesh
      minb = HUGE(minb)
      maxb = TINY(maxb)
      DO iM=1, nMsh
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            xp(:) = x(:,Ac) + lD(:,Ac)
            DO i=1, nsd
               IF (minb(i) .GT. xp(i)) minb(i) = xp(i)
               IF (maxb(i) .LT. xp(i)) maxb(i) = xp(i)
            END DO
         END DO
      END DO

      DO i=1, nsd
         minb(i) = MAX(MINVAL(xpL(i,:)), minb(i))
         maxb(i) = MIN(MAXVAL(xpL(i,:)), maxb(i))
      END DO
      minb(:) = minb(:) - lM%dx
      maxb(:) = maxb(:) + lM%dx

!     Loop over all possible background mesh and find traces of each IFEM
!     integration point onto the elements of the background mesh.
      ALLOCATE(ePtr(2,lM%nG,lM%nEl))
      ePtr = 0
      DO iM=1, nMsh
!        Initialize all the mesh elements whose traces are to be found
         ALLOCATE(incEl(lM%nEl))
         incEl(:) = 1

!        Identify the IFEM elements which are owned by the current
!        process. Search for traces is performed only these elements
         DO e=1, lM%nEl
            b = 0
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               Ac = lM%lN(Ac)
               xp = xpL(:,Ac)
               j = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               j = j + 1
               END DO
               IF (j .EQ. nsd) b = b + 1
            END DO
            IF (b .LT. lM%eNoN) incEl(e) = 0
         END DO

!        Create a master list of elements of the background mesh
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
!           Reset eList if the there is no overlap with iblank field
            b = 0
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               b  = b + iblank(Ac)
            END DO
            IF (b .EQ. 0) eList(e) = 0
         END DO

!        Save to master list
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         ne = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               ne = ne + 1
               masElist(ne) = e
            END IF
         END DO

         IF (.NOT.cm%seq()) THEN
            ALLOCATE(part(lM%nNo))
            CALL IFEM_PARTMSH(lM, msh(iM), lD, eList, part)
            DO e=1, lM%nEl
               b = 0
               DO a=1, lM%eNoN
                  Ac = lM%IEN(a,e)
                  Ac = lM%lN(Ac)
                  IF (part(Ac) .NE. cm%tF()) b = b + 1
               END DO
               IF (b .EQ. lM%eNoN) incEl(e) = 0
            END DO
            DEALLOCATE(part)
         END IF
         DEALLOCATE(eList)
         IF (SUM(incEl) .EQ. 0) THEN
            DEALLOCATE(incEl, masEList)
            CYCLE
         END IF

!        At this stage, only the processes involved in search continue.
!        We have also identified the IFEM nodes local to the process and
!        prepared a master list of background mesh elements involved in
!        the trace search process.
         ALLOCATE(ichk(lM%nEl), rootEl(lM%nEl))
         ichk   = .FALSE.
         rootEl = 0

!        Identify a seed element using master element trace search and
!        the corresponding trace element is chosen as root element.
!        Identify the neighbors of the seed element and form a queue for
!        subsequent search. The root element for the seed element is
!        assigned as a seed search element for the queued neighboring
!        elements
         ALLOCATE(xl(nsd,lM%eNoN))
         Ec = 0
         E_LOOP: DO e=1, lM%nEl
            ichk(e) = .TRUE.
            IF (incEl(e) .EQ. 0) CYCLE
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               Ac = lM%lN(Ac)
               xl(:,a) = xpL(:,Ac)
            END DO
            DO g=1, lM%nG
               xp = 0._RKIND
               DO a=1, lM%eNoN
                  xp = xp + xl(:,a)*lM%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec,
     2            xi)
               IF (Ec .EQ. 0) CYCLE
               ichk(e) = .FALSE.
               EXIT E_LOOP
            END DO
         END DO E_LOOP

         IF (ALL(ichk(:))) THEN
            DEALLOCATE(incEl, masElist, ichk, rootEl, xl)
            CYCLE
         END IF

         ePtr(1,g,e) = Ec
         ePtr(2,g,e) = iM
         rootEl(e)   = Ec
         DO i=lM%eAdj%prow(e), lM%eAdj%prow(e+1)-1
            j = lM%eAdj%pcol(i)
            IF (incEl(j) .EQ. 1) THEN
               CALL ENQUEUE(probeElQ, j)
               rootEl(j) = Ec
            END IF
         END DO

!        The nonlinear grid-grid search begins here.
         DO WHILE (DEQUEUE(probeElQ, e))
            IF (ALL(ichk(:))) EXIT
            IF (ichk(e)) CYCLE
            ichk(e) = .TRUE.

            DO g=1, lM%nG
               DO a=1, lM%eNoN
                  Ac = lM%IEN(a,e)
                  Ac = lM%lN(Ac)
                  xl(:,a) = xpL(:,Ac)
               END DO

               xp = 0._RKIND
               DO a=1, lM%eNoN
                  xp = xp + xl(:,a)*lM%N(a,g)
               END DO

               Ec = rootEl(e)
               ne = msh(iM)%eAdj%prow(Ec+1) - msh(iM)%eAdj%prow(Ec)
               IF (ALLOCATED(eList)) DEALLOCATE(eList)
               ALLOCATE(eList(ne))
               j = 0
               DO i=msh(iM)%eAdj%prow(Ec), msh(iM)%eAdj%prow(Ec+1)-1
                  j = j + 1
                  eList(j) = msh(iM)%eAdj%pcol(i)
               END DO

               CALL FINDE(xp, msh(iM), x, lD, tnNo, ne, eList, Ec, xi)

!              If a trace is not found, then include the neighborhood
!              elements for search. Repeat this twice. If not found yet,
!              search using master element list. If master list search
!              also fails, continue
               IF (Ec .EQ. 0) THEN
                  CALL IFEM_FPSRCH(xp,msh(iM),lD,ne,eList,2,Ec,xi)
                  IF (Ec .EQ. 0) THEN
                     CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist,
     2                  Ec, xi)
                     IF (Ec .EQ. 0) CYCLE
                  END IF
               END IF

!              At this point, the trace is definitely found. Save it.
               ePtr(1,g,e) = Ec
               ePtr(2,g,e) = iM

!              Use the trace as a seed search element for the neigboring
!              IFEM elements in the following trace search
               rootEl(e)  = Ec
               DO i=lM%eAdj%prow(e), lM%eAdj%prow(e+1)-1
                  j = lM%eAdj%pcol(i)
                  IF (incEl(j).EQ.1 .AND. .NOT.ichk(j)) THEN
                     CALL ENQUEUE(probeElQ, j)
                     rootEl(j) = Ec
                  END IF
               END DO
            END DO
         END DO
         IF (ALLOCATED(eList)) DEALLOCATE(eList)
         CALL DESTROY(probeElQ)

!        Perform a brute search on any missed elements
         DO e=1, lM%nEl
            IF (.NOT.ichk(e) .AND. incEl(e).EQ.1) THEN
               DO g=1, lM%nG
                  DO a=1, lM%eNoN
                     Ac = lM%IEN(a,e)
                     Ac = lM%lN(Ac)
                     xl(:,a) = xpL(:,Ac)
                  END DO

                  xp = 0._RKIND
                  DO a=1, lM%eNoN
                     xp = xp + xl(:,a)*lM%N(a,g)
                  END DO

                  CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist,
     2               Ec, xi)
                  IF (Ec .EQ. 0) CYCLE
                  ePtr(1,g,e) = Ec
                  ePtr(2,g,e) = iM
               END DO
            END IF
         END DO

         DEALLOCATE(incEl, masElist, ichk, rootEl, xl)
      END DO

!     Transfer ePtr to trace data structure
      IF (ALLOCATED(lM%trc%gE))   DEALLOCATE(lM%trc%gE)
      IF (ALLOCATED(lM%trc%gptr)) DEALLOCATE(lM%trc%gptr)
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) i = i + 1
         END DO
      END DO
      lM%trc%nG = i
      ALLOCATE(lM%trc%gE(2,i), lM%trc%gptr(2,i))
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) THEN
               i = i + 1
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               lM%trc%gE(1,i)   = e
               lM%trc%gE(2,i)   = g
               lM%trc%gptr(:,i) = ePtr(:,g,e)
               msh(iM)%iGC(Ec)  = 1
            END IF
         END DO
      END DO

!     Check if all traces are found and return if yes
      tt = CPUT()
      i  = cm%reduce(lM%trc%n)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      j = lM%nEl*lM%nG
      IF (i .GE. j) THEN
         DEALLOCATE(xpL, ePtr)
         RETURN
      END IF

!     If all traces are not found, add a fool-proof search on master
!     list. If the search fails here, the element is perhaps distorted
!     and hence, the simulation aborts.
!     First, create a list of all successfully found traces on master
      ALLOCATE(sCount(cm%np()), disps(cm%np()))
      sCount = 0
      disps  = 0
      i  = lM%trc%nG

      tt = CPUT()
      CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
     2   ierr)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      j = SUM(sCount(:))
      sCount = 4*sCount(:)
      DO i=2, cm%np()
         disps(i) = disps(i-1) + sCount(i-1)
      END DO

      ALLOCATE(ptr(4*lM%trc%nG), gptr(4*j))
      DO i=1, lM%trc%nG
         ptr(4*i-3) = lM%trc%gE(1,i)
         ptr(4*i-2) = lM%trc%gE(2,i)
         ptr(4*i-1) = lM%trc%gptr(1,i)
         ptr(4*i)   = lM%trc%gptr(2,i)
      END DO

      i  = 4*lM%trc%nG
      tt = CPUT()
      CALL MPI_GATHERV(ptr, i, mpint, gptr, sCount, disps, mpint,
     2   master, cm%com(), ierr)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      DEALLOCATE(ptr, disps, sCount)

      IF (cm%mas()) THEN
         ALLOCATE(tmpI(2,lM%nG,lM%nEl), ptr(lM%nEl))
         tmpI = 0
         ptr  = 0
         DO i=1, j
            e  = gptr(4*i-3)
            g  = gptr(4*i-2)
            Ec = gptr(4*i-1)
            iM = gptr(4*i)
            tmpI(1,g,e) = Ec
            tmpI(2,g,e) = iM
         END DO

         DO e=1, lM%nEl
            DO g=1, lM%nG
               Ec = tmpI(1,g,e)
               iM = tmpI(2,g,e)
               IF (Ec.EQ.0 .OR. iM.EQ.0) ptr(e) = 1
            END DO
         END DO

         ne = SUM(ptr)
         ALLOCATE(incEl(ne))
         i = 0
         DO e=1, lM%nEl
            IF (ptr(e) .EQ. 1) THEN
               i = i + 1
               incEl(i) = e
            END IF
         END DO
         DEALLOCATE(tmpI, ptr)
      END IF
      DEALLOCATE(gptr)

!     Share the element list to all processes
      tt = CPUT()
      CALL cm%bcast(ne)
      IF (cm%slv()) ALLOCATE(incEl(ne))
      CALL cm%bcast(incEl)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

!     Loop over all the background mesh and find traces of each left out
!     integration point onto the elements of the background mesh
      ALLOCATE(xl(nsd,lM%eNoN))
      DO iM=1, nMsh
!        Create a master list of elements of the background mesh based
!        IFEM bounding box position
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
         END DO
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         i = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               i = i + 1
               masElist(i) = e
            END IF
         END DO

!        Now perform search for each integration point of an element
!        whose trace was not determined earlier
         DO i=1, ne
            e = incEl(i)
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               Ac = lM%lN(Ac)
               xl(:,a) = xpL(:,Ac)
            END DO

            DO g=1, lM%nG
               xp = 0._RKIND
               DO a=1, lM%eNoN
                  xp = xp + xl(:,a)*lM%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec,
     2            xi)
               IF (Ec .NE. 0) THEN
                  ePtr(1,g,e) = Ec
                  ePtr(2,g,e) = iM
               END IF
            END DO
         END DO
         DEALLOCATE(eList, masElist)
      END DO
      DEALLOCATE(xl, incEl, xpL)

!     Transfer ePtr to trace data structure
      IF (ALLOCATED(lM%trc%gE))   DEALLOCATE(lM%trc%gE)
      IF (ALLOCATED(lM%trc%gptr)) DEALLOCATE(lM%trc%gptr)
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) i = i + 1
         END DO
      END DO
      lM%trc%nG = i
      ALLOCATE(lM%trc%gE(2,i), lM%trc%gptr(2,i))
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) THEN
               i = i + 1
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               lM%trc%gE(1,i)   = e
               lM%trc%gE(2,i)   = g
               lM%trc%gptr(:,i) = ePtr(:,g,e)
               msh(iM)%iGC(Ec)  = 1
            END IF
         END DO
      END DO
      DEALLOCATE(ePtr)

!     Abort simulation if all traces are still not found
      tt = CPUT()
      i  = cm%reduce(lM%trc%nG)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      j  = lM%nEl * lM%nG
      IF (i .LT. j) CALL DEBUGIBGPTRCS(lM)

      RETURN
      END SUBROUTINE IFEM_FINDGPTRACES
!--------------------------------------------------------------------
      SUBROUTINE IFEM_PARTMSH(lM, gM, lD, eList, part)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM, gM
      INTEGER(KIND=IKIND), INTENT(IN) :: eList(gM%nEl)
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)
      INTEGER(KIND=IKIND), INTENT(OUT) :: part(lM%nNo)

      INTEGER(KIND=IKIND) :: itr, a, b, e, Ac, Bc, ierr
      REAL(KIND=RKIND) :: tt, f, tol, dS, minS, xp(nsd), xb(nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: tmpI(:), gN(:)

      part = 0

      ALLOCATE(gN(tnNo))
      gN = 0
      DO e=1, gM%nEl
         IF (eList(e) .EQ. 0) CYCLE
         DO a=1, gM%eNoN
            Ac = gM%IEN(a,e)
            gN(Ac) = 1
         END DO
      END DO

      itr = 0
      f   = 0.1_RKIND
      ALLOCATE(tmpI(lM%nNo))
 001  part = 0
      tmpI = 0
      itr  = itr + 1
      f    = 2._RKIND*f
      tol  = (1._RKIND + f)*lM%dx
      DO a=1, lM%nNo
         IF (part(a) .NE. 0) CYCLE
         Ac   = lM%gN(a)
         xp   = ifem%x(:,Ac) + ifem%Ubo(:,Ac)
         minS = HUGE(minS)
         DO b=1, gM%nNo
            Bc = gM%gN(b)
            IF (gN(Bc) .EQ. 0) CYCLE
            xb = x(:,Bc) + lD(:,Bc)
            dS = SQRT( SUM( (xp(:)-xb(:))**2._RKIND ))
            IF (minS .GT. dS) minS = dS
         END DO
         IF (minS .LT. tol) THEN
            part(a) = cm%tF()
         END IF
      END DO

      tt = CPUT()
      CALL MPI_ALLREDUCE(part, tmpI, lM%nNo, mpint, MPI_MAX, cm%com(),
     2   ierr)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      b = 0
      DO a=1, lM%nNo
         IF (tmpI(a) .GT. 0) b = b + 1
      END DO

      IF (b .NE. lM%nNo) THEN
c         wrn = "Found only "//STR(b)//" nodes in pass "//STR(itr)//
c     2      " out of "//STR(lM%nNo)//" nodes"
         IF (itr .GT. 5) err = "Could not partition mesh in "//
     2      STR(itr)//"passes (IFEM_PARTMSH). Mesh could be distorted"//
     3      " or try changing edge size."
         GOTO 001
      END IF

      DEALLOCATE(gN, tmpI)

      RETURN
      END SUBROUTINE IFEM_PARTMSH
!-----------------------------------------------------------------------
      SUBROUTINE IFEM_FPSRCH(xp, lM, lD, ne, eSrch, itMax, Ec, xi)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: xp(nsd), lD(nsd,tnNo)
      INTEGER(KIND=IKIND), INTENT(IN)  :: ne, eSrch(ne), itMax
      INTEGER(KIND=IKIND), INTENT(OUT) :: Ec
      REAL(KIND=RKIND), INTENT(OUT) :: xi(nsd)
      TYPE(mshType), INTENT(IN) :: lM

      LOGICAL flag
      INTEGER(KIND=IKIND) :: a, e, i, El, En, nEl, nEn, iter, is, ie

      LOGICAL, ALLOCATABLE :: eChk(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: eList(:), tmpI(:)

      nEn = ne
      ALLOCATE(eList(nEn), eChk(lM%nEl))
      eList = eSrch
      eChk  = .FALSE.
      iter  = 0

      DO WHILE (iter .LE. itMax)
         ALLOCATE(tmpI(lM%nEl))
         tmpI = 0
         nEl  = nEn
         nEn  = 0
         DO e=1, nEl
            El = eList(e)
            eChk(El) = .TRUE.
         END DO
         DO e=1, nEl
            El = eList(e)
            is = lM%eAdj%prow(El)
            ie = lM%eAdj%prow(El+1) - 1
            DO i=is, ie
               En = lM%eAdj%pcol(i)
               IF (Ec.EQ.En .OR. eChk(En)) CYCLE
               flag = .TRUE.
               DO a=1, nEn
                  IF (En .EQ. tmpI(a)) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (flag) THEN
                  nEn = nEn + 1
                  tmpI(nEn) = En
               END IF
            END DO
         END DO
         DEALLOCATE(eList)
         ALLOCATE(eList(nEn))
         eList(:) = tmpI(1:nEn)
         CALL FINDE(xp, lM, x, lD, tnNo, nEn, eList, Ec, xi)
         DEALLOCATE(tmpI)
         IF (Ec .NE. 0) THEN
            DEALLOCATE(eChk, eList)
            RETURN
         ELSE
            iter = iter + 1
         END IF
      END DO

      DEALLOCATE(eChk, eList)

      RETURN
      END SUBROUTINE IFEM_FPSRCH
!####################################################################
!     Find parametric coordinate with respect to the parent element
!     of the background mesh for each IFEM nodal trace
      SUBROUTINE IFEM_FINDXINDTRC(lM, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      INTEGER a, b, i, Ac, Bc, Ec, iM, eNoN
      REAL(KIND=RKIND) :: xi(nsd), xp(nsd)

      REAL(KIND=RKIND), ALLOCATABLE :: N(:), Nx(:,:), xl(:,:)

      IF (ALLOCATED(lM%trc%xi)) DEALLOCATE(lM%trc%xi)
      ALLOCATE(lM%trc%xi(nsd,lM%trc%n))
      lM%trc%xi = 0._RKIND

      DO i=1, lM%trc%n
         b  = lM%trc%gN(i)
         Bc = lM%gN(b)
         Ec = lM%trc%nptr(1,i)
         iM = lM%trc%nptr(2,i)

!        Transfer to local arrays: background mesh variables
         eNoN = msh(iM)%eNoN
         ALLOCATE(N(eNoN), Nx(nsd,eNoN), xl(nsd,eNoN))
         DO a=1, eNoN
            Ac = msh(iM)%IEN(a,Ec)
            xl(:,a) = x(:,Ac)
            IF (mvMsh) xl(:,a) = xl(:,a) + lD(:,Ac)
         END DO

!        Initialize parametric coordinate
         xi = 0._RKIND
         DO a=1, msh(iM)%nG
            xi = xi + msh(iM)%xi(:,a)
         END DO
         xi = xi / REAL(msh(iM)%nG, KIND=RKIND)

!        Coordinates of the IFEM node
         xp = ifem%x(:,Bc) + ifem%Ubo(:,Bc)

!        Find shape functions and derivatives on the background mesh
!        at the IFEM node
         CALL GETNNX(msh(iM)%eType, eNoN, xl, msh(iM)%xib, msh(iM)%Nb,
     2      xp, xi, N, Nx)

         lM%trc%xi(:,i) = xi(:)

         DEALLOCATE(N, Nx, xl)
      END DO

      RETURN
      END SUBROUTINE IFEM_FINDXINDTRC
!--------------------------------------------------------------------
!     Find parametric coordinate with respect to the parent element
!     of the background mesh for each IFEM Gauss point trace
      SUBROUTINE IFEM_FINDXIGPTRC(lM, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      INTEGER a, b, e, g, i, Ac, Bc, Ec, iM, eNoN
      REAL(KIND=RKIND) :: xi(nsd), xp(nsd)

      REAL(KIND=RKIND), ALLOCATABLE :: Nb(:), N(:), Nx(:,:), xl(:,:)

      IF (ALLOCATED(lM%trc%xiG)) DEALLOCATE(lM%trc%xiG)
      ALLOCATE(lM%trc%xiG(nsd,lM%trc%nG))
      lM%trc%xiG = 0._RKIND

      ALLOCATE(Nb(lM%eNoN))
      DO i=1, lM%trc%nG
         e  = lM%trc%gE(1,i)
         g  = lM%trc%gE(2,i)
         Ec = lM%trc%gptr(1,i)
         iM = lM%trc%gptr(2,i)

         Nb = lM%N(:,g)
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Transfer to local arrays: background mesh variables
         eNoN = msh(iM)%eNoN
         ALLOCATE(N(eNoN), Nx(nsd,eNoN), xl(nsd,eNoN))
         DO a=1, eNoN
            Ac = msh(iM)%IEN(a,Ec)
            xl(:,a) = x(:,Ac)
            IF (mvMsh) xl(:,a) = xl(:,a) + lD(:,Ac)
         END DO

!        Initialize parametric coordinate
         xi = 0._RKIND
         DO a=1, msh(iM)%nG
            xi = xi + msh(iM)%xi(:,a)
         END DO
         xi = xi / REAL(msh(iM)%nG, KIND=RKIND)

!        Coordinates of the IFEM node
         xp = 0._RKIND
         DO b=1, lM%eNoN
            Bc = lM%IEN(b,e)
            xp(:) = xp(:) + Nb(b)*(ifem%x(:,Bc) + ifem%Ubo(:,Bc))
         END DO

!        Find shape functions and derivatives on the background mesh
!        at the IFEM node
         CALL GETNNX(msh(iM)%eType, eNoN, xl, msh(iM)%xib, msh(iM)%Nb,
     2      xp, xi, N, Nx)

         lM%trc%xiG(:,i) = xi(:)

         DEALLOCATE(N, Nx, xl)
      END DO
      DEALLOCATE(Nb)

      RETURN
      END SUBROUTINE IFEM_FINDXIGPTRC
!####################################################################
!     Communication structure for IFEM is initialized here. Here we
!     create a list of traces on master that are local to other
!     processes. The master then gathers all the data, projects flow
!     variables (acceleration, velocity and pressure) and broadcasts to
!     all the processes.
      SUBROUTINE IFEM_SETCOMMU()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

!     Free memory if already allocated
      CALL DESTROY(ifem%cm)

!     Create communication structure for nodal traces
      CALL IFEM_SETCOMMND()

!     Create communication structure for integration point traces
      CALL IFEM_SETCOMMGP()

      RETURN
      END SUBROUTINE IFEM_SETCOMMU
!--------------------------------------------------------------------
      SUBROUTINE IFEM_SETCOMMND()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) i, a, n, Ac, iM, ierr, tag
      REAL(KIND=RKIND) tt

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), ptr(:), rA(:,:),
     2   rReq(:)

!     Map trace pointers local to a process into a nodal vector.
      ALLOCATE(incNd(ifem%tnNo))
      incNd = 0
      DO iM=1, ifem%nMsh
         DO i=1, ifem%msh(iM)%trc%n
            a  = ifem%msh(iM)%trc%gN(i)
            Ac = ifem%msh(iM)%gN(a)
            incNd(Ac) = 1
         END DO
      END DO

!     All the included IFEM nodes are now mapped to a local vector
      n = SUM(incNd)
      ALLOCATE(ptr(n))
      i = 0
      DO a=1, ifem%tnNo
         IF (incNd(a) .NE. 0) THEN
            i = i + 1
            ptr(i) = a
         END IF
      END DO

!     Set IFEM comm data structures for sequential run
      IF (cm%seq()) THEN
         ALLOCATE(ifem%cm%n(1), ifem%cm%gN(n))
         ifem%cm%n(1)  = n
         ifem%cm%gN(:) = ptr
         DEALLOCATE(incNd, ptr)
         RETURN
      END IF

!     Gather no of included nodes from each process on master. Data at
!     these nodes will later be gathered by master from other processes
      ALLOCATE(ifem%cm%n(cm%np()))
      ifem%cm%n = 0
      tt = CPUT()
      CALL MPI_GATHER(n, 1, mpint, ifem%cm%n, 1, mpint, master, 
     2   cm%com(), ierr)

!     The processes that do not have any nodal traces return
      IF (.NOT.cm%mas() .AND. n.EQ.0) THEN
         DEALLOCATE(incNd, ptr)
         ALLOCATE(ifem%cm%gN(0))
         RETURN
      END IF

!     Master receives list of all nodal traces from other processes
      IF (cm%mas()) THEN
         n = SUM(ifem%cm%n)
         a = MAXVAL(ifem%cm%n)
         ALLOCATE(ifem%cm%gN(n), rA(a,cm%np()), rReq(cm%np()))
         ifem%cm%gN = 0
         rA(:,:)  = 0
         DO i=1, cm%np()
            n = ifem%cm%n(i)
            IF (n .EQ. 0) CYCLE
            IF (i .EQ. 1) THEN
               rA(1:n,i) = ptr(:)
            ELSE
               tag = i*100
               CALL MPI_IRECV(rA(1:n,i), n, mpint, i-1, tag, cm%com(),
     2            rReq(i), ierr)
            END IF
         END DO

         DO i=1, cm%np()
            IF (i.EQ.1 .OR. ifem%cm%n(i).EQ.0) CYCLE
            CALL MPI_WAIT(rReq(i), MPI_STATUS_IGNORE, ierr)
         END DO
         DEALLOCATE(rReq)

         a = 0
         DO i=1, cm%np()
            n = ifem%cm%n(i)
            ifem%cm%gN(a+1:a+n) = rA(1:n,i)
            a = a + n
         END DO
      ELSE
         ALLOCATE(ifem%cm%gN(0), rA(0,0))
         tag = cm%tF() * 100
         CALL MPI_SEND(ptr, n, mpint, master, tag, cm%com(), ierr)
      END IF
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      DEALLOCATE(incNd, ptr, rA)

      RETURN
      END SUBROUTINE IFEM_SETCOMMND
!--------------------------------------------------------------------
      SUBROUTINE IFEM_SETCOMMGP()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) i, a, e, nG, Ac, iM, ierr, tag
      REAL(KIND=RKIND) tt

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), ptr(:), rA(:,:),
     2   rReq(:)

!     Map trace pointers local to a process into a nodal vector. Note,
!     however, that the traces point to integration point of an element.
!     Therefore, we set all the nodes of an element that get contribution
!     from a valid trace.
      ALLOCATE(incNd(ifem%tnNo))
      incNd = 0
      DO iM=1, ifem%nMsh
         DO i=1, ifem%msh(iM)%trc%nG
            e = ifem%msh(iM)%trc%gE(1,i)
            DO a=1, ifem%msh(iM)%eNoN
               Ac = ifem%msh(iM)%IEN(a,e)
               incNd(Ac) = 1
            END DO
         END DO
      END DO

!     All the included IFEM nodes are now mapped to a local vector
      nG = SUM(incNd)
      ALLOCATE(ptr(nG))
      i = 0
      DO a=1, ifem%tnNo
         IF (incNd(a) .NE. 0) THEN
            i = i + 1
            ptr(i) = a
         END IF
      END DO

!     Set IFEM comm data structures for sequential run
      IF (cm%seq()) THEN
         ALLOCATE(ifem%cm%nG(1), ifem%cm%gE(nG))
         ifem%cm%nG(1) = nG
         ifem%cm%gE(:) = ptr
         DEALLOCATE(incNd, ptr)
         RETURN
      END IF

!     Gather no of included nodes from each process on master. Data at
!     these nodes will later be gathered by master from other processes
      ALLOCATE(ifem%cm%nG(cm%np()))
      ifem%cm%nG = 0
      tt = CPUT()
      CALL MPI_GATHER(nG, 1, mpint, ifem%cm%nG, 1, mpint, master,
     2   cm%com(), ierr)

!     The processes that do not have any integ. point traces, return
      IF (.NOT.cm%mas() .AND. nG.EQ.0) THEN
         DEALLOCATE(incNd, ptr)
         ALLOCATE(ifem%cm%gE(0))
         RETURN
      END IF

!     Master receives list of all nodal traces from other processes
      IF (cm%mas()) THEN
         nG = SUM(ifem%cm%nG)
         a  = MAXVAL(ifem%cm%nG)
         ALLOCATE(ifem%cm%gE(nG), rA(a,cm%np()), rReq(cm%np()))
         ifem%cm%gE = 0
         rA(:,:)  = 0
         DO i=1, cm%np()
            nG = ifem%cm%nG(i)
            IF (nG .EQ. 0) CYCLE
            IF (i .EQ. 1) THEN
               rA(1:nG,i) = ptr(:)
            ELSE
               tag = i*100
               CALL MPI_IRECV(rA(1:nG,i), nG, mpint, i-1, tag, cm%com(),
     2            rReq(i), ierr)
            END IF
         END DO

         DO i=1, cm%np()
            IF (i.EQ.1 .OR. ifem%cm%nG(i).EQ.0) CYCLE
            CALL MPI_WAIT(rReq(i), MPI_STATUS_IGNORE, ierr)
         END DO
         DEALLOCATE(rReq)

         a = 0
         DO i=1, cm%np()
            nG = ifem%cm%nG(i)
            ifem%cm%gE(a+1:a+nG) = rA(1:nG,i)
            a = a + nG
         END DO
      ELSE
         ALLOCATE(ifem%cm%gE(0), rA(0,0))
         tag = cm%tF() * 100
         CALL MPI_SEND(ptr, nG, mpint, master, tag, cm%com(), ierr)
      END IF
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      DEALLOCATE(incNd, ptr, rA)

      RETURN
      END SUBROUTINE IFEM_SETCOMMGP
!####################################################################
!     Set ighost field
!        ighost is set for both solids and shells
!        For solids, ghost node is a fluid node connected to atleast one
!        iblank=1 solid node. iblank is set to -1 on ghost nodes
!        For shells, ghost nodes are nodes that belong to ghost cells
      SUBROUTINE IFEM_SETIGHOST()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) :: a, b, e, Ac, iM
      REAL(KIND=RKIND) tt

      REAL(KIND=RKIND), ALLOCATABLE :: rG(:)

      ALLOCATE(rG(tnNo))
      rG = 0._RKIND

!     Update mesh and nodal ghost cell pointers based on iblank field
      DO iM=1, nMsh
         DO e=1, msh(iM)%nEl
            b = 0
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               b = b + iblank(Ac)
            END DO
            IF (b.GT.0 .AND. b.LT.msh(iM)%eNoN) THEN
               msh(iM)%iGC(e) = 1
               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%IEN(a,e)
                  b  = iblank(Ac)
                  IF (b .EQ. 0) rG(Ac) = 1._RKIND
               END DO
            END IF
         END DO
      END DO

c!     For Shells, use IFEM traces to identify ghost cells
c      DO iM=1, ifem%nMsh
c         IF (ifem%msh(iM)%lShl) THEN
c            DO b=1, ifem%msh(iM)%trc%n
c               e  = ifem%msh(iM)%trc%ptr(1,b)
c               jM = ifem%msh(iM)%trc%ptr(2,b)
c               DO a=1, msh(jM)%eNoN
c                  Ac = msh(jM)%IEN(a,e)
c                  rG(Ac) = 1._RKIND
c               END DO
c            END DO
c         END IF
c      END DO

      tt = CPUT()
      CALL COMMU(rG)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      DO a=1, tnNo
         IF (rG(a) .GT. 1.E-6_RKIND) iblank(a) = -1
      END DO

      RETURN
      END SUBROUTINE IFEM_SETIGHOST
!####################################################################
!     Write IFEM solution to a vtu file
      SUBROUTINE IFEM_WRITEVTUS(lY, lU)
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lY(nsd,ifem%tnNo)
      REAL(KIND=RKIND), INTENT(IN) :: lU(nsd,ifem%tnNo)

      TYPE(dataType) :: d(ifem%nMsh)
      TYPE(vtkXMLType) :: vtu

      INTEGER(KIND=IKIND) :: iStat, iEq, iOut, iM, a, e, Ac, Ec, nNo,
     2   nEl, s, l, ie, is, nSh, oGrp, outDof, nOut, cOut
      CHARACTER(LEN=stdL) :: fName

      INTEGER(KIND=IKIND), ALLOCATABLE :: outS(:), tmpI(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: outNames(:)

      IF (cm%slv()) THEN
         ifem%savedOnce = .TRUE.
         RETURN
      END IF

      nOut   = 1
      outDof = nsd
      DO iEq=1, nEq
         DO iOut=1, eq(iEq)%nOutIB
            IF (.NOT.eq(iEq)%outIB(iOut)%wtn(1)) CYCLE
            nOut   = nOut + 1
            outDof = outDof + eq(iEq)%outIB(iOut)%l
         END DO
      END DO

      ALLOCATE(outNames(nOut), outS(nOut+1))

!     Prepare all solultions in to dataType d
      nNo = 0
      nEl = 0
      DO iM=1, nMsh
         cOut           = 1
         outS(cOut)     = 1
         outS(cOut+1)   = nsd + 1
         outNames(cOut) = ""

         IF (ifem%msh(iM)%eType .EQ. eType_NRB) err =
     2      " Outputs for NURBS data is under development"

         d(iM)%nNo     = ifem%msh(iM)%nNo
         d(iM)%nEl     = ifem%msh(iM)%nEl
         d(iM)%eNoN    = ifem%msh(iM)%eNoN
         d(iM)%vtkType = ifem%msh(iM)%vtkType

         ALLOCATE(d(iM)%x(outDof,d(iM)%nNo),
     2      d(iM)%IEN(d(iM)%eNoN,d(iM)%nEl))
         DO a=1, ifem%msh(iM)%nNo
            Ac = ifem%msh(iM)%gN(a)
            d(iM)%x(1:nsd,a) = ifem%x(:,Ac)
         END DO

         DO e=1, ifem%msh(iM)%nEl
            d(iM)%IEN(:,e) = ifem%msh(iM)%IEN(:,e)
         END DO

         DO iEq=1, nEq
            DO iOut=1, eq(iEq)%nOutIB
               IF (.NOT.eq(iEq)%outIB(iOut)%wtn(1)) CYCLE
               l  = eq(iEq)%outIB(iOut)%l
               s  = eq(iEq)%s + eq(iEq)%outIB(iOut)%o
               e  = s + l - 1

               cOut = cOut + 1
               is   = outS(cOut)
               ie   = is + l - 1
               outS(cOut+1)   = ie + 1
               outNames(cOut) = "IFEM_"//TRIM(eq(iEq)%outIB(iOut)%name)

               oGrp = eq(iEq)%outIB(iOut)%grp
               SELECT CASE (oGrp)
                  CASE (outGrp_NA)
                  err = "Undefined output grp in VTK"
               CASE (outGrp_Y)
                  DO a=1, ifem%msh(iM)%nNo
                     Ac = ifem%msh(iM)%gN(a)
                     d(iM)%x(is:ie,a) = lY(s:e,Ac)
                  END DO
               CASE (outGrp_D)
                  DO a=1, ifem%msh(iM)%nNo
                     Ac = ifem%msh(iM)%gN(a)
                     d(iM)%x(is:ie,a) = lU(s:e,Ac)
                  END DO
               CASE DEFAULT
                  err = "Undefined output "//
     2               TRIM(eq(iEq)%outIB(iOut)%name)
               END SELECT
            END DO
         END DO

         ALLOCATE(d(iM)%xe(d(iM)%nEl,1))
         IF (.NOT.savedOnce) THEN
            IF (ALLOCATED(ifem%dmnID)) THEN
               d(iM)%xe(:,1) = REAL(ifem%msh(iM)%eId(:), KIND=RKIND)
            ELSE
               d(iM)%xe(:,1) = 1._RKIND
            END IF
         END  IF
         nNo = nNo +  d(iM)%nNo
         nEl = nEl +  d(iM)%nEl
      END DO

      ALLOCATE(tmpV(maxnsd,nNo))

!     Writing to vtu file (master only)
      IF (cTS .GE. 1000) THEN
         fName = STR(cTS)
      ELSE
         WRITE(fName,'(I3.3)') cTS
      END IF

      fName = TRIM(saveName)//"_IFEM_"//TRIM(ADJUSTL(fName))//".vtu"
      dbg = "Writing VTU"

      CALL vtkInitWriter(vtu, TRIM(fName), iStat)
      IF (iStat .LT. 0) err = "VTU file write error (init)"

!     Writing the position data
      iOut = 1
      s    = outS(iOut)
      e    = outS(iOut+1)-1
      nSh  = 0
      tmpV = 0._RKIND
      DO iM=1, ifem%nMsh
         DO a=1, d(iM)%nNo
            tmpV(1:nsd,a+nSh) = d(iM)%x(s:e,a)
         END DO
         nSh = nSh + d(iM)%nNo
      END DO
      CALL putVTK_pointCoords(vtu, tmpV(1:nsd,:), iStat)
      IF (iStat .LT. 0) err = "VTU file write error (coords)"

!     Writing the connectivity data
      nSh = -1
      DO iM=1, ifem%nMsh
         ALLOCATE(tmpI(d(iM)%eNoN,d(iM)%nEl))
         DO e=1, d(iM)%nEl
            tmpI(:,e) = d(iM)%IEN(:,e) + nSh
         END DO
         CALL putVTK_elemIEN(vtu, tmpI, d(iM)%vtkType, iStat)
         IF (iStat .LT. 0) err = "VTU file write error (ien)"
         DEALLOCATE(tmpI)
         nSh = nSh + d(iM)%nNo
      END DO

!     Writing all solutions
      DO iOut=2, nOut
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         s = outS(iOut)
         e = outS(iOut+1) - 1
         l = e - s + 1
         ALLOCATE(tmpV(l, nNo))
         nSh = 0
         DO iM=1, ifem%nMsh
            DO a=1, d(iM)%nNo
               tmpV(:,a+nSh) = d(iM)%x(s:e,a)
            END DO
            nSh = nSh + d(iM)%nNo
         END DO
         CALL putVTK_pointData(vtu, outNames(iOut), tmpV, iStat)
         IF (iStat .LT. 0) err = "VTU file write error (point data)"
      END DO

!     Write element-based variables
      IF (.NOT.savedOnce .OR. mvMsh) THEN
         ifem%savedOnce = .TRUE.
         ALLOCATE(tmpI(1,nEl))
!     Write the domain ID
         IF (ALLOCATED(ifem%dmnID)) THEN
            Ec = 0
            DO iM=1, ifem%nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(1,Ec) = INT(d(iM)%xe(e,1), KIND=IKIND)
               END DO
            END DO
            CALL putVTK_elemData(vtu, 'Domain_ID', tmpI, iStat)
            IF (iStat .LT. 0) err = "VTU file write error (dom id)"
         END IF

!     Write the mesh ID
         IF (ifem%nMsh .GT. 1) THEN
            Ec = 0
            DO iM=1, ifem%nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(1,Ec) = iM
               END DO
            END DO
            CALL putVTK_elemData(vtu, 'Mesh_ID', tmpI, iStat)
            IF (iStat .LT. 0) err = "VTU file write error (mesh id)"
         END IF
         DEALLOCATE(tmpI)
      END IF

      DO iM=1, nMsh
         CALL DESTROY(d(iM))
      END DO

      CALL vtkWriteToFile(vtu, iStat)
      IF (iStat .LT. 0) err = "VTU file write error"

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE IFEM_WRITEVTUS
!####################################################################
!     This is to check/create the txt file
      SUBROUTINE IFEM_CCTXT(lEq, fName, wtn)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eqType), INTENT(IN) :: lEq
      CHARACTER(LEN=stdL), INTENT(IN) :: fName(2)
      LOGICAL, INTENT(IN) :: wtn(2)

      INTEGER(KIND=IKIND), PARAMETER :: prL = 10

      LOGICAL flag
      INTEGER(KIND=IKIND) iM, iFa, fid, iDmn, i
      CHARACTER(LEN=stdL) stmp

      fid = 1
      IF (cm%slv()) RETURN

!     i=1 are the boundary values and i=2 are volume values
      DO i=1, 2
         IF (.NOT.wtn(i)) CYCLE

         INQUIRE(FILE=TRIM(fName(i)), EXIST=flag)
         IF (cTS.NE.0 .AND. flag) THEN
            CALL TRIMFILE(cTS+3,fName(i))
            CYCLE
         END IF

         OPEN(fid, FILE=TRIM(fName(i)))
         IF (i .EQ. 1) THEN
            DO iM=1, ifem%nMsh
               DO iFa=1, ifem%msh(iM)%nFa
                  stmp = ifem%msh(iM)%fa(iFa)%name
                  IF (LEN(TRIM(stmp)) .LE. prL) THEN
                     WRITE(fid,'(A)', ADVANCE='NO')
     2                  ADJUSTR(stmp(1:prL))//" "
                  ELSE
                     WRITE(fid,'(A)',ADVANCE='NO') TRIM(stmp)//" "
                  END IF
               END DO
            END DO
            WRITE(fid,*)
            DO iM=1, ifem%nMsh
               DO iFa=1, ifem%msh(iM)%nFa
                  stmp = STR(ifem%msh(iM)%fa(iFa)%area,prL)
                  WRITE(fid,'(A)', ADVANCE='NO') stmp(1:prL+1)
               END DO
            END DO
         ELSE
            DO iDmn=1, lEq%nDmnIB
               stmp = "DOMAIN-"//STR(lEq%dmnIB(iDmn)%Id)
               IF (lEq%dmnIB(iDmn)%Id .EQ. -1) stmp = "ENTIRE"
               WRITE(fid,'(A)', ADVANCE='NO') ADJUSTR(stmp(1:prL))//" "
            END DO
            WRITE(fid,*)
            DO iDmn=1, lEq%nDmnIB
               stmp = STR(lEq%dmnIB(iDmn)%v,prL)
               WRITE(fid,'(A)', ADVANCE='NO') stmp(1:prL+1)
            END DO
         END IF
         WRITE(fid,*)
         WRITE(fid,*)
         CLOSE(fid)
      END DO

      RETURN
      END SUBROUTINE IFEM_CCTXT
!####################################################################
!     This is to write to txt file
      SUBROUTINE IFEM_WTXT(lEq, m, fName, tmpV, wtn, div)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eqType), INTENT(IN) :: lEq
      INTEGER(KIND=IKIND), INTENT(IN) :: m
      CHARACTER(LEN=stdL), INTENT(IN) :: fName(2)
      REAL(KIND=RKIND), INTENT(IN) :: tmpV(maxnsd,ifem%tnNo)
      LOGICAL, INTENT(IN) :: wtn(2), div

      INTEGER(KIND=IKIND), PARAMETER :: prL = 10

      INTEGER(KIND=IKIND) iM, iFa, fid, i, iDmn
      REAL(KIND=RKIND) tmp

      fid = 1
      DO i=1, 2
         IF (.NOT.wtn(i)) CYCLE

         IF (cm%mas()) OPEN(fid, FILE=TRIM(fName(i)), STATUS='OLD',
     2      POSITION='APPEND')

         IF (i .EQ. 1) THEN
            DO iM=1, ifem%nMsh
               DO iFa=1, ifem%msh(iM)%nFa
                  IF (m .EQ. 1) THEN
                     IF (div) THEN
                        tmp = ifem%msh(iM)%fa(iFa)%area
                        tmp = Integ(ifem%msh(iM)%fa(iFa),tmpV,1)/tmp
                     ELSE
                        tmp = Integ(ifem%msh(iM)%fa(iFa),tmpV,1)
                     END IF
                  ELSE IF (m .EQ. nsd) THEN
                     tmp = Integ(ifem%msh(iM)%fa(iFa),tmpV,1,m)
                  ELSE
                     err = "WTXT only accepts 1 and nsd"
                  END IF
                  IF (cm%mas())
     2               WRITE(fid,'(A)',ADVANCE='NO') STR(tmp,prL)//" "
               END DO
            END DO
         ELSE
            DO iDmn=1, lEq%nDmnIB
               IF (div) THEN
                  tmp = lEq%dmnIB(iDmn)%v
                  tmp = Integ(lEq%dmnIB(iDmn)%Id, tmpV, 1, m)/tmp
               ELSE
                  tmp = Integ(lEq%dmnIB(iDmn)%Id, tmpV, 1, m)
               END IF
               IF (cm%mas())
     2            WRITE(fid,'(A)', ADVANCE='NO') STR(tmp,prL)//" "
            END DO
         END IF
         IF (cm%mas()) THEN
            WRITE(fid,*)
            CLOSE(fid)
         END IF
      END DO

      RETURN
      END SUBROUTINE IFEM_WTXT
!####################################################################
!     Apply Dirichlet boundary conditions on the immersed faces
      SUBROUTINE IFEM_SETBCDIR(Yb, Ub) !**! STILL TO CHECK
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: Yb(nsd+1,ifem%tnNo),
     2   Ub(nsd,ifem%tnNo)

      LOGICAL :: eDir(maxnsd)
      INTEGER(KIND=IKIND) :: iFa, iM, iEq, iBc, a, Ac, nNo, lDof, i

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: nV(:,:), tmpA(:,:), tmpY(:,:)

      iEq = ifem%cEq
      IF (iEq .EQ. 0) RETURN

      DO iBc=1, eq(iEq)%nBcIB
         IF (.NOT.BTEST(eq(iEq)%bcIB(iBc)%bType,bType_Dir)) CYCLE
         eDir = .FALSE.
         lDof = 0
         DO i=1, nsd
            IF (eq(iEq)%bcIB(iBc)%eDrn(i) .NE. 0) THEN
               eDir(i) = .TRUE.
               lDof = lDof + 1
            END IF
         END DO
         IF (lDof .EQ. 0) lDof = nsd

         iFa = eq(iEq)%bcIB(iBc)%iFa
         iM  = eq(iEq)%bcIB(iBc)%iM

!     Prepare a pointer list and normals for a face or a shell
         nNo = ifem%msh(iM)%fa(iFa)%nNo
         ALLOCATE(ptr(nNo), nV(nsd,nNo))
         ptr = 0
         DO a=1, nNo
            ptr(a)  = ifem%msh(iM)%fa(iFa)%gN(a)
            nV(:,a) = ifem%msh(iM)%fa(iFa)%nV(:,a)
         END DO

         ALLOCATE(tmpA(lDof,nNo), tmpY(lDof,nNo))
         CALL IFEM_SETBCDIRL(eq(iEq)%bcIB(iBc),nNo,lDof,nV,tmpA,tmpY)

         IF (ANY(eDir)) THEN
            DO a=1, nNo
               Ac = ptr(a)
               IF (BTEST(eq(iEq)%bcIB(iBc)%bType,bType_impD)) THEN
                  DO i=1, nsd
                     lDof = 0
                     IF (eDir(i)) THEN
                        lDof = lDof + 1
                        Yb(i,Ac) = tmpA(lDof,a)
                        Ub(i,Ac) = tmpY(lDof,a)
                     END IF
                  END DO
               ELSE
                  DO i=1, nsd
                     lDof = 0
                     IF (eDir(i)) THEN
                        lDof = lDof + 1
                        Yb(i,Ac) = tmpY(lDof,a)
                        Ub(i,Ac) = Ub(i,Ac) + Yb(i,Ac)*dt
                     END IF
                  END DO
               END IF
            END DO
         ELSE
            DO a=1, nNo
               Ac = ptr(a)
               IF (BTEST(eq(iEq)%bcIB(iBc)%bType,bType_impD)) THEN
                  Yb(1:nsd,Ac) = tmpA(:,a)
                  Ub(1:nsd,Ac) = tmpY(:,a)
               ELSE
                  Yb(1:nsd,Ac) = tmpY(:,a)
                  Ub(1:nsd,Ac) = Ub(1:nsd,Ac) + Yb(1:nsd,Ac)*dt
               END IF
            END DO
         END IF

         DEALLOCATE(ptr, nV, tmpA, tmpY)
      END DO

      RETURN
      END SUBROUTINE IFEM_SETBCDIR
!--------------------------------------------------------------------
      SUBROUTINE IFEM_SETBCDIRL(lBc, nNo, lDof, nvL, lA, lY)
      USE COMMOD
      IMPLICIT NONE
      TYPE(bcType), INTENT(IN) :: lBc
      INTEGER(KIND=IKIND), INTENT(IN) :: lDof, nNo
      REAL(KIND=RKIND), INTENT(IN) :: nvL(nsd,nNo)
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(lDof,nNo), lY(lDof,nNo)

      INTEGER(KIND=IKIND) :: a, i
      REAL(KIND=RKIND) :: dirY, dirA, nV(nsd)

      IF (BTEST(lBc%bType,bType_gen)) THEN
         IF (lDof .NE. lBc%gm%dof) err = "Inconsistent DOF to apply "//
     2      "Gen BC"
         IF (nNo .NE. SIZE(lBc%gm%d,2)) err = "Inconsistent nodes "//
     2      "to apply Gen BC"
         CALL IGBC(lBc%gm, lY, lA)
         RETURN
      ELSE IF (BTEST(lBc%bType,bType_ustd)) THEN
         CALL IFFT(lBc%gt, dirY, dirA)
      ELSE ! std / cpl
         dirA = 0._RKIND
         dirY = lBc%g
      END IF

      IF (lDof .EQ. nsd) THEN
         DO a=1, nNo
            nV      = nvL(:,a)
            lA(:,a) = dirA*lBc%gx(a)*nV
            lY(:,a) = dirY*lBc%gx(a)*nV
         END DO
      ELSE
         DO a=1, nNo
            DO i=1, lDof
               lA(i,a) = dirA*lBc%gx(a)
               lY(i,a) = dirY*lBc%gx(a)
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE IFEM_SETBCDIRL
!####################################################################
!     Project displacement from IFEM nodes to background mesh (implicit)
      SUBROUTINE IFEM_PRJCTU(Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) :: a, b, e, g, i, Ac, Bc, Ec, iM, jM, eNoN,
     2   eNoNb
      REAL(KIND=RKIND) :: w, Je, Jac, xi(nsd), up(nsd), F(nsd,nsd),tt(2)

      REAL(KIND=RKIND), ALLOCATABLE :: sA(:), Nb(:), Nbx(:,:), xbl(:,:),
     2   ubl(:,:), N(:), Nx(:,:), xl(:,:)

      tt(1) = CPUT()
      ifem%Un = 0._RKIND

      ALLOCATE(sA(tnNo))
      sA = 0._RKIND
      DO iM=1, ifem%nMsh
         eNoNb = ifem%msh(iM)%eNoN
         ALLOCATE(Nb(eNoNb), Nbx(nsd,eNoNb), xbl(nsd,eNoNb),
     2      ubl(nsd,eNoNb))
!        Loop over each trace, as we need to first interpolate flow var
!        at the IFEM integration points based on its trace
         DO i=1, ifem%msh(iM)%trc%nG
            e  = ifem%msh(iM)%trc%gE(1,i)
            g  = ifem%msh(iM)%trc%gE(2,i)
            Ec = ifem%msh(iM)%trc%gptr(1,i)
            jM = ifem%msh(iM)%trc%gptr(2,i)

!           Transfer to local arrays: IFEM mesh variables
            Nb = ifem%msh(iM)%N(:,g)
            IF (ifem%msh(iM)%eType .EQ. eType_NRB)
     2         CALL NRBNNX(ifem%msh(iM), e)
            DO b=1, eNoNb
               Bc = ifem%msh(iM)%IEN(b,e)
               xbl(:,b) = ifem%x(:,Bc)
               ubl(:,b) = ifem%Ubn(:,Bc)
            END DO
            CALL GNN(eNoNb, nsd, ifem%msh(iM)%Nx(:,:,g), xbl, 
     2                                                   Nbx, Je, F)
            IF (ISZERO(Je)) err = " Jac < 0 @ element "//e

            F = 0._RKIND
            IF (nsd .EQ. 3) THEN
               F(1,1) = 1._RKIND
               F(2,2) = 1._RKIND
               F(3,3) = 1._RKIND
               DO b=1, eNoNb
                  F(1,1) = F(1,1) + Nbx(1,b)*ubl(1,b)
                  F(1,2) = F(1,2) + Nbx(2,b)*ubl(1,b)
                  F(1,3) = F(1,3) + Nbx(3,b)*ubl(1,b)

                  F(2,1) = F(2,1) + Nbx(1,b)*ubl(2,b)
                  F(2,2) = F(2,2) + Nbx(2,b)*ubl(2,b)
                  F(2,3) = F(2,3) + Nbx(3,b)*ubl(2,b)

                  F(3,1) = F(3,1) + Nbx(1,b)*ubl(3,b)
                  F(3,2) = F(3,2) + Nbx(2,b)*ubl(3,b)
                  F(3,3) = F(3,3) + Nbx(3,b)*ubl(3,b)
               END DO
            ELSE
               F(1,1) = 1._RKIND
               F(2,2) = 1._RKIND
               DO b=1, eNoNb
                  F(1,1) = F(1,1) + Nbx(1,b)*ubl(1,b)
                  F(1,2) = F(1,2) + Nbx(2,b)*ubl(1,b)
                  F(2,1) = F(2,1) + Nbx(1,b)*ubl(2,b)
                  F(2,2) = F(2,2) + Nbx(2,b)*ubl(2,b)
               END DO
            END IF
            Jac = MAT_DET(F, nsd)
            w   = ifem%msh(iM)%w(g) * Jac * Je

!           Transfer to local arrays: background mesh variables
            eNoN = msh(jM)%eNoN
            ALLOCATE(N(eNoN), Nx(nsd,eNoN), xl(nsd,eNoN))
            DO a=1, eNoN
               Ac = msh(jM)%IEN(a,Ec)
               xl(:,a) = x(:,Ac)
               IF (mvMsh) xl(:,a) = xl(:,a) + Dg(nsd+2:2*nsd+1,Ac)
            END DO

!           Displacement of the integration point
            up = 0._RKIND
            DO b=1, eNoNb
               up = up + Nb(b)*ubl(:,b)
            END DO

!           Find shape functions and derivatives on the background mesh
!           at the integration point.
            xi = ifem%msh(iM)%trc%xiG(:,i)
            CALL GETGNN(nsd, msh(jM)%eType, eNoN, xi, N, Nx)

!           Project flow variables to IFEM nodes
            DO a=1, eNoN
               Ac = msh(jM)%IEN(a,Ec)
               sA(Ac) = sA(Ac) + w*N(a)
               ifem%Un(:,Ac) = ifem%Un(:,Ac) + w*N(a)*up(:)
            END DO

            DEALLOCATE(N, Nx, xl)
         END DO
         DEALLOCATE(Nb, Nbx, xbl, ubl)
      END DO

      tt(2) = CPUT()
      CALL COMMU(ifem%Un)
      CALL COMMU(sA)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt(2)

      DO a=1, tnNo
         IF (.NOT.ISZERO(sA(a))) THEN
            ifem%Un(:,a) = ifem%Un(:,a) / sA(a)
         END IF
      END DO
      DEALLOCATE(sA)

      ifem%callD(1) = ifem%callD(1) + CPUT() - tt(1)

      RETURN
      END SUBROUTINE IFEM_PRJCTU
!####################################################################
!     Predictor step (implicit)
      SUBROUTINE IFEM_PICP()
      USE COMMOD
      IMPLICIT NONE

      INTEGER a
      REAL(KIND=RKIND) coef, tt

      tt = CPUT()
      coef = (eq(ifem%cEq)%gam - 1._RKIND)/eq(ifem%cEq)%gam

      DO a=1, ifem%tnNo
         ifem%Aun(:,a) = ifem%Auo(:,a) * coef
         ifem%Ubn(:,a) = ifem%Ubo(:,a)
      END DO

      DO a=1, tnNo
         ifem%Uo(:,a) = ifem%Un(:,a)
      END DO

      ifem%callD(1) = ifem%callD(1) + CPUT() - tt

      RETURN
      END SUBROUTINE IFEM_PICP
!####################################################################
!     Initiator step (implicit) - compute Au_(n+am) and Ub_(n+af)
      SUBROUTINE IFEM_PICI()
      USE COMMOD
      IMPLICIT NONE

      INTEGER a
      REAL(KIND=8) :: coef(4), tt

      IF (ifem%cpld .NE. ibCpld_I) RETURN
      tt = CPUT()

      coef(1) = 1._RKIND - eq(ifem%cEq)%am
      coef(2) = eq(ifem%cEq)%am
      coef(3) = 1._RKIND - eq(ifem%cEq)%af
      coef(4) = eq(ifem%cEq)%af

      DO a=1, ifem%tnNo
         ifem%Auk(:,a) = ifem%Auo(:,a)*coef(1) + ifem%Aun(:,a)*coef(2)
         ifem%Ubk(:,a) = ifem%Ubo(:,a)*coef(3) + ifem%Ubn(:,a)*coef(4)
      END DO

      DO a=1, tnNo
         ifem%Un(:,a) = ifem%Uo(:,a)*coef(3) + ifem%Un(:,a)*coef(4)
      END DO

      ifem%callD(1) = ifem%callD(1) + CPUT() - tt

      RETURN
      END SUBROUTINE IFEM_PICI
!####################################################################
!     Compute FSI force on the immersed bodies
      SUBROUTINE IFEM_CALCFFSI(Ag, Yg, Dg, Aug, Ubg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo), Aug(nsd,ifem%tnNo), Ubg(nsd,ifem%tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, iFn, nFn, iM, j
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd), tt(2)

      REAL(KIND=RKIND), ALLOCATABLE :: xbl(:,:), aul(:,:), ubl(:,:), 
     2   fN(:,:), N(:), Nx(:,:), lR(:,:)
 

      tt(1) = CPUT()
C     cDmn = DOMAIN(msh(1), ifem%cEq, 1)

!     Initialize residues to 0
      ifem%Rsolid = 0._RKIND
      write(*,*)" Beginning IFEM_CALCFFSI "


!     We compute the fluid-structure interaction force on the solid initial 
!     domain
!     Loop over all IFEM mesh
C       write(*,*)"Beginning loop over nmb of meshes "
      DO iM=1, ifem%nMsh
   !     Loop over all elements of mesh
         eNoN = ifem%msh(iM)%eNoN
         nFn  = MAX(ifem%msh(iM)%nFn, 1)

         ALLOCATE(xbl(nsd,eNoN), aul(nsd,eNoN), ubl(nsd,eNoN), 
     2            fN(nsd,nFn), N(eNoN), Nx(nsd,eNoN), lR(nsd,eNoN))
         
C          write(*,*)"Beginning loop over elemtn mesh "
         DO e=1, ifem%msh(iM)%nEl

            ifem%cDmn = IB_DOMAIN(ifem%msh(iM), ifem%cEq, e)

   !        Update shape functions for NURBS
            IF (ifem%msh(iM)%eType .EQ. eType_NRB) THEN
               CALL NRBNNX(ifem%msh(iM), e)
            END IF

   !        Create local copies
            fN   = 0._RKIND
            DO a=1, eNoN
               Ac = ifem%msh(iM)%IEN(a,e)
               xbl(:,a) = ifem%x(:,Ac)
               ubl(:,a) = Ubg(:,Ac)
               aul(:,a) = Aug(:,Ac)
               IF (ALLOCATED(ifem%msh(iM)%fN)) THEN
                  DO j=1, nFn
                     fN(:,j) = ifem%msh(iM)%fN((j-1)*nsd+1:j*nsd,e)
                  END DO
               END IF
            END DO

C             write(*,*) " Beginning loop over gauss point "
   !        Gauss integration
C             lR = 0._RKIND
            DO g=1, ifem%msh(iM)%nG
               lR = 0._RKIND

               IF (g.EQ.1 .OR. .NOT.ifem%msh(iM)%lShpF) THEN
                  CALL GNN(eNoN, nsd, ifem%msh(iM)%Nx(:,:,g), xbl, Nx, 
     2                              Jac, ksix)
                  IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
               END IF
               w = ifem%msh(iM)%w(g) * Jac
               N = ifem%msh(iM)%N(:,g)

C                write(*,*)"Calling local force assembly"
               IF (nsd .EQ. 3) THEN
                  CALL IFEM_CALCLFFSI3D(eNoN, nFn, w, N, Nx, aul, ubl, 
     2                fN, lR)

               ELSE IF (nsd .EQ. 2) THEN
                  CALL IFEM_CALCLFFSI2D(eNoN, nFn, w, N, Nx, aul, ubl, 
     2                fN, lR)
               END IF

C                write(*,*)"End call local assembly"

!              Assemble to ifem global residue
               DO a=1, eNoN
                  Ac = ifem%msh(iM)%IEN(a,e)
                  DO j=1, nsd ! TODO nsd+1 ??
                     ifem%Rsolid(j,Ac) = ifem%Rsolid(j,Ac) + lR(j,a)
                  END DO
               END DO
            END DO
         
         END DO

         DEALLOCATE(xbl, aul, ubl, fN, N, Nx, lR)

      END DO

      tt(2) = CPUT()
      CALL COMMU(ifem%Rsolid)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt(2)

      ifem%callD(1) = ifem%callD(1) + CPUT() - tt(1)

      RETURN
      END SUBROUTINE IFEM_CALCFFSI
!--------------------------------------------------------------------
!     Compute the 3D FSI force due to IFEM on the background fluid
      SUBROUTINE IFEM_CALCLFFSI3D(eNoN, nFn, w, N, Nx, aul, ubl, fN, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN), 
     2   aul(3,eNoN), ubl(3,eNoN), fN(3,nFn)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(4,eNoN)

      INTEGER(KIND=IKIND) :: a, iEq, iDmn, i, j, k
      REAL(KIND=RKIND) :: Jac, ya_g
      REAL(KIND=RKIND) :: Dm(6,6), F(3,3), Fi(3,3), S(3,3), P(3,3)

      iEq     = ifem%cEq
      iDmn    = ifem%cDmn
      i       = 1
      j       = 2
      k       = 3

      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      F(3,3) = 1._RKIND
      ya_g   = 0._RKIND

      DO a=1, eNoN
         F(1,1)  = F(1,1)  + Nx(1,a)*ubl(i,a)
         F(1,2)  = F(1,2)  + Nx(2,a)*ubl(i,a)
         F(1,3)  = F(1,3)  + Nx(3,a)*ubl(i,a)
         F(2,1)  = F(2,1)  + Nx(1,a)*ubl(j,a)
         F(2,2)  = F(2,2)  + Nx(2,a)*ubl(j,a)
         F(2,3)  = F(2,3)  + Nx(3,a)*ubl(j,a)
         F(3,1)  = F(3,1)  + Nx(1,a)*ubl(k,a)
         F(3,2)  = F(3,2)  + Nx(2,a)*ubl(k,a)
         F(3,3)  = F(3,3)  + Nx(3,a)*ubl(k,a)
      END DO

      Jac = MAT_DET(F, 3)
      Fi  = MAT_INV(F, 3)

!     2nd Piola-Kirchhoff tensor (S) and material stiffness tensor in
!     Voigt notation
      CALL GETPK2CC(eq(iEq)%dmnIB(iDmn), F, nFn, fN, ya_g, S, Dm)

!     1st Piola-Kirchhoff tensor (P)
      P = MATMUL(F, S)

!     Local residue
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(Nx(1,a)*P(1,1) + 
     2      Nx(2,a)*P(1,2) + Nx(3,a)*P(1,3))
         lR(2,a) = lR(2,a) + w*(Nx(1,a)*P(2,1) +
     2      Nx(2,a)*P(2,2) + Nx(3,a)*P(2,3))
         lR(3,a) = lR(3,a) + w*(Nx(1,a)*P(3,1) +
     2      Nx(2,a)*P(3,2) + Nx(3,a)*P(3,3))
      END DO

      RETURN
      END SUBROUTINE IFEM_CALCLFFSI3D
!--------------------------------------------------------------------
!     Compute the 2D FSI force due to IFEM on the background fluid
      SUBROUTINE IFEM_CALCLFFSI2D(eNoN, nFn, w, N, Nx, aul, ubl, fN, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN), 
     2   aul(2,eNoN), ubl(2,eNoN), fN(2,nFn)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(2,eNoN)

      INTEGER(KIND=IKIND) :: a, iEq, iDmn, i, j
      REAL(KIND=RKIND) :: Jac, ya_g
      REAL(KIND=RKIND) :: Dm(3,3), F(2,2), Fi(2,2), S(2,2), P(2,2)

      iEq     = ifem%cEq
      iDmn    = ifem%cDmn

      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      ya_g   = 0._RKIND
      i = 1
      j = 2

      DO a=1, eNoN
         F(1,1)  = F(1,1)  + Nx(1,a)*ubl(i,a)
         F(1,2)  = F(1,2)  + Nx(2,a)*ubl(i,a)
         F(2,1)  = F(2,1)  + Nx(1,a)*ubl(j,a)
         F(2,2)  = F(2,2)  + Nx(2,a)*ubl(j,a)
      END DO

C       write(*,*)"ubl = ", ubl
C       write(*,*)"F = ", F

      Jac = MAT_DET(F, 2)
      Fi  = MAT_INV(F, 2)

C       write(*,*)"Jac = ", Jac

!     2nd Piola-Kirchhoff tensor (S) and material stiffness tensor in
!     Voigt notation
      CALL GETPK2CC(eq(iEq)%dmnIB(iDmn), F, nFn, fN, ya_g, S, Dm)

C       write(*,*)"eq(iEq)%dmnIB(iDmn)%phys", eq(iEq)%dmnIB(iDmn)%phys
C       write(*,*)" S = ", S

C       write(*,*)" GETPK2CC done"
!     1st Piola-Kirchhoff tensor (P)
      P = MATMUL(F, S)

C       write(*,*)"adding to local residue"
!     Local residue
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(Nx(1,a)*P(1,1) + Nx(2,a)*P(1,2))
         lR(2,a) = lR(2,a) + w*(Nx(1,a)*P(2,1) + Nx(2,a)*P(2,2))
      END DO

C       write(*,*)"lR = ", lR

      RETURN
      END SUBROUTINE IFEM_CALCLFFSI2D
!--------------------------------------------------------------------
!####################################################################
      SUBROUTINE IFEM_DOASSEM(eNoN, ptr, lKu, lK)
      USE TYPEMOD
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, ptr(eNoN)
      REAL(KIND=RKIND), INTENT(IN) :: lKu((nsd+1)*nsd,eNoN,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b, rowN, colN, idx

      DO a=1, eNoN
         rowN = ptr(a)
         DO b=1, eNoN
            colN = ptr(b)
            CALL GETCOLPTR()
            ifem%Ku(:,idx) = ifem%Ku(:,idx) + lKu(:,a,b)
            Val(:,idx)   = Val(:,idx)   + lK(:,a,b)
         END DO
      END DO

      RETURN
      CONTAINS
!--------------------------------------------------------------------
         SUBROUTINE GETCOLPTR()
         IMPLICIT NONE

         INTEGER(KIND=IKIND) left, right

         left  = rowPtr(rowN)
         right = rowPtr(rowN+1)
         idx   = (right + left)/2
         DO WHILE (colN .NE. colPtr(idx))
            IF (colN .GT. colPtr(idx)) THEN
               left  = idx
            ELSE
               right = idx
            END IF
            idx = (right + left)/2
         END DO

         RETURN
         END SUBROUTINE GETCOLPTR
!--------------------------------------------------------------------
      END SUBROUTINE IFEM_DOASSEM
!####################################################################
      SUBROUTINE IFEM_RHSUpdate()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER a, i, c
      REAL(KIND=RKIND) ami, tt(2)

      REAL(KIND=RKIND), ALLOCATABLE :: KR(:,:)

      tt(1) = CPUT()

      DO a=1, ifem%tnNo
         DO i=1, nsd
            ifem%Rub(i,a) = ifem%Auk(i,a) - ifem%Yb(i,a)
         END DO
      END DO

      ami = 1._RKIND/eq(ifem%cEq)%am

      IF (nsd .EQ. 3) THEN
         ALLOCATE(KR(4,tnNo))
         KR = 0._RKIND
         DO a=1, tnNo
            DO i=rowPtr(a), rowPtr(a+1)-1
               c = colPtr(i)
               KR(1,a) = KR(1,a) + ifem%Ku(1,i)*ifem%Ru(1,c)
     2          + ifem%Ku(2,i)*ifem%Ru(2,c) + ifem%Ku(3,i)*ifem%Ru(3,c)
               KR(2,a) = KR(2,a) + ifem%Ku(4,i)*ifem%Ru(1,c)
     2          + ifem%Ku(5,i)*ifem%Ru(2,c) + ifem%Ku(6,i)*ifem%Ru(3,c)
               KR(3,a) = KR(3,a) + ifem%Ku(7,i)*ifem%Ru(1,c)
     2          + ifem%Ku(8,i)*ifem%Ru(2,c) + ifem%Ku(9,i)*ifem%Ru(3,c)
               KR(4,a) = KR(4,a) + ifem%Ku(10,i)*ifem%Ru(1,c)
     2         + ifem%Ku(11,i)*ifem%Ru(2,c) + ifem%Ku(12,i)*ifem%Ru(3,c)
            END DO
         END DO

         tt(2) = CPUT()
         CALL COMMU(KR)
         ifem%callD(3) = ifem%callD(3) + CPUT() - tt(2)

         DO a=1, tnNo
            ifem%R(1,a) = ifem%R(1,a) - ami*KR(1,a)
            ifem%R(2,a) = ifem%R(2,a) - ami*KR(2,a)
            ifem%R(3,a) = ifem%R(3,a) - ami*KR(3,a)
            ifem%R(4,a) = ifem%R(4,a) - ami*KR(4,a)
         END DO

      ELSE IF (nsd .EQ. 2) THEN
         ALLOCATE(KR(3,tnNo))
         KR = 0._RKIND
         DO a=1, tnNo
            DO i=rowPtr(a), rowPtr(a+1)-1
               c = colPtr(i)
               KR(1,a) = KR(1,a) + ifem%Ku(1,i)*ifem%Ru(1,c)
     2            + ifem%Ku(2,i)*ifem%Ru(2,c)
               KR(2,a) = KR(1,a) + ifem%Ku(3,i)*ifem%Ru(1,c)
     2            + ifem%Ku(4,i)*ifem%Ru(2,c)
               KR(3,a) = KR(1,a) + ifem%Ku(5,i)*ifem%Ru(1,c)
     2            + ifem%Ku(6,i)*ifem%Ru(2,c)
            END DO
         END DO

         CALL COMMU(KR)

         DO a=1, tnNo
            ifem%R(1,a) = ifem%R(1,a) - ami*KR(1,a)
            ifem%R(2,a) = ifem%R(2,a) - ami*KR(2,a)
            ifem%R(3,a) = ifem%R(3,a) - ami*KR(3,a)
         END DO

      END IF

      ifem%callD(1) = ifem%callD(1) + CPUT() - tt(1)

      RETURN
      END SUBROUTINE IFEM_RHSUpdate
!####################################################################
!     Add contribution from IFEM to the residue (RHS) via MLS
!     We suppose for the moment that the fluid mesh is static
      SUBROUTINE IFEM_CONSTRUCT()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) a, is, iM, mnS, idFCls, nbrSN, idFStc, nd

      REAL(KIND=RKIND), ALLOCATABLE :: Amls(:,:), Bmls(:,:), Wmls(:)
      REAL(KIND=RKIND), ALLOCATABLE :: Pmls(:,:), PPt(:,:), q(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: QMLS(:,:), Ai(:,:), qt(:,:)
      REAL(KIND=RKIND) :: xlf(nsd), xls(nsd), rnorm, wfunc, sum
      REAL(KIND=RKIND) :: Aaux(1,nsd+1)

      ALLOCATE(Amls(nsd+1,nsd+1), Pmls(nsd+1,1), PPt(nsd+1,nsd+1), 
     2         Ai(nsd+1,nsd+1))

      mnS = SIZE(msh(1)%stn%ndStn,2)
      ifem%maxNbrST = mnS
      ALLOCATE(QMLS(mnS,ifem%tnNo))
      QMLS = 0._RKIND

      ifem%Rfluid = 0._RKIND
!     Loop over the nbr of solid meshes 
      DO iM=1, nMsh
C          write(*,*)"inside loop mesh"
!        Compute mesh space discr param if not computed yet         
         CALL GETMESHDIAM(msh(iM))
C          write(*,*)"diam computed"

!        Loop over the solid nodes 
         DO a=1, ifem%tnNo
C             write(*,*)"inside loop solid node"
!           Global id closest fluid node
            idFCls = ifem%clsFNd(a)  
!           Nbr of node in stencil
            nbrSN = msh(iM)%stn%nbrNdStn(idFCls) 
!           Update current solid node coord
            xls = ifem%x(:,a) + ifem%Ubo(:,a)

C             write(*,*)"solid coord:", xls
C             write(*,*)"if closest point :", idFCls
C             write(*,*)"nmb stencil node", nbrSN

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
!              Extract fluid coordinates xlf
               idFStc = msh(iM)%stn%ndStn(idFCls,is)
               xlf = x(:,idFStc)

C                write(*,*)"idFStc = ", idFStc
C                write(*,*)"xlf = ", xlf

               Pmls(1,1) = 1._RKIND
               Pmls(2,1) = xlf(1)
               Pmls(3,1) = xlf(2)

C                write(*,*)" Pmls ", Pmls
               
!              Compute cubic spline value
               rnorm = msh(iM)%diam
               Wmls(is) = WHTFUNC(xls,xlf,rnorm)
!              Check that rnorm is good 
               rnorm = DIST(xls,xlf)
               rnorm = rnorm / (msh(iM)%diam*1.15)
C                write(*,*)"rnorm is ", rnorm
C                write(*,*)"Wmls is ", Wmls(is)


               PPt = Wmls(is)*MATMUL(Pmls, TRANSPOSE(Pmls))
C                write(*,*)"PPt = ",PPt

               Amls = Amls + PPt
C                write(*,*)"Matrix A inside is ",Amls

               Bmls(:,is) = Wmls(is)*Pmls(:,1)

               sum = sum + Wmls(is)
            END DO

C             write(*,*)"Sum W = ", sum
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

C             write(*,*) ""//""
C             write(*,*) ""//""
C             write(*,*) ""//""

         END DO
      END DO

C       write(*,*)"QMLS = ", QMLS
C       write(*,*)"ifem%Rsolid = ", ifem%Rsolid
C       write(*,*)"ifem%Rfluid = ", ifem%Rfluid

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

            ifem%Rfluid(:,idFStc) = ifem%Rfluid(:,idFStc) + 
     2                                     QMLS(is,a)*ifem%Rsolid(:,a)
         END DO
      END DO

      write(*,*)"Rfluid built done"

      IF(.NOT.ALLOCATED(ifem%QMLS)) ALLOCATE(ifem%QMLS(mnS,ifem%tnNo))
      ifem%QMLS = QMLS



      DEALLOCATE(Amls, Pmls, PPt, Ai, QMLS)

      RETURN
      END SUBROUTINE IFEM_CONSTRUCT
!####################################################################
      SUBROUTINE IFEM_RASSEMBLY()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) a

!     Final residual assembly 
      DO a=1, tnNo
         R(1:nsd,a) = R(1:nsd,a) + ifem%Rfluid(:,a)
      END DO


      RETURN
      END SUBROUTINE IFEM_RASSEMBLY
!####################################################################
!     Interpolate flow velocity and pressure at IFEM nodes from background
!     mesh. Use IFEM velocity to compute IFEM displacement (for explicit
!     coupling)
      SUBROUTINE IFEM_INTERPVEL(Yg, Dg, iT)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)
      INTEGER(KIND=IKIND) :: iT

      REAL(KIND=RKIND) :: tt
      INTEGER(KIND=IKIND) :: a, is, ie, i
      REAL(KIND=RKIND), ALLOCATABLE :: Yl(:,:)

!     For the moment we are using only Yn, the displacement in case 
!     of ALE will be consider in the future

      write(*,*)"Inside IFEM_INTERPVEL"

      tt = CPUT()
      is = eq(ifem%cEq)%s
      ie = eq(ifem%cEq)%e

C       write(*,*)"is = ", is, ", ie = ", ie

      ALLOCATE(Yl(nsd+1,tnNo))
      DO a=1, tnNo
         Yl(:,a) = Yg(is:ie,a) ! keeping only the fluid vel and pressure 
      END DO

C       write(*,*)"calling IFEM_FINDSOLVEL"
!      CALL IFEM_FINDSOLVEL(nsd+1, Yl, Dg, ifem%Yb)

      CALL IFEM_FINDSOLVEL_MLS(nsd+1, Yl, Dg, ifem%Yb)

      DEALLOCATE(Yl)
C       write(*,*)"done calling IFEM_FINDSOLVEL"


!     Computing new solid location using A-B scheme given the new interpolated 
!     fluid velocity   

C       write(*,*)"ifem%xCuo ", ifem%xCuo
C       write(*,*)""//""
C       write(*,*)"ifem%xCu", ifem%xCu
C       write(*,*)""//""
      
!     Solid update location via displacement   
      IF(iT .GE. 3) THEN 
         DO a=1, ifem%tnNo
            DO i=1, nsd
               ifem%Ubo(i,a) = ifem%Ubo(i,a)
     2                 + dt*0.5_RKIND*( 3._RKIND * ifem%Auo(i,a) 
     3                 - ifem%Auoo(i,a) )

               ifem%xCu(i,a) = ifem%x(i,a) + ifem%Ubo(i,a)
               IF (ABS(ifem%xCu(i,a)).GT.1._RKIND) THEN 
                  write(*,*)"OHHH NOOOO, out of fluid !"
                  CALL EXIT(1)
               END IF
            END DO
         END DO
      ELSE 
         DO a=1, ifem%tnNo
            DO i=1, nsd
               ifem%Ubo(i,a) = ifem%Ubo(i,a) + dt*ifem%Yb(i,a)

               ifem%xCu(i,a) = ifem%x(i,a) + ifem%Ubo(i,a)
               IF (ABS(ifem%xCu(i,a)).GT.1._RKIND) THEN 
                  write(*,*)"OHHH NOOOO, out of fluid !"
                  CALL EXIT(1)
               END IF
            END DO
         END DO
      END IF
      write(*,*)"Done AB scheme" 
      ifem%Auoo = ifem%Auo
      ifem%Auo = ifem%Yb









      ifem%callD(1) = ifem%callD(1) + CPUT() - tt

      RETURN
      END SUBROUTINE IFEM_INTERPVEL
!####################################################################
!     Interpolate data at IFEM nodes from background mesh using fem
!     nodal function 
      SUBROUTINE IFEM_FINDSOLVEL(m, Ug, Dg, Ub)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m
      REAL(KIND=RKIND), INTENT(IN) :: Ug(m,tnNo), Dg(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(OUT) :: Ub(nsd,ifem%tnNo) !new solid vel

      REAL(KIND=RKIND) :: xi(nsd), up(nsd)
      REAL(KIND=RKIND), ALLOCATABLE :: N(:), Nx(:,:), xl(:,:), ul(:,:)
      INTEGER(KIND=IKIND) :: a, jM, Ac, e, eNoN, idFEl

      jM = 1

      Ub = 0._RKIND

      DO e=1, ifem%tnNo ! maybe just nNo
!        Find fluid elem in which is it in the def config 
         idFEl = ifem%clsFElm(e)
         
         eNoN = msh(jM)%eNoN
         ALLOCATE(N(eNoN), Nx(nsd,eNoN), xl(nsd,eNoN), ul(nsd,eNoN))
         
         ul = 0._RKIND
         DO a=1, eNoN
            Ac = msh(jM)%IEN(a,idFEl)
            ul(:,a) = Ug(1:nsd,Ac)
            xl(:,a) = x(:,Ac)
C             IF (mvMsh) xl(:,a) = xl(:,a) + Dg(nsd+2:2*nsd+1,Ac)
         END DO

!        Find shape functions and derivatives on the background mesh
!        at the IFEM node
         xi = ifem%x(:,e)
         CALL GETGNN(nsd, msh(jM)%eType, eNoN, xi, N, Nx)

!           Use computed shape functions to interpolate flow variables
         up = 0._RKIND
         DO a=1, eNoN
            up = up + N(a)*ul(:,a)
         END DO

         Ub(:,e) = Ub(:,e) + up(:)

         DEALLOCATE(N, Nx, xl, ul)
      END DO

      RETURN
      END SUBROUTINE IFEM_FINDSOLVEL
!####################################################################
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
!--------------------------------------------------------------------
      SUBROUTINE IFEM_INTERPND(m, Ug, Dg, Ub)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m
      REAL(KIND=RKIND), INTENT(IN) :: Ug(m,tnNo), Dg(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(OUT) :: Ub(m,ifem%tnNo)

      INTEGER(KIND=IKIND) :: a, b, i, Ac, Bc, Ec, iM, jM, eNoN
      REAL(KIND=RKIND) :: xi(nsd), up(m), tt

      REAL(KIND=RKIND), ALLOCATABLE :: wgt(:), N(:), Nx(:,:), xl(:,:),
     2   ul(:,:)

!     Copy ifem%Ubn
      ALLOCATE(wgt(ifem%tnNo))
      wgt = 0._RKIND
      Ub  = 0._RKIND
      DO iM=1, ifem%nMsh
!        Loop over each trace, as we need to first interpolate flow var
!        at the IFEM nodes based on its trace
         DO i=1, ifem%msh(iM)%trc%n
            b  = ifem%msh(iM)%trc%gN(i)
            Bc = ifem%msh(iM)%gN(b)
            Ec = ifem%msh(iM)%trc%nptr(1,i)
            jM = ifem%msh(iM)%trc%nptr(2,i)

!           Transfer to local arrays: background mesh variables
            eNoN = msh(jM)%eNoN
            ALLOCATE(N(eNoN), Nx(nsd,eNoN), xl(nsd,eNoN), ul(m,eNoN))
            ul = 0._RKIND
            DO a=1, eNoN
               Ac = msh(jM)%IEN(a,Ec)
               ul(:,a) = Ug(:,Ac)
               xl(:,a) = x(:,Ac)
               IF (mvMsh) xl(:,a) = xl(:,a) + Dg(nsd+2:2*nsd+1,Ac)
            END DO

!           Find shape functions and derivatives on the background mesh
!           at the IFEM node
            xi = ifem%msh(iM)%trc%xi(:,i)
            CALL GETGNN(nsd, msh(jM)%eType, eNoN, xi, N, Nx)

!           Use computed shape functions to interpolate flow variables
            up = 0._RKIND
            DO a=1, eNoN
               up = up + N(a)*ul(:,a)
            END DO

            Ub(:,Bc) = Ub(:,Bc) + up(:)
            wgt(Bc)  = wgt(Bc)  + 1._RKIND

            DEALLOCATE(N, Nx, xl, ul)
         END DO
      END DO

!     Synchronize Yb across all the processes
      tt = CPUT()
      CALL IB_SYNCN(Ub)
      CALL IB_SYNCN(wgt)
!     TODO
!     CALL IFEM_SYNCN(Ub)
!     CALL IFEM_SYNCN(wgt)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      DO a=1, ifem%tnNo
         IF (.NOT.ISZERO(wgt(a))) THEN
            Ub(:,a) = Ub(:,a) / wgt(a)
         END IF
      END DO

      DEALLOCATE(wgt)

      RETURN
      END SUBROUTINE IFEM_INTERPND
!--------------------------------------------------------------------
      SUBROUTINE IFEM_INTERPGP(m, Ug, Dg, Ub)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m
      REAL(KIND=RKIND), INTENT(IN) :: Ug(m,tnNo), Dg(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(INOUT) :: Ub(m,ifem%tnNo)

      INTEGER(KIND=IKIND) :: a, b, e, g, i, Ac, Bc, Ec, iM, jM, eNoN,
     2   eNoNb
      REAL(KIND=RKIND) :: w, Jac, xi(nsd), up(m), Gmat(nsd,nsd), tt

      REAL(KIND=RKIND), ALLOCATABLE :: Nb(:), Nbx(:,:), xbl(:,:), N(:),
     2   Nx(:,:), xl(:,:), ul(:,:), sA(:)

      ALLOCATE(sA(ifem%tnNo))
      sA = 0._RKIND
      Ub = 0._RKIND

!     Use L2 projection with mass lumping to project flow variables from
!     background fluid mesh to IFEM
      DO iM=1, ifem%nMsh
         eNoNb = ifem%msh(iM)%eNoN
         ALLOCATE(Nb(eNoNb), Nbx(nsd,eNoNb), xbl(nsd,eNoNb))
!        Loop over each trace, as we need to first interpolate flow var
!        at the IFEM integration points based on its trace
         DO i=1, ifem%msh(iM)%trc%nG
            e  = ifem%msh(iM)%trc%gE(1,i)
            g  = ifem%msh(iM)%trc%gE(2,i)
            Ec = ifem%msh(iM)%trc%gptr(1,i)
            jM = ifem%msh(iM)%trc%gptr(2,i)

!           Transfer to local arrays: IFEM mesh variables
            Nb = ifem%msh(iM)%N(:,g)
            IF (ifem%msh(iM)%eType .EQ. eType_NRB)
     2         CALL NRBNNX(ifem%msh(iM), e)
            DO b=1, eNoNb
               Bc = ifem%msh(iM)%IEN(b,e)
               xbl(:,b) = ifem%x(:,Bc) + ifem%Ubo(:,Bc)
            END DO
            CALL GNN(eNoNb, nsd, ifem%msh(iM)%Nx(:,:,g), xbl, Nbx, Jac,
     2         Gmat)
            IF (ISZERO(Jac)) err = " Jac < 0 @ element "//e
            w = ifem%msh(iM)%w(g) * Jac

!           Transfer to local arrays: background mesh variables
            eNoN = msh(jM)%eNoN
            ALLOCATE(N(eNoN), Nx(nsd,eNoN), xl(nsd,eNoN), ul(m,eNoN))
            ul = 0._RKIND
            DO a=1, eNoN
               Ac = msh(jM)%IEN(a,Ec)
               ul(:,a) = Ug(:,Ac)
               xl(:,a) = x(:,Ac)
               IF (mvMsh) xl(:,a) = xl(:,a) + Dg(nsd+2:2*nsd+1,Ac)
            END DO

!           Find shape functions and derivatives on the background mesh
!           at the integration point.
            xi = ifem%msh(iM)%trc%xiG(:,i)
            CALL GETGNN(nsd, msh(jM)%eType, eNoN, xi, N, Nx)

!           Use the computed shape functions to interpolate flow var at
!           the IFEM integration point
            up = 0._RKIND
            DO a=1, eNoN
               up = up + N(a)*ul(:,a)
            END DO

!           Project flow variables to IFEM nodes
            DO b=1, eNoNb
               Bc = ifem%msh(iM)%IEN(b,e)
               Ub(:,Bc) = Ub(:,Bc) + w*Nb(b)*up(:)
               sA(Bc)   = sA(Bc)   + w*Nb(b)
            END DO

            DEALLOCATE(N, Nx, xl, ul)
         END DO
         DEALLOCATE(Nb, Nbx, xbl)
      END DO

!     Synchronize Yb across all the processes
      tt = CPUT()
!      CALL IFEM_SYNCG(Ub)
      CALL IB_SYNCG(Ub)
!     CALL IFEM_SYNCG(sA)
      CALL IB_SYNCG(sA)
      ifem%callD(3) = ifem%callD(3) + CPUT() - tt

      DO a=1, ifem%tnNo
         IF (.NOT.ISZERO(sA(a))) THEN
            Ub(:,a) = Ub(:,a) / sA(a)
         END IF
      END DO

      DEALLOCATE(sA)

      RETURN
      END SUBROUTINE IFEM_INTERPGP
!####################################################################
!     Write IFEM call duration
      SUBROUTINE IFEM_OUTCPUT()
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=RKIND) rtmp
      CHARACTER(LEN=stdL) sOut

      IF (ifem%cpld .EQ. ibCpld_I) THEN
         ifem%callD(1) = ifem%callD(1)/REAL(eq(ifem%cEq)%itr,KIND=RKIND)
         ifem%callD(3) = ifem%callD(3)/REAL(eq(ifem%cEq)%itr,KIND=RKIND)
      END IF
      ifem%callD(3) = ifem%callD(3) + ifem%callD(4)
      ifem%callD(1) = ifem%callD(1) + ifem%callD(2) + ifem%callD(3)

      std = REPEAT("-",69)
      WRITE(sOut,'(F6.2)') ifem%callD(1)
      WRITE(sOut,'(A)') " IFEM call duration: "//TRIM(sOut)//' sec'
      rtmp = 100._RKIND*ifem%callD(3)/ifem%callD(1)
      WRITE(sOut,'(A)') TRIM(sOut)//" (comm."//
     2   STR(NINT(rtmp, KIND=IKIND),3)//"%)"
      rtmp = 100._RKIND*ifem%callD(2)/ifem%callD(1)
      WRITE(sOut,'(A)') TRIM(sOut)//", (updt."//
     2   STR(NINT(rtmp, KIND=IKIND),3)//"%)"
      std = sOut
      std = REPEAT("-",69)

      RETURN
      END SUBROUTINE IFEM_OUTCPUT
!####################################################################
C !     Debugs FSI force on the IFEM
C       SUBROUTINE DEBUGIBR(incNd)
C       USE COMMOD
C       USE ALLFUN
C       IMPLICIT NONE
C       INTEGER(KIND=IKIND), INTENT(IN) :: incNd(ifem%tnNo)

C       INTEGER(KIND=IKIND) :: a, i, Ac, iM, fid
C       REAL(KIND=RKIND) :: s, lo, hi, av
C       CHARACTER(LEN=stdL) :: fName

C       INTEGER(KIND=IKIND), ALLOCATABLE :: lI(:), gI(:)
C       REAL(KIND=RKIND), ALLOCATABLE :: lR(:,:), gR(:,:)

C !     DEBUG ifem%R
C       fid = 1289
C       IF (cm%mas()) THEN
C          WRITE(fName,'(A)') TRIM(appPath)//"dbg_ibR_hist.dat"
C          IF (cTS .EQ. 1) THEN
C             OPEN(fid,FILE=TRIM(fName))
C             CLOSE(fid,STATUS='DELETE')
C          END IF
C          OPEN(fid,FILE=TRIM(fName),POSITION='APPEND')
C          WRITE(fid,'(A)',ADVANCE='NO') STR(cTS)
C       END IF

C       DO iM=1, nMsh
C          ALLOCATE(lR(nsd+1,msh(iM)%nNo), lI(msh(iM)%nNo))
C          IF (cm%mas()) THEN
C             ALLOCATE(gR(nsd+1,msh(iM)%gnNo), gI(msh(iM)%gnNo))
C          ELSE
C             ALLOCATE(gR(0,0), gI(0))
C          END IF
C          lR = 0._RKIND
C          lI = 0
C          DO a=1, msh(iM)%nNo
C             Ac = msh(iM)%gN(a)
C             lR(:,a) = ifem%R(:,Ac)
C             lI(a)   = incNd(Ac)
C          END DO
C          gR = GLOBAL(msh(iM), lR)
C          gI = GLOBAL(msh(iM), lI)
C          DEALLOCATE(lR, lI)

C          IF (cm%mas()) THEN
C !           first, momentum residue
C             av =  0._RKIND
C             lo =  1.E+6_RKIND
C             hi = -1.E+6_RKIND
C             DO a=1, msh(iM)%gnNo
C                IF (gI(a) .EQ. 0) CYCLE
C                s = 0._RKIND
C                DO i=1, nsd
C                   s = s + gR(i,a)**2._RKIND
C                END DO
C                s = SQRT(s)
C                av = av + s
C                IF (s .LT. lo) lo = s
C                IF (s .GT. hi) hi = s
C             END DO
C             av = av / REAL(SUM(gI(:)), KIND=RKIND)
C             WRITE(fid,'(A)',ADVANCE='NO') " "//STR(av)//" "//STR(hi-lo)

C !           Pressure/continuity
C             av =  0._RKIND
C             lo =  1.E+6_RKIND
C             hi = -1.E+6_RKIND
C             DO a=1, msh(iM)%gnNo
C                IF (gI(a) .EQ. 0) CYCLE
C                s = gR(nsd+1,a)
C                av = av + s
C                IF (s .LT. lo) lo = s
C                IF (s .GT. hi) hi = s
C             END DO
C             av = av / REAL(SUM(gI(:)), KIND=RKIND)
C             WRITE(fid,'(A)') " "//STR(av)//" "//STR(hi-lo)
C             CLOSE(fid)
C          END IF
C          DEALLOCATE(gR, gI)
C       END DO

C       iM = cm%reduce(iM)

C       RETURN
C       END SUBROUTINE DEBUGIBR
C !--------------------------------------------------------------------
C !     Debugs IFEM nodal traces
C       SUBROUTINE DEBUGIBNDTRCS(lM)
C       USE COMMOD
C       USE ALLFUN
C       IMPLICIT NONE
C       TYPE(mshType), INTENT(IN) :: lM

C       INTEGER(KIND=IKIND) :: i, a, n, Ac, Ec, iM, fid, ierr
C       CHARACTER(LEN=stdL) :: fName
C       REAL(KIND=RKIND) xp(nsd)

C       INTEGER(KIND=IKIND), ALLOCATABLE :: sCount(:), disps(:), lN(:),
C      2   gN(:), nPtr(:,:)
C       REAL(KIND=RKIND), ALLOCATABLE :: xpL(:,:)

C       ALLOCATE(xpL(nsd,lM%nNo))
C       DO Ac=1, ifem%tnNo
C          a = lM%lN(Ac)
C          xpL(:,a) = ifem%x(:,Ac) + ifem%Ubo(:,Ac)
C       END DO

C       ALLOCATE(sCount(cm%np()), disps(cm%np()))
C       sCount = 0
C       disps  = 0
C       i = lM%trc%n
C       CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
C      2   ierr)

C       n = SUM(sCount(:))
C       sCount = 3*sCount(:)
C       DO i=2, cm%np()
C          disps(i) = disps(i-1) + sCount(i-1)
C       END DO

C       ALLOCATE(lN(3*lM%trc%n), gN(3*n))
C       DO i=1, lM%trc%n
C          lN(3*i-2) = lM%trc%gN(i)
C          lN(3*i-1) = lM%trc%nptr(1,i)
C          lN(3*i)   = lM%trc%nptr(2,i)
C       END DO

C       i = 3*lM%trc%n
C       CALL MPI_GATHERV(lN, i, mpint, gN, sCount, disps, mpint, master,
C      2   cm%com(), ierr)

C       DEALLOCATE(lN, disps, sCount)

C       WRITE(fName,'(A)') TRIM(appPath)//"dbg_IFEM_nd_trc_"//STR(cTS)//
C      2   ".dat"
C       IF (cm%mas()) THEN
C          ALLOCATE(nPtr(2,lM%nNo))
C          nPtr = 0
C          DO i=1, n
C             a  = gN(3*i-2)
C             Ec = gN(3*i-1)
C             iM = gN(3*i)
C             nPtr(1,a) = Ec
C             nPtr(2,a) = iM
C          END DO

C          fid = 1289
C          OPEN(fid,FILE=TRIM(fName))
C          WRITE(fid,'(A)') "List of failed traces on mesh: "//
C      2      TRIM(lM%name)
C          DO a=1, lM%nNo
C             Ac = lM%gN(a)
C             Ec = nPtr(1,a)
C             iM = nPtr(2,a)
C             IF (Ec.EQ.0 .OR. iM.EQ.0) THEN
C                xp = xpL(:,a)
C                WRITE(fid,'(2X,A)',ADVANCE='NO') STR(a)//" "//STR(Ac)
C                DO i=1, nsd
C                   WRITE(fid,'(A)',ADVANCE='NO') " "//STR(xp(i))
C                END DO
C                WRITE(fid,'(A)')
C             END IF
C          END DO
C          CLOSE(fid)
C       ELSE
C          ALLOCATE(nPtr(0,0))
C       END IF
C       DEALLOCATE(xpL, gN, nPtr)

C       i = cm%reduce(i)
C       err = "ERROR: Failed to detect all the nodal traces on "//
C      2   TRIM(lM%name)//". See "//TRIM(fName)//" for more information."

C       RETURN
C       END SUBROUTINE DEBUGIBNDTRCS
C !--------------------------------------------------------------------
C !     Debugs IFEM integration point traces
C       SUBROUTINE DEBUGIBGPTRCS(lM)
C       USE COMMOD
C       USE ALLFUN
C       IMPLICIT NONE
C       TYPE(mshType), INTENT(IN) :: lM

C       INTEGER(KIND=IKIND) :: i, a, e, g, n, Ac, Ec, iM, fid, ierr
C       CHARACTER(LEN=stdL) :: fName
C       REAL(KIND=RKIND) xp(nsd)

C       INTEGER(KIND=IKIND), ALLOCATABLE :: sCount(:), disps(:), lE(:),
C      2   gE(:), eptr(:,:,:)
C       REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:)

C       ALLOCATE(sCount(cm%np()), disps(cm%np()))
C       sCount = 0
C       disps  = 0
C       i = lM%trc%nG
C       CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
C      2   ierr)

C       n = SUM(sCount(:))
C       sCount = 4*sCount(:)
C       DO i=2, cm%np()
C          disps(i) = disps(i-1) + sCount(i-1)
C       END DO

C       ALLOCATE(lE(4*lM%trc%nG), gE(4*n))
C       DO i=1, lM%trc%nG
C          lE(4*i-3) = lM%trc%gE(1,i)
C          lE(4*i-2) = lM%trc%gE(2,i)
C          lE(4*i-1) = lM%trc%gptr(1,i)
C          lE(4*i)   = lM%trc%gptr(2,i)
C       END DO

C       i = 4*lM%trc%nG
C       CALL MPI_GATHERV(lE, i, mpint, gE, sCount, disps, mpint, master,
C      2   cm%com(), ierr)

C       DEALLOCATE(lE, disps, sCount)

C       WRITE(fName,'(A)') TRIM(appPath)//"dbg_IFEM_gp_trc_"//STR(cTS)//
C      2   ".dat"
C       IF (cm%mas()) THEN
C          ALLOCATE(eptr(2,lM%nG,lM%nEl))
C          eptr = 0
C          DO i=1, n
C             e  = gE(4*i-3)
C             g  = gE(4*i-2)
C             Ec = gE(4*i-1)
C             iM = gE(4*i)
C             ePtr(1,g,e) = Ec
C             ePtr(2,g,e) = iM
C          END DO

C          fid = 1289
C          OPEN(fid,FILE=TRIM(fName))
C          WRITE(fid,'(A)') "List of failed traces on mesh: "//
C      2      TRIM(lM%name)
C          ALLOCATE(xl(nsd,lM%eNoN))
C          DO e=1, lM%nEl
C             DO a=1, lM%eNoN
C                Ac = lM%IEN(a,e)
C                xl(:,a) = ifem%x(:,Ac) + ifem%Ubo(:,Ac)
C             END DO
C             DO g=1, lM%nG
C                Ec = ePtr(1,g,e)
C                iM = ePtr(2,g,e)
C                IF (Ec.EQ.0 .OR. iM.EQ.0) THEN
C                   xp = 0._RKIND
C                   DO a=1, lM%eNoN
C                      xp(:) = xp(:) + lM%N(a,g)*xl(:,a)
C                   END DO
C                   WRITE(fid,'(2X,A)',ADVANCE='NO') STR(g)//" "//STR(e)
C                   DO i=1, nsd
C                      WRITE(fid,'(A)',ADVANCE='NO') " "//STR(xp(i))
C                   END DO
C                   WRITE(fid,'(A)')
C                END IF
C             END DO
C          END DO
C          DEALLOCATE(xl)
C          CLOSE(fid)
C       ELSE
C          ALLOCATE(ePtr(0,0,0))
C       END IF
C       DEALLOCATE(gE, ePtr)

C       i = cm%reduce(i)
C       err = "ERROR: Failed to detect all the traces on "//
C      2   TRIM(lM%name)//". See "//TRIM(fName)//" for more information."

C       RETURN
C       END SUBROUTINE DEBUGIBGPTRCS
!####################################################################
