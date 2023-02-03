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
!     This routine applies penalty-based contact model for multiple 
!     thick structures with large displacement 
!
!-----------------------------------------------------------------------

      SUBROUTINE EVAL_CNT_PATTERN !(Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
C       REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      REAL(KIND=RKIND) :: Dg(tDof,tnNo)
      INTEGER(KIND=IKIND) :: iMS, iMM, iBndS, eNoNb, eNoN, elmS, a, Ac, 
     2                       elmM, g, iBndM, Int, nbrCntElm
      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:), ptrM(:) 
      REAL(KIND=RKIND), ALLOCATABLE :: xbl(:,:), dbl(:,:), sN(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xblM(:,:), dblM(:,:), mN(:)
      REAL(KIND=RKIND) :: xqp(nsd), Jac, IntP(nsd)
      TYPE(faceType) :: iFaS, iFaM

      REAL :: quadPntCurr(3,nsd)

!     Let's consider for the moment that the first solid is the slave
!     and the second is the master. Only make solid to master detection first 

!     Define id mesh master and slave
      iMS = 1
      iMM = 2

      Dg = 0._RKIND

!     Allocate here temporary quantities, this should be done elsewhere in 
!     the future. We are considering initially that only 2 faces can 
!     get in contact, and that it's only the first face in each mesh. 
!     We can send matching meshes ids from file later

!     For the moment we start with only slave contact. One master and one slave 

      nbrCntElm = msh(iMS)%fa(1)%nEl * msh(iMS)%fa(1)%nG !+ msh(iMM)%fa(1)%nEl
      ALLOCATE( faceContElm(nbrCntElm), volContElm(nbrCntElm,2), 
     2                                      intPointCnt(nbrCntElm,nsd) )

      faceContElm = -1 
      volContElm = -1 
      intPointCnt = -1 

!     Just to plot the quan pont 

!     Loop over slave boundary element 
      DO iBndS = 1, msh(iMS)%nFa

         IF( iBndS .GE. 2 ) CYCLE

         iFaS = msh(iMS)%fa(iBndS)
         eNoNb = iFaS%eNoN

         eNoN = msh(iMS)%eNoN

C          write(*,*)" Face info: "
C          write(*,*)" eNoN face ", eNoNb
C          write(*,*)" nG ", iFaS%nG
C          write(*,*)" nEl ", iFaS%nEl
C          write(*,*)" mesh eNoN ", eNoN

!        Loop on element and compute the normal 
         DO elmS = 1, iFaS%nEl

            ALLOCATE(ptr(eNoN), xbl(nsd,eNoN), dbl(nsd,eNoN), sN(nsd), 
     2              ptrM(eNoN), xblM(nsd,eNoN), dblM(nsd,eNoN), mN(nsd))

            DO a=1, eNoNb
               Ac      = iFaS%IEN(a,elmS)
               ptr(a)  = Ac
               xbl(:,a) = x(:,Ac)
               dbl(:,a) = Dg(:,Ac)
            END DO

!           Compute the coordinates and the normal vector in the current 
!           configuration
            xbl = xbl + dbl

           CALL NORMAL(xbl,eNoNb,sN)

            !        Gauss integration 1
            DO g=1, iFaS%nG
!               This is not working here since we need to compute the normal 
!               in the current configuration, while GNNB computes it in the initial config       
C                CALL GNNB(iFaS,elmS, g, nsd-1, eNoNb, iFaS%Nx(:,:,g), sN)
C                Jac = SQRT(NORM(sN))
C                sN  = sN/Jac

!              For each quadrature point find the coordinates in the current 
!              configuration 
               xqp = 0._RKIND
               DO a=1, eNoNb
                  xqp = xqp + xbl(:,a)*iFaS%N(a,g)
               END DO

               IF(elmS .EQ. 1 ) THEN 
                  quadPntCurr(g,:) = xqp
               END IF 

C                write(*,*)" Cooridnates gp current slav: ", xqp

!              Loop over master boundary element 
               DO iBndM = 1, msh(iMM)%nFa

C                   write(*,*)" Inside master face ", iBndM

                  IF( iBndM .GE. 2 ) CYCLE
                  iFaM = msh(iMM)%fa(iBndM)

                  DO elmM = 1, iFaM%nEl

C                      write(*,*)" Looking at master element face ", elmM
!                    Get master element coordinates    
                     DO a=1, eNoNb
                        Ac        = iFaM%IEN(a,elmM)
                        ptrM(a)   = Ac
                        xblM(:,a) = x(:,Ac)
                        dblM(:,a) = Dg(:,Ac)
                     END DO

!                    Find element coordinates in the current configuration 
!                    and search for possible contact, at the same time compute the 
!                    gap. The final contact is the closest detected 
                     xblM = xblM + dblM

!                    Same here: 
C                      CALL GNNB(iFaM, elmM, g, nsd-1, eNoNb, 
C      2                                               iFaM%Nx(:,:,g), mN)
C                      Jac = SQRT(NORM(mN))
C                      mN  = mN/Jac

                     CALL NORMAL(xblM,eNoNb,mN)

C                      write(*,*)" Normal master: ", mN

C                      IF( nsd .EQ. 3) THEN 
                        IntP = 0._RKIND
                        CALL RAYTOTRIINTER( xblM, eNoN, sN, xqp, Int, 
     2                                                            IntP )
C                      ELSE 
C                         CALL RAYTOSEGINTER()
C                      END IF

C                      write(*,*)" Int = ", Int
                     IF( Int .EQ. 1) THEN 
C                         write(*,*)" Found possible contact: " 
C                         write(*,*)" Master coordinates element: "
C                         write(*,*)" ", xblM(:,1)
C                         write(*,*)" ", xblM(:,2)
C                         write(*,*)" ", xblM(:,3)
C                         write(*,*)""   
C                         write(*,*)" xqp = ", xqp
C                         write(*,*)" direction = ", sN
C                         write(*,*)""
C                         write(*,*)" Insersection point ", IntP
C                         write(*,*)""

                        faceContElm((elmS-1)*iFaS%nG+g) = elmM  
                        volContElm((elmS-1)*iFaS%nG+g,1) = iFaS%gE(elmS)    
                        volContElm((elmS-1)*iFaS%nG+g,2) = iFaM%gE(elmM)    
                        intPointCnt((elmS-1)*iFaS%nG+g,:) = IntP
                     END IF

                  END DO 
               END DO
            END DO

            DEALLOCATE(ptr, xbl, dbl, sN)
            DEALLOCATE(ptrM, xblM, dblM, mN)
         END DO
      END DO

C       write(*,*)" faceContElm: " , faceContElm
C       write(*,*)" volContElm slave: " , volContElm(:,1)
C       write(*,*)" volContElm master: " , volContElm(:,2)
C       write(*,*)" First 3 int point intPointCnt: " 
C       write(*,*)" Quap point " , quadPntCurr(1,:)
C       write(*,*)" " , intPointCnt(1,:)
C       write(*,*)" Quap point " , quadPntCurr(2,:)
C       write(*,*)" " , intPointCnt(2,:)
C       write(*,*)" Quap point " , quadPntCurr(3,:)
C       write(*,*)" " , intPointCnt(3,:)

      RETURN
      END SUBROUTINE EVAL_CNT_PATTERN

!#######################################################################
!     Routine to compute the normal of a face given the coordinates
      SUBROUTINE NORMAL(xFace,eNoN,NCur)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: xFace(nsd,eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: NCur(nsd)

      REAL(KIND=RKIND) n1, n2, n3, nrm, covBas(nsd-1,nsd)


      IF( nsd .EQ. 2) THEN 
         n1 =   xFace(2,2) - xFace(2,1) 
         n2 = - (xFace(1,2) - xFace(1,1))
         nrm = SQRT(n1*n1 + n2*n2)
         NCur(1) = n1/nrm
         NCur(2) = n2/nrm

C          write(*,*)" The normal of face ", xFace(:,1), " ",xFace(:,2)
C          write(*,*)" is : ", NCur
      ELSE 
         covBas(1,:) = xFace(:,2) - xFace(:,1) 
         covBas(2,:) = xFace(:,3) - xFace(:,1) 
         n1 = covBas(1,2)*covBas(2,3)-covBas(1,3)*covBas(2,2);
         n2 = covBas(1,3)*covBas(2,1)-covBas(1,1)*covBas(2,3);
         n3 = covBas(1,1)*covBas(2,2)-covBas(1,2)*covBas(2,1);
         nrm = SQRT(n1*n1 + n2*n2 + n3*n3)
         NCur(1) = n1/nrm
         NCur(2) = n2/nrm
         NCur(3) = n3/nrm

C          write(*,*)" The normal of face ", xFace(:,1)
C          write(*,*)" ",xFace(:,2)
C          write(*,*) " ",xFace(:,3)
C          write(*,*) " is : ", NCur
      END IF

      RETURN
      END SUBROUTINE NORMAL

!#######################################################################
!     Routine to compute the 3-dimensional intersection of a ray with a triangle 
      SUBROUTINE RAYTOTRIINTER( xl, eNoN, dir, P0, Int, IntP )
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: xl(nsd,eNoN), dir(nsd), P0(nsd)
      REAL(KIND=RKIND), INTENT(OUT) :: IntP(nsd)
      INTEGER(KIND=IKIND), INTENT(OUT) :: Int

      REAL(KIND=RKIND) u(nsd), v(nsd), w(nsd), n(nsd), w0(nsd), 
     2     uxv(nsd, 2), nrm, a, b, rdir, uu, uv, vv, wu, wv, D, s, t

      Int = 1

      u = xl(:,2) - xl(:,1)
      v = xl(:,3) - xl(:,1)
      w0 = P0 - xl(:,1)

      uxv(:,1) = u 
      uxv(:,2) = v 

      n = CROSS(uxv)

      nrm = NORM(n,n)
!     Check if the triangle is degenerate      
      IF ( SQRT(nrm) .LT. 1.E-6 ) THEN 
         Int = -1 
      END IF
        
      a = -NORM(n,w0);
      b = NORM(n,dir);
!     Check if the ray is  parallel to the triangle plane      
      IF (ABS(b) .LT. 1.E-6) THEN 
!        The ray lies in triangle plane         
         IF (a .EQ. 0) THEN 
            Int = 2     
         ELSE  
!           The ray is disjoint from the plane     
            Int = 0
         END IF
      END IF            
   
      ! Get intersect point of ray with triangle plane       
      rdir = a / b

      ! Get intersect point between the ray and the plane described from the 3 vertices 
      IntP = P0 + rdir*dir

      ! Now let understand whether the point is inside the triangle or not 
       
      uu = NORM(u,u)
      uv = NORM(u,v)
      vv = NORM(v,v)

      w = IntP - xl(:,1)

      wu = NORM(w,u)
      wv = NORM(w,v)
      D = uv * uv - uu * vv

!     Get and test parametric coords
      s = (uv * wv - vv * wu) / D
      IF ((s .LT. 0._RKIND) .OR. (s .GT. 1._RKIND)) THEN 
!        IntP is outside the triangle
         Int = 0
      END IF
      
      t = (uv * wu - uu * wv) / D
      IF ((t .LT. 0._RKIND) .OR. ((s + t) .GT. 1._RKIND)) THEN 
!        IntP is outside the triangle
         Int = 0
      END IF

      RETURN
      END SUBROUTINE RAYTOTRIINTER
!#######################################################################

      SUBROUTINE NL_CONTACTFORCES(Yg,Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) :: a, iMS, iMM, iFaC, Ac, Ec, eNoN, eNoNb, g, 
     2          iFaS, i, j, k, b

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptrM(:), ptrS(:)

      REAL(KIND=RKIND) :: w, Jac, yp(nsd), xp(nsd), xpRef(nsd), 
     2          xpPar(nsd), xi0(nsd), Nrm(nsd), sN(nsd), mN(nsd), 
     3          ypRef(nsd), Ks(nsd,nsd), gap0, gap, gamma, afu, aux, 
     4          nxny, Fx(nsd,nsd), Tnx(nsd,nsd), Fxi(nsd,nsd), 
     5          Ox(nsd,nsd), nnT(nsd,nsd), sNsN(nsd,1), ElastM

      REAL(KIND=RKIND), ALLOCATABLE :: xblM(:,:), 
     2          dblM(:,:), xbl(:,:), dbl(:,:), xl(:,:), dl(:,:),
     3          xlM(:,:),  yl(:,:), ylM(:,:), s_lR(:,:), sm_lK(:,:,:), 
     4          s_Nw(:), s_Nwx(:,:), s_Nwxi(:,:), m_Nw(:), m_Nwx(:,:), 
     5          m_Nwxi(:,:), m_lR(:,:), ss_lK(:,:,:), ms_lK(:,:,:), 
     6          mm_lK(:,:,:), dlM(:,:) 

      TYPE(fsType) :: fs(2)
      TYPE(faceType) :: FacS, FacM
      LOGICAL :: cntActive

      write(*,*)" Inside NL_CONTACTFORCES: "

      iMS  = 1  
      iMM  = 2  
      iFaC = 1 

      ElastM = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)

      gap0 = 0.05_RKIND
C       gamma = 1._RKIND*ElastM/(0.13**2_RKIND)
      gamma = 1.0E6
      cntActive = .FALSE.

      afu     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i       = eq(cEq)%s
      j       = i + 1
      k       = j + 1

!     In the future we will have a loop over all the meshes, or at least 
!     the coupled (with possible contact) meshes, and a loop on their faces 

!     Before the loop over the quadrature point for each slave element 
!     we need to prepare the coupling, unfortunately this is changing, 
!     so we can't do it before, or we can recompute it if it is different. 

!     We can get all the info of the slave element before the loop over the 
!     quadrature point, while for the master we have to compute everything 
!     inside the loop 

!     Loop over slave contact face element 

      FacS  = msh(iMS)%fa(iFaC) 
      eNoNb = FacS%eNoN
      eNoN  = msh(iMS)%eNoN
      FacM  = msh(iMM)%fa(iFaC) 

!     Set function spaces for velocity and pressure on mesh
!     We suppose the slave and master func spaces to be equal 
      CALL GETTHOODFS(fs, msh(iMS), .TRUE., 1)

      ALLOCATE( xbl(nsd,eNoNb), dbl(nsd,eNoNb), ptrM(eNoN),
     2    xblM(nsd,eNoNb), dblM(nsd,eNoNb), xlM(nsd,eNoN), 
     3    ptrS(eNoN), xl(nsd,eNoN), yl(tDof,eNoN), s_lR(dof,eNoN),
     4    ss_lK(dof*dof,eNoN,eNoN), s_Nw(fs(1)%eNoN), dlM(nsd,eNoN), 
     5    ylM(tDof,eNoN), s_Nwxi(nsd,fs(1)%eNoN), m_Nw(fs(1)%eNoN), 
     6    m_Nwx(nsd,fs(1)%eNoN), m_Nwxi(nsd,fs(1)%eNoN), dl(tDof,eNoN),
     7    sm_lK(dof*dof,eNoN,eNoN), s_Nwx(nsd,fs(1)%eNoN) )


      DO iFaS = 1, FacS%nEl

!        Get element coordinates in current and reference configuration

         DO a=1, eNoNb
            Ac       = FacS%IEN(a,iFaS)
            xbl(:,a) = x(:,Ac)
            dbl(:,a) = Dg(:,Ac)
         END DO

!        Create local copies of fluid velocity and position vector
         Ec = FacS%gE(iFaS)
         DO a=1, eNoN
            Ac = msh(iMS)%IEN(a,Ec)
            ptrS(a)  = Ac
            xl(:,a) = x(:,Ac)
            !! check this, do we need x at the ref or current configuration?
C             IF (mvMsh) xl(:,a) = xl(:,a) + Dg(nsd+2:2*nsd+1,Ac)
!           xl must be in the current configuration 
C             xl(:,a) = xl(:,a) + Dg(1:nsd,Ac)
            yl(:,a) = Yg(:,Ac)
            dl(:,a) = Dg(1:nsd,Ac)
         END DO

!        Compute the boundary coordinates and the normal vector in the current 
!        configuration
         xbl = xbl + dbl
         CALL NORMAL(xbl,eNoNb,sN) 
         xbl = xbl - dbl

!        Define initial guess for the Newton's iterations, called during the  
!        localization of a point in the parent element, from the 
!        current configuration 
         xi0 = 0._RKIND
         DO g=1, fs(1)%nG
            xi0 = xi0 + fs(1)%xi(:,g)
         END DO
         xi0 = xi0 / REAL(fs(1)%nG, KIND=RKIND)

!        Local residue initialization 
         s_lR = 0._RKIND

!        Loop over quadrature points 
         DO g=1, FacS%nG

!           Local tanget contrib initialization 
            ss_lK = 0._RKIND
            sm_lK = 0._RKIND
            cntActive = .FALSE.


!-----      Slave element dependencies 

C             write(*,*)" Beginning Slave " 

!           Get element normal, weight
            CALL GNNB(FacS, iFaS, g, nsd-1, eNoNb, FacS%Nx(:,:,g), Nrm)
            Jac = SQRT(NORM(Nrm))
            Nrm  = Nrm/Jac
            w   = FacS%w(g) * Jac

!           For each quadrature point find the coordinates in the current 
!           configuration. This is x. 
            xbl = xbl + dbl
            xp = 0._RKIND
            DO a=1, eNoNb
               xp = xp + xbl(:,a)*FacS%N(a,g)
            END DO

            xbl = xbl - dbl
            xpRef = 0._RKIND
            DO a=1, eNoNb
               xpRef = xpRef + xbl(:,a)*FacS%N(a,g)
            END DO


!           Compute Nw and Nwxi of the mesh at the face integration
!           point using Newton method. Then calculate Nwx.
            xpPar = xi0
            CALL GETNNX(fs(1)%eType, fs(1)%eNoN, xl, fs(1)%xib,
     2         fs(1)%Nb, xpRef, xpPar, s_Nw, s_Nwxi)
            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF)
     2         CALL GNN(fs(1)%eNoN, nsd, s_Nwxi, xl, s_Nwx, Jac, Ks)


C !---  PRINT SLAVE INFOS
C             write(*,*)""
C             write(*,*)" Slave info print "
C             write(*,*)" x is ", xp 
C             write(*,*)" normal is ", sN
C             write(*,*)" element coord bnd : "
C             write(*,*)" ", xbl(:,1) 
C             write(*,*)" ", xbl(:,2) 
C             write(*,*)" ", xbl(:,3) 
C             write(*,*)" volume element coord : "
C             write(*,*)" ", xl(:,1) 
C             write(*,*)" ", xl(:,2) 
C             write(*,*)" ", xl(:,3) 
C             write(*,*)" ", xl(:,4) 
C             write(*,*)" quad point parent elem ", xpPar
C             write(*,*)" value Nw in parent elem ", s_Nw
C             write(*,*)" value Nwxi,1 in parent elem ", s_Nwxi(1,:)
C !---  END PRINT SLAVE 


!-----      Master element dependencies 

C             write(*,*)" Beginning Master "

!           Get y in the current configuration 
            yp = intPointCnt((iFaS-1)*FacS%nG + g, :) 

!           Get the master id face and element coordinates   
            DO a=1, eNoNb
               Ac       = FacM%IEN(a,faceContElm((iFaS-1)*FacS%nG + g))
               xblM(:,a) = x(:,Ac)
               dblM(:,a) = Dg(:,Ac)
            END DO

!           Compute normal in the current configuration 
            xblM = xblM + dblM
            CALL NORMAL(xblM,eNoNb,mN) 

!           Create local copies of fluid velocity and position vector
            Ec = volContElm((iFaS-1)*FacS%nG + g, 2)
            DO a=1, eNoN
               Ac = msh(iMM)%IEN(a,Ec)
               ptrM(a)  = Ac
               xlM(:,a) = x(:,Ac)
               !! check this, do we need x at the ref or current configuration?
               xlM(:,a) = xlM(:,a) + Dg(1:nsd,Ac)
               ylM(:,a) = Yg(:,Ac)
               dlM(:,a) = Dg(1:nsd,Ac)
            END DO

!           Compute Nw and Nwxi of the master mesh at the face integration
!           point using Newton method, starting from the current configuration
            ypRef = xi0
            CALL GETNNX(fs(1)%eType, fs(1)%eNoN, xlM, fs(1)%xib,
     2         fs(1)%Nb, yp, ypRef, m_Nw, m_Nwxi)

!           Evaluate test gradient test function in the 
!           parent element, starting from reference configuration
            xlM = xlM - dlM
            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF)
     2         CALL GNN(fs(1)%eNoN, nsd, m_Nwxi, xlM, m_Nwx, Jac, Ks)


!---  PRINT MASTER INFOS
C             write(*,*)""
C             write(*,*)" Master info print "
C             write(*,*)" y is ", yp
C             write(*,*)" normal is ", mN
C             write(*,*)" element coord bnd : "
C             write(*,*)" ", xblM(:,1) 
C             write(*,*)" ", xblM(:,2) 
C             write(*,*)" ", xblM(:,3) 
C             write(*,*)" volume element coord : "
C             write(*,*)" ", xlM(:,1) 
C             write(*,*)" ", xlM(:,2) 
C             write(*,*)" ", xlM(:,3) 
C             write(*,*)" ", xlM(:,4) 
C             write(*,*)" quad point parent elem ", ypRef
C             write(*,*)" value Nw in parent elem ", m_Nw
C             write(*,*)" value Nwxi,1 in parent elem ", m_Nwxi(1,:)
!---  END PRINT SLAVE 

!           We should now have all the necessary information. 
!           We can proceed with the assembly               

! ------------------------ Beginning assembly --------------------------
! ----------------------------------------------------------------------

! ------------------------ Gap computation 
!           gap = ( y  - x ) . nx 
            gap = 0._RKIND
            DO a=1, nsd
               gap = gap + (yp(a) - xp(a))*sN(a)
            END DO

            IF( gap .LT. gap0 ) THEN 
               cntActive = .TRUE.
               write(*,*)" ****  Contact active **** "
               write(*,*)" gap = ", gap
            ELSE
               CYCLE
            END IF

!           If we are here, it's because contact is active at this quad point
! ------------------------ Residue assembly 

            DO a=1, eNoN
               s_lR(1,a) = s_lR(1,a) + gamma*w*s_Nw(a)*(gap0-gap)*sN(1)
               s_lR(2,a) = s_lR(2,a) + gamma*w*s_Nw(a)*(gap0-gap)*sN(2)
               s_lR(3,a) = s_lR(3,a) + gamma*w*s_Nw(a)*(gap0-gap)*sN(3)
            END DO

! ------------------------ T1_1
! \int Heaviside(g0-g)/(nx.ny)   ( u(X) - u(Y) ) . ny  * ( v(X) - v(Y) ) . nx 

C !           Compute the displacement at X and Y            
C             DO a=1, eNoN
C                dX(1) = s_Nw(a)*dl(i,a)
C                dX(2) = s_Nw(a)*dl(j,a)
C                dX(3) = s_Nw(a)*dl(k,a)

C                dY(1) = m_Nw(a)*dlM(i,a)
C                dY(2) = m_Nw(a)*dlM(j,a)
C                dY(3) = m_Nw(a)*dlM(k,a)
C             END DO

!           Compute nx.ny
            nxny = NORMS(sN, mN)
          
            DO a = 1, eNoN ! row
               DO b = 1, eNoN ! col

!---              v1 equation
!                 u1(X), v1(X)
                  aux = (gamma*s_Nw(b)*mN(1)*s_Nw(a)*sN(1)) / nxny
                  ss_lK(1,a,b) = ss_lK(1,a,b) + w * afu * aux 
!                 u1(Y), v1(X)
                  aux = - (gamma*m_Nw(b)*mN(1)*s_Nw(a)*sN(1)) / nxny
                  sm_lK(1,a,b) = sm_lK(1,a,b) + w * afu * aux 

!                 u2(X), v1(X)
                  aux = (gamma*s_Nw(b)*mN(2)*s_Nw(a)*sN(1)) / nxny
                  ss_lK(2,a,b) = ss_lK(2,a,b) + w * afu * aux 
!                 u2(Y), v1(X)
                  aux = - (gamma*m_Nw(b)*mN(2)*s_Nw(a)*sN(1)) / nxny
                  sm_lK(2,a,b) = sm_lK(2,a,b) + w * afu * aux 

!                 u3(X), v1(X)
                  aux = (gamma*s_Nw(b)*mN(3)*s_Nw(a)*sN(1)) / nxny
                  ss_lK(3,a,b) = ss_lK(3,a,b) + w * afu * aux 
!                 u3(Y), v1(X)
                  aux = - (gamma*m_Nw(b)*mN(3)*s_Nw(a)*sN(1)) / nxny
                  sm_lK(3,a,b) = sm_lK(3,a,b) + w * afu * aux 

!---              v2 equation
!                 u1(X), v2(X)
                  aux = (gamma*s_Nw(b)*mN(1)*s_Nw(a)*sN(2)) / nxny
                  ss_lK(4,a,b) = ss_lK(4,a,b) + w * afu * aux 
!                 u1(Y), v2(X)
                  aux = - (gamma*m_Nw(b)*mN(1)*s_Nw(a)*sN(2)) / nxny
                  sm_lK(4,a,b) = sm_lK(4,a,b) + w * afu * aux 

!                 u2(X), v2(X)
                  aux = (gamma*s_Nw(b)*mN(2)*s_Nw(a)*sN(2)) / nxny
                  ss_lK(5,a,b) = ss_lK(5,a,b) + w * afu * aux 
!                 u2(Y), v2(X)
                  aux = - (gamma*m_Nw(b)*mN(2)*s_Nw(a)*sN(2)) / nxny
                  sm_lK(5,a,b) = sm_lK(5,a,b) + w * afu * aux 

!                 u3(X), v2(X)
                  aux = (gamma*s_Nw(b)*mN(3)*s_Nw(a)*sN(2)) / nxny
                  ss_lK(6,a,b) = ss_lK(6,a,b) + w * afu * aux 
!                 u3(Y), v2(X)
                  aux = - (gamma*m_Nw(b)*mN(3)*s_Nw(a)*sN(2)) / nxny
                  sm_lK(6,a,b) = sm_lK(6,a,b) + w * afu * aux 

!---              v3 equation
!                 u1(X), v3(X)
                  aux = (gamma*s_Nw(b)*mN(1)*s_Nw(a)*sN(3)) / nxny
                  ss_lK(7,a,b) = ss_lK(7,a,b) + w * afu * aux 
!                 u1(Y), v3(X)
                  aux = - (gamma*m_Nw(b)*mN(1)*s_Nw(a)*sN(3)) / nxny
                  sm_lK(7,a,b) = sm_lK(7,a,b) + w * afu * aux 

!                 u2(X), v3(X)
                  aux = (gamma*s_Nw(b)*mN(2)*s_Nw(a)*sN(3)) / nxny
                  ss_lK(8,a,b) = ss_lK(8,a,b) + w * afu * aux 
!                 u2(Y), v3(X)
                  aux = - (gamma*m_Nw(b)*mN(2)*s_Nw(a)*sN(3)) / nxny
                  sm_lK(8,a,b) = sm_lK(8,a,b) + w * afu * aux 

!                 u3(X), v3(X)
                  aux = (gamma*s_Nw(b)*mN(3)*s_Nw(a)*sN(3)) / nxny
                  ss_lK(9,a,b) = ss_lK(9,a,b) + w * afu * aux 
!                 u3(Y), v3(X)
                  aux = - (gamma*m_Nw(b)*mN(3)*s_Nw(a)*sN(3)) / nxny
                  sm_lK(9,a,b) = sm_lK(9,a,b) + w * afu * aux 
               END DO
            END DO 

! ------------------------ T1_2
! coef / h^2 alfa \int H(g0-g) g / ( (nx.ny) | a_3(x) |  ) ny. Tnx . (Orth2 . u_r - Oth1.ur )   * ( v(X) - V(Y) ) . nx  

            Fx      = 0._RKIND
            Fx(1,1) = 1._RKIND
            Fx(2,2) = 1._RKIND
            Fx(3,3) = 1._RKIND
            Tnx = Fx
            DO a = 1, eNoN
               Fx(1,1) = Fx(1,1) + s_Nwx(1,a)*dl(i,a)
               Fx(1,2) = Fx(1,2) + s_Nwx(2,a)*dl(i,a)
               Fx(1,3) = Fx(1,3) + s_Nwx(3,a)*dl(i,a)
               Fx(2,1) = Fx(2,1) + s_Nwx(1,a)*dl(j,a)
               Fx(2,2) = Fx(2,2) + s_Nwx(2,a)*dl(j,a)
               Fx(2,3) = Fx(2,3) + s_Nwx(3,a)*dl(j,a)
               Fx(3,1) = Fx(3,1) + s_Nwx(1,a)*dl(k,a)
               Fx(3,2) = Fx(3,2) + s_Nwx(2,a)*dl(k,a)
               Fx(3,3) = Fx(3,3) + s_Nwx(3,a)*dl(k,a)
            END DO
            Fxi = MAT_INV(Fx, 3)

            sNsN(:,1) = sN
            nnT = MATMUL(sNsN,TRANSPOSE(sNsN))
            Tnx = Tnx - nnT

            Ox = MATMUL(Tnx, TRANSPOSE(Fxi))

            DO a = 1, eNoN ! row
               DO b = 1, eNoN ! col

!---              v1 equation
!                 u1(X), v1(X)
                  aux = (mN(1)*Ox(1,1) + mN(2)*Ox(2,1) + 
     2               mN(3)*Ox(3,1)) * s_Nwx(1,b)*sN(1) * s_Nw(a)*sN(1) 
                  aux = gamma * gap * aux / nxny 
                  ss_lK(1,a,b) = ss_lK(1,a,b) + w * afu * aux            

!                 u2(X), v1(X)
                  aux = (mN(1)*Ox(1,1) + mN(2)*Ox(2,1) + 
     2               mN(3)*Ox(3,1)) * s_Nwx(1,b)*sN(2) * s_Nw(a)*sN(1) 
                  aux = gamma * gap * aux / nxny 
                  ss_lK(2,a,b) = ss_lK(2,a,b) + w * afu * aux 

!                 u3(X), v1(X)
                  aux = (mN(1)*Ox(1,1) + mN(2)*Ox(2,1) + 
     2               mN(3)*Ox(3,1)) * s_Nwx(1,b)*sN(3) * s_Nw(a)*sN(1) 
                  aux = gamma * gap * aux / nxny 
                  ss_lK(3,a,b) = ss_lK(3,a,b) + w * afu * aux 

!---              v2 equation
!                 u1(X), v2(X)
                  aux = (mN(1)*Ox(1,1) + mN(2)*Ox(2,1) + 
     2               mN(3)*Ox(3,1)) * s_Nwx(2,b)*sN(1) * s_Nw(a)*sN(2) 
                  aux = gamma * gap * aux / nxny 
                  ss_lK(4,a,b) = ss_lK(4,a,b) + w * afu * aux 

!                 u2(X), v2(X)
                  aux = (mN(1)*Ox(1,1) + mN(2)*Ox(2,1) + 
     2               mN(3)*Ox(3,1)) * s_Nwx(2,b)*sN(2) * s_Nw(a)*sN(2) 
                  aux = gamma * gap * aux / nxny 
                  ss_lK(5,a,b) = ss_lK(5,a,b) + w * afu * aux 

!                 u3(X), v2(X)
                  aux = (mN(1)*Ox(1,1) + mN(2)*Ox(2,1) + 
     2               mN(3)*Ox(3,1)) * s_Nwx(2,b)*sN(3) * s_Nw(a)*sN(2) 
                  aux = gamma * gap * aux / nxny 
                  ss_lK(6,a,b) = ss_lK(6,a,b) + w * afu * aux 

!---              v3 equation
!                 u1(X), v3(X)
                  aux = (mN(1)*Ox(1,1) + mN(2)*Ox(2,1) + 
     2               mN(3)*Ox(3,1)) * s_Nwx(3,b)*sN(1) * s_Nw(a)*sN(3) 
                  aux = gamma * gap * aux / nxny 
                  ss_lK(7,a,b) = ss_lK(7,a,b) + w * afu * aux 

!                 u2(X), v3(X)
                  aux = (mN(1)*Ox(1,1) + mN(2)*Ox(2,1) + 
     2               mN(3)*Ox(3,1)) * s_Nwx(3,b)*sN(2) * s_Nw(a)*sN(3) 
                  aux = gamma * gap * aux / nxny
                  ss_lK(8,a,b) = ss_lK(8,a,b) + w * afu * aux 

!                 u3(X), v3(X)
                  aux = (mN(1)*Ox(1,1) + mN(2)*Ox(2,1) + 
     2               mN(3)*Ox(3,1)) * s_Nwx(3,b)*sN(3) * s_Nw(a)*sN(3) 
                  aux = gamma * gap * aux / nxny
                  ss_lK(9,a,b) = ss_lK(9,a,b) + w * afu * aux 

               END DO
            END DO 

! ------------------------ T2
! coef / h^2 alfa \int ( [ -g]_+ / | a_3(x) |  )  ( v(X) - V(Y) ) . Tnx. ( Orth2.u_r- Orth1.u_s  )  

            DO a = 1, eNoN ! row
               DO b = 1, eNoN ! col

!---              v1 equation
!                 u1(X), v1(X)
                  aux = (Ox(1,1)*s_Nwx(1,b)*sN(1) + 
     2                   Ox(1,2)*s_Nwx(2,b)*sN(1) + 
     3                   Ox(1,3)*s_Nwx(3,b)*sN(1)) * s_Nw(a) 
                  aux = gamma * (gap0 - gap) * aux  
                  ss_lK(1,a,b) = ss_lK(1,a,b) - w * afu * aux            

!                 u2(X), v1(X)
                  aux = (Ox(1,1)*s_Nwx(1,b)*sN(2) + 
     2                   Ox(1,2)*s_Nwx(2,b)*sN(2) + 
     3                   Ox(1,3)*s_Nwx(3,b)*sN(2)) * s_Nw(a) 
                  aux = gamma * (gap0 - gap) * aux  
                  ss_lK(2,a,b) = ss_lK(2,a,b) - w * afu * aux 

!                 u3(X), v1(X)
                  aux = (Ox(1,1)*s_Nwx(1,b)*sN(3) + 
     2                   Ox(1,2)*s_Nwx(2,b)*sN(3) + 
     3                   Ox(1,3)*s_Nwx(3,b)*sN(3)) * s_Nw(a) 
                  aux = gamma * (gap0 - gap) * aux  
                  ss_lK(3,a,b) = ss_lK(3,a,b) - w * afu * aux 

!---              v2 equation
!                 u1(X), v2(X)
                  aux = (Ox(2,1)*s_Nwx(1,b)*sN(1) + 
     2                   Ox(2,2)*s_Nwx(2,b)*sN(1) + 
     3                   Ox(2,3)*s_Nwx(3,b)*sN(1)) * s_Nw(a) 
                  aux = gamma * (gap0 - gap) * aux  
                  ss_lK(4,a,b) = ss_lK(4,a,b) - w * afu * aux 

!                 u2(X), v2(X)
                  aux = (Ox(2,1)*s_Nwx(1,b)*sN(2) + 
     2                   Ox(2,2)*s_Nwx(2,b)*sN(2) + 
     3                   Ox(2,3)*s_Nwx(3,b)*sN(2)) * s_Nw(a) 
                  aux = gamma * (gap0 - gap) * aux  
                  ss_lK(5,a,b) = ss_lK(5,a,b) - w * afu * aux 

!                 u3(X), v2(X)
                  aux = (Ox(2,1)*s_Nwx(1,b)*sN(3) + 
     2                   Ox(2,2)*s_Nwx(2,b)*sN(3) + 
     3                   Ox(2,3)*s_Nwx(3,b)*sN(3)) * s_Nw(a) 
                  aux = gamma * (gap0 - gap) * aux  
                  ss_lK(6,a,b) = ss_lK(6,a,b) - w * afu * aux 

!---              v3 equation
!                 u1(X), v3(X)
                  aux = (Ox(3,1)*s_Nwx(1,b)*sN(1) + 
     2                   Ox(3,2)*s_Nwx(2,b)*sN(1) + 
     3                   Ox(3,3)*s_Nwx(3,b)*sN(1)) * s_Nw(a) 
                  aux = gamma * (gap0 - gap) * aux  
                  ss_lK(7,a,b) = ss_lK(7,a,b) - w * afu * aux 

!                 u2(X), v3(X)
                  aux = (Ox(3,1)*s_Nwx(1,b)*sN(2) + 
     2                   Ox(3,2)*s_Nwx(2,b)*sN(2) + 
     3                   Ox(3,3)*s_Nwx(3,b)*sN(2)) * s_Nw(a) 
                  aux = gamma * (gap0 - gap) * aux 
                  ss_lK(8,a,b) = ss_lK(8,a,b) - w * afu * aux 

!                 u3(X), v3(X)
                  aux = (Ox(3,1)*s_Nwx(1,b)*sN(3) + 
     2                   Ox(3,2)*s_Nwx(2,b)*sN(3) + 
     3                   Ox(3,3)*s_Nwx(3,b)*sN(3)) * s_Nw(a) 
                  aux = gamma * (gap0 - gap) * aux 
                  ss_lK(9,a,b) = ss_lK(9,a,b) - w * afu * aux 

               END DO
            END DO 




!           The assembly as to be done for each quadrature point separately 
            CALL DOASSEM_MTR(eNoN, ptrS, ptrS, ss_lK)  
            CALL DOASSEM_MTR(eNoN, ptrS, ptrM, sm_lK) 

         END DO

         CALL DOASSEM_RES(eNoN, ptrS, s_lR)

C          ss_lK = 0._RKIND
C          CALL DOASSEM(eNoN, ptrS, ss_lK, s_lR)


      END DO

      DEALLOCATE( xbl, dbl, ptrM, xblM, dblM, xlM, ptrS, xl, yl, s_lR,
     2   ss_lK, sm_lK, s_Nw, s_Nwx, s_Nwxi, m_Nw, m_Nwx, m_Nwxi, ylM, 
     3   dlM, dl )

      RETURN
      END SUBROUTINE NL_CONTACTFORCES


!####################################################################
!     This subroutine assembels the element stiffness matrix into the
!     global stiffness matrix (Val sparse matrix formatted as a vector)
      SUBROUTINE DOASSEM_RES (d, eqN, lR)
      USE TYPEMOD
      USE COMMOD, ONLY: dof, rowPtr, colPtr, R, Val
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: d, eqN(d)
      REAL(KIND=RKIND), INTENT(IN) :: lR(dof,d)

      INTEGER(KIND=IKIND) a, rowN

      DO a=1, d
         rowN = eqN(a)
         IF (rowN .EQ. 0) CYCLE
         R(:,rowN) = R(:,rowN) + lR(:,a)
      END DO

      RETURN
      END SUBROUTINE DOASSEM_RES
!####################################################################

!####################################################################
!     This subroutine assembels the element stiffness matrix into the
!     global stiffness matrix (Val sparse matrix formatted as a vector)
      SUBROUTINE DOASSEM_MTR (d, eqRow, eqCol, lK)
      USE TYPEMOD
      USE COMMOD, ONLY: dof, rowPtr, colPtr, Val
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: d, eqRow(d), eqCol(d)
      REAL(KIND=RKIND), INTENT(IN) :: lK(dof*dof,d,d)

      INTEGER(KIND=IKIND) a, b, ptr, rowN, colN, left, right

      DO a=1, d
         rowN = eqRow(a)
         IF (rowN .EQ. 0) CYCLE
         DO b=1, d
            colN = eqCol(b)
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
      END SUBROUTINE DOASSEM_MTR


!#######################################################################
!#######################################################################