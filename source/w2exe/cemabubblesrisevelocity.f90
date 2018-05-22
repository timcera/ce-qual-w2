!*==cemabubblesrisevelocity.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine CEMABUBBLESRISEVELOCITY(Radius, Rhog, Risevelocity)
 
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(8) :: Radius, Rhog, Risevelocity
     intent (in) Rhog
     intent (out) Risevelocity
     intent (inout) Radius
!
! Local variables
!
     real(8) :: dynvisc, eo, h, j, m, nd, radius1, reynolds, rhow,             &
              & risevelocity1, risevelocity2, sigma, w
!
!*** End of declarations rewritten by SPAG
!
 
 
    !Calculate Rise Velocity
     rhow = 1000            !kg/m³
     dynvisc = 0.001002     !kg/m/s
     sigma = 0.0725         !N/m
 
     if(Radius*1000<=1)then  !<= 1 mm
         nd = 4.*rhow*(rhow - Rhog)*9.8*Radius**3/(3.*dynvisc**2)
         w = DLOG10(nd)
         if(nd<=73.)reynolds = nd/24. - 1.7569D-4*nd**2 + 6.9252D-7*nd**3 -    &
                             & 2.3027D-10*nd**4
         if(nd>73.0 .AND. nd<=580.0)                                           &
          & reynolds = 10**( - 1.7095 + 1.33438*w - 0.11591*w**2)
         if(nd>580.)reynolds = 10**( - 1.81391 + 1.34671*w - 0.12427*w**2 +    &
                             & 0.006344*w**3)
         Risevelocity = reynolds*dynvisc/(rhow*Radius)
     endif
 
     if(Radius*1000.0<=15.0 .AND. Radius*1000.0>1.0)then    !<= 15 mm
         m = 9.8*dynvisc**4*(rhow - Rhog)/(rhow**2*sigma**3)
         eo = 9.8*(rhow - Rhog)*Radius**2/sigma
         h = (4./3.)*eo*m**( - 0.149)*(dynvisc/dynvisc)**( - 0.14)
         if(h<59.3)then
             j = 0.94*h**0.757
         else
             j = 3.42*h**0.441
         endif
         Risevelocity = dynvisc/(rhow*Radius)*m**( - 0.149)*(j - 0.857)
     endif
 
     if(Radius*1000.0>15.0 .AND. Radius*1000.0<=18.0)then    !<= 15 mm to 18 mm
 
         radius1 = Radius
         Radius = 0.015
         m = 9.8*dynvisc**4*(rhow - Rhog)/(rhow**2*sigma**3)
         eo = 9.8*(rhow - Rhog)*Radius**2/sigma
         h = (4./3.)*eo*m**( - 0.149)*(dynvisc/dynvisc)**( - 0.14)
         if(h<59.3)then
             j = 0.94*h**0.757
         else
             j = 3.42*h**0.441
         endif
         risevelocity1 = dynvisc/(rhow*Radius)*m**( - 0.149)*(j - 0.857)
 
         Radius = 0.018
         risevelocity2 = 0.711*SQRT(9.8*Radius*(rhow - Rhog)/rhow)
 
         Radius = radius1
 
         Risevelocity = ((Radius - 0.015)*risevelocity2 +                      &
                      & risevelocity1*(0.018 - Radius))/(0.018 - 0.015)
 
     endif
 
     if(Radius*1000.0>18.0)Risevelocity = 0.711*SQRT(9.8*Radius*(rhow - Rhog)  &
      & /rhow)                   !> 18 mm
 
 
     end subroutine CEMABUBBLESRISEVELOCITY
