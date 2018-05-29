!*==cemadisgasphasedistribution.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
     subroutine CEMADISGASPHASEDISTRIBUTION(Ctotal, Molwt, Gastemp, Henryconst,&
       & Vwtr, Cgasph, Cliqph)
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(8) :: Cgasph, Cliqph, Ctotal, Gastemp, Henryconst, Molwt, Vwtr
     intent (in) Ctotal, Gastemp, Henryconst
     intent (out) Cgasph, Cliqph
!
! Local variables
!
!
!*** End of declarations rewritten by SPAG
!
 
 
     Cgasph = Ctotal/(1 + gasconst_r*Gastemp/Henryconst)
     Cliqph = Ctotal*(gasconst_r*Gastemp/Henryconst)                           &
            & /(1 + gasconst_r*Gastemp/Henryconst)
 
     end subroutine CEMADISGASPHASEDISTRIBUTION
