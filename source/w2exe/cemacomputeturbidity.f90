!*==cemacomputeturbidity.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine CEMACOMPUTETURBIDITY
 
 
! Type declarations
     use GLOBAL
     use SCREENC
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(8) :: celltssvalue
     integer(2) :: k
!
!*** End of declarations rewritten by SPAG
!
 
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 do k = kt, KB(segnumi)
                    !CellTSSValue = C1(K,SegNumI,6)
                     celltssvalue = C1(k, segnumi, nturb)   ! cb 2/18/13
                    !C1(K,SegNumI,6) = exp(CoeffA_Turb*log(CellTSSValue) + CoeffB_Turb)
                    !C2(K,SegNumI,6) = exp(CoeffA_Turb*log(CellTSSValue) + CoeffB_Turb)
                     C1(k, segnumi, nturb)                                     &
                       & = EXP(coeffa_turb*LOG(celltssvalue) + coeffb_turb)                      ! cb 2/18/13
                     C2(k, segnumi, nturb)                                     &
                       & = EXP(coeffa_turb*LOG(celltssvalue) + coeffb_turb)
                 enddo !K
             enddo
         enddo
     enddo
 
 
     end subroutine CEMACOMPUTETURBIDITY
