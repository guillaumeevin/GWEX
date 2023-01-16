! compile with the command "R CMD SHLIB GWex_disag.f90" in a command prompt

!=========================== indexx ==========================
! Use the Heapsort algorithm to index an array arrin of length n.
! Output the array indx such that arrin(indx(j)) is in ascending
! order for j = 1,2,...,n.  The input quantities n and arrin are
! not changed.
!
! This is a Numerical Recipes routine, but modified by one
! line to work if n equals 1.
subroutine indexx(arr,n,indx)

  parameter (M=7,NSTACK=50)
  integer n,indx(n)
  double precision arr(n)
  integer i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
  double precision a

  do j=1,n
     indx(j)=j
  enddo

  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.M) then
     do j=l+1,ir
        indxt=indx(j)
        a=arr(indxt)
        do i=j-1,1,-1
           if(arr(indx(i)).le.a) goto 2
           indx(i+1)=indx(i)
        enddo
        i=0
2       indx(i+1)=indxt
     enddo
     if(jstack.eq.0) return

     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2

  else
     k=(l+ir)/2
     itemp=indx(k)
     indx(k)=indx(l+1)
     indx(l+1)=itemp

     if(arr(indx(l+1)).gt.arr(indx(ir))) then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
     endif

     if(arr(indx(l)).gt.arr(indx(ir))) then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
     endif

     if(arr(indx(l+1)).gt.arr(indx(l))) then
        itemp=indx(l+1)
        indx(l+1)=indx(l)
        indx(l)=itemp
     endif

     i=l+1
     j=ir
     indxt=indx(l)
     a=arr(indxt)
3    continue
     i=i+1
     if(arr(indx(i)).lt.a) goto 3

4    continue
     j=j-1
     if(arr(indx(j)).gt.a) goto 4
     if(j.lt.i) goto 5
     itemp=indx(i)
     indx(i)=indx(j)
     indx(j)=itemp
     goto 3

5    indx(l)=indx(j)
     indx(j)=indxt
     jstack=jstack+2
     if(jstack.gt.NSTACK) return
     if(ir-i+1.ge.j-l)then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     endif
  endif
  goto 1

end subroutine indexx


!=========================== disag3DayGWexPrec_F ==========================
! For each 3-day simulations at n gauges, we find the closest field in observations for which
! 3-day precipitation structures are available.
subroutine disag3DayGWexPrec_F(Yobs, Y3obs, mObs, cObs, Y3sim, mSim, cSim, nTobs,&
    &nStat, nTsim, nLagScore, Ysim, codeDisag)
    implicit none
    !  Input/Output
    integer, intent(in) :: nTobs, nStat, nTsim, nLagScore
    integer, intent(in) :: mObs(nTobs) ! index of the month (1..12) for the matrix of 3-day observations
    integer, intent(in) :: mSim(nTsim) ! index of the month (1..12) for the matrix of 3-day simulations
    integer, intent(in) :: cObs(nTobs) ! class (=1,2,3,4) of precip for observations
    integer, intent(in) :: cSim(nTsim) ! class (=1,2,3,4) of precip for simulations
    double precision, intent(in), dimension(nTobs*3,nStat) :: Yobs    ! Matrix of observations: nTobs*3 [days] x nStat [stations]
    double precision, intent(in), dimension(nTobs,nStat) :: Y3obs  ! Matrix of observations amounts for 3-day periods: nTobs [3-days] x nStat [stations]
    double precision, intent(in), dimension(nTsim,nStat) :: Y3sim   ! Matrix of simulated amounts for 3-day periods: nTsim [3-days] x nStat [stations]
    double precision, intent(out), dimension(nTsim*3,nStat) :: Ysim ! Matrix of disag. simulated amounts: nTsim*3 [days] x nStat [stations]
    double precision, intent(out), dimension(nTsim,nStat) :: codeDisag ! Matrix indicating how it has been disaggregated

    !  Locals
    integer :: nBestField,i,j,k,iDay,jDay,j3Day,iLag
    double precision :: naVal,rmseIJ, r
    PARAMETER(nBestField=10,naVal=-9999.)
    double precision :: adimObs(nStat), adimSim(nStat)
    double precision, dimension(3) :: Yobs3D
    double precision, dimension(nTobs) :: rmseI
    integer, dimension(nTobs) :: indBestRMSEI, indBestrmseDayI
    integer, dimension(nBestField) :: indBestFieldI
    logical :: notSameClass

    ! For each simulated 3-day period
    do i = 1, nTsim
       ! index of the beginning of the 3-day period corresponding to this day, in daily matrices
       iDay = (i-1)*3

       !======= Step 1: Minimize Score for this selection
       ! sum of square differences for the closest fields on two time steps

       ! for the first time step
       if(i <= nLagScore) then
           do j = 1, nTobs
              ! same month and class
              notSameClass = (mSim(i)/=mObs(j)) .OR. (cSim(i)/=cObs(j))
              ! if any observed value is missing in the observed prec. field or if the months do not correspond
              ! we discard this observed day
              if((any(Y3obs(j,:) == naVal)) .OR. notSameClass) then
                rmseI(j) = 1E30
              else
                ! absolute differences between adimensioned precipitation for this day
                if(sum(Y3sim(i,:))==0) then
                    adimSim = 0
                else
                    adimSim = Y3sim(i,:)/sum(Y3sim(i,:))
                end if

                 if(sum(Y3obs(j,:))==0) then
                    adimObs = 0
                else
                    adimObs = Y3obs(j,:)/sum(Y3obs(j,:))
                end if
                rmseI(j) = sum(abs(adimSim-adimObs))
              end if
           enddo
       else
         ! discard first elements
         do j = 1, nLagScore
            rmseI(j) = 1E30
         end do

         ! for the next days, compute score
         loopStationScore: do j = (nLagScore+1), nTobs
              ! same month and class
              notSameClass = (mSim(i)/=mObs(j)) .OR. (cSim(i)/=cObs(j))
              ! if any observed value is missing in the observed prec. field or if the months do not correspond
              ! we discard this observed day
              if(any(Y3obs(j,:) == naVal) .OR. notSameClass) then
                rmseI(j) = 1E30
              else
                jDay = (j-1)*3
                ! absolute differences between adimensioned precipitation for this day
                if(sum(Y3sim(i,:))==0) then
                    adimSim = 0
                else
                    adimSim = Y3sim(i,:)/sum(Y3sim(i,:))
                end if

                 if(sum(Y3obs(j,:))==0) then
                    adimObs = 0
                else
                    adimObs = Y3obs(j,:)/sum(Y3obs(j,:))
                end if
                rmseIJ = sum(abs(adimSim-adimObs))
                ! add differences for the previous days, just non na values
                do iLag = 1, nLagScore
                    if(any(Yobs(jDay+1-iLag,:)==naVal)) then
                       rmseI(j) = 1E30
                       cycle loopStationScore
                    else
                        if(sum(Ysim(iDay+1-iLag,:))==0) then
                            adimSim = 0
                        else
                            adimSim = Ysim(iDay+1-iLag,:)/sum(Ysim(iDay+1-iLag,:))
                        end if

                         if(sum(Yobs(jDay+1-iLag,:))==0) then
                            adimObs = 0
                        else
                            adimObs = Yobs(jDay+1-iLag,:)/sum(Yobs(jDay+1-iLag,:))
                        end if
                        rmseIJ = rmseIJ + sum(abs(adimSim-adimObs))
                    end if
                end do
                rmseI(j) = rmseIJ
              end if
           enddo loopStationScore
       end if

       call indexx(rmseI,nTobs,indBestRMSEI)
       indBestFieldI = indBestRMSEI(1:nBestField)

       !======= Step 3: Look at the different case and disaggregate =====
       loopStationDisag: do k = 1, nStat
         ! initialise code
         codeDisag(i,k) = naVal

         !!!!! case 1: no occurrence for this period and this station, nothing to disaggregate.
         if(Y3sim(i,k)==0.) then
            codeDisag(i,k) = 0.
            Ysim((iDay+1):(iDay+3),k) = 0.

         !!!!! case 2: loop at the closest fields, if there is the same number of occurrences, we take the observed
         ! temporal structure (less restrictive than the exact same structure)
         else
            do j = 1, nBestField
              jDay = (indBestFieldI(j)-1)*3

              ! check if there is an observed precitation for the selected observed field and for this station
              if(Y3obs(indBestFieldI(j),k)>0) then
                ! update code of disaggregation
                codeDisag(i,k) = dble(j)
                ! simulated precipitation for these 3 days are observed precipitation for the close field, rescaled
                ! by the simulated cumulated value for these 3 days
                Yobs3D = Yobs((jDay+1):(jDay+3),k)
                Ysim((iDay+1):(iDay+3),k) = Yobs3D*Y3sim(i,k)/sum(Yobs3D)
                ! codes to analyse how large 3-day intensities are disaggregated
                if(any(Ysim((iDay+1):(iDay+3),k)>maxval(Yobs(:,k)))) then
                    codeDisag(i,k) = codeDisag(i,k) + 10000.
                end if
                ! get out of the loop for this station
                cycle loopStationDisag
              end if
            end do
         end if

         !!!!! case 3: if we did not find similar structure then we find, for this station, days with similar amounts
         if(codeDisag(i,k)==naVal) then
               ! for the first time step
               if(i <= nLagScore) then
                   do j = 1, nTobs
                      ! same month and class
                      notSameClass = (mSim(i)/=mObs(j)) .OR. (cSim(i)/=cObs(j))
                      ! if any observed value is missing in the observed prec. field or if the months do not correspond
                      ! we discard this observed day
                      if((Y3obs(j,k) == naVal) .OR. notSameClass) then
                        rmseI(j) = 1E30
                      else
                        rmseI(j) = abs(Y3sim(i,k)-Y3obs(j,k))
                      end if
                   enddo
               else
                 ! discard first elements
                 do j = 1, nLagScore
                    rmseI(j) = 1E30
                 end do

                 ! for the next days, compute score
                 do j = (nLagScore+1), nTobs
                      ! same month and class
                      notSameClass = (mSim(i)/=mObs(j)) .OR. (cSim(i)/=cObs(j))
                      ! if any observed value is missing in the observed prec. field or if the months do not correspond
                      ! or if there is no observed precip for this day
                      ! we discard this observed day
                      if(((Y3obs(j,k) == naVal) .OR. (Y3obs(j,k) == 0)) .OR. notSameClass) then
                        rmseI(j) = 1E30
                      else
                        rmseIJ = abs(Y3sim(i,k)-Y3obs(j,k))
                        ! add differences for the previous days, just non na values
                        jDay = (j-1)*3
                        if(any(Yobs((jDay+1-nLagScore):jDay,k)==naVal)) then
                            rmseIJ = 1E30
                        else
                            do iLag = 1, nLagScore
                               rmseIJ = rmseIJ + abs(Ysim(iDay+1-iLag,k)-Yobs(jDay+1-iLag,k))
                            end do
                        end if
                        rmseI(j) = rmseIJ
                      end if
                   enddo
               end if

            ! pick one of 10 days with similar occ structures and intensities at this station
            call indexx(rmseI,nTobs,indBestrmseDayI)
            codeDisag(i,k) = 2000
            call RANDOM_NUMBER(r)
            j3Day = indBestrmseDayI(int(r*10)+1)
            jDay = (j3Day-1)*3
            Yobs3D = Yobs((jDay+1):(jDay+3),k)
            Ysim((iDay+1):(iDay+3),k) = Yobs3D*Y3sim(i,k)/sum(Yobs3D)
            ! codes to analyse how large 3-day intensities are disaggregated
            if(any(Ysim((iDay+1):(iDay+3),k)>maxval(Yobs(:,k)))) then
               codeDisag(i,k) = codeDisag(i,k) + 10000.
            end if
        end if
    enddo loopStationDisag
  enddo
end subroutine disag3DayGWexPrec_F
