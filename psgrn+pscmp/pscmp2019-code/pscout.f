      subroutine pscout(ierr)
      use pscalloc
      implicit none
      integer*4 ierr
c
c     Last modified: Potsdam, Feb, 2019, by R. Wang
c
      integer*4 i,iadd,nadd,it,ieq,ieq1,ieq2,isnap
      integer*4 itr,itr1,itr2,j,inx,irec,lend,m
      real*8 dt,t,t1,t2,wl,wr,swap
      real*8 s0xx,s0yy,s0zz,s0xy,s0yz,s0zx,p0
      real*8 sxx,syy,szz,sxy,syz,szx,p
      real*8 strikeopt1,dipopt1,rakeopt1,strikeopt2,dipopt2,rakeopt2
      real*8 cfs1,cfs2,cfs0,cfsmaxdif,cfsopt,cfsfix
      real*8 sig1,sig2,sigfix,sigopt,sigopt1,sigopt2,rakeopt
      real*8 obsout(14),obsadd(15)
      character*13 cmptxt(14),cmptxtadd(15)
      character*6 shortcmptxt(14)
c
c     DATA OUTPUT
c     ===========
c
      allocate(rtxt(nrec),stat=ierr)
      if(ierr.ne.0)stop ' Error in pscout: rtxt not allocated!'
c
      do lend=80,1,-1
        if(outdir(lend:lend).ne.' ')goto 100
      enddo
100   continue
c
      if(lend.lt.1)then
        stop 'Error in edcmain: wrong for output dirtory!'
      endif
c
      shortcmptxt(1)='   Un_'
      shortcmptxt(2)='   Ue_'
      shortcmptxt(3)='   Ud_'
      shortcmptxt(4)='  Snn_'
      shortcmptxt(5)='  See_'
      shortcmptxt(6)='  Sdd_'
      shortcmptxt(7)='  Sne_'
      shortcmptxt(8)='  Sed_'
      shortcmptxt(9)='  Sdn_'
      shortcmptxt(10)='   Tn_'
      shortcmptxt(11)='   Te_'
      shortcmptxt(12)='  Rot_'
      shortcmptxt(13)='   Gd_'
      shortcmptxt(14)='   Gr_'
c
      cmptxt(1)='   Disp_north'
      cmptxt(2)='    Disp_east'
      cmptxt(3)='    Disp_down'
      cmptxt(4)='    Stress_nn'
      cmptxt(5)='    Stress_ee'
      cmptxt(6)='    Stress_dd'
      cmptxt(7)='    Stress_ne'
      cmptxt(8)='    Stress_ed'
      cmptxt(9)='    Stress_dn'
      cmptxt(10)='       Tilt_n'
      cmptxt(11)='       Tilt_e'
      cmptxt(12)='     Rotation'
      cmptxt(13)='        Geoid'
      cmptxt(14)='      Gravity'
      iadd=0
      if(insar.eq.1)then
        iadd=iadd+1
        cmptxtadd(iadd)='     Disp_LOS'
      endif
      if(icfs.eq.1)then
        iadd=iadd+1
        cmptxtadd(iadd)='      CFS_Max'
        iadd=iadd+1
        cmptxtadd(iadd)='      CFS_Mas'
        iadd=iadd+1
        cmptxtadd(iadd)='  CFS_Mas_Opt'
        iadd=iadd+1
        cmptxtadd(iadd)='    Sigma_Mas'
        iadd=iadd+1
        cmptxtadd(iadd)=' Rake_Mas_Opt'
        iadd=iadd+1
        cmptxtadd(iadd)='      CFS_Opt'
        iadd=iadd+1
        cmptxtadd(iadd)='  Sigma_Opt_1'
        iadd=iadd+1
        cmptxtadd(iadd)=' Strike_Opt_1'
        iadd=iadd+1
        cmptxtadd(iadd)='    Dip_Opt_1'
        iadd=iadd+1
        cmptxtadd(iadd)='   Rake_Opt_1'
        iadd=iadd+1
        cmptxtadd(iadd)='  Sigma_Opt_2'
        iadd=iadd+1
        cmptxtadd(iadd)=' Strike_Opt_2'
        iadd=iadd+1
        cmptxtadd(iadd)='    Dip_Opt_2'
        iadd=iadd+1
        cmptxtadd(iadd)='   Rake_Opt_2'
c
c       p0 = -skempton*(sigma0(1)+sigma0(2)+sigma0(3))/3.d0
c       for undrained conditions
c       p0 = 0 for drained conditions (during preseismic period)
c
        p0=0.d0
        call prestress(sigma0(1),sigma0(2),sigma0(3),strike0,dip0,rake0,
     &                 p0,friction,cfs0,s0xx,s0yy,s0zz,s0xy,s0yz,s0zx)
      endif
      nadd=iadd
      inx=int(alog10(0.1+float(nrec)))+1
      do irec=1,nrec
        m=irec
        do j=1,inx
          i=m/10**(inx-j)
          rtxt(irec)(j:j)=char(ichar('0')+i)
          m=m-i*10**(inx-j)
        enddo
      enddo
c
      if(nt.gt.1)then
        dt=twindow/dble(nt-1)
      else
        dt=twindow
      endif
c
      do i=1,14
        if(itout(i).eq.1)then
          outputfile=outdir(1:lend)//toutfile(i)
          open(30,file=outputfile,status='unknown')
          write(30,'(a13,$)')'    Time_day'
          do irec=1,nrec-1
            write(30,'(a13,$)')shortcmptxt(i)//rtxt(irec)(1:inx)
          enddo
          write(30,'(a13)')shortcmptxt(i)//rtxt(nrec)(1:inx)
          do it=1,nt
            t=dble(it-1)*dt
            write(30,'(E13.4,$)')(t+torigin)/DAY2SEC
            do irec=1,nrec-1
              obsout(i)=obs(it,irec,i)
              do ieq=1,neq
                if(t.ge.eqtime(ieq))then
                  obsout(i)=obsout(i)+coobs(ieq,irec,i)
                endif
              enddo
              write(30,'(E13.4,$)')obsout(i)
            enddo
            obsout(i)=obs(it,nrec,i)
            do ieq=1,neq
              if(t.ge.eqtime(ieq))then
                obsout(i)=obsout(i)+coobs(ieq,nrec,i)
              endif
            enddo
            write(30,'(E13.4)')obsout(i)
          enddo
          close(30)
        endif
      enddo
c
c     scenario outputs
c
      do isnap=1,nsnap
        if(tsnap(isnap).gt.twindow)then
          nwarn=nwarn+1
          print *,' Warning in pecout: scenario outside time window,'
          print *,' no output for the ',isnap,'. scenario!'
          goto 500
        endif
        outputfile=outdir(1:lend)//scoutfile(isnap)
        open(30,file=outputfile,status='unknown')
        write(30,'(a26,$)')'     Lat[deg]     Lon[deg]'
        do i=1,14
          write(30,'(a13,$)')cmptxt(i)
        enddo
        write(30,'(14a13)')(cmptxtadd(i),i=1,nadd)
        if(onlysc)then
          it=min0(1+idint(tsnap(isnap)/dt),nt)
          itr1=0
          do itr=1,ntr
            if(it.eq.itsnap(itr))itr1=itr
          enddo
          if(itr1.eq.0)then
            stop 'Error in pscout: snapshot time(-) not found!'
          endif
c
          it=min0(it+1,nt)
          itr2=0
          do itr=1,ntr
            if(it.eq.itsnap(itr))itr2=itr
          enddo
          if(itr2.eq.0)then
            stop 'Error in pscout: snapshot time(+) not found!'
          endif
          t1=dble(itsnap(itr1)-1)*dt
          t2=dble(itsnap(itr2)-1)*dt
        else
          itr1=min0(1+idint(tsnap(isnap)/dt),nt)
          itr2=min0(itr1+1,nt)
          t1=dble(itr1-1)*dt
          t2=dble(itr2-1)*dt
        endif
        ieq1=0
        do ieq=1,neq
          if(t1.lt.eqtime(ieq).and.eqtime(ieq).le.tsnap(isnap))then
            ieq1=ieq
            t1=eqtime(ieq)
          endif
        enddo
        ieq2=0
        do ieq=neq,1,-1
          if(t2.gt.eqtime(ieq).and.eqtime(ieq).ge.tsnap(isnap))then
            ieq2=ieq
            t2=eqtime(ieq)
          endif
        enddo
        if(ieq1.eq.0)then
          do irec=1,nrec
            do i=1,14
              obs1(irec,i)=obs(itr1,irec,i)
            enddo
          enddo
        else
          do irec=1,nrec
            do i=1,14
              obs1(irec,i)=poobs(ieq1,irec,i)
            enddo
          enddo
        endif
        if(ieq2.eq.0)then
          do irec=1,nrec
            do i=1,14
              obs2(irec,i)=obs(itr2,irec,i)
            enddo
          enddo
        else
          do irec=1,nrec
            do i=1,14
              obs2(irec,i)=poobs(ieq2,irec,i)
            enddo
          enddo
        endif
        if(t2-t1.gt.0.1d-03*dt)then
          wr=(tsnap(isnap)-t1)/(t2-t1)
        else
          wr=0.d0
        endif
        wl=1.d0-wr
        do irec=1,nrec
          iadd=0
          if(insar.eq.1)then
            iadd=iadd+1
            obsadd(iadd)=xlos*(wl*obs1(irec,1)+wr*obs2(irec,1))
     &                  +ylos*(wl*obs1(irec,2)+wr*obs2(irec,2))
     &                  +zlos*(wl*obs1(irec,3)+wr*obs2(irec,3))
c
            do ieq=1,neq
              if(tsnap(isnap).ge.eqtime(ieq))then
                obsadd(iadd)=obsadd(iadd)+xlos*coobs(ieq,irec,1)
     &                 +ylos*coobs(ieq,irec,2)+zlos*coobs(ieq,irec,3)
              endif
            enddo
            obsadd(iadd)=obsadd(iadd)
          endif
          if(icfs.eq.1)then
            sxx=wl*obs1(irec,4)+wr*obs2(irec,4)
            syy=wl*obs1(irec,5)+wr*obs2(irec,5)
            szz=wl*obs1(irec,6)+wr*obs2(irec,6)
            sxy=wl*obs1(irec,7)+wr*obs2(irec,7)
            syz=wl*obs1(irec,8)+wr*obs2(irec,8)
            szx=wl*obs1(irec,9)+wr*obs2(irec,9)
            do ieq=1,neq
              if(tsnap(isnap).ge.eqtime(ieq))then
                sxx=sxx+coobs(ieq,irec,4)
                syy=syy+coobs(ieq,irec,5)
                szz=szz+coobs(ieq,irec,6)
                sxy=sxy+coobs(ieq,irec,7)
                syz=syz+coobs(ieq,irec,8)
                szx=szx+coobs(ieq,irec,9)
              endif
            enddo
c
c           p = excess pore pressure under undrained conditions
c
            p=-skempton*(sxx+syy+szz)/3.d0
c
            call cfs3dopt(s0xx+sxx,s0yy+syy,s0zz+szz,s0xy+sxy,
     &                    s0yz+syz,s0zx+szx,p,friction,1,
     &                    strike0,dip0,rake0,
     &                    cfsopt,sigopt,strikeopt1,dipopt1,rakeopt1,
     &                    strikeopt2,dipopt2,rakeopt2)
c
c           add CFS_Max_Dif
c
            iadd=iadd+1
            obsadd(iadd)=cfsopt-cfs0
c
            call cfsmas(sxx,syy,szz,sxy,syz,szx,
     &                  p,friction,cfsfix,sigfix,strike0,dip0,rake0)
c
c           add CFS_Mas_fix
c
            iadd=iadd+1
            obsadd(iadd)=cfsfix
c
            call cfsmasopt(s0xx+sxx,s0yy+syy,s0zz+szz,s0xy+sxy,
     &                     s0yz+syz,s0zx+szx,
     &                     p,friction,cfs1,sig1,
     &                     strike0,dip0,rakeopt)
c
            call cfsmas(s0xx,s0yy,s0zz,s0xy,s0yz,s0zx,
     &                  p0,friction,cfs2,sig2,
     &                  strike0,dip0,rakeopt)
c
c           CFS_Mas_Opt
c
            iadd=iadd+1
            obsadd(iadd)=cfs1-cfs2
c
c           Sig_Mas
c
            iadd=iadd+1
            obsadd(iadd)=sigfix
c
c           Rake_Mas_Opt
c
            iadd=iadd+1
            obsadd(iadd)=rakeopt
c
            call cfsmas(s0xx,s0yy,s0zz,s0xy,s0yz,s0zx,
     &                  p0,friction,cfs1,sig1,
     &                  strikeopt1,dipopt1,rakeopt1)
c
c           Rake_3D_Opt
c
            iadd=iadd+1
            obsadd(iadd)=cfsopt-cfs1
c
c           (Sigma, Strke, Dip, Rake)_Opt_1
c
            iadd=iadd+1
            obsadd(iadd)=sigopt-sig1
            iadd=iadd+1
            obsadd(iadd)=strikeopt1
            iadd=iadd+1
            obsadd(iadd)=dipopt1
            iadd=iadd+1
            obsadd(iadd)=rakeopt1
c
            call cfsmas(s0xx,s0yy,s0zz,s0xy,s0yz,s0zx,
     &                  p0,friction,cfs2,sig2,
     &                  strikeopt2,dipopt2,rakeopt2)
c
c           (Sigma, Strke, Dip, Rake)_Opt_2
c
            iadd=iadd+1
            obsadd(iadd)=sigopt-sig2
            iadd=iadd+1
            obsadd(iadd)=strikeopt2
            iadd=iadd+1
            obsadd(iadd)=dipopt2
            iadd=iadd+1
            obsadd(iadd)=rakeopt2
          endif
          write(30,'(2f13.6,$)')latrec(irec),lonrec(irec)
          do i=1,14
            obsout(i)=wl*obs1(irec,i)+wr*obs2(irec,i)
            do ieq=1,neq
              if(tsnap(isnap).ge.eqtime(ieq))then
                obsout(i)=obsout(i)+coobs(ieq,irec,i)
              endif
            enddo
          enddo
          write(30,'(29E13.5)')(obsout(i),i=1,14),(obsadd(i),i=1,nadd)
        enddo
        close(30)
500     continue
      enddo
      if(icfs.eq.1)then
        write(*,*)'------------------------------------------------'
     &          //'----------------'
        print *,'     Pre-stress Tensor (e/n/d = east/north/down) [MPa]'
        write(*,'(a)')'       Snn       See       Sdd'
     &              //'       Sne       Sed       Sdn'
        write(*,'(6f10.4)')s0xx*1.0d-06,s0yy*1.0d-06,s0zz*1.0d-06,
     &                     s0xy*1.0d-06,s0yz*1.0d-06,s0zx*1.0d-06
        write(*,'(a)')'      Principal pre-stresses [MPa]'
        write(*,'(3f10.4)')(sigma0(i)*1.0d-06,i=1,3)
        write(*,'(a)')'      Maximum pre- Coulomb stress [MPa]'
        write(*,'(f10.4)')cfs0*1.0d-06
      endif
c
      write(*,*)'------------------------------------------------'
     &        //'----------------'
      write(*,'(a,E12.4,$)')' Seismic moment = ',mw
      if(mw.gt.0.d0)then
        mw=(dlog10(mw)-9.1d0)/1.5d0
      endif
      write(*,'(a,f4.2,a)')' Nm (Mw = ',mw,') by tensor summation '
      write(*,'(a,E12.4,$)')' or ',mwsca
      if(mwsca.gt.0.d0)then
        mwsca=(dlog10(mwsca)-9.1d0)/1.5d0
      endif
      write(*,'(a,f4.2,a)')' Nm (Mw = ',mwsca,
     &                     ') by scalar summation used for'
      write(*,'(a)')' calculating the deformation field.'
      write(*,*)'------------------------------------------------'
     &        //'----------------'
c
      return
      end
