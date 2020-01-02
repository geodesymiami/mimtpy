      program pscmain
      use pscalloc
      implicit none
c
c     work space
c
      integer*4 ierr
      integer*4 runtime
      integer*4 time
c
      nwarn=0
c
      print *,'########################################################'
      print *,'#                                                      #'
      print *,'#                  Welcome to                          #'
      print *,'#                                                      #'
      print *,'#       PPPP     SSSS    CCCC   M   M    PPPP          #'
      print *,'#       P   P   S       C       MM MM    P   P         #'
      print *,'#       PPPP     SSS    C       M M M    PPPP          #'
      print *,'#       P           S   C       M   M    P             #'
      print *,'#       P       SSSS     CCCC   M   M    P             #'
      print *,'#                                                      #'
      print *,'#                    Version 2019                      #'
      print *,'#             (update of version 2008b)                #'
      print *,'# (new input format, add. output, dynamic memory, etc) #'
      print *,'#                                                      #'
      print *,'#                        by                            #'
      print *,'#                   Rongjiang Wang                     #'
      print *,'#                (wang@gfz-potsdam.de)                 #'
      print *,'#                                                      #'
      print *,'#               Helmholtz Centre Potsdam               #'
      print *,'#      GFZ German Research Centre for Geosciences      #'
      print *,'#              Last modified: March 2019               #'
      print *,'########################################################'
      print *,'                                                        '
      write(*,'(a,$)')' Please type the file name of input data: '
      read(*,'(a)')inputfile
      runtime=time()
c
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN INPUT PARAMETERS
c     ========================
c00000000000000000000000000000000000000000000000000000000000000000000000
      print *,'... read input data ...'
      call pscgetinp(ierr)
c00000000000000000000000000000000000000000000000000000000000000000000000
c      MAIN PROCESSING
c      ===============
c00000000000000000000000000000000000000000000000000000000000000000000000
      nwarn=0
      if(iesmodel.eq.0)then
        print *,'... using the analytical Okada solutions ...'
      else
        print *,'... using the numerical Green function approach ...'
      endif
      call pscgrn(ierr)
      print *,'... outputs ...'
      call pscout(ierr)
c00000000000000000000000000000000000000000000000000000000000000000000000
c      END OF STANDARD PROCESSING
c      ==========================
c00000000000000000000000000000000000000000000000000000000000000000000000
      runtime=time()-runtime
      if(nwarn.eq.0)then
        write(*,'(a)')'################################################'
        write(*,'(a)')'#                                              #'
        write(*,'(a)')'#        End of computations with PSCMP        #'
        write(*,'(a)')'#                                              #'
        write(*,'(a,i10,a)')'#        Run time: ',runtime,
     &                                             ' sec              #'
        write(*,'(a)')'################################################'
      else if(nwarn.lt.10)then
        write(*,'(a)')'################################################'
        write(*,'(a,i10,a)')'#        Run time: ',runtime,
     &                                             ' sec              #'
        write(*,'(a,i2,a)')'   Sorry, there have been',nwarn,
     &                     ' warnings.    '
        write(*,'(a)')'           Results may be inaccurate!           '
        write(*,'(a)')'################################################'
      else
        write(*,'(a)')'################################################'
        write(*,'(a,i10,a)')'#        Run time: ',runtime,
     &                                             ' sec              #'
        write(*,'(a)')'  Sorry, there have been more than 10 warnings  '
        write(*,'(a)')'           Results may be inaccurate!           '
        write(*,'(a)')'################################################'
      endif
c
      stop
      end
