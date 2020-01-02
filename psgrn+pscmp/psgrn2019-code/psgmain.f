      program psgmain
      use psgalloc
      implicit none
c
c     work space
c
      integer*4 ierr
      integer*4 runtime
      integer*4 time
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#                  Welcome to                        #'
      print *,'#                                                    #'
      print *,'#     PPPP     SSSS    GGGG   RRRR    N   N          #'
      print *,'#     P   P   S       G       R   R   NN  N          #'
      print *,'#     PPPP     SSS    G GGG   RRRR    N N N          #'
      print *,'#     P           S   G   G   R R     N  NN          #'
      print *,'#     P       SSSS     GGGG   R  R    N   N          #'
      print *,'#                                                    #'
      print *,'#                  Version 2019                      #'
      print *,'# (update of version 2008b -> use of dynamic memory) #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#             Helmholtz Centre Potsdam               #'
      print *,'#    GFZ German Research Centre for Geosciences      #'
      print *,'#             Last modified: Feb. 2019               #'
      print *,'######################################################'
      print *,'                                                      '
c
      nwarn=0
      write(*,'(a,$)')' Please type the file name of input data: '
      read(*,'(a)')inputfile
      runtime=time()
c
c     read input file
c
      call psggetinp(ierr)
c
c     construction of sublayers
c
      call psgsublay(ierr)
c
c     main computation procedure
c
      call psgprocess(ierr)
c
      runtime=time()-runtime
      if(nwarn.eq.0)then
        write(*,'(a)')'################################################'
        write(*,'(a)')'#                                              #'
        write(*,'(a)')'#        End of computations with PSGRN        #'
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
