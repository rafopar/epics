      subroutine read_status(ich, status)
      implicit none
      character*100 pvname
      integer*4 ich, iret_code, lenpv
      integer*2 status
      call get_pvname(ich, pvname, lenpv)
      call g_value(pvname(1:lenpv)//'.STI', 'short', status, 1,
     &  iret_code)
      if (iret_code .ne. 0) print *, ' return code from get', iret_code
      return
      end
