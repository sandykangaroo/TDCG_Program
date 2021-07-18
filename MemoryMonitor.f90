module ModMemoryMonitor
    implicit none
     
    contains    

subroutine memorymonitor(memoryused)
    use Kernel32
    use ISO_C_Binding 
    Implicit None
    type MY_MEMORYSTATUSEX
    sequence
        integer(DWORD) dwlength
        integer(DWORD) dwMemeoryLoad
        integer(DWORDLONG) ullTotalPhys     
        integer(DWORDLONG) ullAvailPhys       
        integer(DWORDLONG) ullTotalPageFile       
        integer(DWORDLONG) ullAvailPageFile   
        integer(DWORDLONG) ullTotalVirtual  
        integer(DWORDLONG) ullAvailVirtual   
        integer(DWORDLONG) ullAvailExtendedVirtual
    end type MY_MEMORYSTATUSEX
  
    type(T_MEMORYSTATUSEX) stMemStat
    type(MY_MEMORYSTATUSEX) stMyMemStat
    integer i
    integer(DWORDLONG) memoryused
  
    stMemStat%dwLength = sizeof(stMemStat)   ! 必须有这一句，否则函数出错
    i = GlobalMemoryStatusEx ( stMemStat )
    stMyMemStat = transfer(stMemStat, stMyMemStat)
  !write(*,*) "There is  ",stMemStat%dwMemoryLoad,"% percent of memory in use."                         !内存使用百分比
  !write(*,*) "There are ",stMyMemStat%ullTotalPhys/1024/1024," total MB of physical memory."      !物理内存总量
  write(*,*) "There are ",stMyMemStat%ullAvailPhys/1024/1024," free  MB of physical memory."       !物理内存可用量
  !write(*,*) "There are ",stMyMemStat%ullTotalPageFile/1024," total KB of paging file."                     !页面文件总量
  !write(*,*) "There are ",stMyMemStat%ullAvailPageFile/1024," free  KB of paging file."                      !页面文件可用量
  !write(*,*) "There are ",stMyMemStat%ullTotalVirtual/1024," total KB of virtual memory."                 !虚拟内存总量
  !write(*,*) "There are ",stMyMemStat%ullAvailVirtual/1024," free  KB of virtual memory."                 !虚拟内存可用量
  !write(*,*) "There are ",stMyMemStat%ullAvailExtendedVirtual/1024," free  KB of extended memory."
  !memoryused = stMyMemStat%ullTotalPhys/1024/1024-stMyMemStat%ullAvailPhys/1024/1024

end subroutine memorymonitor  

endmodule ModMemoryMonitor