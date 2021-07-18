!======================================================================
!  
!  TDCG Program ---- Three Demision Cartesian Grid Program
!  
!                          ����������������
!                        ����  ��������������
!                        ��������������������
!                        ��������������������
!                        ��������������������
!                        ��������
!                        ����������������
!  ��                  ��������
!  ��              ������������
!  ����        ��������������������
!  ������    ������������������  ��
!  ����������������������������
!  ����������������������������
!    ������������������������
!      ����������������������
!        ������������������
!          ��������������
!            ������  ����
!            ����      ��
!            ��        ��
!            ����      ����
!  
!  1...
!  2...
!  3...
!  Owuuuuuu~~~~~~~~
!  
!  An unreliable CFD solver based on the Cartesian grid
!  
!======================================================================
!  Author:
!  Xueliang Li
!  Central South University, Changsha, China
!  lixueliang@csu.edu.cn
!======================================================================
!======================================================================
    program TDCGmain
    use ModPrecision
    implicit none
    real(R8):: tStart
    real(R8):: tEnd

    print*,'Welcome TDCGprogram'
    call CPU_TIME(tStart)
    call TDCGRead
    call TDCGPerporcessing
    call TDCGMesh
    call TDCGInitAll
    call TDCGSolver
    call TDCGOutput('OK')
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') 'Program running time: ', tEnd-tStart

    end program TDCGmain
!======================================================================
    subroutine TDCGRead
    implicit none

    print*, 'Reading input file NameList.inp......'
        call ReadInp
    print*,'Done'
    print*,'Reading geometry file......'
        call ReadGeometry
    print*,'Done'

    endsubroutine TDCGRead
!======================================================================
    subroutine TDCGPerporcessing
    !call BuildGeoBBOX
    endsubroutine TDCGPerporcessing
!======================================================================
    subroutine TDCGMesh
    use ModPrecision
    use ModInpGlobal
    use ModInpMesh
    use ModMemoryMonitor
    use Kernel32
    use ISO_C_Binding 
    implicit none
    integer::i,a,b
    integer(DWORDLONG)  ::ii
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time
    real(R8):: tStartG   ! Start time
    real(R8):: tEndG     ! End time
    print*,'Generating mesh......'
    
    call memorymonitor(ii)

    call CPU_TIME(tStart)
    if (Restart) return
    call CPU_TIME(tStartG)
    call GenerateBGMesh
    call CPU_TIME(tEndG)
    call initFindNeighbor
    call initSurfaceAdapt
    call CPU_TIME(tEnd)

    !call GetMinDistance
    b=0
     !inquire( file = fort , exist = file_exist )
 !    open(1,file='fort.102')
 !    do i=1,6990126
 !    read(1,*) a
 !    b=b+a
 !enddo
 !write(*,*)b
 !close(1)
     

    write(*,'(1X,A,F10.5)') "BackG Mesh generation time: ", tEndG-tStartG
    write(*,'(1X,A,F10.2)') "Total Mesh generation time: ", tEnd-tStart
    print*,'Done'

    end subroutine TDCGMesh
!======================================================================
    subroutine TDCGInitAll
    end subroutine TDCGInitAll
!======================================================================
    subroutine TDCGSolver
    use ModSolve
    use ModInpGlobal
    use ModMesh
    use ModInpMesh
    implicit none
    TimeStep=CFL*(BGCellSize(1)/2**InitRefineLVL)
    end subroutine TDCGSolver
!======================================================================
    subroutine TDCGOutput(TimeStepStr)
    use ModPrecision
    use ModInpGlobal
    use ModMemoryMonitor
    use Kernel32
    use ISO_C_Binding 
    implicit none
    character(*),INTENT(IN) :: TimeStepStr
    integer(DWORDLONG)  ::i 
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time
    
    call memorymonitor(i)
    !write(*,*)"There is  ",i,"MB of memory in use." 

    call CPU_TIME(tStart)
    if (OutputFormat=='.plt') then
        CALL OutputFlowFieldBinary(TimeStepStr,0)
    elseif (OutputFormat=='.szplt') then
        CALL OutputFlowFieldBinary(TimeStepStr,1)
    else
        CALL OutputFlowFieldASCII(TimeStepStr)
    endif
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "Subroutine-Output time: ", tEnd-tStart

    endsubroutine TDCGOutput
!======================================================================
!======================================================================