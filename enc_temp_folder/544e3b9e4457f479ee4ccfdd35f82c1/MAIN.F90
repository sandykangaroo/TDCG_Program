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
    implicit none
    integer::i
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    print*,'Generating mesh......'

    call CPU_TIME(tStart)
    if (Restart) return
    call GenerateBGMesh
    call initFindNeighbor
    call initSurfaceAdapt
    ! call initSmoothMesh
    call CPU_TIME(tEnd)
    ! call GetMinDistance

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
    !call GetMinDistance
    TimeStep=CFL*(BGCellSize(1)/2**InitRefineLVL)
    end subroutine TDCGSolver
!======================================================================
    subroutine AABBtime
    use ModPrecision
    use ModMesh
    use ModInpGlobal
    use ModMeshTools
    use ModInpMesh
    implicit none
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time
    integer :: i, j, k
    type(octCell),pointer :: t
    
    call CPU_TIME(tStart)
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        t=>Cell(i, j, k)
        call initCellCross(t)
    enddo
    enddo
    enddo
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "AABB time: ", tEnd-tStart
    
    end subroutine AABBtime
    
!======================================================================
    subroutine TDCGOutput(TimeStepStr)
    use ModPrecision
    use ModInpGlobal
    implicit none
    character(*),INTENT(IN) :: TimeStepStr
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

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