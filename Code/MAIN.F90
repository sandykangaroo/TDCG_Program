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
    use ModPrecision
    use ModInpGlobal,only: GeometryFormat
    implicit none
    real(R8):: tStart
    real(R8):: tEnd

    call CPU_TIME(tStart)

    call ReadInp
    if (GeometryFormat=='stl') then
        CALL ReadStl
    elseif (GeometryFormat=='facet') then
        CALL ReadFacet
    else
        stop 'error GeometryFormat'
    endif

    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "Total read time: ", tEnd-tStart

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
    real(R8):: tStartG   ! Start time
    real(R8):: tEndG     ! End time
    print*,'Generating mesh......'

    call CPU_TIME(tStart)
    if (Restart) return
    call GenerateBGMesh
    call initSurfaceAdapt
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
    subroutine TDCGOutput(TimeStepStr)
    use ModPrecision
    use ModInpGlobal
    implicit none
    character(*),INTENT(IN) :: TimeStepStr
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    if (OutputFormat=='plt') then
        CALL OutputFlowFieldBinary(TimeStepStr,0)
    elseif (OutputFormat=='szplt') then
        CALL OutputFlowFieldBinary(TimeStepStr,1)
    elseif (OutputFormat=='dat') then
        CALL OutputFlowFieldASCII(TimeStepStr)
    else
        stop 'error OutputFormat'
    endif
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "Subroutine-Output time: ", tEnd-tStart

    endsubroutine TDCGOutput
!======================================================================
!======================================================================