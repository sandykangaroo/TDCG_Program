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
    use ModTime
    use ModPrecision
    implicit none
    real(R8):: tProgramStart
    real(R8):: tProgramEnd

    print*,'Welcome TDCGprogram'
    call CPU_TIME(tProgramStart)
    call TDCGRead
    call TDCGPerporcessing
    call TDCGMesh
    call TDCGInitAll
    call TDCGSolver
    call TDCGOutput('OK')
    call CPU_TIME(tProgramEnd)
    print*,'Program running time: ', tProgramEnd-tProgramStart

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
    use ModInpGlobal
    use ModTime
    implicit none

    print*,'Generating mesh......'
    call CPU_TIME(tStart)
    if (Restart) return
        call GenerateBGMesh
    call CPU_TIME(tEnd)
    print*,"Subroutine-BGMeshCross time: ", tEnd-tStart

    call CPU_TIME(tStart)
        call initSurfaceAdapt
    call CPU_TIME(tEnd)
    print*,"Subroutine-SurfaceAdapt time: ", tEnd-tStart
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
    TimeStep=CFL*(BGStep(1)/2**InitRefineLVL)
    end subroutine TDCGSolver
!======================================================================
    subroutine TDCGOutput(TimeStepStr)
    use ModInpGlobal
    use ModTime
    implicit none
    character(*),INTENT(IN) :: TimeStepStr

    call CPU_TIME(tStart)
    if (OutputFormat=='.plt') then
        CALL OutputFlowFieldBinary(TimeStepStr,0)
    elseif (OutputFormat=='.szplt') then
        CALL OutputFlowFieldBinary(TimeStepStr,1)
    else
        CALL OutputFlowFieldASCII(TimeStepStr)
    endif
    call CPU_TIME(tEnd)
    print*,"Subroutine-Output time: ", tEnd-tStart

    endsubroutine TDCGOutput
!======================================================================
!======================================================================