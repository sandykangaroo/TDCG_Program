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
    implicit none

    print*,'Welcome TDCGprogram'
    call TDCGReadInp
    call TDCGMesh
    call TDCGInitAll
    call TDCGSolver
    call TDCGOutput('OK')
    print*,' '

    end program TDCGmain
!======================================================================
    subroutine TDCGMesh
    use ModInpGlobal
    implicit none

    if (Restart) return
    print*,'Reading geometry......'
        call ReadGeometry
    print*,'Done'
    print*,'Generating mesh......'
        call GenerateBGMesh
        call SurfaceAdapt
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
    implicit none
    character(*),INTENT(IN) :: TimeStepStr

    if (OutputFormat=='.plt') then
        CALL OutputFlowFieldBinary(TimeStepStr,0)
    elseif (OutputFormat=='.szplt') then
        CALL OutputFlowFieldBinary(TimeStepStr,1)
    else
        CALL OutputFlowFieldASCII(TimeStepStr)
    endif

    endsubroutine TDCGOutput
!======================================================================
!======================================================================