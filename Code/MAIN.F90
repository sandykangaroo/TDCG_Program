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
    call dellocateTri
    call TDCGOutput('OK')
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') 'Program running time: ', tEnd-tStart

    end program TDCGmain
!======================================================================
    subroutine TDCGRead
    use ModTypDef
    use ModInpGlobal,only: GeometryFormat
    implicit none
    real(R8)                    :: tStart
    real(R8)                    :: tEnd

    call CPU_TIME(tStart)

        call ReadInp
    if (GeometryFormat=='stl') then
        CALL ReadStl
    elseif (GeometryFormat=='facet') then
        CALL ReadFacet
    else
        stop 'error GeometryFormat'
    endif
    ! call initCheckVector

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
    use ModWalldistance
    implicit none
    integer::i,a,b
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

    call GetWallDistance

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
    implicit none
    character(*),INTENT(IN) :: TimeStepStr
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    select case(OutputFormat)
    case ('plt')
        CALL OutputFlowFieldBinary(TimeStepStr,0)
    case ('szplt')
        CALL OutputFlowFieldBinary(TimeStepStr,1)
    case ('dat')
        CALL OutputFlowFieldASCII(TimeStepStr)
    case default
        stop 'error OutputFormat'
    end select
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "Subroutine-Output time: ", tEnd-tStart

    endsubroutine TDCGOutput
!======================================================================
!======================================================================
    subroutine dellocateTri
    use ModGeometry
    use ModTypDef
    use ModKDTree

    call aaa(KDTree(1)%root)
    DEALLOCATE(KDTree)
    ! DEALLOCATE(body)
    contains
        recursive subroutine aaa(tree)
        implicit none
        type(KDT_node),pointer::tree
        if(ASSOCIATED(tree%left))then
            call aaa(tree%left)
        endif
        if(ASSOCIATED(tree%right))then
            call aaa(tree%right)
        endif
        DEALLOCATE(tree%the_data)
        NULLIFY(tree%parent)
        DEALLOCATE(tree)

        endsubroutine aaa
    endsubroutine