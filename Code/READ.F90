!======================================================================
    subroutine ReadInp
    use ModInpGlobal
    use ModInpInflow
    use ModInpMesh
    use ModInpNRRset
    implicit none
    logical :: exist
    integer :: ios

    NAMELIST /Global/   Aniso, nStep, nSave, nAdaptStep, Debug, CFL,  &
                        NRR, Chimera, Limiter, Restart, nGeometry,    &
                        OutputName, GeometryName, OutputFormat,    &
                        MeshOnly, GeometryFormat
    NAMELIST /Mesh/     nCell, DomainMin, DomainMax, InitRefineLVL,       &
                        AdaptRefineLVL, cIntersectMethod, useKDT
    NAMELIST /NRRset/   NRRSeed, NRRLength, NRRTheta
    NAMELIST /Inflow/   Alpha, Beta, Re, T00, Ma00, Gama00, Rgas

    INQUIRE(file='NameList.inp', exist=exist)
        if (.NOT.exist) stop 'Error====> Not found file NameList.inp'
    open(unit=21, file='NameList.inp', iostat=ios, status="old", action="read")
        if ( ios /= 0 ) stop "Error====> Error opening file NameList.inp"
        read (21, NML=Global);          rewind (21)
        if (.NOT.Restart) then
            read (21, NML=Mesh);        rewind (21)
            read (21, NML=Inflow);      rewind (21)
        endif
        if (NRR) read (21, NML=NRRset); rewind (21)
    close(21)
    end subroutine ReadInp
!======================================================================
    subroutine ReadPLT

    endsubroutine ReadPLT
!======================================================================
!======================================================================
!======================================================================