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
                        OutputNameStr, GeometryName, OutputFormat
    NAMELIST /Mesh/     nCell, DomainMin, DomainMax, InitRefineLVL,       &
                        AdaptRefineLVL, cIntersectMethod,             &
                        PaintingAlgorithmMethod
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
    ! subroutine ReadGeometry
    ! use ModInpGlobal
    ! use ModMesh
    ! implicit none
    ! integer :: ios, i, j 
    ! character(10)::FileForm

    ! print*,'File name: ', GeometryName
    ! open(unit=31, file=GeometryName, iostat=ios, status="old", action="read")
    ! if ( ios /= 0 ) stop ("Error opening file: "//GeometryName)
    ! read(31, fmt="(G10.1)", iostat=ios) FileForm
    ! if ( ios /= 0 ) stop ("Error reading file: "//GeometryName)
    ! if (.NOT.FileForm=="FACET FILE")              &
    !     stop ("Error file header: "//GeometryName)
    ! read(31,"(///I9)") nGeoPoints
    ! print*,'Geometry points count: ', nGeoPoints
    ! ALLOCATE(Geometry(nGeoPoints,3))
    ! read(31,"(3(E23.15,1X))",err=10) ((Geometry(i,j),j=1,3),i=1,nGeoPoints)
    ! read(31,"(/)")
    ! read(31,*) nGeoFaces
    ! print*,'Geometry faces count: ', nGeoFaces
    ! ALLOCATE(GeoFace(nGeoFaces,3))
    ! read(31,"(3(I7,1X))") ((GeoFace(i,j),j=1,3),i=1,nGeoFaces)
    ! close(31)
    ! RETURN
    ! 10  continue
    !     backspace(31)
    !     read(31,"(3(E14.6,1X))") ((Geometry(i,j),j=1,3),i=1,nGeoPoints)
    !     read(31,"(/)")
    !     read(31,*) nGeoFaces
    !     print*,'Geometry faces count: ', nGeoFaces
    !     ALLOCATE(GeoFace(nGeoFaces,3))
    !     read(31,"(3(I7,1X))") ((GeoFace(i,j),j=1,3),i=1,nGeoFaces)
    !     close(31)
    !     RETURN
    ! endsubroutine ReadGeometry
!======================================================================
    subroutine ReadGeometryBinary
    endsubroutine ReadGeometryBinary
!======================================================================
!======================================================================