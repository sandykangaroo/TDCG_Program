!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------!--------------------------------------------------------------------------------------------------
! Define Precision 
    module DefinePrecisionMod
        implicit none
! Double Precision 
        integer,parameter:: DoublePrec=kind(1.0d0)                
        integer,parameter::       Prec=DoublePrec
    end module  DefinePrecisionMod
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
! Define dynamic data structure
    module DefineDataTypeMod
        use DefinePrecisionMod
        implicit none
! One dimension dynamic real array
        type Real1DArray
            real(prec) R1D(:)
        end type Real1DArray
    
! Two dimension dynamic real array
        type Real2DArray
            real(prec) R2D(:,:)
        end type Real2DArray
    
! One dimension dynamic integer array
        type Integer1DArray
            integer I1D(:)
        end type Integer1DArray
    
! Two dimension dynamic integer array
        type Integer2DArray
            integer I2D(:,:)
        end type Integer2DArray
    
! Cartesian Cell's point coordinate
        type point
            real(prec) x, y, z
        end type point
    
! Triangle's three vertexs£¬odered by right-hand rule
        type triangle
            type(point) tripoint1, tripoint2, tripoint3
        end type triangle
    
! Cartesian grid data structure
! nroot      = number of the back-ground grid (root node) 
! levelx     = x-direct level
! levely     = y-direct level
! levelz     = z-direct level
! cross      = 0 the cartesian cell is outside the object surface
!            = 1 the cartesian cell intersects the object surface, and cell center is outside
!            = 2 the cartesian cell intersects the object surface, and cell center is inside
!            = 3 the cartesian cell is inside the object surface
! curve      = 0 the cartesian cell need not to be refined 
!            = 1 the cartesian cell need to be refined 
! fsplittype = 0,the cartesian cell is isotropical 
!            = 1,the cartesian cell is obtained by refined in x-direction
!            = 2,the cartesian cell is obtained by refined in y-direction
!            = 3,the cartesian cell is obtained by refined in z-direction
!            = 4,the cartesian cell is obtained by refined in xy-direction
!            = 5,the cartesian cell is obtained by refined in xz-direction
!            = 6,the cartesian cell is obtained by refined in yz-direction
! location   = 0 the cartesian cell does not have a son 
!            = 1,2,3,...8 the cartesian cell's location among siblings
! node       number of the eight vertexs
    
        type grid
            integer :: nroot,ngrid
            integer :: levelx, levely, levelz
            integer :: cross
            integer :: fsplittype,location
            integer :: node(8)
            real(prec) :: U(5)
            type(grid),pointer :: father,son(8),neighbor(6)     
        end type grid
    end module DefineDataTypeMod
    
        
! Define globally share constants
    module GlobalConstantsMod
        use DefinePrecisionMod 
        implicit none
    
! constants used in sutherlan'law, =110.3 in NSMB5.0
        real(prec),parameter:: C00=100.4

        real(prec),parameter:: Rgas=287    ![J/(kg.K)]
        real(prec),parameter:: gama=1.4, gama1=gama-1.0
        real(prec),parameter:: pr=0.72, prt=0.9
        real(prec),parameter:: PI=3.1415926535897932d0   

        integer,parameter:: SchemeNND2=1,SchemeWENO3=3          
        integer,parameter:: FluxRoe=1,Fluxcentral=2                                                      
        integer,parameter:: BCWall=2, BCSymmetry=3, BCFarfield=4                                        
        integer,parameter:: TimeRK3=1,TimeLUSGS=0                                                       
        integer,parameter:: TurSA=1,TurSST=2,TurKW=3  

    end module GlobalConstantsMod
!--------------------------------------------------------------------------------------------------
! Define globally share variables
    module GlobalvariablesMod
        use DefinePrecisionMod 
        implicit none

        integer,save:: nStep,nSave,nAdaptStep,nthreads
        integer,save:: Flagturblencemodel,FlagScheme,Flagflux,TimeMethod,nEqutions,Flagviscous
! Initial cell number
        integer,save:: nCellx,nCelly,nCellz,Total
        integer,save:: InitialRefineLevel,SurfaceAdaptRefineLevel,CurveAdaptRefineLevel,FlowAdaptRefineLevel
        integer,save:: NRRSeedNumber
        real(prec),save:: NRRLength,NRRTheta
        real(prec),save:: Ma,Re,gamma,CFL,Alfa,Beta,Twall,Kt00,Wt00,vt00
        real(prec),save:: DomainX,DomainY,DomainZ,h
    
        real(prec),save:: CpuTimeBegin, CpuTimeEnd
        integer,save::  CpuTimeUse
    
    
        
    ! primary variables of inflow
        real(prec),save:: rou00, u00, v00, w00, p00, T00
            
    end module GlobalvariablesMod

!--------------------------------------------------------------------------------------------------
