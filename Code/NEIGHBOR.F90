    module ModNeighbor
        use ModMesh
        use ModTypDef
        use ModInpMesh
        implicit none
    contains

    function Neighbor(c,dirct)
    implicit none
    type(FTTCell),pointer :: Neighbor
    type(FTTCell),pointer :: c
    integer               :: dirct

    select case (dirct)
    case (1)
        Neighbor => NeighborX1(c)
    case (2)
        Neighbor => NeighborX2(c)
    case (3)
        Neighbor => NeighborY1(c)
    case (4)
        Neighbor => NeighborY2(c)
    case (5)
        Neighbor => NeighborZ1(c)
    case (6)
        Neighbor => NeighborZ2(c)
    end select
    endfunction Neighbor

! X- neighbor
    recursive function NeighborX1(c) RESULT(Nei)
        implicit none
        type(FTTCell),pointer :: c,cf,cinitial,cfn,Nei
        !cf,c'father;cinitial,c'initial X- Nei;cfn,c'father's X- Nei

        cf=>c%father
        if(.not.ASSOCIATED(cf))then
            if(c%nBGCell(1)==1)then
                Nei=>null()
            else
                Nei=>OctCell(c%nBGCell(1)-1,c%nBGCell(2),c%nBGCell(3))
            endif
            return
        endif

        cfn=>cf%Octson%Neighbor1
        if (.not.ASSOCIATED(cfn)) then
            select case (c%fSplitTyp)
            case(0)
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==4)then
                    Nei=>null()
                elseif(c%Location==5)then
                    Nei=>null()
                elseif(c%Location==6)then
                    Nei=>cf%Octson%son(5)
                elseif(c%Location==7)then
                    Nei=>cf%Octson%son(8)
                elseif(c%Location==8)then
                    Nei=>null()
                endif
            case(1)
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                endif
            case(2)
                    Nei=>null()
            case(3)
                    Nei=>null()
            case(4)
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==4)then
                    Nei=>null()
                endif
            case(5)
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==4)then
                    Nei=>null()
                endif
            case(6)
                    Nei=>null()
            end select
        else
            select case (c%fSplitTyp)
            case(0)
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0 .or. &
                            cfn%Octson%son(1)%fSplitTyp==1 .or. &
                            cfn%Octson%son(1)%fSplitTyp==4 .or. &
                            cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(2)
                        else !cfn%Octson%son(1)%fSplitTyp==2  3  6
                            Nei=>cfn%Octson%son(1)
                        endif
                    endif
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if( cfn%Octson%son(1)%fSplitTyp==0 .or. &
                            cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(3)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else !cfn%Octson%son(1)%fSplitTyp==1  2  5  6
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                elseif(c%Location==5)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(6)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or. &
                                cfn%Octson%son(1)%fSplitTyp==3.or. &
                                cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(3)
                        else
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                elseif(c%Location==6)then
                    Nei=>cf%Octson%son(5)
                elseif(c%Location==7)then
                    Nei=>cf%Octson%son(8)
                elseif(c%Location==8)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(7)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or. &
                                cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        else !cfn%Octson%son(1)%fSplitTyp==4  5  6
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                endif
            case(1)
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(2)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                endif
            case(2)
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                                Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                                Nei=>cfn%Octson%son(1)
                        else
                                Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==2)then
                                Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4)then
                                Nei=>cfn%Octson%son(3)
                        else
                                Nei=>cfn
                        endif
                    endif
                endif
            case(3)
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(3)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            case(4)
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(3)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            case(5)
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(3)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            case (6)
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or. &
                            cfn%Octson%son(1)%fSplitTyp==3.or. &
                            cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(1)
                        else !cfn%Octson%son(1)%fSplitTyp==0  1  4  5
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    else ! if(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0.or. &
                            cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(3)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else !cfn%Octson%son(1)%fSplitTyp==1  2  5  6
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(7)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or.&
                               cfn%Octson%son(1)%fSplitTyp==2.or. &
                               cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4.or.&
                               cfn%Octson%son(1)%fSplitTyp==5.or. &
                               cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(6)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or.&
                               cfn%Octson%son(1)%fSplitTyp==3.or. &
                               cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(3)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                endif
            end select
        endif
    end function NeighborX1

! X+ Nei
    recursive function NeighborX2(c) RESULT(Nei)
        implicit none
        type(FTTCell),pointer :: c,cf,cinitial,cfn,Nei
        !cf,c'father;cinitial,c'initial X+ Nei;cfn,c'father's X+ Nei

        cf=>c%father
        if(.not.ASSOCIATED(cf))then
            if(c%nBGCell(1)==nCell(1))then
                Nei=>null()
            else
                Nei=>OctCell(c%nBGCell(1)+1,c%nBGCell(2),c%nBGCell(3))
            endif
            return
        endif

        cfn=>cf%Octson%Neighbor2
        if(ASSOCIATED(cfn))then
            if(c%fSplitTyp==0)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        Nei=>cfn%Octson%son(1)
                    endif
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or.&
                               cfn%Octson%son(1)%fSplitTyp==3.or. &
                               cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2.or. &
                               cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==5)then
                    Nei=>cf%Octson%son(6)
                elseif(c%Location==6)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if (cfn%Octson%son(1)%fSplitTyp==1.or.&
                            cfn%Octson%son(1)%fSplitTyp==2.or. &
                            cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(5)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5.or. &
                               cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                elseif(c%Location==7)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(8)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2.or. &
                               cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4.or. &
                               cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                elseif(c%Location==8)then
                    Nei=>cf%Octson%son(7)
                endif
            elseif(c%fSplitTyp==1)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==2)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if (cfn%Octson%son(1)%fSplitTyp==1.or.&
                            cfn%Octson%son(1)%fSplitTyp==2.or. &
                            cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==3)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if (cfn%Octson%son(1)%fSplitTyp==1.or.&
                            cfn%Octson%son(1)%fSplitTyp==3.or. &
                            cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(4)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==4)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or.&
                           cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(4)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(3)
                endif
            elseif(c%fSplitTyp==5)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or.&
                           cfn%Octson%son(1)%fSplitTyp==3.or. &
                           cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(4)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(3)
                endif
            elseif(c%fSplitTyp==6)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        Nei=>cfn%Octson%son(1)      !may be bug
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==0.or. &
                               cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(4)
                        else
                            Nei=>cfn%Octson%son(1)
                        endif
                    endif
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(8)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2.or. &
                               cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4.or. &
                               cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(5)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or.&
                               cfn%Octson%son(1)%fSplitTyp==2.or. &
                               cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5.or. &
                               cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                endif
            endif
        else ! .not.ASSOCIATED(cfn)
            if(c%fSplitTyp==0)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    Nei=>null()
                elseif(c%Location==3)then
                    Nei=>null()
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==5)then
                    Nei=>cf%Octson%son(6)
                elseif(c%Location==6)then
                    Nei=>null()
                elseif(c%Location==7)then
                    Nei=>null()
                elseif(c%Location==8)then
                    Nei=>cf%Octson%son(7)
                endif
            elseif(c%fSplitTyp==1)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    Nei=>null()
                endif
            elseif(c%fSplitTyp==2)then
                    Nei=>null()
            elseif(c%fSplitTyp==3)then
                    Nei=>null()
            elseif(c%fSplitTyp==4)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    Nei=>null()
                elseif(c%Location==3)then
                    Nei=>null()
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(3)
                endif
            elseif(c%fSplitTyp==5)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    Nei=>null()
                elseif(c%Location==3)then
                    Nei=>null()
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(3)
                endif
            elseif(c%fSplitTyp==6)then
                    Nei=>null()
            endif
        endif
        return
    end function NeighborX2

! Y- Nei
    recursive function NeighborY1(c) RESULT(Nei)
        implicit none
            type(FTTCell),pointer :: c,cf,cinitial,cfn,Nei
            !cf,c'father;cinitial,c'initial Y- Nei;cfn,c'father's Y- Nei

        cf=>c%father
        if(.not.ASSOCIATED(cf))then
            if(c%nBGCell(2)==1)then
                Nei=>null()
            else
                Nei=>OctCell(c%nBGCell(1),c%nBGCell(2)-1,c%nBGCell(3))
            endif
            return
        endif

        cfn=>cf%Octson%Neighbor3
        if(ASSOCIATED(cfn))then
            if(c%fSplitTyp==0)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or.&
                           cfn%Octson%son(1)%fSplitTyp==3.or. &
                           cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==0.or. &
                               cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2.or. &
                               cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==5.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==0.or. &
                               cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(3)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        endif
                    endif
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==5)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(8)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2.or. &
                               cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4.or. &
                               cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                elseif(c%Location==6)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(7)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or.&
                               cfn%Octson%son(1)%fSplitTyp==2.or. &
                               cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif (cfn%Octson%son(1)%fSplitTyp==4.or.&
                                cfn%Octson%son(1)%fSplitTyp==5.or. &
                                cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                elseif(c%Location==7)then
                    Nei=>cf%Octson%son(6)
                elseif(c%Location==8)then
                    Nei=>cf%Octson%son(5)
                endif
            elseif(c%fSplitTyp==1)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(2)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(3)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==2)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(2)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                endif
            elseif(c%fSplitTyp==3)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==4)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(4)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(3)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(1)
                endif
            elseif(c%fSplitTyp==5)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(4)
                        elseif (cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==3.or. &
                                cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(1)
                        elseif (cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(3)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else                         !  1  2  5  6
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(7)
                        elseif (cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif (cfn%Octson%son(1)%fSplitTyp==4.or.&
                                cfn%Octson%son(1)%fSplitTyp==5.or. &
                                cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(8)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif (cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif (cfn%Octson%son(1)%fSplitTyp==4.or. &
                                cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==6)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            endif
        else ! .not.ASSOCIATED(cfn)
            if(c%fSplitTyp==0)then
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>null()
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==5)then
                    Nei=>null()
                elseif(c%Location==6)then
                    Nei=>null()
                elseif(c%Location==7)then
                    Nei=>cf%Octson%son(6)
                elseif(c%Location==8)then
                    Nei=>cf%Octson%son(5)
                endif
            elseif(c%fSplitTyp==1)then
                    Nei=>null()
            elseif(c%fSplitTyp==2)then
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                endif
            elseif(c%fSplitTyp==3)then
                    Nei=>null()
            elseif(c%fSplitTyp==4)then
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>null()
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(1)
                endif
            elseif(c%fSplitTyp==5)then
                    Nei=>null()
            elseif(c%fSplitTyp==6)then
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==4)then
                    Nei=>null()
                endif
            endif
        endif
        return
    end function NeighborY1

! Y+ Nei
    recursive function NeighborY2(c) RESULT(Nei)
        implicit none
            type(FTTCell),pointer :: c,cf,cinitial,cfn,Nei
            !cf,c'father;cinitial,c'initial Y+ Nei;cfn,c'father's Y+ Nei

        cf=>c%father
        if(.not.ASSOCIATED(cf))then
            if(c%nBGCell(2)==nCell(2))then
                Nei=>null()
            else
                Nei=>OctCell(c%nBGCell(1),c%nBGCell(2)+1,c%nBGCell(3))
            endif
            return
        endif

        cfn=>cf%Octson%Neighbor4
        if(ASSOCIATED(cfn))then
            if(c%fSplitTyp==0)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==4.or. &
                           cfn%Octson%son(1)%fSplitTyp==5.or. &
                           cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(2)
                        elseif (cfn%Octson%son(1)%fSplitTyp==2.or.&
                                cfn%Octson%son(1)%fSplitTyp==3.or. &
                                cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(1)
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        Nei=>cfn%Octson%son(1)
                    endif
                elseif(c%Location==5)then
                    Nei=>cf%Octson%son(8)
                elseif(c%Location==6)then
                    Nei=>cf%Octson%son(7)
                elseif(c%Location==7)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(6)
                        elseif (cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==3.or. &
                                cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(3)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                elseif(c%Location==8)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(5)
                        elseif (cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif (cfn%Octson%son(1)%fSplitTyp==5.or. &
                                cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==1)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or.&
                           cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==2)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==3)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or.&
                           cfn%Octson%son(1)%fSplitTyp==3.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==4)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or.&
                           cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==5)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        Nei=>cfn%Octson%son(1)
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or.&
                           cfn%Octson%son(1)%fSplitTyp==3.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(2)
                        else
                            Nei=>cfn%Octson%son(1)
                        endif
                    endif
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(6)
                        elseif (cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==3.or. &
                                cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(3)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(5)
                        elseif (cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif (cfn%Octson%son(1)%fSplitTyp==5.or. &
                                cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==6)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or.&
                           cfn%Octson%son(1)%fSplitTyp==3.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(3)
                endif
            endif
        else ! .not.ASSOCIATED(cfn)
            if(c%fSplitTyp==0)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==3)then
                    Nei=>null()
                elseif(c%Location==4)then
                    Nei=>null()
                elseif(c%Location==5)then
                    Nei=>cf%Octson%son(8)
                elseif(c%Location==6)then
                    Nei=>cf%Octson%son(7)
                elseif(c%Location==7)then
                    Nei=>null()
                elseif(c%Location==8)then
                    Nei=>null()
                endif
            elseif(c%fSplitTyp==1)then
                    Nei=>null()
            elseif(c%fSplitTyp==2)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    Nei=>null()
                endif
            elseif(c%fSplitTyp==3)then
                    Nei=>null()
            elseif(c%fSplitTyp==4)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==3)then
                    Nei=>null()
                elseif(c%Location==4)then
                    Nei=>null()
                endif
            elseif(c%fSplitTyp==5)then
                    Nei=>null()
            elseif(c%fSplitTyp==6)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    Nei=>null()
                elseif(c%Location==3)then
                    Nei=>null()
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(3)
                endif
            endif
        endif
        return
    end function NeighborY2

! Z- Nei
    recursive function NeighborZ1(c) RESULT(Nei)
        implicit none
            type(FTTCell),pointer :: c,cf,cinitial,cfn,Nei
            !cf,c'father;cinitial,c'initial Z- Nei;cfn,c'father's Z- Nei

        cf=>c%father
        if(.not.ASSOCIATED(cf))then
            if(c%nBGCell(3)==1)then
                Nei=>null()
            else
                Nei=>OctCell(c%nBGCell(1),c%nBGCell(2),c%nBGCell(3)-1)
            endif
            return
        endif

        cfn=>cf%Octson%Neighbor5
        if(ASSOCIATED(cfn))then
            if(c%fSplitTyp==0)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(5)
                        elseif  (cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5.or. &
                               cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(6)
                        elseif (cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==3.or. &
                                cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(3)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(7)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        else                        !4   5   6
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(8)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4.or. &
                                cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                elseif(c%Location==5)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==6)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==7)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==8)then
                    Nei=>cf%Octson%son(4)
                endif
            elseif(c%fSplitTyp==1)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(4)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(3)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==2)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==3)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                endif
            elseif(c%fSplitTyp==4)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(5)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5.or. &
                                cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(6)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==3.or. &
                                cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(3)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        endif
                    endif
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(7)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        else
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0)then
                            Nei=>cfn%Octson%son(8)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2.or. &
                                cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==4.or. &
                                cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==5)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(4)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(3)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(1)
                endif
            elseif(c%fSplitTyp==6)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2)then
                            Nei=>cfn%Octson%son(1)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(4)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(3)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(1)
                endif
            endif
        else ! .not.ASSOCIATED(cfn)
            if(c%fSplitTyp==0)then
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>null()
                elseif(c%Location==3)then
                    Nei=>null()
                elseif(c%Location==4)then
                    Nei=>null()
                elseif(c%Location==5)then
                    Nei=>cf%Octson%son(1)
                elseif(c%Location==6)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==7)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==8)then
                    Nei=>cf%Octson%son(4)
                endif
            elseif(c%fSplitTyp==1)then
                    Nei=>null()
            elseif(c%fSplitTyp==2)then
                    Nei=>null()
            elseif(c%fSplitTyp==3)then
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(1)
                endif
            elseif(c%fSplitTyp==4)then
                    Nei=>null()
            elseif(c%fSplitTyp==5)then
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>null()
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(1)
                endif
            elseif(c%fSplitTyp==6)then
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>null()
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(1)
                endif
            endif
        endif
        return
    end function NeighborZ1

! Z+ Nei
    recursive function NeighborZ2(c) RESULT(Nei)
        implicit none
            type(FTTCell),pointer :: c,cf,cinitial,cfn,Nei
            !cf,c'father;cinitial,c'initial Z+ Nei;cfn,c'father's Z+ Nei

        cf=>c%father
        if(.not.ASSOCIATED(cf))then
            if(c%nBGCell(3)==nCell(3))then
                Nei=>null()
            else
                Nei=>OctCell(c%nBGCell(1),c%nBGCell(2),c%nBGCell(3)+1)
            endif
            return
        endif

        cfn=>cf%Octson%Neighbor6
        if(ASSOCIATED(cfn))then
            if(c%fSplitTyp==0)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(5)
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(6)
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(7)
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(8)
                elseif(c%Location==5)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        Nei=>cfn%Octson%son(1)
                    endif
                elseif(c%Location==6)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or.&
                           cfn%Octson%son(1)%fSplitTyp==3.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                elseif(c%Location==7)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(3)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                elseif(c%Location==8)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==1.or.&
                                cfn%Octson%son(1)%fSplitTyp==3.or. &
                                cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==1)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or.&
                           cfn%Octson%son(1)%fSplitTyp==3.or. &
                           cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==2)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or.&
                           cfn%Octson%son(1)%fSplitTyp==3.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==3)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    Nei=>cfn
                endif
            elseif(c%fSplitTyp==4)then
                if(c%Location==1)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        Nei=>cfn%Octson%son(1)
                    endif
                elseif(c%Location==2)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or.&
                           cfn%Octson%son(1)%fSplitTyp==3.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(3)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn%Octson%son(2)
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==0.or. &
                           cfn%Octson%son(1)%fSplitTyp==4)then
                            Nei=>cfn%Octson%son(4)
                        elseif(cfn%Octson%son(1)%fSplitTyp==2.or. &
                               cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(2)
                        else
                            Nei=>cfn%Octson%son(1)
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==5)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or. &
                           cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==1.or.&
                           cfn%Octson%son(1)%fSplitTyp==3.or. &
                           cfn%Octson%son(1)%fSplitTyp==5)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitTyp==6)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==3)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(2)
                        elseif(cfn%Octson%son(1)%fSplitTyp==3)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.ASSOCIATED(cfn%Octson))then
                        Nei=>cfn
                    elseif(ASSOCIATED(cfn%Octson))then
                        if(cfn%Octson%son(1)%fSplitTyp==2.or.&
                           cfn%Octson%son(1)%fSplitTyp==3.or. &
                           cfn%Octson%son(1)%fSplitTyp==6)then
                            Nei=>cfn%Octson%son(1)
                        else
                            Nei=>cfn
                        endif
                    endif
                endif
            endif
        else ! .not.ASSOCIATED(cfn)
            if(c%fSplitTyp==0)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(5)
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(6)
                elseif(c%Location==3)then
                    Nei=>cf%Octson%son(7)
                elseif(c%Location==4)then
                    Nei=>cf%Octson%son(8)
                elseif(c%Location==5)then
                    Nei=>null()
                elseif(c%Location==6)then
                    Nei=>null()
                elseif(c%Location==7)then
                    Nei=>null()
                elseif(c%Location==8)then
                    Nei=>null()
                endif
            elseif(c%fSplitTyp==1)then
                if(c%Location==1)then
                    Nei=>null()
                elseif(c%Location==2)then
                    Nei=>null()
                endif
            elseif(c%fSplitTyp==2)then
                    Nei=>null()
            elseif(c%fSplitTyp==3)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(2)
                elseif(c%Location==2)then
                    Nei=>null()
                endif
            elseif(c%fSplitTyp==4)then
                    Nei=>null()
            elseif(c%fSplitTyp==5)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==3)then
                    Nei=>null()
                elseif(c%Location==4)then
                    Nei=>null()
                endif
            elseif(c%fSplitTyp==6)then
                if(c%Location==1)then
                    Nei=>cf%Octson%son(4)
                elseif(c%Location==2)then
                    Nei=>cf%Octson%son(3)
                elseif(c%Location==3)then
                    Nei=>null()
                elseif(c%Location==4)then
                    Nei=>null()
                endif
            endif
        endif
        return
    end function NeighborZ2

    end module ModNeighbor
