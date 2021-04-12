    module ModNeighbor
        use ModMesh
        use ModTypDef 
        use ModInpMesh
        implicit none

    contains

! Initial X- neighbor    
    function InitialNeighborX1(c)
        implicit none
        type(OctCell),pointer :: InitialNeighborX1, c
        integer               :: i,m,n,l,j,k
        
        m=nCell(1)
        n=nCell(2)
        l=nCell(3)
        
        if(mod(c%nCell-1,m)==0)then
            InitialNeighborX1=>null()     
        else
            InitialNeighborX1=>Cell(i-1,j,k)
        endif
        return    
    end function InitialNeighborX1

! Initial X+ neighbor     
    function InitialNeighborX2(c)
        implicit none
        type(OctCell),pointer :: InitialNeighborX2, c
        integer               :: i,m,n,l,j,k
        
        m=nCell(1)
        n=nCell(2)
        l=nCell(3)
        
        if(mod(c%nCell,m)==0)then
            InitialNeighborX2=>null()     
        else
            InitialNeighborX2=>Cell(i+1,j,k)
        endif
        return    
    end function InitialNeighborX2
    
! Initial Y- neighbor     
    function InitialNeighborY1(c)
        implicit none
        type(OctCell),pointer :: InitialNeighborY1, c
        integer               :: i,m,n,l,j,k
        
        m=nCell(1)
        n=nCell(2)
        l=nCell(3)
        
        if(mod(c%nCell,m*n)<=m .and. mod(c%nCell,m*n)/=0)then
            InitialNeighborY1=>null()     
        else
            InitialNeighborY1=>Cell(i,j-1,k)
        endif
        return    
    end function InitialNeighborY1 
    
! Initial Y+ neighbor     
    function InitialNeighborY2(c)
        implicit none
        type(OctCell),pointer :: InitialNeighborY2, c
        integer               :: i,m,n,l,j,k
        
        m=nCell(1)
        n=nCell(2)
        l=nCell(3)
        
        if(mod(c%nCell,m*n)>(n-1)*m .or. mod(c%nCell,m*n)==0)then
            InitialNeighborY2=>null()     
        else
            InitialNeighborY2=>Cell(i,j+1,k)
        endif
        return    
    end function InitialNeighborY2 
    
! Initial Z- neighbor     
    function InitialNeighborZ1(c)
        implicit none
        type(OctCell),pointer :: InitialNeighborZ1, c
        integer               :: i,m,n,l,j,k
        
        m=nCell(1)
        n=nCell(2)
        l=nCell(3)
        
        if(c%nCell<=m*n)then
            InitialNeighborZ1=>null()     
        else
            InitialNeighborZ1=>Cell(i,j,k-1)
        endif
        return    
    end function InitialNeighborZ1
    
! Initial Z+ neighbor     
    function InitialNeighborZ2(c)
        implicit none
        type(OctCell),pointer :: InitialNeighborZ2, c
        integer               :: i,m,n,l,j,k
        
        m=nCell(1)
        n=nCell(2)
        l=nCell(3)
        
        if(c%nCell>m*n*(l-1))then
            InitialNeighborZ2=>null()     
        else
            InitialNeighborZ2=>Cell(i,j,k+1)
        endif
        return    
    end function InitialNeighborZ2
    
! X- neighbor 
    recursive function NeighborX1(c)
        implicit none
        type(OctCell),pointer :: c,cf,cinitial,cfn,NeighborX1
        !cf,c'father;cinitial,c'initial X- neighbor;cfn,c'father's X- neighbor
        
        cinitial=>InitialNeighborX1(c)
        NeighborX1=>null()
        cf=>c%father
        
        if(.not.associated(cf))then
            NeighborX1=>cinitial
            return
        endif
        
        cfn=>NeighborX1(cf)
        if(associated(cfn))then
            if(c%fSplitType==0)then
                if(c%Location==1)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==1 .or. &
                           cfn%son1%fSplitType==4 .or. &
                           cfn%son1%fSplitType==5)then
                            NeighborX1=>cfn%son2
                        else !cfn%son1%fSplitType==0  2  3  6
                            NeighborX1=>cfn%son1   
                        endif
                    endif
                elseif(c%Location==2)then
                    NeighborX1=>cf%son1
                elseif(c%Location==3)then
                    NeighborX1=>cf%son4
                elseif(c%Location==4)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
                            NeighborX1=>cfn%son3
                        elseif(cfn%son1%fSplitType==3)then                         
                            NeighborX1=>cfn%son1
                        else
                            NeighborX1=>cfn%son2                !cfn%son1%fSplitType==1  2  5  6
                        endif
                    endif
                elseif(c%Location==5)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==0)then
                            NeighborX1=>cfn%son6
                        elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==4)then                         
                            NeighborX1=>cfn%son2
                        elseif(cfn%son1%fSplitType==2)then      
                            NeighborX1=>cfn%son1
                        elseif(cfn%son1%fSplitType==5)then      
                            NeighborX1=>cfn%son3
                        else      
                            NeighborX1=>cfn%son4                      
                        endif
                    endif
                elseif(c%Location==6)then
                    NeighborX1=>cf%son5
                elseif(c%Location==7)then
                    NeighborX1=>cf%son8
                elseif(c%Location==8)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==0)then
                            NeighborX1=>cfn%son7
                        elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then                         
                            NeighborX1=>cfn%son2
                        else                         !cfn%son1%fSplitType==4  5  6
                            NeighborX1=>cfn%son3                      
                        endif
                    endif
                endif
            elseif(c%fSplitType==1)then
                if(c%Location==1)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==1)then
                            NeighborX1=>cfn%son2
                        else
                            NeighborX1=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    NeighborX1=>cf%son1
                endif
            elseif(c%fSplitType==2)then
                if(c%Location==1)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==4)then
                             NeighborX1=>cfn%son2
                        elseif(cfn%son1%fSplitType==2)then
                             NeighborX1=>cfn%son1
                        else
                             NeighborX1=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2)then
                             NeighborX1=>cfn%son2
                        elseif(cfn%son1%fSplitType==4)then
                             NeighborX1=>cfn%son3
                        else
                             NeighborX1=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitType==3)then
                if(c%Location==1)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==5)then
                            NeighborX1=>cfn%son2
                        elseif(cfn%son1%fSplitType==3)then
                            NeighborX1=>cfn%son1
                        else
                            NeighborX1=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3)then
                            NeighborX1=>cfn%son2
                        elseif(cfn%son1%fSplitType==5)then
                            NeighborX1=>cfn%son3
                        else
                            NeighborX1=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitType==4)then
                if(c%Location==1)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==4)then
                            NeighborX1=>cfn%son2
                        elseif(cfn%son1%fSplitType==2)then
                            NeighborX1=>cfn%son1
                        else
                            NeighborX1=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    NeighborX1=>cf%son1
                elseif(c%Location==3)then
                    NeighborX1=>cf%son4
                elseif(c%Location==4)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2)then
                            NeighborX1=>cfn%son2
                        elseif(cfn%son1%fSplitType==4)then
                            NeighborX1=>cfn%son3
                        else
                            NeighborX1=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitType==5)then
                if(c%Location==1)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==5)then
                            NeighborX1=>cfn%son2
                        elseif(cfn%son1%fSplitType==3)then
                            NeighborX1=>cfn%son1
                        else
                            NeighborX1=>cfn
                        endif
                    endif
                elseif(c%Location==2)then
                    NeighborX1=>cf%son1
                elseif(c%Location==3)then
                    NeighborX1=>cf%son4
                elseif(c%Location==4)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3)then
                            NeighborX1=>cfn%son2
                        elseif(cfn%son1%fSplitType==5)then
                            NeighborX1=>cfn%son3
                        else
                            NeighborX1=>cfn
                        endif
                    endif
                endif
            elseif(c%fSplitType==6)then
                if(c%Location==1)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==6)then
                            NeighborX1=>cfn%son1
                        else
                            NeighborX1=>cfn%son2                !cfn%son1%fSplitType==0  1  4  5  
                        endif
                    endif
                elseif(c%Location==2)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
                            NeighborX1=>cfn%son3
                        elseif(cfn%son1%fSplitType==3)then
                            NeighborX1=>cfn%son1
                        else                              !cfn%son1%fSplitType==1  2  5  6
                            NeighborX1=>cfn%son2               
                        endif
                    endif
                elseif(c%Location==3)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==0)then
                            NeighborX1=>cfn%son7
                        elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
                            NeighborX1=>cfn%son2
                        elseif(cfn%son1%fSplitType==4.or.cfn%son1%fSplitType==5.or.cfn%son1%fSplitType==6)then     
                            NeighborX1=>cfn%son3               
                        endif
                    endif
                elseif(c%Location==4)then
                    if(.not.associated(cfn%son1))then
                        NeighborX1=>cfn
                    elseif(associated(cfn%son1))then
                        if(cfn%son1%fSplitType==0)then
                            NeighborX1=>cfn%son6
                        elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==4)then
                            NeighborX1=>cfn%son2
                        elseif(cfn%son1%fSplitType==2)then     
                            NeighborX1=>cfn%son1
                        elseif(cfn%son1%fSplitType==5)then     
                            NeighborX1=>cfn%son3 
                        elseif(cfn%son1%fSplitType==6)then     
                            NeighborX1=>cfn%son4     
                        endif
                    endif
                endif
            endif
        endif
        return
    end function NeighborX1
    
! ! X+ neighbor 
!     recursive function NeighborX2(c)
!         implicit none
!         type(OctCell),pointer :: c,cf,cinitial,cfn,NeighborX2 
!         !cf,c'father;cinitial,c'initial X+ neighbor;cfn,c'father's X+ neighbor
        
!         cinitial=>InitialNeighborX2(c)
!         NeighborX2=>null()
!         cf=>c%father
        
!         if(.not.associated(cf))then
!             NeighborX2=>cinitial
!             return
!         endif
        
!         cfn=>NeighborX2(cf)
!         if(associated(cfn))then
!             if(c%fSplitType==0)then
!                 if(c%Location==1)then
!                     NeighborX2=>cf%son2
!                 elseif(c%Location==2)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         NeighborX2=>cfn%son1
!                     endif
!                 elseif(c%Location==3)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
!                             NeighborX2=>cfn%son4
!                         elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==5)then
!                             NeighborX2=>cfn%son1
!                         elseif(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==6)then
!                             NeighborX2=>cfn%son2
!                         endif
!                     endif
!                 elseif(c%Location==4)then
!                     NeighborX2=>cf%son3
!                 elseif(c%Location==5)then    
!                     NeighborX2=>cf%son6
!                 elseif(c%Location==6)then  
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==4)then
!                             NeighborX2=>cfn%son1
!                         elseif(cfn%son1%fSplitType==3)then
!                             NeighborX2=>cfn%son2
!                         elseif(cfn%son1%fSplitType==0)then
!                             NeighborX2=>cfn%son5
!                         elseif(cfn%son1%fSplitType==5.or.cfn%son1%fSplitType==6)then
!                             NeighborX2=>cfn%son4
!                         endif
!                     endif
!                 elseif(c%Location==7)then    
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==1)then
!                             NeighborX2=>cfn%son1
!                         elseif(cfn%son1%fSplitType==0)then
!                             NeighborX2=>cfn%son8
!                         elseif(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                             NeighborX2=>cfn%son2
!                         elseif(cfn%son1%fSplitType==4.or.cfn%son1%fSplitType==5)then
!                             NeighborX2=>cfn%son4
!                         elseif(cfn%son1%fSplitType==6)then
!                             NeighborX2=>cfn%son3    
!                         endif
!                     endif
!                 elseif(c%Location==8)then
!                     NeighborX2=>cf%son7
!                 endif                           
!             elseif(c%fSplitType==1)then
!                 if(c%Location==1)then
!                     NeighborX2=>cf%son2
!                 elseif(c%Location==2)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==1)then
!                             NeighborX2=>cfn%son1
!                         else
!                             NeighborX2=>cfn
!                         endif           
!                     endif
!                 endif 
!             elseif(c%fSplitType==2)then
!                 if(c%Location==1)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==4)then
!                             NeighborX2=>cfn%son1
!                         else
!                             NeighborX2=>cfn
!                         endif
!                     endif        
!                 elseif(c%Location==2)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==2)then
!                             NeighborX2=>cfn%son2
!                         elseif(cfn%son1%fSplitType==4)then
!                             NeighborX2=>cfn%son4
!                         elseif(cfn%son1%fSplitType==1)then
!                             NeighborX2=>cfn%son1
!                         else
!                             NeighborX2=>cfn
!                         endif
!                     endif     
!                 endif 
!             elseif(c%fSplitType==3)then
!                 if(c%Location==1)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==5)then ! ==1 need to be noticed
!                             NeighborX2=>cfn%son1
!                         else
!                             NeighborX2=>cfn
!                         endif
!                     endif
!                 elseif(c%Location==2)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==1)then ! ==1 need to be noticed
!                             NeighborX2=>cfn%son1
!                         elseif(cfn%son1%fSplitType==3)then
!                             NeighborX2=>cfn%son2
!                         elseif(cfn%son1%fSplitType==5)then
!                             NeighborX2=>cfn%son4
!                         else
!                             NeighborX2=>cfn
!                         endif
!                     endif
!                 endif 
!             elseif(c%fSplitType==4)then
!                 if(c%Location==1)then
!                     NeighborX2=>cf%son2
!                 elseif(c%Location==2)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==4)then 
!                             NeighborX2=>cfn%son1
!                         else
!                             NeighborX2=>cfn
!                         endif
!                     endif
!                 elseif(c%Location==3)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==1)then 
!                             NeighborX2=>cfn%son1
!                         elseif(cfn%son1%fSplitType==2)then
!                             NeighborX2=>cfn%son2
!                         elseif(cfn%son1%fSplitType==4)then
!                             NeighborX2=>cfn%son4
!                         else
!                             NeighborX2=>cfn
!                         endif
!                     endif
!                 elseif(c%Location==4)then    
!                     NeighborX2=>cf%son3 
!                 endif 
!             elseif(c%fSplitType==5)then
!                 if(c%Location==1)then
!                     NeighborX2=>cf%son2 
!                 elseif(c%Location==2)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==5)then 
!                             NeighborX2=>cfn%son1
!                         else
!                             NeighborX2=>cfn
!                         endif
!                     endif
!                 elseif(c%Location==3)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==1)then 
!                             NeighborX2=>cfn%son1
!                         elseif(cfn%son1%fSplitType==3)then
!                             NeighborX2=>cfn%son2
!                         elseif(cfn%son1%fSplitType==5)then
!                             NeighborX2=>cfn%son4
!                         else
!                             NeighborX2=>cfn
!                         endif
!                     endif
!                 elseif(c%Location==4)then    
!                     NeighborX2=>cf%son3
!                 endif 
!             elseif(c%fSplitType==6)then
!                 if(c%Location==1)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         NeighborX2=>cfn%son1                    !may be bug   
!                     endif
!                 elseif(c%Location==2)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==6)then
!                             NeighborX2=>cfn%son2
!                         elseif(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
!                             NeighborX2=>cfn%son4
!                         else
!                             NeighborX2=>cfn%son1
!                         endif                                        
!                     endif
!                 elseif(c%Location==3)then
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==0)then
!                             NeighborX2=>cfn%son8
!                         elseif(cfn%son1%fSplitType==1)then
!                             NeighborX2=>cfn%son1
!                         elseif(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                             NeighborX2=>cfn%son2
!                         elseif(cfn%son1%fSplitType==4.or.cfn%son1%fSplitType==5)then
!                             NeighborX2=>cfn%son4 
!                         elseif(cfn%son1%fSplitType==6)then
!                             NeighborX2=>cfn%son3
!                         endif                                        
!                     endif
!                 elseif(c%Location==4)then    
!                     if(.not.associated(cfn%son1))then
!                         NeighborX2=>cfn
!                     elseif(associated(cfn%son1))then
!                         if(cfn%son1%fSplitType==0)then
!                             NeighborX2=>cfn%son5
!                         elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==4)then
!                             NeighborX2=>cfn%son1
!                         elseif(cfn%son1%fSplitType==3)then
!                             NeighborX2=>cfn%son2
!                         elseif(cfn%son1%fSplitType==5.or.cfn%son1%fSplitType==6)then
!                             NeighborX2=>cfn%son4 
!                         endif                                        
!                     endif
!                 endif 
!             endif
!         endif
!         return
!     end function NeighborX2   
    
! ! Y- neighbor
!     recursive function NeighborY1(c)
!         implicit none
!             type(OctCell),pointer :: c,cf,cinitial,cfn,NeighborY1 
!             !cf,c'father;cinitial,c'initial Y- neighbor;cfn,c'father's Y- neighbor
        
!             cinitial=>InitialNeighborY1(c)
!             NeighborY1=>null()
!             cf=>c%father
        
!             if(.not.associated(cf))then
!                 NeighborY1=>cinitial
!                 return
!             endif
        
!             cfn=>NeighborY1(cf)
!             if(associated(cfn))then
!                 if(c%fSplitType==0)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==5)then
!                                 NeighborY1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
!                                 NeighborY1=>cfn%son4
!                             elseif(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==6)then
!                                 NeighborY1=>cfn%son2
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==5.or.cfn%son1%fSplitType==6)then
!                                 NeighborY1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
!                                 NeighborY1=>cfn%son3
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborY1=>cfn%son1
!                             endif
!                         endif
!                     elseif(c%Location==3)then
!                         NeighborY1=>cf%son2
!                     elseif(c%Location==4)then
!                         NeighborY1=>cf%son1
!                     elseif(c%Location==5)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborY1=>cfn%son8
!                             elseif(cfn%son1%fSplitType==1)then
!                                 NeighborY1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                                 NeighborY1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==4.or.cfn%son1%fSplitType==5)then
!                                 NeighborY1=>cfn%son4
!                             elseif(cfn%son1%fSplitType==6)then
!                                 NeighborY1=>cfn%son3    
!                             endif
!                         endif
!                     elseif(c%Location==6)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborY1=>cfn%son7
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                                 NeighborY1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==4.or.cfn%son1%fSplitType==5.or.cfn%son1%fSplitType==6)then
!                                 NeighborY1=>cfn%son3
!                             endif
!                         endif
!                     elseif(c%Location==7)then
!                         NeighborY1=>cf%son6
!                     elseif(c%Location==8)then
!                         NeighborY1=>cf%son5
!                     endif
!                 elseif(c%fSplitType==1)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1)then
!                                 NeighborY1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==4)then
!                                 NeighborY1=>cfn%son4
!                             elseif(cfn%son1%fSplitType==2)then
!                                 NeighborY1=>cfn%son2    
!                             else
!                                 NeighborY1=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2)then
!                                 NeighborY1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==4)then
!                                 NeighborY1=>cfn%son3    
!                             else
!                                 NeighborY1=>cfn
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==2)then 
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2)then
!                                 NeighborY1=>cfn%son2
!                             else
!                                 NeighborY1=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         NeighborY1=>cf%son1
!                     endif
!                 elseif(c%fSplitType==3)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2)then
!                                 NeighborY1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborY1=>cfn%son1
!                             else
!                                 NeighborY1=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                                 NeighborY1=>cfn%son2
!                             else
!                                 NeighborY1=>cfn
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==4)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1)then
!                                 NeighborY1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==2)then
!                                 NeighborY1=>cfn%son2    
!                             elseif(cfn%son1%fSplitType==4)then
!                                 NeighborY1=>cfn%son4
!                             else
!                                 NeighborY1=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2)then
!                                 NeighborY1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==4)then
!                                 NeighborY1=>cfn%son3
!                             else
!                                 NeighborY1=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==3)then
!                         NeighborY1=>cf%son2
!                     elseif(c%Location==4)then
!                         NeighborY1=>cf%son1
!                     endif
!                 elseif(c%fSplitType==5)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
!                                 NeighborY1=>cfn%son4
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==5)then
!                                 NeighborY1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==6)then
!                                 NeighborY1=>cfn%son2
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
!                                 NeighborY1=>cfn%son3
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborY1=>cfn%son1
!                             else                         !  1  2  5  6
!                                 NeighborY1=>cfn%son2
!                             endif
!                         endif
!                     elseif(c%Location==3)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborY1=>cfn%son7
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                                 NeighborY1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==4.or.cfn%son1%fSplitType==5.or.cfn%son1%fSplitType==6)then         
!                                 NeighborY1=>cfn%son3
!                             endif
!                         endif
!                     elseif(c%Location==4)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborY1=>cfn%son8
!                             elseif(cfn%son1%fSplitType==1)then
!                                 NeighborY1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then         
!                                 NeighborY1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==4.or.cfn%son1%fSplitType==5)then         
!                                 NeighborY1=>cfn%son4
!                             elseif(cfn%son1%fSplitType==6)then         
!                                 NeighborY1=>cfn%son3
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==6)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==6)then
!                                 NeighborY1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborY1=>cfn%son1
!                             else         
!                                 NeighborY1=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         NeighborY1=>cf%son1
!                     elseif(c%Location==3)then
!                         NeighborY1=>cf%son4
!                     elseif(c%Location==4)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                                 NeighborY1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==6)then         
!                                 NeighborY1=>cfn%son3
!                             else         
!                                 NeighborY1=>cfn
!                             endif
!                         endif
!                     endif
!                 endif
!             endif
!             return
!     end function NeighborY1
    
! ! Y+ neighbor
!     recursive function NeighborY2(c)
!         implicit none
!             type(OctCell),pointer :: c,cf,cinitial,cfn,NeighborY2 
!             !cf,c'father;cinitial,c'initial Y+ neighbor;cfn,c'father's Y+ neighbor
        
!             cinitial=>InitialNeighborY2(c)
!             NeighborY2=>null()
!             cf=>c%father
        
!             if(.not.associated(cf))then
!                 NeighborY2=>cinitial
!                 return
!             endif
        
!             cfn=>NeighborY2(cf)
!             if(associated(cfn))then
!                 if(c%fSplitType==0)then
!                     if(c%Location==1)then
!                         NeighborY2=>cf%son4
!                     elseif(c%Location==2)then
!                         NeighborY2=>cf%son3
!                     elseif(c%Location==3)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==4.or.cfn%son1%fSplitType==5.or.cfn%son1%fSplitType==0)then
!                                 NeighborY2=>cfn%son2
!                             elseif(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==6)then         
!                                 NeighborY2=>cfn%son1
!                             endif
!                         endif
!                     elseif(c%Location==4)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             NeighborY2=>cfn%son1
!                         endif
!                     elseif(c%Location==5)then
!                         NeighborY2=>cf%son8
!                     elseif(c%Location==6)then
!                         NeighborY2=>cf%son7
!                     elseif(c%Location==7)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborY2=>cfn%son6
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==4)then         
!                                 NeighborY2=>cfn%son2
!                             elseif(cfn%son1%fSplitType==2)then         
!                                 NeighborY2=>cfn%son1
!                             elseif(cfn%son1%fSplitType==5)then         
!                                 NeighborY2=>cfn%son3
!                             elseif(cfn%son1%fSplitType==6)then         
!                                 NeighborY2=>cfn%son4
!                             endif
!                         endif
!                     elseif(c%Location==8)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborY2=>cfn%son5
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==4)then         
!                                 NeighborY2=>cfn%son1
!                             elseif(cfn%son1%fSplitType==3)then         
!                                 NeighborY2=>cfn%son2
!                             elseif(cfn%son1%fSplitType==5.or.cfn%son1%fSplitType==6)then         
!                                 NeighborY2=>cfn%son4
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==1)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==4)then         
!                                 NeighborY2=>cfn%son1
!                             else         
!                                 NeighborY2=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==4)then         
!                                 NeighborY2=>cfn%son2
!                             elseif(cfn%son1%fSplitType==2)then         
!                                 NeighborY2=>cfn%son1
!                             else
!                                 NeighborY2=>cfn
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==2)then 
!                     if(c%Location==1)then
!                         NeighborY2=>cf%son2
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2)then         
!                                 NeighborY2=>cfn%son1
!                             else
!                                 NeighborY2=>cfn
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==3)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==6)then         
!                                 NeighborY2=>cfn%son1
!                             else
!                                 NeighborY2=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2)then         
!                                 NeighborY2=>cfn%son1
!                             elseif(cfn%son1%fSplitType==3)then         
!                                 NeighborY2=>cfn%son2
!                             elseif(cfn%son1%fSplitType==6)then         
!                                 NeighborY2=>cfn%son4
!                             else
!                                 NeighborY2=>cfn
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==4)then
!                     if(c%Location==1)then
!                         NeighborY2=>cf%son4
!                     elseif(c%Location==2)then
!                         NeighborY2=>cf%son3
!                     elseif(c%Location==3)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==4)then         
!                                 NeighborY2=>cfn%son2
!                             elseif(cfn%son1%fSplitType==2)then         
!                                 NeighborY2=>cfn%son1
!                             else
!                                 NeighborY2=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==4)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==4)then         
!                                 NeighborY2=>cfn%son1
!                             else
!                                 NeighborY2=>cfn
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==5)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             NeighborY2=>cfn%son1
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==6)then         
!                                 NeighborY2=>cfn%son2
!                             else      
!                                 NeighborY2=>cfn%son1
!                             endif
!                         endif
!                     elseif(c%Location==3)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then         
!                                 NeighborY2=>cfn%son6
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==4)then       
!                                 NeighborY2=>cfn%son2
!                             elseif(cfn%son1%fSplitType==2)then       
!                                 NeighborY2=>cfn%son1 
!                             elseif(cfn%son1%fSplitType==5)then       
!                                 NeighborY2=>cfn%son3
!                             elseif(cfn%son1%fSplitType==6)then       
!                                 NeighborY2=>cfn%son4
!                             endif
!                         endif
!                     elseif(c%Location==4)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then         
!                                 NeighborY2=>cfn%son5
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==4)then       
!                                 NeighborY2=>cfn%son1
!                             elseif(cfn%son1%fSplitType==3)then       
!                                 NeighborY2=>cfn%son2 
!                             elseif(cfn%son1%fSplitType==5.or.cfn%son1%fSplitType==6)then       
!                                 NeighborY2=>cfn%son4
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==6)then
!                     if(c%Location==1)then
!                         NeighborY2=>cf%son2
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==6)then           
!                                 NeighborY2=>cfn%son1
!                             else       
!                                 NeighborY2=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==3)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborY2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2)then           
!                                 NeighborY2=>cfn%son1
!                             elseif(cfn%son1%fSplitType==3)then           
!                                 NeighborY2=>cfn%son2 
!                             elseif(cfn%son1%fSplitType==6)then           
!                                 NeighborY2=>cfn%son4
!                             else
!                                 NeighborY2=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==4)then
!                         NeighborY2=>cf%son3
!                     endif
!                 endif
!             endif
!             return
!     end function NeighborY2    
 
! ! Z- neighbor
!     recursive function NeighborZ1(c)
!         implicit none
!             type(OctCell),pointer :: c,cf,cinitial,cfn,NeighborZ1 
!             !cf,c'father;cinitial,c'initial Z- neighbor;cfn,c'father's Z- neighbor
        
!             cinitial=>InitialNeighborZ1(c)
!             NeighborZ1=>null()
!             cf=>c%father
        
!             if(.not.associated(cf))then
!                 NeighborZ1=>cinitial
!                 return
!             endif
        
!             cfn=>NeighborZ1(cf)
!             if(associated(cfn))then
!                 if(c%fSplitType==0)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborZ1=>cfn%son5
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==4)then
!                                 NeighborZ1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==5.or.cfn%son1%fSplitType==6)then
!                                 NeighborZ1=>cfn%son4
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborZ1=>cfn%son6
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==4)then
!                                 NeighborZ1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==2)then
!                                 NeighborZ1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==5)then
!                                 NeighborZ1=>cfn%son3
!                             elseif(cfn%son1%fSplitType==6)then
!                                 NeighborZ1=>cfn%son4   
!                             endif
!                         endif
!                     elseif(c%Location==3)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborZ1=>cfn%son7
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             else                        !4   5   6
!                                 NeighborZ1=>cfn%son3   
!                             endif
!                         endif
!                     elseif(c%Location==4)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborZ1=>cfn%son8
!                             elseif(cfn%son1%fSplitType==1)then
!                                 NeighborZ1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==4.or.cfn%son1%fSplitType==5)then
!                                 NeighborZ1=>cfn%son4
!                             elseif(cfn%son1%fSplitType==6)then
!                                 NeighborZ1=>cfn%son3  
!                             endif
!                         endif
!                     elseif(c%Location==5)then
!                         NeighborZ1=>cf%son1  
!                     elseif(c%Location==6)then
!                         NeighborZ1=>cf%son2  
!                     elseif(c%Location==7)then
!                         NeighborZ1=>cf%son3  
!                     elseif(c%Location==8)then
!                         NeighborZ1=>cf%son4  
!                     endif
!                 elseif(c%fSplitType==1)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1)then
!                                 NeighborZ1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==5)then
!                                 NeighborZ1=>cfn%son4
!                             else
!                                 NeighborZ1=>cfn  
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==5)then
!                                 NeighborZ1=>cfn%son3
!                             else
!                                 NeighborZ1=>cfn  
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==2)then 
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2)then
!                                 NeighborZ1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==6)then
!                                 NeighborZ1=>cfn%son4
!                             else    
!                                 NeighborZ1=>cfn  
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==6)then
!                                 NeighborZ1=>cfn%son3
!                             else    
!                                 NeighborZ1=>cfn  
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==3)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             else
!                                 NeighborZ1=>cfn  
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         NeighborZ1=>cf%son1
!                     endif
!                 elseif(c%fSplitType==4)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborZ1=>cfn%son5
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==4)then
!                                 NeighborZ1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2 
!                             elseif(cfn%son1%fSplitType==5.or.cfn%son1%fSplitType==6)then
!                                 NeighborZ1=>cfn%son4
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborZ1=>cfn%son6
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==4)then
!                                 NeighborZ1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==2)then
!                                 NeighborZ1=>cfn%son1 
!                             elseif(cfn%son1%fSplitType==5)then
!                                 NeighborZ1=>cfn%son3
!                             elseif(cfn%son1%fSplitType==6)then
!                                 NeighborZ1=>cfn%son4
!                             endif
!                         endif
!                     elseif(c%Location==3)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborZ1=>cfn%son7
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             else
!                                 NeighborZ1=>cfn%son3
!                             endif
!                         endif
!                     elseif(c%Location==4)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0)then
!                                 NeighborZ1=>cfn%son8
!                             elseif(cfn%son1%fSplitType==1)then
!                                 NeighborZ1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2 
!                             elseif(cfn%son1%fSplitType==4.or.cfn%son1%fSplitType==5)then
!                                 NeighborZ1=>cfn%son4
!                             elseif(cfn%son1%fSplitType==6)then
!                                 NeighborZ1=>cfn%son3
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==5)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1)then
!                                 NeighborZ1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2 
!                             elseif(cfn%son1%fSplitType==5)then
!                                 NeighborZ1=>cfn%son4
!                             else
!                                 NeighborZ1=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==5)then
!                                 NeighborZ1=>cfn%son3
!                             else
!                                 NeighborZ1=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==3)then
!                         NeighborZ1=>cf%son2
!                     elseif(c%Location==4)then
!                         NeighborZ1=>cf%son1
!                     endif
!                 elseif(c%fSplitType==6)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2)then
!                                 NeighborZ1=>cfn%son1
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==6)then
!                                 NeighborZ1=>cfn%son4
!                             else
!                                 NeighborZ1=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ1=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3)then
!                                 NeighborZ1=>cfn%son2
!                             elseif(cfn%son1%fSplitType==6)then
!                                 NeighborZ1=>cfn%son3
!                             else
!                                 NeighborZ1=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==3)then
!                         NeighborZ1=>cf%son2
!                     elseif(c%Location==4)then
!                         NeighborZ1=>cf%son1
!                     endif
!                 endif
!             endif
!             return
!     end function NeighborZ1        

! ! Z+ neighbor
!     recursive function NeighborZ2(c)
!         implicit none
!             type(OctCell),pointer :: c,cf,cinitial,cfn,NeighborZ2 
!             !cf,c'father;cinitial,c'initial Z+ neighbor;cfn,c'father's Z+ neighbor
        
!             cinitial=>InitialNeighborZ2(c)
!             NeighborZ2=>null()
!             cf=>c%father
        
!             if(.not.associated(cf))then
!                 NeighborZ2=>cinitial
!                 return
!             endif
        
!             cfn=>NeighborZ2(cf)
!             if(associated(cfn))then
!                 if(c%fSplitType==0)then
!                     if(c%Location==1)then
!                         NeighborZ2=>cf%son5
!                     elseif(c%Location==2)then
!                         NeighborZ2=>cf%son6
!                     elseif(c%Location==3)then
!                         NeighborZ2=>cf%son7
!                     elseif(c%Location==4)then
!                         NeighborZ2=>cf%son8
!                     elseif(c%Location==5)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             NeighborZ2=>cfn%son1
!                         endif
!                     elseif(c%Location==6)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==6)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn%son2
!                             endif
!                         endif
!                     elseif(c%Location==7)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
!                                 NeighborZ2=>cfn%son3
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn%son2
!                             endif
!                         endif
!                     elseif(c%Location==8)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
!                                 NeighborZ2=>cfn%son4
!                             elseif(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==5)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn%son2
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==1)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==5)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==5)then
!                                 NeighborZ2=>cfn%son2
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==2)then 
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==6)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==6)then
!                                 NeighborZ2=>cfn%son2
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==3)then
!                     if(c%Location==1)then
!                         NeighborZ2=>cf%son2
!                     elseif(c%Location==2)then
!                         NeighborZ2=>cfn
!                     endif
!                 elseif(c%fSplitType==4)then
!                     if(c%Location==1)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             NeighborZ2=>cfn%son1
!                         endif
!                     elseif(c%Location==2)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==6)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn%son2
!                             endif
!                         endif
!                     elseif(c%Location==3)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
!                                 NeighborZ2=>cfn%son3
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn%son2
!                             endif
!                         endif
!                     elseif(c%Location==4)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==0.or.cfn%son1%fSplitType==4)then
!                                 NeighborZ2=>cfn%son4
!                             elseif(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==6)then
!                                 NeighborZ2=>cfn%son2
!                             else
!                                 NeighborZ2=>cfn%son1
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==5)then
!                     if(c%Location==1)then
!                         NeighborZ2=>cf%son4
!                     elseif(c%Location==2)then
!                         NeighborZ2=>cf%son3
!                     elseif(c%Location==3)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==5)then
!                                 NeighborZ2=>cfn%son2
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn
!                             endif
!                         endif
!                     elseif(c%Location==4)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==1.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==5)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn
!                             endif
!                         endif
!                     endif
!                 elseif(c%fSplitType==6)then
!                     if(c%Location==1)then
!                         NeighborZ2=>cf%son4
!                     elseif(c%Location==2)then
!                         NeighborZ2=>cf%son3
!                     elseif(c%Location==3)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==6)then
!                                 NeighborZ2=>cfn%son2
!                             elseif(cfn%son1%fSplitType==3)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn   
!                             endif
!                         endif
!                     elseif(c%Location==4)then
!                         if(.not.associated(cfn%son1))then
!                             NeighborZ2=>cfn
!                         elseif(associated(cfn%son1))then
!                             if(cfn%son1%fSplitType==2.or.cfn%son1%fSplitType==3.or.cfn%son1%fSplitType==6)then
!                                 NeighborZ2=>cfn%son1
!                             else
!                                 NeighborZ2=>cfn   
!                             endif
!                         endif
!                     endif
!                 endif
!             endif
!             return
!     end function NeighborZ2        

    end module ModNeighbor
    
    