program sinctwodim
    implicit none

    integer:: xnodes,ynodes, memstatus, io, iostatus,i,eigen
    double precision:: r0xx, r0yy, rminx, rminy, rmaxx, rmaxy, Fxx, Fxy, Fyy, mass, Lx,Ly
    double precision, allocatable:: Xmatrix(:), Ymatrix(:), Kmatrix(:), Vmatrix(:),&
    Tmatrix(:,:), Hmatrix(:,:), eigenvalues(:), pointmatrix(:,:), Umatrix(:,:)

    open(newunit=io, file="input.txt", status="old", action="read", iostat = iostatus)
    if(iostatus /=0) then
        stop "***File open error***"
    end if
    read(io,*)
    read(io,*) mass, xnodes, ynodes
    read(io,*)
    read(io,*) rminx,rmaxx
    read(io,*)
    read(io,*) rminy,rmaxy
    read(io,*)
    read(io,*) r0xx,r0yy
    read(io,*)
    read(io,*) Fxx,Fyy,Fxy
    read(io,*)
    read(io,*) eigen
    close(io)

    allocate(Xmatrix(xnodes),Ymatrix(ynodes),Kmatrix(xnodes*ynodes),Vmatrix(xnodes*ynodes),&
    Tmatrix(xnodes*ynodes,xnodes*ynodes), Hmatrix(xnodes*ynodes,xnodes*ynodes),&
    eigenvalues(xnodes*ynodes), pointmatrix(3,(5*xnodes+1)*(5*ynodes+1)), Umatrix(3,xnodes*ynodes), STAT = memstatus)
    if(memstatus /= 0) then
    stop "*** Not enough memory for allocation***"
    end if
    Lx = rmaxx - rminx
    Ly= rmaxy - rminy

    call discretegrid(rminx,rmaxx,xnodes,Xmatrix)
    call discretegrid(rminy,rmaxy,ynodes,Ymatrix)
    call twodimkinetic(Kmatrix,xnodes,ynodes,mass,Lx,Ly)
    call twodimpotential(Xmatrix, Ymatrix, Vmatrix, xnodes, ynodes,&
        r0xx, r0yy, Fxx*0.5, Fxy, Fyy*0.5)
    call transformmatrix(xnodes, ynodes, Tmatrix)
    call FBRHamiltonian(Hmatrix,xnodes,ynodes,Vmatrix,Kmatrix,Tmatrix)
    call DIAG2(xnodes*ynodes,xnodes*ynodes,eigenvalues,Hmatrix)
    open(newunit=io, file="output.txt")
    if(iostatus /=0) then
        stop "***Output file error***"
    end if
    write(io,*) "Computed first 20 eigenvalues of Hamiltonian", "Nodes", xnodes*ynodes
    write(io,*) "n", "energy, hartree "
    do i=1,min(xnodes*ynodes,20)
        write(io,*) i, eigenvalues(i)
        print*, i, eigenvalues(i)
    end do
    write(io,*) "Hamiltonian eigenvectors"
    do i=1,min(xnodes*ynodes,20)
        write(io,*) "Stationary State No. ", i, "Energy", eigenvalues(i)
        write(io,*) Hmatrix(:,i)
    end do

    !Print the scatter plot of nth stationary state.
    close(io)
    print*, "Output saved."
    open(newunit=io, file="plots.txt")
    if(iostatus /=0) then
        stop "***Output file error***"
    end if
    write(io,*) "Wavefunction of state n=",eigen
    call FunctionPlot3D(Hmatrix,xnodes,ynodes,pointmatrix,eigen,rminx,&
    rminy,rmaxx,rmaxy)
        do i=1,(5*xnodes+1)*(5*ynodes+1)
            write(io,*) pointmatrix(:,i)
        end do
    print*,"Plot saved"
    close(io)

    print*, "Output saved."
    open(newunit=io, file="Uplots.txt")
    if(iostatus /=0) then
        stop "***Output file error***"
    end if
    write(io,*) "U matrix plot"
    call Uplot3D(Vmatrix,xnodes,ynodes,Umatrix,rminx,rminy,rmaxx,rmaxy)
        do i=1,xnodes*ynodes
            write(io,*) Umatrix(:,i)
        end do
    print*,"Plot saved"
    close(io)

    deallocate(Xmatrix,Ymatrix,Kmatrix,Vmatrix,Tmatrix, Hmatrix, &
    eigenvalues, pointmatrix, Umatrix, STAT = memstatus)

    if(memstatus /= 0) then
    stop "*** Deallocation fail***"
    end if

end program



!generate a discrete grid with node points
subroutine discretegrid(r_min,r_max,nodes,xmatrix)
    integer:: nodes, i
    double precision:: r_min, r_max, xmatrix(nodes), step
    step = (r_max-r_min)/(nodes+1)
    do i=1,nodes
    xmatrix(i) = r_min + step*i
    end do
return
end subroutine

!Generate the 2D kinetic energy operator matrix for particle in a box repn
!Because U=0 inside of PIB T_ij = H0_ij for the basis set
subroutine twodimkinetic(Kmat,xnodes,ynodes,mass,Lx,Ly)
    integer i,j,k,l, xnodes, ynodes,a,temp
    double precision Kmat(xnodes*ynodes)
    double precision mass, xconst, Lx,Ly, yconst
    xconst =  (acos(-1.0d0)/Lx)**2/(2*mass)
    yconst = (acos(-1.0d0)/Ly)**2/(2*mass)
    !K_ab = a^2+b^2*pi^2/(2mL^2)
    do i=1,xnodes
        do j=1,ynodes
            Kmat((i-1)*ynodes+j) = xconst*i*i + yconst*j*j
        end do
    end do

return
end subroutine

!Generate 2D harmonic quadratic form as a 1D matrix
subroutine twodimpotential(Xmatrix, Ymatrix, Vmatrix, xnodes, ynodes,&
    r0xx, r0yy, Fxx, Fxy, Fyy)
    integer xnodes, ynodes, i,j
    double precision r0xx, r0yy, Fxx, Fyy, Fxy
    double precision Xmatrix(xnodes), Ymatrix(ynodes), Vmatrix(xnodes*ynodes)
    do i=1, xnodes
        do j=1, ynodes
            Vmatrix(ynodes*(i-1)+j) = Fxx*(Xmatrix(i)-r0xx)**2 +Fyy*(Ymatrix(j)-r0yy)**2+&
            Fxy*(Xmatrix(i)-r0xx)*(Ymatrix(j)-r0yy)
            !Vmatrix(ynodes*(i-1)+j) = 2.0d-2
            !Vmatrix(ynodes*(i-1)+j) = Fxx*cosh(sqrt((Xmatrix(i)-r0xx)**2+(Ymatrix(j)-r0yy)**2))

        end do
    end do
return
end subroutine

!Generate transform matrix from fbr to dvr
subroutine transformmatrix(xnodes,ynodes,tmatrix)
    integer i,j,k,l, xnodes, ynodes
    double precision tmatrix(xnodes*ynodes,xnodes*ynodes)
    double precision temp, xtempconst, xtempconst2, ytempconst, ytempconst2,a

    xtempconst = sqrt(2.0/(xnodes+1))
    xtempconst2= acos(-1.0d0)/(xnodes+1)
    ytempconst = sqrt(2.0/(ynodes+1))
    ytempconst2= acos(-1.0d0)/(ynodes+1)

!For 2D sinc basis, the transform product is a direct product:
! T_xy = T_x (X) T_y
!First compute the multiple copies of T_y matrix, then compute direct
!product.
    do i=1,ynodes
        do j=1,i
            temp = sin(i*j*ytempconst2)*ytempconst
            tmatrix(j,i) = temp
            tmatrix(i,j) = temp
        end do
    end do
    do i=1,ynodes*xnodes
        do j=1,xnodes*ynodes
            k = mod(i,ynodes)
            l = mod(j,ynodes)
            if (k==0) then
            k=ynodes
            end if
            if (l==0) then
            l=ynodes
            end if
            a = tmatrix(k,l)
            tmatrix(i,j)=a
        end do
    end do
    do i=1, xnodes
        do j=1, xnodes
            temp = sin(i*j*xtempconst2)*xtempconst
            do k=1, ynodes
                do l=1, ynodes
                    a = tmatrix((i-1)*ynodes+k,(j-1)*ynodes+l)
                    tmatrix((i-1)*ynodes+k,(j-1)*ynodes+l) = a*temp

                end do
            end do
        end do
    end do

return
end subroutine

!Generate Hamiltonian matrix for FBR representation
subroutine FBRHamiltonian(Hmatrix, xnodes, ynodes, Vmatrix,Kmatrix,Tmatrix)
    integer i,j, xnodes, ynodes
    double precision Hmatrix(xnodes*ynodes,xnodes*ynodes), Vmatrix(xnodes*ynodes),&
    Kmatrix(xnodes*ynodes),Tmatrix(xnodes*ynodes,xnodes*ynodes), temp
    do i=1, xnodes*ynodes
        Hmatrix(i,i) = Kmatrix(i)
        do j=1, i
            temp=0
            do k=1,xnodes*ynodes
                temp = temp + Tmatrix(k,i)*Vmatrix(k)*Tmatrix(j,k)
            end do
            Hmatrix(i,j) = temp + Hmatrix(i,j)
            Hmatrix(j,i) = temp + Hmatrix(j,i)
        end do
    end do
return
end subroutine

!Generate 2D wavefunction plot for 3D visualisation
subroutine FunctionPlot3D(Hmatrix,xnodes,ynodes,pointmatrix, eigen,&
    rminx,rminy,rmaxx,rmaxy)
    integer i,j, xnodes, ynodes, eigen,k,l
    double precision Hmatrix(xnodes*ynodes,xnodes*ynodes), pointmatrix(3,(5*xnodes+1)*(5*ynodes+1))
    double precision rminx,rminy,rmaxx,rmaxy, stepx, stepy, const,const1,const2, Lx,Ly, temp
    Lx = rmaxx-rminx
    Ly = rmaxy-rminy
    stepx = Lx/(5*xnodes)
    stepy = Ly/(5*ynodes)
    !phi_ij = 2/sqrt(LxLy) * aij * sin(ipix/Lx)* sin(jpiy/Ly)
    const = 2.0/sqrt(Lx*Ly)
    const1 = acos(-1.0)/Lx
    const2 = acos(-1.0)/Ly
    !Generate coordinate points
    do i=1,(5*xnodes+1)
        do j=1,(5*ynodes+1)
            pointmatrix(1,(5*ynodes+1)*(i-1)+j)= rminx+stepx*(i-1)
            pointmatrix(2,(5*ynodes+1)*(i-1)+j)= rminy+stepy*(j-1)
        end do
    end do
    !Compute function value
    do i=1,(5*xnodes+1)
        do j=1,(5*ynodes+1)
        temp = 0.0
            do k=1,xnodes
                do l=1,ynodes
                temp = temp + const* Hmatrix(xnodes*(k-1)+l,eigen)* &
                sin(const1*pointmatrix(1,(5*ynodes+1)*(i-1)+j)*k)*&
                sin(const2*pointmatrix(2,(5*ynodes+1)*(i-1)+j)*l)
                end do
            end do
        pointmatrix(3,(5*ynodes+1)*(i-1)+j)= temp
        end do
    end do

return
end subroutine

subroutine Uplot3D(Vmat,xnodes,ynodes,Upointmatrix,&
    rminx,rminy,rmaxx,rmaxy)
    integer i,j, xnodes, ynodes
    double precision Vmat(xnodes*ynodes), Upointmatrix(3,xnodes*ynodes)
    double precision rminx,rminy,rmaxx,rmaxy,xstep, ystep, Lx,Ly
    Lx = rmaxx-rminx
    Ly = rmaxy-rminy
    xstep = Lx/xnodes
    ystep = Ly/ynodes
    !Generate coordinate points
    do i=1,xnodes
        do j=1,ynodes
            Upointmatrix(1,ynodes*(i-1)+j)= rminx+xstep*(i-1)
            Upointmatrix(2,ynodes*(i-1)+j)= rminy+ystep*(j-1)
        end do
    end do
    !Compute function value
    do i=1,xnodes
        do j=1,ynodes
            Upointmatrix(3,ynodes*(i-1)+j)= Vmat(ynodes*(i-1)+j)
        end do
    end do
return
end subroutine
