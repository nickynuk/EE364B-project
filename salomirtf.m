function [ H ] = salomirtf( dimension,dx,dy,dz,orientation )
%creates the salomir transfer function in k-space given a volume size and
%resolution.  Volume size is the number of voxels in a direction, while the 
%resolution is the dimension of the voxels in mm. scanner coordinate system
%used
%Dimension is the dimensions of the (desired image) in terms of [rows,
%columns, slices, echos]
%%orientation refers to patient orientation with respect to the magnetic field orientation. 
%choices are 'axial', 'coronal', 'sagittal', or and arbitrary orientation with respect to the main field [theta, phi] (spherical
%coordinates)
%it is a string that specifies the orientation of the patient's body
%with respect to the magnetic field.  Default is axial, choices are axial,
%sagittal and coronal

    %check for odd dimensions;
    
if mod(dimension(3),2)==1
    dimension(3)=dimension(3)-1;
end
%y direction
if mod(dimension(1),2)==1
    dimension(1)=dimension(1)-1;
end
%x direction
if mod(dimension(2),2)==1
    dimension(2)=dimension(2)-1;
end

if strcmp(orientation,'sagittal')==1
    Nz=dimension(1);
    Ny=dimension(2);
    Nx=dimension(3);
    %calculate sample width in k space
    dKx=2*pi/Nx/dx;
    dKy=2*pi/Ny/dy;
    dKz=2*pi/Nz/dz;
    %-----------------------------------------------------------------------------
    %construct Salomir Filter Transfer Function
    Kx=[0:1:(floor(Nx/2)-1) -(floor(Nx/2)):1:-1]*dKx;
    temp=ones(dimension(1:3));
    for n=1:Nx
        temp(:,:,n)=Kx(n);
    end
    Kx=temp;
    Ky=[0:1:(floor(Ny/2)-1) -(floor(Ny/2)):1:-1]*dKy;
    Ky=Ky';
    for n=1:Ny
        temp(n,:,:)=Ky(n);
    end
    Ky=temp;
    Kz=[0:1:(floor(Nz/2)-1) -(floor(Nz/2)):1:-1]*dKz;
    for n=1:Nz
        temp(:,n,:)=Kz(n);
    end
    Kz=temp;
    H=ones(dimension);
    A=Kz.^2;
    B=(Kx.^2+Ky.^2+Kz.^2);
    B(find(B==0))=1*10^-10;
    H=1/3-A./B;
    
elseif strcmp(orientation,'coronal')==1
    Nz=dimension(1);
    Nx=dimension(2);
    Ny=dimension(3);
    %calculate sample width in k space
    dKx=2*pi/Nx/dx;
    dKy=2*pi/Ny/dy;
    dKz=2*pi/Nz/dz;

    %-----------------------------------------------------------------------------
    %construct Salomir Filter Transfer Function
    Kx=[0:1:(floor(Nx/2)-1) -(floor(Nx/2)):1:-1]*dKx;
    temp=ones(dimension(1:3));
    for n=1:dimension(2)
        temp(:,n,:)=Kx(n);
    end
    Kx=temp;
    Ky=[0:1:(floor(Ny/2)-1) -(floor(Ny/2)):1:-1]*dKy;
    Ky=Ky';
    for n=1:dimension(3);
        temp(:,:,n)=Ky(n);
    end
    Ky=temp;
    Kz=[0:1:(floor(Nz/2)-1) -(floor(Nz/2)):1:-1]*dKz;
    for n=1:dimension(1)
        temp(n,:,:)=Kz(n);
    end
    Kz=temp;
    H=ones(dimension);
    A=Kz.^2;
    B=(Kx.^2+Ky.^2+Kz.^2);
    B(find(B==0))=1*10^-10;
    H=1/3-A./B;
elseif isfloat(orientation)&&length(orientation)==2
    theta=orientation(1);
    phi=orientation(2);
    Nz=dimension(1);
    Nx=dimension(2);
    Ny=dimension(3);
    %calculate sample width in k space
    dKx=2*pi/Nx/dx;
    dKy=2*pi/Ny/dy;
    dKz=2*pi/Nz/dz;

    %-----------------------------------------------------------------------------
    %construct Salomir Filter Transfer Function
    Kx=[0:1:(floor(Nx/2)-1) -(floor(Nx/2)):1:-1]*dKx;
    temp=ones(dimension(1:3));
    for n=1:dimension(2)
        temp(:,n,:)=Kx(n);
    end
    Kx=temp;
    Ky=[0:1:(floor(Ny/2)-1) -(floor(Ny/2)):1:-1]*dKy;
    Ky=Ky';
    for n=1:dimension(3);
        temp(:,:,n)=Ky(n);
    end
    Ky=temp;
    Kz=[0:1:(floor(Nz/2)-1) -(floor(Nz/2)):1:-1]*dKz;
    for n=1:dimension(1)
        temp(n,:,:)=Kz(n);
    end
    Kz=temp;
    H=ones(dimension);
    A=(Kz*cos(theta)*cos(phi)-Ky*sin(theta)*cos(phi)+Kx*sin(phi)).^2;
    B=(Kx.^2+Ky.^2+Kz.^2);
    B(find(B==0))=1*10^-10;
    H=1/3-A./B;
else
    Ny=dimension(1);
    Nx=dimension(2);
    Nz=dimension(3);
    %calculate sample width in k space
    dKx=2*pi/Nx/dx;
    dKy=2*pi/Ny/dy;
    dKz=2*pi/Nz/dz;
    %-----------------------------------------------------------------------------
    %construct Salomir Filter Transfer Function
    Kx=[0:1:(floor(Nx/2)-1) -(floor(Nx/2)):1:-1]*dKx;
    temp=ones(dimension(1:3));
    for n=1:dimension(2)
        temp(:,n,:)=Kx(n);
    end
    Kx=temp;
    Ky=[0:1:(floor(Ny/2)-1) -(floor(Ny/2)):1:-1]*dKy;
    Ky=Ky';
    for n=1:dimension(1);
        temp(n,:,:)=Ky(n);
    end
    Ky=temp;
    Kz=[0:1:(floor(Nz/2)-1) -(floor(Nz/2)):1:-1]*dKz;
    for n=1:dimension(3)
        temp(:,:,n)=Kz(n);
    end
    Kz=temp;
    H=ones(dimension);
    A=Kz.^2;
    B=(Kx.^2+Ky.^2+Kz.^2);
    B(find(B==0))=1*10^-10;
    H=1/3-A./B;
end
end