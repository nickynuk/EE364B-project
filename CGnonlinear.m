%% CG for nonlinear objective with l1 norm regularization
clear
clc
close all
%% open phase nifti file
file=open_nii('/Users/nickynuk/Dropbox/Stanford class/EE364b/EE364Project/gky_7t_oreint_SWI_neutral_SENSE_3_1_1_unwrap_pdfphase31MASK.nii');

%ensure even number of slices
dimension=size(file.img);

%check for odd dimensions;
%z direction
if mod(dimension(3),2)==1
    file.img=file.img(:,:,1:(dimension(3)-1),:);
    file.hdr.dime.dim(4)=file.hdr.dime.dim(4)-1;
end
%y direction
if mod(dimension(1),2)==1
    file.img=file.img(1:(dimension(1)-1),:,:,:);
    file.hdr.dime.dim(2)=file.hdr.dime.dim(2)-1;
end
%x direction
if mod(dimension(2),2)==1
    file.img=file.img(1:(dimension(2))-1,:,:,:);
    file.hdr.dime.dim(3)=file.hdr.dime.dim(3)-1;
end
dimension=size(file.img);%update dimension
%check for data sets with only one echo
if size(dimension)<4
    dimension=[dimension,1];
end

% resolution in m
dy=file.hdr.dime.pixdim(2)/1000;
dx=file.hdr.dime.pixdim(3)/1000;
dz=file.hdr.dime.pixdim(4)/1000;

phaseims = file.img;
Mask = zeros(dimension);
Mask(phaseims~=0)=1;
%% open magnitude nifti file
file=open_nii('/Users/nickynuk/Dropbox/Stanford class/EE364b/EE364Project/gky_7t_oreint_SWI_neutral_SENSE_3_1_0.nii');
file.img = file.img(:,:,:,1);
%ensure even number of slices
dimension=size(file.img);

%check for odd dimensions;
%z direction
if mod(dimension(3),2)==1
    file.img=file.img(:,:,1:(dimension(3)-1),:);
    file.hdr.dime.dim(4)=file.hdr.dime.dim(4)-1;
end
%y direction
if mod(dimension(1),2)==1
    file.img=file.img(1:(dimension(1)-1),:,:,:);
    file.hdr.dime.dim(2)=file.hdr.dime.dim(2)-1;
end
%x direction
if mod(dimension(2),2)==1
    file.img=file.img(1:(dimension(2))-1,:,:,:);
    file.hdr.dime.dim(3)=file.hdr.dime.dim(3)-1;
end
dimension=size(file.img);%update dimension
%check for data sets with only one echo
if size(dimension)<4
    dimension=[dimension,1];
end

% resolution in m
dy=file.hdr.dime.pixdim(2)/1000;
dx=file.hdr.dime.pixdim(3)/1000;
dz=file.hdr.dime.pixdim(4)/1000;

M = Mask.*file.img;
M = M/max(M(:));
%% define transfer function
D=salomirtf(dimension,dx,dy,dz,'axial');

%% define gradient masks
Wx = zeros(dimension);
for i = 1:dimension(3)
    tmp = imfilter(phaseims(:,:,i),fspecial('log',[1,5],.1));
    Wx(:,:,i) = 1-im2bw(tmp);
end
%Wx = Wx(:,1:end-1,:);

Wy = zeros(dimension);
for i = 1:dimension(3)
    tmp = imfilter(phaseims(:,:,i),fspecial('log',[5,1],.1));
    Wy(:,:,i) = 1-im2bw(tmp);
end
%Wy = Wy(1:end-1,:,:);

Wz = zeros(dimension);
for i = 1:dimension(2)
    tmp = imfilter(squeeze(phaseims(:,i,:)),fspecial('log',[1,5],.1));
    Wz(:,i,:) = 1-im2bw(tmp);
end
%Wz = Wz(:,:,1:end-1);

%% define gradient matrices
[Ex,Ey,Ez] = gradienttf(dimension,dx,dy,dz,'axial');

%% define quadratic constraint function
PIMS = fftn(phaseims);

%% define nonlinear objective and norm 1 regularization term
z = zeros(1,1,2);
z(1,1,:) = [-1,1];
lambda =1;
nln_obj = @(fftnx) lambda*sum(sum(sum(abs(M.*(exp(1i*ifftn(D.*fftnx))-exp(1i*phaseims)).^2))));
n1_obj  = @(x) sum(sum(sum(abs(Wx.*convn(x,[-1,1],'same'))+abs(Wy.*convn(x,[-1,1]','same'))+abs(Wz.*convn(x,z,'same')))));
%% gradient of nonlinear objective
cosphase = Mask.^2.*cos(phaseims);
sinphase = Mask.^2.*sin(phaseims);
dfo = @(fftnx) -2*(ifftn(D.*fftn(-cosphase.*sin(ifftn(D.*fftnx))+sinphase.*cos(ifftn(D.*fftnx)))));
%% Initialization of CG based algorithm
x0 = zeros(dimension);
fvec = [];
fvec2 = [];
tvec = [];
x = x0;
X = fftn(x);
maxiter = 250;
scale = 0.2;
for repeat = 1:8;
   
    %% compute gradient of norm1 regularization term
    Cx = Wx(:,1:end-1,:).*convn(x,[-1,1],'valid');
    Cy = Wy(1:end-1,:,:).*convn(x,[-1,1]','valid');
    Cz = Wz(:,:,1:end-1).*convn(x,z,'valid');
    norm1Cx = zeros(dimension);
    norm1Cx(:,1:end-1,:) = abs(Cx);
    norm1Cx(1:end-1,:,:) = norm1Cx(1:end-1,:,:) + abs(Cy);
    norm1Cx(:,:,1:end-1) = norm1Cx(:,:,1:end-1) + abs(Cz);
    norm1Cx(norm1Cx == 0) = 1e-2;
    norm1Cxall = sum(sum(sum(norm1Cx)));
    
    % (A_1 + A_2)deltax = B_1 + B_2
    B1 = - (1/norm1Cxall) * ...
           ( Wx.*convn(Cx,[1,-1],'full') ...
          + Wy.*convn(Cy,[1,-1]','full') ...
          + Wz.*convn(Cz,-z,'full') );
    B2 = -2*lambda*real(1i * ifftn(D.* fftn((M.^2).*exp(-1i*ifftn(D.*X) + 1i*phaseims)) ));
    B = B1 + B2;
    A_1 = @(deltax) (1/norm1Cxall) * ( ...
          ( Wx.*convn( Wx(:,1:end-1,:).*convn(deltax,[-1,1],'valid'),[1,-1],'full') ) ...
         +( Wy.*convn( Wy(1:end-1,:,:).*convn(deltax,[-1,1]','valid'), [1,-1]', 'full')) ...
         +( Wz.*convn( Wz(:,:,1:end-1).*convn(x,z,'valid'), -z, 'full')) ...
         );
    A_2 = @(deltax) 2*lambda* ifftn(D.* fftn((M.^2) .*ifftn(D.*fftn(deltax))));
    Afun = @(deltax) A_1(deltax) + A_2(deltax);
    
    % Solve for deltax using conjugate gradient 
    [ deltax ] = CG( Afun, B , maxiter,x);
    
    % Line search
    t = 1;
    alpha = 1/4; 
    beta = 1/2;
    while (nln_obj(fftn(x + t*deltax)) > nln_obj(X) +  alpha*t*sum(sum(sum(dfo(X).*deltax))) )
       t = beta*t;    
       if (t <= 10^-4)
           break;
       end
    end
    
    % update x
    x = x + t*deltax; % Compute fourier transform of x
    X = fftn(x);

    fval = nln_obj(X); % save nonlinear objective
    fvec2 = [fvec2, n1_obj(x)]; % save norm1 regularizatoion value
    fvec = [fvec, fval]; 
    
end
file.img = x;
file.hdr.dime.dim(5) = 1;
save_untouch_nii(file,[file.fileprefix,'_nln_norm1'])
