%% SCP for nonlinear, nonconvex QSM formulation Using TV Norm
clear
clc
close all
%% open phase nifti file
file=open_nii('/Users/grant/Documents/MATLAB/EE364bProjectData/GKY_7T_Parrec/gky_7t_oreint_SWI_neutral_SENSE_3_1_1_unwrap_pdfphase31MASK.nii');

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
file=open_nii('/Users/grant/Documents/MATLAB/EE364bProjectData/GKY_7T_Parrec/gky_7t_oreint_SWI_neutral_SENSE_3_1_0.nii');
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


Wy = zeros(dimension);
for i = 1:dimension(3)
    tmp = imfilter(phaseims(:,:,i),fspecial('log',[5,1],.1));
    Wy(:,:,i) = 1-im2bw(tmp);
end


Wz = zeros(dimension);
for i = 1:dimension(2)
    tmp = imfilter(squeeze(phaseims(:,i,:)),fspecial('log',[1,5],.1));
    Wz(:,i,:) = 1-im2bw(tmp);
end

%% define nonlinear objective function
fo = @(fftnx) sum(sum(sum(abs(M.*(exp(1i*ifftn(D.*fftnx))-exp(1i*phaseims)).^2)))); 
cosphase = Mask.^2.*cos(phaseims);
sinphase = Mask.^2.*sin(phaseims);
dfo = @(fftnx,x) -2*(ifftn(D.*fftn(-cosphase.*sin(ifftn(D.*fftnx))+sinphase.*cos(ifftn(D.*fftnx)))));

%% define smoothness constraint function
z = zeros(1,1,2);
z(1,1,:) = [-1,1];
epsilon2 = 2.5e6;
%TV-norm Version
c2 = @(x)epsilon2 - sum(sum(sum(sqrt((Wx.*convn(x,[-1,1],'same')).^2+(Wy.*convn(x,[-1,1]','same')).^2+(Wz.*convn(x,z,'same')).^2))));

%% initial guess
x = zeros(dimension);
%% SQP

fprintf('Initializing SCP \n');

% compute gradient of the objectives and constraints
%TV norm
Cx = Wx(:,1:end-1,:).*convn(x,[-1,1],'valid');
Cy = Wy(1:end-1,:,:).*convn(x,[-1,1]','valid');
Cz = Wz(:,:,1:end-1).*convn(x,z,'valid');
normCx = zeros(dimension);
XX = Cx.^2;
YY = Cy.^2;
ZZ = Cz.^2;
normCx(:,1:end-1,:) = XX;
normCx(1:end-1,:,:) = normCx(1:end-1,:,:) + YY;
normCx(:,:,1:end-1) = normCx(:,:,1:end-1) + ZZ;
normCx = sqrt(normCx);
normCx(normCx == 0) = 1e-2;
DC2 = - (convn(Wx(:,1:end-1,:).*Cx./normCx(:,1:end-1,:),[1,-1],'full') + convn(Wy(1:end-1,:,:).*Cy./normCx(1:end-1,:,:),[1,-1]','full') + convn(Wz(:,:,1:end-1).*Cz./normCx(:,:,1:end-1),-z,'full'));
%compute constraints and gradient of constraints
X= fftn(x);
C2 = c2(x);
lambda2 = 1/(1000*C2);%Initialize Lagrange Multipliers
dFo = dfo(X);
fplot = fo(X); %initialize plot of objective function values
%% initialize arrays to store gradient and step history for Hessian approximation using L-DFP
fprintf('Begin SCP \n');
for repeat = 1:40;%set maximum iterations of outer loop   
%% Solve QP
%minimize (dF'*p + .5*p'*B*p)
%subject to 
%lambda : c+dc'p>=0
%cvx_end
fprintf('Begin QP\n');
tic
if DC2(:)'*DC2(:) <= 1e-6
    lambda2 = 1/(1000*C2);
    p = -dFo;
    fprintf('!!!!')
else
    lambda2 = (DC2(:)'*dFo(:)-C2)/(DC2(:)'*DC2(:));
    if lambda2 < 0;
        lambda2 = 0;
    end
    p =  lambda2*DC2-dFo;
end
toc
%% update SQP step
%line search
alpha = 1;
beta = .01;
while c2(x+alpha*p)<0
    alpha = .5*alpha;
    if alpha<1e-6
        break;
    end
end
alpha
merit = @(x)fo(fftn(x)) - lambda2*c2(x);
meritx = merit(x);
dmerit = dFo(:)-lambda2*DC2(:);
while merit(x+alpha*p) > meritx + beta*alpha*dmerit'*p(:);
    alpha = .5*alpha;
    if alpha<1e-6
        %fprintf('SCP stuck \n');
        break;
    end
end
alpha
%save last step
x = x + alpha*p;

%gradient of TV norm
Cx = Wx(:,1:end-1,:).*convn(x,[-1,1],'valid');
Cy = Wy(1:end-1,:,:).*convn(x,[-1,1]','valid');
Cz = Wz(:,:,1:end-1).*convn(x,z,'valid');
normCx = zeros(dimension);
XX = Cx.^2;
YY = Cy.^2;
ZZ = Cz.^2;
normCx(:,1:end-1,:) = XX;
normCx(1:end-1,:,:) = normCx(1:end-1,:,:) + YY;
normCx(:,:,1:end-1) = normCx(:,:,1:end-1) + ZZ;
normCx = sqrt(normCx);
normCx(normCx == 0) = 1e-2;
DC2 = - (convn(Wx(:,1:end-1,:).*Cx./normCx(:,1:end-1,:),[1,-1],'full') + convn(Wy(1:end-1,:,:).*Cy./normCx(1:end-1,:,:),[1,-1]','full') + convn(Wz(:,:,1:end-1).*Cz./normCx(:,:,1:end-1),-z,'full'));

%compute constraints and gradient of constraints
X = fftn(x);
C2 = c2(x);
dFo = dfo(X);

%update x
fplot = [fplot fo(X)]

if abs((fplot(end-1)-fplot(end))/fplot(end-1))<=1e-4
    fprintf('SCP Converged \n');
    break;
end

end
%%save file
% file.img = x;
% file.hdr.dime.dim(5) = 1;
% save_untouch_nii(file,[file.fileprefix,'_SCP2_5e6'])


