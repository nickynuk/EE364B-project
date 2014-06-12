function [ x ] = CG( Hx,Gx,maxiter,x )
%% Function description
%Computes inv(H)*Gx using Conjugate Gradients
%Hx is a function handle
%Gx is the gradient

%% Initialize vectors and parameters
i=0; %counts number of iterations
epsq=0.01;%tolerance
dimension = size(Gx);
b = Hx(Gx);

%Compute residual (r=b-Ax)
r = b - Hx(Hx(x));
rsnew = 1000;
p=r;
rsold=sum(sum(sum(r.^2)));
gaphist = [];

while (i<maxiter) & (rsnew > epsq) % stop early when residual is increasing
    Ap=Hx(Hx(p));
    alpha=rsold/sum(sum(sum(double(p).*double(Ap))));
    x=x+alpha*p;
    r=r-alpha*Ap;
    rsnew=sum(sum(sum(r.^2)))
    if (sqrt(rsnew)<epsq | rsnew > rsold)
        break;
    end
    p=r+rsnew/rsold*p;
    rsold=rsnew;
    i=i+1;
    gap = Hx(x);
    gaphist = [gaphist norm(gap(:)-Gx(:))/norm(Gx(:))];
end
figure ()
plot(gaphist);
end