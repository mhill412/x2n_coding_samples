function convect(Ra) 

    % Originally Written by Paul Tackley October 2012
    % modified for ESS311, 2022

    % Infinite Prandtl number convection (e.g. mantle convection)
    % using a finite volume discretization (staggered grid)
    %  - velocities at cell faces
    %  - temperature at cell centres
    %  - streamfunction and vorticity (when used) at cell corners

    % Equations are nondimensionalised!

    % Boundary conditions are: 
    %  - free-slip on all boundaries
    %  - dT/dx=0 on sides, T=1 at bottom and T=0 at top

    % If variable_viscosity=true then it uses a DIRECT solver
    % If variable_viscosity=false then it uses an iterative MULTIGRID solver
    %  to solve the streamfunction-vorticity version of the equations
    % 

    % ---------SET INPUT PARAMETERS HERE------------

    % number of grid points; nz=vertical. 
    % For variable_viscosity=false use a POWER OF 2 for each

    %Ra=1e+5;   % Rayleigh number
    H=0;      % Internal heating rate
    initial_temperature=0.5;   % Initial temperature (0 to 1)

    variable_viscosity=true;   % true or false; switch on/off variable viscosity
    viscosity_contrast_temperature=1e6;  % factor variation with T
    viscosity_contrast_depth=1;      % factor variation with depth

    nx=64; nz=32;   % #grid points nz=vertical
    nsteps=5000; % 10000   % number of time steps to take

    draw_colorbar=false;    % If graphics doesn't work properly, try setting this to false
    
    %disp('Iteration: ');  % addedd for ESS311    
    % --------DO NOT EDIT ANYTHING BELOW THIS LINE------

    figure

    % ---------------------initialise fields

    A_viscosity=log(viscosity_contrast_temperature);
    B_viscosity=log(viscosity_contrast_depth);

    d=1/nz; dsq=d^2;

    s=zeros(nx+1,nz+1); w=zeros(size(s)); 
    T=zeros(nx,nz);   
    vx=zeros(size(s)); vz=zeros(size(s));
    eta=zeros(nx,nz,2);

    T=initial_temperature + 0.1*rand(nx,nz);


    for step=1:nsteps

    % ---------------------take timestep

        if variable_viscosity

            % calculate viscosity field
            for j=1:nz
                depth=(0.5+nz-j)*d;
                viscosity(:,j) = exp(-A_viscosity*(T(:,j)-0.5)+B_viscosity*(depth-0.5));
            end

            % calculate rhs at vz points (=-Ra*T)
            rhs(1:nx,2:nz) = -Ra * 0.5*(T(1:nx,1:nz-1)+T(1:nx,2:nz));

            % call direct solver        
            [vx vz p] = direct_solve_Stokes_sparse(viscosity,rhs,d);

        else  % not variable viscosity: use stream-function
            % calculate streamfunction from dT/dx

            %  first, right-hand side = Ra*(dT/dx)
            %    calculated at cell corners

            rhs=zeros(nx+1,nz+1);
            for i=2:nx
                for j=2:nz
                    dT_dx=0.5*(T(i,j)+T(i,j-1)-T(i-1,j)-T(i-1,j-1))/d;
                    rhs(i,j) = -Ra*dT_dx;
                end
            end

            %    next, call poisson solver twice to find stream-function

            w = poisson_multigrid (rhs,w);
            s = poisson_multigrid (  w,s);

            %   calculate velocities from streamfunction

            for i=1:nx
                ip=1+mod(i,nx);
                for j=1:nz
                    vx(i,j)=-(s(i,j+1)-s(i,j))/d;
                    vz(i,j)= (s(ip,j) -s(i,j))/d;
                end
            end

        end

        % determine stable time-step
        dt_diff=0.1*dsq; % diffusive
        vmax=max(max(max(abs(vx))),max(max(abs(vz))));
        dt_adv=0.5*d/vmax;
        dt=min(dt_diff,dt_adv);

      
	%disp(step);  % comment out for ESS311
	%fprintf('%i.',step);  % added for ESS311
	%disp(dt);  % comment out for ESS311

        % time-step temperature
        %    - diffusion & internal heating
        T = T + dt*( del2T(T) + H); 

        %    - advection    
        T = donor_cell_advection(T,vx,vz,dt);

        % finally, plot
        %subplot(2,1,1);   %  commented out for ESS311 
        str1=['(yellow = hot, blue = cold)'];
	str2=['Step: ',num2str(step),', Ra = ',num2str(Ra,'%2.1e'),')'];
	contourf(T'); title(["Relative Temperature",str1,str2]); axis equal;  hold on
	xlabel('Non-dimensional horizontal distance');
        ylabel('Non-dimensional vertical distance');

        if draw_colorbar
            colorbar; 
        end
        quiver(vx(1:nx,1:nz)', vz(1:nx,1:nz)'); hold off

	% Commented out for ESS311, to simplify the plot
        %if variable_viscosity
        %    subplot(2,1,2); contourf(log10(viscosity')); title('log10(Viscosity)')
        %else
        %    subplot(2,1,2); contourf(s'); title('Stream-function')
        %end
        %if draw_colorbar
        %    colorbar; 
        %end
        pause(0.001)

    end
end


%-----------------------------------------------------------------------

function [delsqT] = del2T (T)
%   calculates del2 of temperature field assuming:
%    - zero flux side boundaries 
%    - T=0 at top
%    - T=1 at bottom

    [nx nz] = size(T);
    delsqT = zeros(size(T));
    dsq = 1./(nz-1)^2;
    
    for i=1:nx
        im=max(i-1,1); ip=min(i+1,nx);
        
        for j=1:nz
                    
            T_xm=T(im,j); T_xp=T(ip,j); 
            if j==1   % near lower boundary where T=1
                T_zm=2-T(i,j);
            else
                T_zm=T(i,j-1);
            end
            if j==nz % near upper boundary where T=0
                T_zp=-T(i,j);
            else
                T_zp=T(i,j+1);
            end
            
            delsqT(i,j)=(T_xm+T_xp+T_zm+T_zp-4*T(i,j))/dsq;

        end
    end
    
end

%-----------------------------------------------------------------------

function [Tnew] = donor_cell_advection(T,vx,vz,dt)

    [nx nz] = size(T);
    Tnew = zeros(size(T));
    d = 1./(nz-1);

    for i=1:nx    
        for j=1:nz
        
            % Donor cell advection (div(VT))
            
            if vx(i,j)>0
                flux_xm = T(i-1,j)*vx(i,j);  
            else
                flux_xm = T(i  ,j)*vx(i,j);  
            end
            if i<nx
                if vx(i+1,j)>0
                    flux_xp = T(i  ,j)*vx(i+1,j);  
                else
                    flux_xp = T(i+1,j)*vx(i+1,j);  
                end
            else
                flux_xp = 0;
            end
    
            if vz(i,j)>0
                flux_zm = T(i,j-1)*vz(i,j);  
            else
                flux_zm = T(i  ,j)*vz(i,j);  
            end
            if j<nz
                if vz(i,j+1)>=0
                    flux_zp = T(i,j)*vz(i,j+1);  
                else
                    flux_zp = T(i,j+1)*vz(i,j+1);  
                end
            else
                flux_zp = 0;
            end

            dT_dt = (flux_xm-flux_xp+flux_zm-flux_zp)/d;
            Tnew(i,j) = T(i,j) + dt*dT_dt;
        
        end
    end
 
end
    
%-----------------------------------------------------------------------

function [vx vz p] = direct_solve_Stokes_sparse(eta,rhs,h)

% direct solver for 2D variable viscosity Stokes flow
% sxx,x + szx,z - p,x = 0
% szz,z + sxz,x - p,z = rhs        
% (dimensional rhs=density*g; nondimensionally, rhs=-Ra*T)

% 2*(eta_xp*(u_xp-u0)/h-eta_xm*(u0-u_xm)/h)/h +
%   (eta_zp*(u_zp-u0-w_xpzp-w_xmzp)/h-eta_zm*(u0-u_zm-w_xpzm-w_xmzm)/h)/h=0

% Boundary conditions:
%  shear stress=0 -> put etaxz=0 on the boundaries
%  impermeable -> no action needed

% rhs should have dimensions of the number of cells

[nx nz]=size(rhs);
oh=1/h;
ohsq=1/h^2;

f=zeros(nx*nz*3,1);
idx=3; idz=nx*3;

n_non0=(11+11+4)*3; % estimated number of non-zero matrix entries
r=zeros(1,n_non0); s=r; s=r;  % faster to set these to ~correct size
n=0; % number of non-zero point

for iz=1:nz
    for ix=1:nx
        icell=ix-1+(iz-1)*nx;
        ieqx=icell*3+1; vxc=ieqx;
        ieqz=ieqx+1;    vzc=ieqz;
        ieqc=ieqx+2;     pc=ieqc;
        
        % viscosity points
        etaii_c  = eta(ix,iz);
        if ix>1
            etaii_xm = eta(ix-1,iz);
        end
        if iz>1
            etaii_zm = eta(ix,iz-1);
        end
        if (ix>1 && iz>1)
            etaxz_c  = (eta(ix,iz)*eta(ix-1,iz)*eta(ix,iz-1)*eta(ix-1,iz-1))^0.25;
        else
            etaxz_c = 0;
        end
        if (ix>1 && iz<nz) 
            etaxz_zp = (eta(ix,iz+1)*eta(ix-1,iz+1)*eta(ix,iz)*eta(ix-1,iz))^0.25;
        else
            etaxz_zp = 0;
        end
        if (ix<nx && iz>1)   
            etaxz_xp  = (eta(ix+1,iz)*eta(ix,iz)*eta(ix+1,iz-1)*eta(ix,iz-1))^0.25;
        else
            etaxz_xp = 0;
        end
        
        xmom_zero_eta = (etaii_c==0 && etaii_xm==0 && etaxz_c==0 && etaxz_zp==0);
        zmom_zero_eta = (etaii_c==0 && etaii_zm==0 && etazx_c==0 && etazx_xp==0);
        
        % x-momentum
       
        if ix>1 &&  ~xmom_zero_eta
            n=n+1; r(n)=ieqx; c(n)=vxc    ; s(n)=-ohsq*(2*etaii_c+2*etaii_xm+etaxz_c+etaxz_zp);
            n=n+1; r(n)=ieqx; c(n)=vxc-idx; s(n)= 2*ohsq*etaii_xm;
            n=n+1; r(n)=ieqx; c(n)=vzc    ; s(n)=-ohsq*etaxz_c;
            n=n+1; r(n)=ieqx; c(n)=vzc-idx; s(n)= ohsq*etaxz_c;
            n=n+1; r(n)=ieqx; c(n)=pc     ; s(n)=-oh;
            n=n+1; r(n)=ieqx; c(n)=pc-idx ; s(n)= oh;
            if (ix+1<=nx)
                n=n+1; r(n)=ieqx; c(n)=vxc+idx; s(n)=2*ohsq*etaii_c;
            end
            if (iz+1<=nz)
                n=n+1; r(n)=ieqx; c(n)=vxc+idz    ; s(n)= ohsq*etaxz_zp;
                n=n+1; r(n)=ieqx; c(n)=vzc+idz    ; s(n)= ohsq*etaxz_zp;
                n=n+1; r(n)=ieqx; c(n)=vzc+idz-idx; s(n)=-ohsq*etaxz_zp;
            end
            if (iz>1)
                n=n+1; r(n)=ieqx; c(n)=vxc-idz; s(n)=ohsq*etaxz_c;
            end
            f(ieqx)             =   0;
        else
            n=n+1; r(n)=ieqx; c(n)=vxc; s(n)=1;
            f(ieqx)             =    0;
        end
        
        % z-momentum
        if iz>1 && ~zmom_zero_eta
            n=n+1; r(n)=ieqz; c(n)=vzc    ; s(n)=-ohsq*(2*etaii_c+2*etaii_zm+etaxz_c+etaxz_xp);
            n=n+1; r(n)=ieqz; c(n)=vzc-idz; s(n)= 2*ohsq*etaii_zm;
            n=n+1; r(n)=ieqz; c(n)=vxc    ; s(n)=-ohsq*etaxz_c;
            n=n+1; r(n)=ieqz; c(n)=vxc-idz; s(n)= ohsq*etaxz_c;
            n=n+1; r(n)=ieqz; c(n)=pc     ; s(n)=-oh;
            n=n+1; r(n)=ieqz; c(n)=pc-idz ; s(n)= oh;
            if iz+1<=nz
                n=n+1; r(n)=ieqz; c(n)=vzc+idz ; s(n)= 2*ohsq*etaii_c;
            end
            if ix+1<=nx
                n=n+1; r(n)=ieqz; c(n)=vzc+idx    ; s(n)=  ohsq*etaxz_xp;
                n=n+1; r(n)=ieqz; c(n)=vxc+idx    ; s(n)=  ohsq*etaxz_xp;
                n=n+1; r(n)=ieqz; c(n)=vxc+idx-idz; s(n)= -ohsq*etaxz_xp;
            end
            if ix>1
                n=n+1; r(n)=ieqz; c(n)=vzc-idx ; s(n)= ohsq*etaxz_c;
            end
            f(ieqz)             =   rhs(ix,iz);
        else
            n=n+1; r(n)=ieqz; c(n)=vzc ; s(n)=1;
            f(ieqz)             = 0;
        end
        
        % continuity
        if (ix==1 && iz==1) || (xmom_zero_eta && zmom_zero_eta)
            n=n+1; r(n)=ieqc; c(n)=pc; s(n)=1;
        else
            n=n+1; r(n)=ieqc; c(n)=vxc; s(n)=-oh;
            n=n+1; r(n)=ieqc; c(n)=vzc; s(n)=-oh;
            if ix+1<=nx
                n=n+1; r(n)=ieqc; c(n)=vxc+idx; s(n)=oh;
            end
            if iz+1<=nz
                n=n+1; r(n)=ieqc; c(n)=vzc+idz; s(n)=oh;
            end
        end
        f(ieqc)             =   0;
        
    end
end
size(s);   % semicolan added for ESS 311
m=sparse(r,c,s);
solution = m\f;
for iz=1:nz
    for ix=1:nx
        i0=(ix-1)+(iz-1)*nx;
        vx(ix,iz)=solution(i0*3+1); % x-velocity
        vz(ix,iz)=solution(i0*3+2); % z-velocity
        p (ix,iz)=solution(i0*3+3); % pressure
    end
end
vx(nx+1,:)=0  ; %impermeable right boundary
vz(:,nz+1)=0  ; %impermeable top boundary
p(:,:) = p(:,:) - mean(mean(p(1:nx,1:nz))); % zero mean pressure

end

%-----------------------------------------------------------------------
 
function [f] = poisson_multigrid (R,fstart)
%function [f] = poisson_multigrid (R,fguess)
% calls poisson_vee repeatedly until desired convergence is reached
% R and f should be a power-of-two in horizontal direction and a 
%  power-of-two plus one in the vertical direction (second dimension)

f=fstart;
Rmean=mean(mean(abs(R)));

residue = delsqzp(f)-R; resmean=mean(mean(abs(residue)));

while resmean>1e-3*Rmean
    
    f=poisson_vee(R,f);

    residue = delsqzp(f)-R; resmean=mean(mean(abs(residue)));

end

end

%--------------------------------------------------------------------

function [f] = poisson_vee (R,fguess)
% function [f] = poisson_vee (R)
% solve poisson equation in 2D using multigrid approach: del2(f)=R
% Boundary conditions:
% - sides are periodic
% - top & bottom are zero
% This routine calls itself recursively to find the correction to the next
% coarsest level

s=size(R); nx=s(1);nz=s(2); %surf(R');title('Rhs');pause
h=1/(nz-1);     % grid spacing, assuming vertical dimension=1

f=fguess;
alpha=0.9;

if nz>5

    for iter=1:2               %  smoothing steps
        f=relaxdelsq(f,R,alpha);
    end
    
    residue = delsqzp(f)-R;    % residue to coarse grid
    rescoarse = residue(1:2:nx,1:2:nz); coarseguess=zeros(size(rescoarse));
    
    coarsecorrection = poisson_vee(rescoarse,coarseguess); % recursive call of this routine
    
    finecorrection = zeros(size(R));             % linear interpolation to fine grid
    finecorrection(1:2:nx,1:2:nz) = coarsecorrection;
    for j=2:2:nz-1
        finecorrection(:,j)=0.5*(finecorrection(:,j-1)+finecorrection(:,j+1));
    end
    for i=2:2:nx-1
        finecorrection(i,:)=0.5*(finecorrection(i-1,:)+finecorrection(i+1,:));
    end
    
    f = f - finecorrection/2;    % add fine correction to solution
    
    for iter=1:2               % more smoothing steps
       f=relaxdelsq(f,R,alpha);
    end
    
else
    
    for iter=1:50               % reach solution using multiple iters
       f=relaxdelsq(f,R,alpha);
    end

end

end

%--------------------------------------------------------------------

function [f] = relaxdelsq(a,rhs,alpha)
% Gauss-Seidel relax delsquared assuming sides, top&bottom are zero
s=size(a); nx=s(1); nz=s(2); ohsq=(nz-1)^2; afac=alpha/(4*ohsq);
d=zeros(s); rmean=0;
for irb = 0:1 %red-black
    for i=2:nx-1
        for j = 2+mod(i+irb,2):2:nz-1   % don't relax top & bottom boundaries
            res = (a(i+1,j)+a(i-1,j)+a(i,j-1)+a(i,j+1)-4*a(i,j))*ohsq-rhs(i,j);
            a(i,j) = a(i,j) + afac*res;
            rmean = rmean+abs(res);
        end
    end
end
rmean=rmean/(nx*(nz-1));
f=a;
end

%--------------------------------------------------------------------

function [d] = delsqzp(a)
% delsquared assuming sides, top&bottom are zero
s=size(a); nx=s(1); nz=s(2); ohsq=(nz-1)^2;
d=zeros(s);
for j = 2:nz-1
    for i = 2:nx-1
        d(i,j)=(a(i+1,j)+a(i-1,j)+a(i,j-1)+a(i,j+1)-4*a(i,j))*ohsq;
    end
end

end

