% This is a 2-D CG-FFT VIE solver for modeling communication systems in 2-D inhomogeneous environments.

% Shuai S. A. Yuan
% E-mail: shuaiyuan1997@zju.edu.cn

function field
clc;clear

tic
wavelength=1;  %--wavelength
k0=2*pi/wavelength; %--wave number in free space

z=120*pi;      %--wave impedance in free space
epr=1.00000000e+000-1.044000000e+008*1j;         % ---permittivity  ****
%epr=5;
deltax=wavelength/2 /sqrt(1);  %--grid size
r=100;  %length
a=sqrt(deltax*deltax/pi);        %--effective area
const=-j*z*epr/(k0*(epr-1));

N_c=1;
gap=1*r;


% x1=0;
% y1=4.5;


x_source=230;
y_source=350;


% M=round((gap*(N_c-1)+2*r)/(deltax))+10;
% N=M;
area_length=1000;
M=round(area_length/(deltax))+10;
N=M;

Ein=zeros(M,N);
model=zeros(M,N);
model2=zeros(M,N);
Z=zeros(M,N);
xc=-M/2*deltax+deltax/2;
yc=N/2*deltax-deltax/2;
epr_arr=const*ones(M,N);
load model
for f_n=1:N  %  ÁÐ
    for f_m=1:M  %  ÐÐ
        x=xc+(f_n-1)*deltax;
        y=yc-(f_m-1)*deltax;
        flag=p(f_m,f_n);
        RR=sqrt((x-x_source)^2+(y-y_source)^2);
        Ein_total(f_m,f_n)=(-1j/4)*besselh(0,2,k0*RR);
        %Ein_total(f_m,f_n)=exp(-j*k0*x);
        if (flag)
            model(f_m,f_n)=1;

            % Ein(f_m,f_n)=exp(-j*k0*x);
            Ein(f_m,f_n)=(-1j/4)*besselh(0,2,k0*RR);
        end

        R=sqrt((x-xc)^2+(y-yc)^2);
        Z(f_m,f_n)=z*pi*a/2*besselj(1,k0*a)*besselh(0,2,k0*R); % impedance matrix  %
    end
end
Z(1,1)=z*pi*a/2*besselh(1,2,k0*a);

% pcolor(abs(Ein_total))
% shading interp
% axis equal
% colorbar

figure(1)
pcolor(model)
shading interp
axis equal
colorbar
toc


Zf1=Z(:,end:-1:2);
Zf2=Z(end:-1:2,:);
Zf3=Zf2(:,end:-1:2);
Z=[Z,Zf1;Zf2,Zf3];


J=zeros(M,N);  % initial current % set
%load initial_result_e3.mat
%J=J;
r0=-Ein+fft_mv(Z,J,M,N,model,const,epr_arr);%
test=norm(r0);
p1=-fft_mv_trans(Z,r0,M,N,model,const,epr_arr);
q1=-p1;
change3=q1.*conj(q1);
con3=sum(sum(change3));
changex=Ein.*conj(Ein);
conx=sum(sum(changex));

tic
for n=1:10^8

    %  A^{*}r_{n-1}
    change1=change3;
    con1=con3;
    %  Ap_{n}
    Ap=fft_mv(Z,p1,M,N,model,const,epr_arr);
    change2=Ap.*conj(Ap);
    con2=sum(sum(change2));
    %  alpha
    alpha=con1/con2;
    %  Update J
    J=J+alpha*p1;
    %  Update r
    r0=r0+alpha*Ap;
    change4=r0.*conj(r0);

    %  Beta
    q1=fft_mv_trans(Z,r0,M,N,model,const,epr_arr);
    change3=q1.*conj(q1);
    con3=sum(sum(change3));
    beta=con3/con1;
    %  p1
    p1=-q1+beta*p1;

    con4=sum(sum(change4))/conx;


    if (mod(n,10)==0)
        con4;
    end

    %  error truncation
    if con4<1*10^(-3)
        n
        break;
    end

end
toc
% figure(1)
% pcolor(abs(J));
% shading interp
% axis equal
% colorbar




A_field=zeros(M,N);
A_field=fft_field(Z,J,M,N,model,const,epr_arr);


%E_field=-1j*k0*z*A_field+Ein_total;%
E_field=-A_field+Ein_total;%
E_field=E_field/max(max(E_field));
E_field=E_field.*abs(model-1);
figure(2);
pcolor(10*log10(abs(E_field)));
shading interp
axis equal
colorbar
hold on
% axis equal
% colorbar




% field calculation
function current=fft_field(Z,J,M,N,model,const,epr_arr)
current=ifft2(fft2(Z).*fft2(J,2*M-1,2*N-1));   
current=current([1:M],[1:N]);                 
current=epr_arr.*J+current;                    
%current=current.*model;                       

%  Z*j
function current=fft_mv(Z,J,M,N,model,const,epr_arr)
current=ifft2(fft2(Z).*fft2(J,2*M-1,2*N-1));   
current=current([1:M],[1:N]);               
current=epr_arr.*J+current;                    
current=current.*model;                       

%  Z'*j
function current=fft_mv_trans(Z,J,M,N,model,const,epr_arr)
current=ifft2(fft2(conj(Z)).*fft2(J,2*M-1,2*N-1));
current=current([1:M],[1:N]);                      
% current=current.';                               
current=current+conj(epr_arr).*J;                 
current=current.*model;                         





