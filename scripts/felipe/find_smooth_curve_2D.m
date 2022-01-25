function [X_sm,Y_sm,h]=find_smooth_curve_2D(M,n_sk,n_sm,h,t_sk,w_sk,t_sm)

%   Here I plot the msh skeleton
    lambda=0.2;
    Mau=[M;M(1,:)];
    [x0,y0,sgma0,Nx_vert,Ny_vert]=centros_y_valores(Mau,lambda);
%    if length(t_sk)~=n_sk
%        disp('n_sk')
        [t_sk,w_sk]=nodos_GaussLegendre(0,1,n_sk);
%    end
%    [X_sk,Y_sk,W_sk]=setup_integration_nodes(Mau,t,w);
    [X_sk,Y_sk,W_sk,Nx_sk,Ny_sk]=setup_integration_nodes(Mau,t_sk,w_sk);
%    if length(t_sm)~=n_sm
%        disp('n_sm')
        [t_sm,~]=nodos_GaussLegendre(0,1,n_sm);
%    end
    [Xb_sm,Yb_sm,psNx_sm,psNy_sm]=setup_base_points(Mau,t_sm,Nx_vert,Ny_vert);
    [X_sm,Y_sm,h]=mi_interp_newton(h,Xb_sm,Yb_sm,psNx_sm,psNy_sm,X_sk,Y_sk,W_sk,Nx_sk,Ny_sk,x0,y0,sgma0);
%     toc
%     figure
%     plot([X_sm; X_sm(1)],[Y_sm; Y_sm(1)],'LINEWIDTH',3)
%     hold on
%     M=M(1:npts_contour,:);
%     plot([M(:,1); M(1,1)],[M(:,2); M(1,2)])
%     hold on
%     plot([M(:,1); M(1,1)],[M(:,2); M(1,2)],'o')
%     axis square
%     grid
end


function [x0,y0,sgma,Nx_vert,Ny_vert]=centros_y_valores(Mau,lambda)
    [n,~]=size(Mau);
    n=n-1;
    PP=(Mau(1:n,:)+Mau(2:n+1,:))/2;
    x0=PP(:,1);
    y0=PP(:,2);
    d=sqrt((Mau(1:n,1)-Mau(2:n+1,1)).^2+(Mau(1:n,2)-Mau(2:n+1,2)).^2);
    sgma=d*lambda;
    
    P1=Mau(1:n,:);
    P2=Mau(2:n+1,:);
    Pdif=P2-P1;
    Norm_x=Pdif(:,2);
    Norm_y=-Pdif(:,1);

    Norm_x=[Norm_x(end); Norm_x];
    Norm_y=[Norm_y(end); Norm_y];
    d_Norm=sqrt(Norm_x.^2+Norm_y.^2);
    Norm_x=Norm_x./d_Norm;
    Norm_y=Norm_y./d_Norm;
    
    Nx_vert=(Norm_x(1:n)+Norm_x(2:n+1))/2;
    Ny_vert=(Norm_y(1:n)+Norm_y(2:n+1))/2;

%     figure
%     size(Mau)
%     size(Nx_vert)
%     size(Ny_vert)
%  plot([Mau(:,1)],[Mau(:,2)])
%  hold on
%     quiver(Mau(1:n,1),Mau(1:n,2),Nx_vert,Ny_vert)
%     return
end





function [X_sk,Y_sk,W_sk,Nx_sk,Ny_sk]=setup_integration_nodes(Mau,t,w)
    [n,~]=size(Mau);
    n=n-1;
    P1=Mau(1:n,:);
    P2=Mau(2:n+1,:);
    Pdif=P2-P1;
    Norm_x=Pdif(:,2);
    Norm_y=-Pdif(:,1);

    dl=sqrt((Mau(1:n,1)-Mau(2:n+1,1)).^2+(Mau(1:n,2)-Mau(2:n+1,2)).^2);

    Norm_x=Norm_x./dl;
    Norm_y=Norm_y./dl;
    
    ones_aux=ones(size(t));
    N_x=Norm_x*ones_aux;
    N_y=Norm_y*ones_aux;
    
    X_sk=P1(:,1)+(P2(:,1)-P1(:,1))*t;
    Y_sk=P1(:,2)+(P2(:,2)-P1(:,2))*t;

    W_sk=dl*w;

    X_sk=reshape(X_sk',[],1);
    Y_sk=reshape(Y_sk',[],1);
    W_sk=reshape(W_sk',[],1);
    
    Nx_sk=reshape(N_x',[],1);
    Ny_sk=reshape(N_y',[],1);

% figure
%     plot(X_sk,Y_sk,'.')
%     hold on
%     quiver(X_sk,Y_sk,Nx_sk,Ny_sk)    
%     sum(W_sk)
end


function [Xb_sm,Yb_sm,psNx_sm,psNy_sm]=setup_base_points(Mau,t,Nx_sk,Ny_sk)
    [n,~]=size(Mau);
    n=n-1;
    P1=Mau(1:n,:);
    P2=Mau(2:n+1,:);
    
    Xb_sm=P1(:,1)+(P2(:,1)-P1(:,1))*t;
    Yb_sm=P1(:,2)+(P2(:,2)-P1(:,2))*t;
    

    Xb_sm=reshape(Xb_sm',[],1);
    Yb_sm=reshape(Yb_sm',[],1);


    Nx_sk_au=[Nx_sk;Nx_sk(1)];
    Ny_sk_au=[Ny_sk;Ny_sk(1)];
    
    
    psNx_sm=Nx_sk_au(1:n)*(1-t)+Nx_sk_au(2:n+1)*(t);
    psNy_sm=Ny_sk_au(1:n)*(1-t)+Ny_sk_au(2:n+1)*(t);

    psNx_sm=reshape(psNx_sm',[],1);
    psNy_sm=reshape(psNy_sm',[],1);
%figure
%    quiver(X_sm,Y_sm,psNx_sm,psNy_sm)    
   
end



function [SGMA,dxSGMA,dySGMA]=sgma_eval(X,Y,x0,y0,sgma0)
    lambda=0.2;
    di=(x0-X').^2+(y0-Y').^2;
    
    alpha=1/max(sgma0)^2/5;

%    alpha=1./(2*(max(sgma0)/lambda/10).^2);
    DDExp=exp(-alpha.*di);
    F=sgma0.'*DDExp;
    D=sum(DDExp);
    SGMA=(F./D)';
    dxDDExp=-alpha*2*(x0-X').*DDExp;
    dyDDExp=-alpha*2*(y0-Y').*DDExp;    
    dFdx=sgma0.'*dxDDExp;
    dDdx=sum(dxDDExp);
    dFdy=sgma0.'*dyDDExp;
    dDdy=sum(dyDDExp);
    dxSGMA=(-(dFdx.*D-F.*dDdx)./D.^2)';
    dySGMA=(-(dFdy.*D-F.*dDdy)./D.^2)';
end


function [F,gradF_x,gradF_y]=eval_F(X,Y,X_sk,Y_sk,W_sk,Nx_sk,Ny_sk,x0,y0,sgma0)

    [sgma,dxSGMA,dySGMA]=sgma_eval(X,Y,x0,y0,sgma0);

    d2=(X_sk-X').^2+(Y_sk-Y').^2;
    [DDD,DDDp,DDDsgma]=mi_phi_der(d2,sgma);    
    
    F=(((Nx_sk.*W_sk)')*((DDD).*(X_sk-X'))+((Ny_sk.*W_sk)')*((DDD).*(Y_sk-Y')));
    F=F';
    a=length(dxSGMA);
    TTT_dxSGMA=sparse(1:a,1:a,dxSGMA);
    TTT_dySGMA=sparse(1:a,1:a,dySGMA);

    gradF_x_aux=-((DDDp).*(X_sk-X'))+DDDsgma*TTT_dxSGMA;
    gradF_x_1=(((Nx_sk.*W_sk)')*((gradF_x_aux).*(X_sk-X'))+((Ny_sk.*W_sk)')*((gradF_x_aux).*(Y_sk-Y')));
    gradF_x_2=-((Nx_sk.*W_sk)')*((DDD));
    gradF_x=(gradF_x_1+gradF_x_2)';
    
    gradF_y_aux=-((DDDp).*(Y_sk-Y'))+DDDsgma*TTT_dySGMA;
    gradF_y_1=(((Nx_sk.*W_sk)')*((gradF_y_aux).*(X_sk-X'))+((Ny_sk.*W_sk)')*((gradF_y_aux).*(Y_sk-Y')));
    gradF_y_2=-((Ny_sk.*W_sk)')*((DDD));   
    gradF_y=(gradF_y_1+gradF_y_2)';   

    F=F-.5;
        
end

function [F,Fp,Fsgma]=mi_phi_der_backup(r2,sgma)
    %This is (1/r)*(d Phi/dr)
    a=length(sgma);
    TTT=sparse(1:a,1:a,1./(2*sgma.^2));
    TTT3=sparse(1:a,1:a,1./(sgma.^3));
    Mexp=exp(-r2*TTT);
    F=-((2.*Mexp)./r2 - 2./r2)/(4*pi);
    Fp=(- (Mexp - 1)./(r2.^2) - (Mexp./r2)*TTT)/pi;
    Fp=-Fp;
    Fsgma=Mexp*TTT3/(2*pi);
    Fsgma=-Fsgma;
    F=F;
    
    
%     x=r2*TTT;
%     y=expm1ox(x);
%     F=-y*TTT/(2*pi);
%     y2=expm1pxox2(-x);
%     Fp=-(Mexp.*y2)*(TTT*TTT)/pi;

end


function [F,Fp,Fsgma]=mi_phi_der(r2,sgma)
    %This is (1/r)*(d Phi/dr)
    a=length(sgma);
    TTT=sparse(1:a,1:a,1./(2*sgma.^2));
    TTT3=sparse(1:a,1:a,1./(sgma.^3));
    Mexp=exp(-r2*TTT);
    x=r2*TTT;
    y=expm1ox(x);
    F=-y*TTT/(2*pi);
    
    y2=expm1pxox2_v2(x);
    Fp=(y2)*(TTT*TTT)/pi;

%    F=-((Mexp-1)./r2)/(2*pi);
%    F=-((2.*Mexp)./r2 - 2./r2)/(4*pi);
%    Fp=(- (Mexp - 1)./(r2.^2) - (Mexp./r2)*TTT)/pi;
%    Fp=-Fp;
    
    
    
    Fsgma=Mexp*TTT3/(2*pi);
    Fsgma=-Fsgma;
%    F=F;
end



function y=expm1ox(x)
    y=zeros(size(x));
    index1=find(abs(x)<.05);
    index2=find(abs(x)>=.05);
    y(index2)=(exp(-x(index2))-1)./x(index2);
    y(index1)=- x(index1).^10/39916800 + x(index1).^9/3628800 - x(index1).^8/362880 + x(index1).^7/40320 - x(index1).^6/5040 + x(index1).^5/720 - x(index1).^4/120 + x(index1).^3/24 - x(index1).^2/6 + x(index1)/2 - 1;
end

function y=expm1pxox2(x)
    y=zeros(size(x));
    index1=find(abs(x)<.1);
    index2=find(abs(x)>=.1);
    y(index2)=(exp(-x(index2))-1+x(index2))./x(index2).^2;
    y(index1)=- x(index1).^9/39916800 + x(index1).^8/3628800 - x(index1).^7/362880 + x(index1).^6/40320 - x(index1).^5/5040 + x(index1).^4/720 - x(index1).^3/120 + x(index1).^2/24 - x(index1)/6 + .5;
end

%x^11/518918400 - x^10/43545600 + x^9/3991680 - x^8/403200 + x^7/45360 - x^6/5760 + x^5/840 - x^4/144 + x^3/30 - x^2/8 + x/3 - 1/2

function y=expm1pxox2_v2(x)
    y=zeros(size(x));
    index1=find(abs(x)<.1);
    index2=find(abs(x)>=.1);
    y(index2)=(exp(-x(index2))-1+x(index2).*exp(-x(index2)))./x(index2).^2;
    y(index1)= x(index1).^9/3991680 - x(index1).^8/403200 + x(index1).^7/45360 - x(index1).^6/5760 + x(index1).^5/840 - x(index1).^4/144 + x(index1).^3/30 - x(index1).^2/8 + x(index1)/3 - .5;
%    y(index1)=- x(index1).^9/39916800 + x(index1).^8/3628800 - x(index1).^7/362880 + x(index1).^6/40320 - x(index1).^5/5040 + x(index1).^4/720 - x(index1).^3/120 + x(index1).^2/24 - x(index1)/6 + .5;
end






function [F,Fp]=mifyfp(h,Xb_sm,Yb_sm,psNx_sm,psNy_sm,X_sk,Y_sk,W_sk,Nx_sk,Ny_sk,x0,y0,sgma0,indices)
    X=Xb_sm(indices)+psNx_sm(indices).*h(indices);
    Y=Yb_sm(indices)+psNy_sm(indices).*h(indices);
    [F,gradF_x,gradF_y]=eval_F(X,Y,X_sk,Y_sk,W_sk,Nx_sk,Ny_sk,x0,y0,sgma0);
    Fp=gradF_x.*psNx_sm(indices)+gradF_y.*psNy_sm(indices);
end

function [X_sm,Y_sm,h]=mi_interp_newton(h_in,Xb_sm,Yb_sm,psNx_sm,psNy_sm,X_sk,Y_sk,W_sk,Nx_sk,Ny_sk,x0,y0,sgma0)
    tol=1e-10;
    maxiter=30;
    err=tol+1;
    niter=0;
    if sum(~isfinite(h_in))>0        
        h=zeros(size(Xb_sm));
    elseif length(h_in)~=length(Xb_sm)
        h=zeros(size(Xb_sm));
    else
        h=h_in;
    end
%    h=zeros(size(Xb_sm));
    f=zeros(size(Xb_sm));
    fp=zeros(size(Xb_sm));
    indices=(1:length(h))';
    while (niter<=maxiter)&(err>tol)
        [f_out,fp_out]=mifyfp(h,Xb_sm,Yb_sm,psNx_sm,psNy_sm,X_sk,Y_sk,W_sk,Nx_sk,Ny_sk,x0,y0,sgma0,indices);
        f(indices)=f_out;
        fp(indices)=fp_out;
        delt=f./fp;
        h=h-delt;
        niter=niter+1;
        err_vect=(abs(delt));
        indices=find(err_vect>tol);
        err=max(max(err_vect));
    end
    X_sm=Xb_sm+psNx_sm.*h;
    Y_sm=Yb_sm+psNy_sm.*h;    
end


function [X_sm,Y_sm,h]=mi_interp_newton_backup(Xb_sm,Yb_sm,psNx_sm,psNy_sm,X_sk,Y_sk,W_sk,Nx_sk,Ny_sk,x0,y0,sgma0)
    tol=1e-10;
    maxiter=30;
    err=tol+1;
    niter=0;
    h=zeros(size(Xb_sm));
    while (niter<=maxiter)&(err>tol)
        [f,fp]=mifyfp(h,Xb_sm,Yb_sm,psNx_sm,psNy_sm,X_sk,Y_sk,W_sk,Nx_sk,Ny_sk,x0,y0,sgma0);
        delt=f./fp;
        err=max(max(abs(delt)));
        h=h-delt;
        niter=niter+1;
    end
    niter;
    err;
    X_sm=Xb_sm+psNx_sm.*h;
    Y_sm=Yb_sm+psNy_sm.*h;
%figure
%plot(X_sm,Y_sm,'LINEWIDTH',4)
    
end

