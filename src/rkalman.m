function [x, G, V, P, th] = rkalman(sys, c, tau, y, x_0, V_0)
    % This version is a bit different from the original version in
    % https://mathworks.com/matlabcentral/fileexchange/54308-robust-kalman-filtering-package?focused=5751770&tab=function
    % The difference is only in the inputs of the function and removing the exitence of 
    % gain condition from the code.
    % We are very thankfull of Prof. Zorzi for making the code available.
    %
    %
    % RKALMAN  robust Kalman estimator.
    %
    %   [x_e,G,V] = RKALMAN(A,B,C,D,y,c,tau) compute the "delayed" robust 
    %   Kalman estimator for the nominal discrete-time model
    %      
    %      x[n+1] = Ax[n] + Bv[n]         {State equation}
    %        y[n] = Cx[n] + Dv[n]         {Measurements}
    %
    %   with disturbance and measurement noise v with variance I. The robust 
    %   "delayed" estimator uses only past measurements up to y[n-1] to 
    %   generate the optimal robust estimate x_e[n] of x[n] and is easier to 
    %   embed in digital control loops. The equation of the robust "delayed" 
    %   estimator:
    %
    %      x_e[n+1|n] = Ax_e[n|n-1] + G (y[n] - Cx_e[n|n-1])
    %
    %   RKALMAN returns the delayed estimate x_e, the estimator gain G, the 
    %   least feavorable covariance matrix V of the estimation error 
    %   x[n+1]-x_e[n+1|n].
    %
    %   The robust estimator is designed knowing that the true model belongs to 
    %   a ball about the nominal one. The models in that ball are such that 
    %   the Tau-divergence between them and the nominal model is less than a 
    %   certain tolerance. To design the robust Kalman estimator, it is 
    %   necessary to specify: 
    %    
    %    - the tolerance c (striclty positive)
    %    - the parameter tau (in the interval [0,1]) of the Tau-divergence
    %    - the measurements y
    %
    %   For more details see: "Robust Kalman filtering under incremental 
    %   model perturbations" by M. Zorzi
    %
    %   See also RKITERATION, MAXTOL.


    %   Author(s): Mattia Zorzi 20-8-2015



    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;


    % parameters 
    n=size(A,1);
    p=size(C,1);
    m=size(B,2);
    T=size(y,1);

    % transform the model
    A=A-B*D'*(D*D')^-1*C;
    B=[(B*(eye(m)-D'*(D*D')^-1*D)*B')^0.5 zeros(n,m-n)];
    D=[zeros(p,n) (D*D')^0.5];
    Q=B*B';
    R=D*D';

    % check conditions
    assert(c>0,'Tolerance c must be positive')
    assert(rank(ctrb(A,Q))>=n,'The model must be reachable')
    assert(rank(ctrb(A',C'))>=n,'The model must be observable')
    cN = maxtol(A,B,C,D,tau,2*n);
%     if c>cN
%         warning('Tolerance c is too large: the filter gain may not exist')
%     end

    % init
    x = zeros(T,n);
    V = zeros(n,n,T+1);
    P = zeros(n,n,T+1);
    G = zeros(n,p,T+1);
    th = zeros(T,1);
    x(1,:) = x_0';
    V(:,:,1) = V_0;


    % iterative part
    for k=1:T
        [x(k+1,:), V(:,:,k+1), G(:,:,k+1), P(:,:,k+1), th(k)]=rkiteration(A,B,C,D,V(:,:,k),tau,c,x(k,:),y(k,:));
    end

    % resize
    x=x(2:T+1,:);
    P=P(:,:,2:T+1);
    V=V(:,:,2:T+1);
    G=G(:,:,2:T+1);
end
function cN=maxtol(A,B,C,D,tau,N)
    %MAXTOL maximum tolerance for which the robust Kalman estimator is 
    %       guaranteed to converge.
    %
    %   [c] = MAXTOL(A,B,C,D,tau) returns the maximum tolerance for which the 
    %   robust Kalman estimator is guaranteed to converge. 
    %
    %   The inputs of MAXTOL are the nominal discrete-time model associated to
    %   the robust Kalman filter
    %      
    %      x[n+1] = Ax[n] + Bv[n]         {State equation}
    %        y[n] = Cx[n] + Dv[n]         {Measurements}
    %
    %   with disturbance and measurement noise v with variance I (independent), 
    %   the parameter tau (in [0,1]) of the Tau-divergence.
    %
    %   For more details see: "Robust Kalman filtering under incremental 
    %   model perturbations" by M. Zorzi 
    %
    %   See also RKALMAN, RKITERATION.

    %   Author(s): Mattia Zorzi 20-8-2015


    % check inputs
    if nargin==5
        N=2*n;
    end

    n=size(A,1);
    m=size(B,2);
    p=size(C,1);

    % construction of the Gramians
    DR=kron(eye(N),D);
    Re=[];
    Ob=[];
    ObR=[];
    H=[];
    L=[];
    for k=1:N
        Re=[Re A^(k-1)*B]; 
        Ob=[C*A^(k-1); Ob];
        ObR=[A^(k-1); ObR];
        T=[];
        for l=1:N
            if l<=N-(k-1)
                T=[T zeros(p,m)];
            else
                T=[T C*A^(l-N+(k-1)-1)*B];
            end
        end
        H=[T; H];
        T=[];
        for l=1:N
            if l<=N-(k-1)
                T=[T zeros(n,m)];
            else
                T=[T A^(l-N+(k-1)-1)*B];
            end
        end
        L=[T; L];
    end

    Om = Ob'*inv(DR*DR'+ H*H')*Ob;
    J = ObR -L*H'*inv(DR*DR'+H*H')*Ob;
    M= L*inv(eye(N*m)+H'*inv(DR*DR')*H)*L';

    % computation of phiN tilde
    phiNtilde = 1/max(eig(M));

    % computation of phiN 
    value=1;
    t1=0;
    t2=(1-10^-10)*phiNtilde;
    while abs(value)>=10^-9
           theta=0.5*(t1+t2);
           Sth = -eye(N*n)+theta*M;
           Omth = Om +J'*theta*(Sth)^-1*J;
           value=-min(eig(Omth));
           if value>0
                t2=theta;
           else
                t1=theta;
           end
    end
    phiN=min(phiNtilde,theta);

    % computation of the initial condition (that is P_q bar) 
    Pq=dare(A',C',B*B',D*D');
    lambda=min(eig(Pq));

    % computation of theta_N
    if tau==1
        thN=-log(1-phiN*lambda)/lambda;
    else
        thN=(1-(1-lambda*phiN)^(1-tau))/((1-tau)*lambda);
    end
    if tau>=0 & tau<1
        thNmax=((1-tau)*max(eig(Pq)))^-1;
    else
        thNmax=+inf;
    end
    % computation of cN

    if thN>=thNmax
        cN=+inf;
    else
        if tau==0
            cN=(log(det(eye(n)-thN*Pq))+trace((eye(n)-thN*Pq)^-1)-n);
        end
        if tau>0 & tau<1
            Lq=chol(Pq)';
            cN=trace(-1/(tau*(1-tau))*(eye(n)-thN*(1-tau)*Lq'*Lq)^(tau/(tau-1))+1/(1-tau)*(eye(n)-thN*(1-tau)*Lq'*Lq)^(1/(tau-1))+1/tau*eye(n));
        end
        if tau==1
            Lq=chol(Pq)';
            cN=trace(expm(thN*Lq'*Lq)*(thN*Lq'*Lq-eye(n))+eye(n));
        end
    end
end
function [x_pred, V, G, P, th]=rkiteration(A,B,C,D,V,tau,c,x,y)
    %RKIteration robust Kalman iteration 
    %
    %   [x_pred, V_next, G]=RKIteration(A,B,C,D,V,tau,c,x,y) performs one 
    %   iteration of the "delayed" robust Kalman estimator for the nominal 
    %    discrete-time model
    %      
    %      x[n+1] = Ax[n] + Bw[n]         {State equation}
    %        y[n] = Cx[n] + Dw[n]         {Measurements}
    %
    %   with disturbance w with variance I.
    %   The equation of the robust "delayed" estimator:
    %
    %      x[n+1|n] = Ax[n|n-1] + G (y[n] - Cx[n|n-1])
    %
    %   where y[n] past measurement, x[n|n-1] past estimate, and G robust
    %   Kalman gain. 
    %
    %   RKALMAN returns the new estimate x[n+1|n], the least feavorable
    %   covariance matrix V_next of the estimation error x[n+1]-x[n+1|n], and 
    %   the estimator gain G.
    %
    %   The robust estimator is designed knowing that the true model belongs to 
    %   a ball about the nominal one. The models in that ball are such that 
    %   the Tau-divergence between them and the nominal model is less than a 
    %   certain tolerance. To design the estimator gain G, it is necessary to 
    %   specify: 
    %    
    %    - the tolerance c (striclty positive)
    %    - the parameter tau (in the interval [0,1]) of the Tau-divergence
    %    - the least favorable covariance matrix V of the previous estimation  
    %      error x[n]-x[n|n-1].
    %
    %   For more details see: "Robust Kalman filtering under incremental 
    %   model perturbations" by M. Zorzi
    %
    %   See also RKALMAN, MAXTOL.


    %   Author(s): Mattia Zorzi 20-8-2015


    % parameters 

    n=size(A,1);
    p=size(C,1);
    m=size(B,2);

    Q=B*B';
    R=D*D';

    % robust Kalman gain
    G = A*V*C'*(C*V*C'+R)^-1;

    % compute the the prediction
    x_pred=(A*x'+G*(y-C*x'))';


    % update the least favorable covariance matrix V
    P = (A-G*C)*V*(A-G*C)'+(B-G*D)*(B-G*D)';
    L=chol(P)';
    value=1;
    t1=0;
    if tau==1
        t2=10/max(eig(P)); % because exp(700) is (more or less) the maximum number that matlab can store
    else
        e = eig(P);
        r = max(abs(e));
        t2=(1-10^-5)*((1-tau)*r)^-1;
    end
    while abs(value)>=10^-9
       th=0.5*(t1+t2);
       if tau==0
           value=trace(inv(eye(n)-th*P)-eye(n)) + log(det(eye(n)-th*P))-c;
       end
       if tau>0 & tau<1
          value=trace(-1/(tau*(1-tau))*(eye(n)-(1-tau)*th*L'*L)^(tau/(tau-1))+1/(1-tau)*(eye(n)-(1-tau)*th*L'*L)^(1/(tau-1))+1/tau*eye(n))-c;
       end
       if tau==1
           value=trace(th*L'*L*expm(th*L'*L)-expm(th*L'*L)+eye(n))-c;
       end
       if value>0
            t2=th;
       else
            t1=th;
       end
    end
    Vold=V;
    if tau==1
        V = L*expm(th*L'*L)*L';
    else
        V = L*(eye(n)-(1-tau)*th*L'*L)^(1/(tau-1))*L';
    end
end