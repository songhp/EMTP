function [s_hat,A_index,Count,loglikelihood,delta2]=DORE(H,Ht,y,r,varargin)
% Double Overrelaxation ECME (DORE) Thresholding routine: 
% Coded by Kun Qiu (kqiu@iastate.edu)
% updated Feb 14, 2010
% 改动： 迭代停止改为相对误差 gap=norm(y - H(s_hat))/norm(y);
% 
% Function usage
% =======================================================
% INPUT(compulsory):
% H:                              the sensing matrix (H can be either an matrix or function handle that computes H*x)
% Ht:                             the adjoint function handle that computes H^T*y (if H is a matrix, just input [])
% y:                               the measurement column vector
% r:                                the sparsity level of the signal to be recovered
% 
% INPUT(optional):
% 'IsHOrthonormal':     whether H is orthonormal: 
%                                   0: no
%                                   1: yes
%                                   (default=1)
% 'InitialSig':                 the mx1 initial signal estimate
%                                   (default=zeros(m,1))
% 'Thresh':                   tolerance threshold for stopping the algorithm
%                                   (default=1e-14)
% 'Visibility':                 Option to view the reconstruction process (valued in {0,1})
%                                   0: do not view
%                                   1: do view
%                                   (default=0)
% ========================================================
% OUTPUT:
% s_hat:              the signal estimate
% A_index:          the estimated support set
% Count:             Count of the number of EM iterations
% loglikelihood:  the value of the log-likelihood function of the final estimate
% delta2:             the variance component estimate
% ========================================================

if (nargin-length(varargin))~=4
    error('Missing required inputs!');
end

RunMode=1;

if ~isa(H, 'function_handle')
    H_temp=H;
    Ht=@(x) H'*x;
    H=@(x) H*x;
else
    RunMode=RunMode+10;
end

N=length(y);
m=length(Ht(y));

%Setting default values for the optional inputs
z_init=zeros(m,1);
thresh=1e-14;
Visibility=0;

%Read the optional inputs
if (rem(length(varargin),2)==1)
    error('Optional inputs must go by pairs!');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case upper('IsHOrthonormal')
                RunMode=RunMode-1+varargin{i+1};
            case upper('Inverse')
                RunMode=101;
                HHt_inv = varargin{i+1};
                if isa(HHt_inv,'float')
                   HHt_inv =@(z) HHt_inv*z; 
                elseif      ~isa(HHt_inv,'function_handle') 
                    error('Inverse needs to be a functrion handle or a matrix.')
                end
            case upper('InitialSig')
                z_init=varargin{i+1};
            case upper('Thresh')
                thresh=varargin{i+1};
            case upper('Visibility')
                Visibility=varargin{i+1};     
            otherwise
                error(['Unrecognized optional input: ''' varargin{i} '''']);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RunMode==0
    HHt_inv=inv(H_temp*H_temp');
    HHt_inv_y=HHt_inv*y;
    clear H_temp;
        
    %Initialization
    z=z_init;
    [z_sort,z_index]=sort(abs(z),'descend');
    A_index=z_index(1:r)';
    s_A=z(A_index);
    s_hat_pre2=zeros(m,1);
    s_hat_pre2(A_index)=s_A;
    Hs_hat_pre2=H(s_hat_pre2);
    HHt_invHs_hat_pre2=HHt_inv*Hs_hat_pre2;

    z=s_hat_pre2+Ht(y-Hs_hat_pre2);
    [z_sort,z_index]=sort(abs(z),'descend');
    A_index=z_index(1:r)';
    s_A=z(A_index);
    s_hat_pre=zeros(m,1);
    s_hat_pre(A_index)=s_A;
    Hs_hat_pre=H(s_hat_pre);
    HHt_invHs_hat_pre=HHt_inv*Hs_hat_pre;
    delta2_pre=(y-Hs_hat_pre)'*(HHt_inv_y-HHt_invHs_hat_pre)/N;
    loglikelihood_pre=-N/2*log(2*pi*delta2_pre)-N/2;

    exit_flag=0;
    p=0;
    while ~exit_flag
        p=p+1;
        %STEP 1
        z=s_hat_pre+Ht(HHt_inv_y-HHt_invHs_hat_pre);
        [z_sort,z_index]=sort(abs(z),'descend');
        A_index=z_index(1:r)';
        s_A=z(A_index);
        s_bar=zeros(m,1);
        s_bar(A_index)=s_A;
        A_index_bar=A_index;
        Hs_bar=H(s_bar);
        HHt_invHs_bar=HHt_inv*Hs_bar;
        delta2_bar=(y-Hs_bar)'*(HHt_inv_y-HHt_invHs_bar)/N;
        loglikelihood_bar=-N/2*log(2*pi*delta2_bar)-N/2;
        %STEP 2
        alpha1=(Hs_bar-Hs_hat_pre)'*(HHt_inv_y-HHt_invHs_bar)/((Hs_bar-Hs_hat_pre)'*(HHt_invHs_bar-HHt_invHs_hat_pre));
        z_prime=s_bar+alpha1*(s_bar-s_hat_pre);
        %STEP 3
        Hz_prime=(1+alpha1)*Hs_bar-alpha1*Hs_hat_pre;
        HHt_invHz_prime=(1+alpha1)*HHt_invHs_bar-alpha1*HHt_invHs_hat_pre;
        alpha2=(Hz_prime-Hs_hat_pre2)'*(HHt_inv_y-HHt_invHz_prime)/((Hz_prime-Hs_hat_pre2)'*(HHt_invHz_prime-HHt_invHs_hat_pre2));
        z_prime=z_prime+alpha2*(z_prime-s_hat_pre2);
        %STEP 4
        [z_sort,z_index]=sort(abs(z_prime),'descend');
        A_index=z_index(1:r)';
        s_A=z_prime(A_index);
        s_hat_prime=zeros(m,1);
        s_hat_prime(A_index)=s_A;
        Hs_hat_prime=H(s_hat_prime);
        HHt_invHs_hat_prime=HHt_inv*Hs_hat_prime;
        delta2_prime=(y-Hs_hat_prime)'*(HHt_inv_y-HHt_invHs_hat_prime)/N;
        loglikelihood_prime=-N/2*log(2*pi*delta2_prime)-N/2;

        %STEP 5
        if loglikelihood_prime>loglikelihood_bar
            s_hat=s_hat_prime;
            delta2=delta2_prime;
            Hs_hat=Hs_hat_prime;
            HHt_invHs_hat=HHt_invHs_hat_prime;
            loglikelihood=loglikelihood_prime;
        else
            s_hat=s_bar;
            A_index=A_index_bar;
            delta2=delta2_bar;
            Hs_hat=Hs_bar;
            HHt_invHs_hat=HHt_invHs_bar;
            loglikelihood=loglikelihood_bar;
        end

%         gap=norm(s_hat-s_hat_pre)^2/m; %修改收敛条件
gap=norm(y - H(s_hat))/norm(y);

        if Visibility
            clc;
            display(['Iteration=',num2str(p),', gap=',num2str(gap),' (target=',num2str(thresh),')']);
        end
        if gap<thresh  || p>100 %修改收敛条件
            exit_flag=1;
        else
            s_hat_pre2=s_hat_pre;
            s_hat_pre=s_hat;
            Hs_hat_pre2=Hs_hat_pre;
            Hs_hat_pre=Hs_hat;
            HHt_invHs_hat_pre2=HHt_invHs_hat_pre;
            HHt_invHs_hat_pre=HHt_invHs_hat;
        end
    end

    Count=p;
    return;
end

if RunMode==100|RunMode==101
    HHt_inv_y=HHt_inv(y);
    clear H_temp;
        
    %Initialization
    z=z_init;
    [z_sort,z_index]=sort(abs(z),'descend');
    A_index=z_index(1:r)';
    s_A=z(A_index);
    s_hat_pre2=zeros(m,1);
    s_hat_pre2(A_index)=s_A;
    Hs_hat_pre2=H(s_hat_pre2);
    HHt_invHs_hat_pre2=HHt_inv(Hs_hat_pre2);

    z=s_hat_pre2+Ht(y-Hs_hat_pre2);
    [z_sort,z_index]=sort(abs(z),'descend');
    A_index=z_index(1:r)';
    s_A=z(A_index);
    s_hat_pre=zeros(m,1);
    s_hat_pre(A_index)=s_A;
    Hs_hat_pre=H(s_hat_pre);
    HHt_invHs_hat_pre=HHt_inv(Hs_hat_pre);
    delta2_pre=(y-Hs_hat_pre)'*(HHt_inv_y-HHt_invHs_hat_pre)/N;
    loglikelihood_pre=-N/2*log(2*pi*delta2_pre)-N/2;

    exit_flag=0;
    p=0;
    while ~exit_flag
        p=p+1;
        %STEP 1
        z=s_hat_pre+Ht(HHt_inv_y-HHt_invHs_hat_pre);
        [z_sort,z_index]=sort(abs(z),'descend');
        A_index=z_index(1:r)';
        s_A=z(A_index);
        s_bar=zeros(m,1);
        s_bar(A_index)=s_A;
        A_index_bar=A_index;
        Hs_bar=H(s_bar);
        HHt_invHs_bar=HHt_inv(Hs_bar);
        delta2_bar=(y-Hs_bar)'*(HHt_inv_y-HHt_invHs_bar)/N;
        loglikelihood_bar=-N/2*log(2*pi*delta2_bar)-N/2;
        %STEP 2
        alpha1=(Hs_bar-Hs_hat_pre)'*(HHt_inv_y-HHt_invHs_bar)/((Hs_bar-Hs_hat_pre)'*(HHt_invHs_bar-HHt_invHs_hat_pre));
        z_prime=s_bar+alpha1*(s_bar-s_hat_pre);
        %STEP 3
        Hz_prime=(1+alpha1)*Hs_bar-alpha1*Hs_hat_pre;
        HHt_invHz_prime=(1+alpha1)*HHt_invHs_bar-alpha1*HHt_invHs_hat_pre;
        alpha2=(Hz_prime-Hs_hat_pre2)'*(HHt_inv_y-HHt_invHz_prime)/((Hz_prime-Hs_hat_pre2)'*(HHt_invHz_prime-HHt_invHs_hat_pre2));
        z_prime=z_prime+alpha2*(z_prime-s_hat_pre2);
        %STEP 4
        [z_sort,z_index]=sort(abs(z_prime),'descend');
        A_index=z_index(1:r)';
        s_A=z_prime(A_index);
        s_hat_prime=zeros(m,1);
        s_hat_prime(A_index)=s_A;
        Hs_hat_prime=H(s_hat_prime);
        HHt_invHs_hat_prime=HHt_inv(Hs_hat_prime);
        delta2_prime=(y-Hs_hat_prime)'*(HHt_inv_y-HHt_invHs_hat_prime)/N;
        loglikelihood_prime=-N/2*log(2*pi*delta2_prime)-N/2;

        %STEP 5
        if loglikelihood_prime>loglikelihood_bar
            s_hat=s_hat_prime;
            delta2=delta2_prime;
            Hs_hat=Hs_hat_prime;
            HHt_invHs_hat=HHt_invHs_hat_prime;
            loglikelihood=loglikelihood_prime;
        else
            s_hat=s_bar;
            A_index=A_index_bar;
            delta2=delta2_bar;
            Hs_hat=Hs_bar;
            HHt_invHs_hat=HHt_invHs_bar;
            loglikelihood=loglikelihood_bar;
        end

%         gap=norm(s_hat-s_hat_pre)^2/m;
gap=norm(y - H(s_hat))/norm(y);

        if Visibility
            clc;
            display(['Iteration=',num2str(p),', gap=',num2str(gap),' (target=',num2str(thresh),')']);
        end
        if gap<thresh || p>100 %修改收敛条件
            exit_flag=1;
        else
            s_hat_pre2=s_hat_pre;
            s_hat_pre=s_hat;
            Hs_hat_pre2=Hs_hat_pre;
            Hs_hat_pre=Hs_hat;
            HHt_invHs_hat_pre2=HHt_invHs_hat_pre;
            HHt_invHs_hat_pre=HHt_invHs_hat;
        end
    end

    Count=p;
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RunMode==10
       
    clear H_temp;
    ee=zeros(N,1);
    ind=randperm(N);
    ee(ind(1))=1;
    H_row=Ht(ee);
    HHt_inv=ones(N,1)/(norm(H_row)^2);
    HHt_inv_y=HHt_inv.*y;
    clear H_temp;
        
    %Initialization
    z=z_init;
    [z_sort,z_index]=sort(abs(z),'descend');
    A_index=z_index(1:r)';
    s_A=z(A_index);
    s_hat_pre2=zeros(m,1);
    s_hat_pre2(A_index)=s_A;
    Hs_hat_pre2=H(s_hat_pre2);
    HHt_invHs_hat_pre2=HHt_inv.*Hs_hat_pre2;

    z=s_hat_pre2+Ht(y-Hs_hat_pre2);
    [z_sort,z_index]=sort(abs(z),'descend');
    A_index=z_index(1:r)';
    s_A=z(A_index);
    s_hat_pre=zeros(m,1);
    s_hat_pre(A_index)=s_A;
    Hs_hat_pre=H(s_hat_pre);
    HHt_invHs_hat_pre=HHt_inv.*Hs_hat_pre;
    delta2_pre=(y-Hs_hat_pre)'*(HHt_inv_y-HHt_invHs_hat_pre)/N;
    loglikelihood_pre=-N/2*log(2*pi*delta2_pre)-N/2;

    exit_flag=0;
    p=0;
    while ~exit_flag
        p=p+1;
        %STEP 1
        z=s_hat_pre+Ht(HHt_inv_y-HHt_invHs_hat_pre);
        [z_sort,z_index]=sort(abs(z),'descend');
        A_index=z_index(1:r)';
        s_A=z(A_index);
        s_bar=zeros(m,1);
        s_bar(A_index)=s_A;
        A_index_bar=A_index;
        Hs_bar=H(s_bar);
        HHt_invHs_bar=HHt_inv.*Hs_bar;
        delta2_bar=(y-Hs_bar)'*(HHt_inv_y-HHt_invHs_bar)/N;
        loglikelihood_bar=-N/2*log(2*pi*delta2_bar)-N/2;
        %STEP 2
        alpha1=(Hs_bar-Hs_hat_pre)'*(HHt_inv_y-HHt_invHs_bar)/((Hs_bar-Hs_hat_pre)'*(HHt_invHs_bar-HHt_invHs_hat_pre));
        z_prime=s_bar+alpha1*(s_bar-s_hat_pre);
        %STEP 3
        Hz_prime=(1+alpha1)*Hs_bar-alpha1*Hs_hat_pre;
        HHt_invHz_prime=(1+alpha1)*HHt_invHs_bar-alpha1*HHt_invHs_hat_pre;
        alpha2=(Hz_prime-Hs_hat_pre2)'*(HHt_inv_y-HHt_invHz_prime)/((Hz_prime-Hs_hat_pre2)'*(HHt_invHz_prime-HHt_invHs_hat_pre2));
        z_prime=z_prime+alpha2*(z_prime-s_hat_pre2);
        %STEP 4
        [z_sort,z_index]=sort(abs(z_prime),'descend');
        A_index=z_index(1:r)';
        s_A=z_prime(A_index);
        s_hat_prime=zeros(m,1);
        s_hat_prime(A_index)=s_A;
        Hs_hat_prime=H(s_hat_prime);
        HHt_invHs_hat_prime=HHt_inv.*Hs_hat_prime;
        delta2_prime=(y-Hs_hat_prime)'*(HHt_inv_y-HHt_invHs_hat_prime)/N;
        loglikelihood_prime=-N/2*log(2*pi*delta2_prime)-N/2;

        %STEP 5
        if loglikelihood_prime>loglikelihood_bar
            s_hat=s_hat_prime;
            delta2=delta2_prime;
            Hs_hat=Hs_hat_prime;
            HHt_invHs_hat=HHt_invHs_hat_prime;
            loglikelihood=loglikelihood_prime;
        else
            s_hat=s_bar;
            A_index=A_index_bar;
            delta2=delta2_bar;
            Hs_hat=Hs_bar;
            HHt_invHs_hat=HHt_invHs_bar;
            loglikelihood=loglikelihood_bar;
        end

%         gap=norm(s_hat-s_hat_pre)^2/m;
gap=norm(y - H(s_hat))/norm(y);


        if Visibility
            clc;
            display(['Iteration=',num2str(p),', gap=',num2str(gap),' (target=',num2str(thresh),')']);
        end
        if gap<thresh || p>100 %修改收敛条件
            exit_flag=1;
        else
            s_hat_pre2=s_hat_pre;
            s_hat_pre=s_hat;
            Hs_hat_pre2=Hs_hat_pre;
            Hs_hat_pre=Hs_hat;
            HHt_invHs_hat_pre2=HHt_invHs_hat_pre;
            HHt_invHs_hat_pre=HHt_invHs_hat;
        end
    end

    Count=p;
    return;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RunMode==11|RunMode==1
    clear H_temp;
    %Initialization
    z=z_init;
    [z_sort,z_index]=sort(abs(z),'descend');
    A_index=z_index(1:r)';
    s_A=z(A_index);
    s_hat_pre2=zeros(m,1);
    s_hat_pre2(A_index)=s_A;
    Hs_hat_pre2=H(s_hat_pre2);

    z=s_hat_pre2+Ht(y-Hs_hat_pre2);
    [z_sort,z_index]=sort(abs(z),'descend');
    A_index=z_index(1:r)';
    s_A=z(A_index);
    s_hat_pre=zeros(m,1);
    s_hat_pre(A_index)=s_A;
    Hs_hat_pre=H(s_hat_pre);
    delta2_pre=norm(y-Hs_hat_pre)^2/N;
    loglikelihood_pre=-N/2*log(2*pi*delta2_pre)-N/2;

    exit_flag=0;
    p=0;
    while ~exit_flag
        p=p+1;
        %STEP 1
        z=s_hat_pre+Ht(y-Hs_hat_pre);           %%%%%%%% Costly
        [z_sort,z_index]=sort(abs(z),'descend');
        A_index=z_index(1:r)';
        s_A=z(A_index);
        s_bar=zeros(m,1);
        s_bar(A_index)=s_A;
        A_index_bar=A_index;
        Hs_bar=H(s_bar);                        %%%%%%%% Costly
        y_Hs_bar=y-Hs_bar;
        delta2_bar=norm(y_Hs_bar)^2/N;
        loglikelihood_bar=-N/2*log(2*pi*delta2_bar)-N/2;
        %STEP 2
        H_s_bar_s_hat_pre=Hs_bar-Hs_hat_pre;
        alpha1=H_s_bar_s_hat_pre'*y_Hs_bar/norm(H_s_bar_s_hat_pre)^2;
        z_prime=s_bar+alpha1*(s_bar-s_hat_pre);
        %STEP 3
        Hz_prime=(1+alpha1)*Hs_bar-alpha1*Hs_hat_pre;
        H_z_prime_s_hat_pre2=Hz_prime-Hs_hat_pre2;
        alpha2=(H_z_prime_s_hat_pre2)'*(y-Hz_prime)/norm(H_z_prime_s_hat_pre2)^2;
        z_prime=z_prime+alpha2*(z_prime-s_hat_pre2);
        %STEP 4
        [z_sort,z_index]=sort(abs(z_prime),'descend');
        A_index=z_index(1:r)';
        s_A=z_prime(A_index);
        s_hat_prime=zeros(m,1);
        s_hat_prime(A_index)=s_A;
        Hs_hat_prime=H(s_hat_prime);            %%%%%%%% Costly
        delta2_prime=norm(y-Hs_hat_prime)^2/N;
        loglikelihood_prime=-N/2*log(2*pi*delta2_prime)-N/2;

        %STEP 5
        if loglikelihood_prime>loglikelihood_bar
            s_hat=s_hat_prime;
            delta2=delta2_prime;
            Hs_hat=Hs_hat_prime;
            loglikelihood=loglikelihood_prime;
        else
            s_hat=s_bar;
            A_index=A_index_bar;
            delta2=delta2_bar;
            Hs_hat=Hs_bar;
            loglikelihood=loglikelihood_bar;
        end

%         gap=norm(s_hat-s_hat_pre)^2/m;
gap=norm(y - H(s_hat))/norm(y);


        if Visibility
            clc;
            display(['Iteration=',num2str(p),', gap=',num2str(gap),' (target=',num2str(thresh),')']);
        end
        if gap<thresh  || p>100 %修改收敛条件
            exit_flag=1;
        else
            s_hat_pre2=s_hat_pre;
            s_hat_pre=s_hat;
            Hs_hat_pre2=Hs_hat_pre;
            Hs_hat_pre=Hs_hat;
        end
    end

    Count=p;
    return;
end
 
 
 

