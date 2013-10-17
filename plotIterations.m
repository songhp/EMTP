% plotIterations

clear all
clc
close all

mytime = datestr(now);
disp(sprintf('start at [%s]',mytime));


path(path, './subfunctions');
% path(path,'./Algorithms');

AlgName = {'DORE','EMTP-k', 'EMTP-\beta=0.8',};
% delta = 0.4;
% rho = 0.2;
delta = 0.3;
rho = 0.3;
NN = 200:400:2200;
NumITER = 100;
thresh  = 1e-6;           %Convergence tolerance for Hard thresholding methods


ERR_DORE    = zeros(length(NN),NumITER);
t_DORE      = zeros(length(NN),NumITER);

ERR_EMTP    = zeros(length(NN),NumITER);
t_EMTP      = zeros(length(NN),NumITER);

ERR_EMTPbeta    = zeros(length(NN),NumITER);
t_EMTPbeta      = zeros(length(NN),NumITER);

for ITER = 1:NumITER
    knum = 0;
    for n = NN
        knum = knum+1;
        display(['ITERATION: ' num2str(ITER) '; n: ' num2str(n)])
        
        m = ceil(delta * n);
        k = ceil(rho * m);
        P =randn(m,n);
        %P =round(rand(M,N))-0.5;
        for nn=1:n
            P(:,nn)=P(:,nn)/norm(P(:,nn));
        end
        
        s_TRUE      = zeros(n,1);
        s_TRUE(1:k) = randn(k,1);
        %s_TRUE(1:k) = ones(k,1);
        
        y           = P*s_TRUE;
        
        invPPt      = inv(P*P');
        H           = @(z) P*z;
        Ht          = @(z) P'*z;
        invHHt      = @(z) invPPt*z;

        %Reconstruction
        
                
        %display('solveing by DORE')
        tic;
        [s_DORE,A_index_DORE,Count_DORE]       = DORE(H,Ht,y,k,'Thresh',thresh,'Inverse',invHHt,'visibility',0,'IsHOrthonormal',0);
%         [s_DORE,A_index_DORE,Count_DORE] =  DORE(P, [], y, k, 'Thresh',thresh,'visibility',0,'IsHOrthonormal',0);
%         t_DORE(knum,ITER)                      = toc;
        toc
        t_DORE(knum,ITER)                      = Count_DORE;
        ERR_DORE(knum,ITER)                    = norm(s_DORE-s_TRUE)/norm(s_TRUE);

        
   
        tic;        
        [s_EMTP, Out]                       = EMTP(P, y, k, 1, 'Tolerance',thresh);
        toc
%         t_EMTP(knum,ITER)                   = toc;
        t_EMTP(knum,ITER)                   = Out.iter;
        ERR_EMTP(knum,ITER)                 = norm(s_EMTP-s_TRUE)/norm(s_TRUE);
        
        
        tic;        
        [s_EMTPbeta, Out]                       = EMTPbeta(P, y, 0.8, 'Tolerance',thresh);
        toc
%         t_EMTP(knum,ITER)                   = toc;
        t_EMTPbeta(knum,ITER)                   = Out.iter;
        ERR_EMTPbeta(knum,ITER)                 = norm(s_EMTPbeta-s_TRUE)/norm(s_TRUE);
        
        
    end
end

tt = datevec(now);
save(['plotIterations', '_',num2str(tt(6)),'.mat']);



figure
plot(NN,(mean(t_DORE,2)),'r:', 'LineWidth',2,'MarkerSize',10)
hold on
plot(NN,(mean(t_EMTP,2)),'b--', 'LineWidth',2,'MarkerSize',10)
plot(NN,(mean(t_EMTPbeta,2)),'k-', 'LineWidth',2,'MarkerSize',10)
legend(AlgName{1}, AlgName{2}, AlgName{3}, 2)
% title('m/n=0.4, k/m=0.2, n=200:400:2200', 'fontsize',14);
title('m/n=0.3, k/m=0.3, n=200:400:2200', 'fontsize',14);
ylabel('Number of Iterations')
%title('Speed')
% ylabel('Computation time in seconds')
xlabel('n')
grid on
hold off
% fn = strcat('Fig_SNR_fix_M=400_Iter', '.fig'); 
% saveas(gcf, fn) 
% %
% set(gcf,'color','none');
% set(gca,'color','none');
% set(gcf,'InvertHardCopy','off');
% print -depsc2 Fig_SNR_fix_M=400_Iter.eps
%
% fn = strcat('Fig_SNR_fix_M=400', '.fig'); 
% saveas(gcf, fn) 
%
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
% print -depsc2 plotIterations42.eps
print -depsc2 plotIterations33.eps

disp(sprintf('start at [%s]',mytime));
disp(sprintf('over at  [%s]',datestr(now)));