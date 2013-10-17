function  PhaseTrans


clear all

mytime = datestr(now);
disp(sprintf('start at [%s]',mytime));


dimensions = 400;  % signal length
numtrials = 100;   
numgradations = 16;
numiterations = 100;	% one stopping criterion
tolerance = 1e-6;	% another stopping criterion           

delta = linspace(0.05, 0.5,numgradations);  % indeterminacy m/n

rho = linspace(0.05, 0.5,numgradations); % sparsity k/m


numalgorithms = 4; % 算法个数

successrecord = zeros(numalgorithms, numgradations,numgradations,numtrials);
computationtimes = zeros(numalgorithms, numgradations,numgradations,numtrials);

Phiseedcounter = 100;
Xseedcounter = 1000;


% % uncomment the above if you wish to generate new simulations
% load data_sparsity.mat

% for kk = 1:numalgorithms
    
measurementcounter = 1;    
for nummeasurements = ceil(delta*dimensions)


SPfail = 0;
DOREfail = 0;
EMTPfail = 0;
EMTPbetafail = 0;

		sparsitycounter = 1;
		for sparsity = ceil(rho*nummeasurements)
			% sample a sensing matrix from spherical ensemble
            Phi = MatrixEnsemble(nummeasurements, dimensions, 'USE');
            % ensemble:
            % 'USE' - Uniform spherical ensemble
            % 'RSE' - Random signs ensemble
            % 'RST' - Partial RST (Real Fourier) ensemble
            % 'Hadamard' - Partial Hadamard ensemble			

			% test recovery algorithms for numtrials
            tstart = tic;            
			for ii=1:numtrials
				x = buildSparseSignal(dimensions,sparsity, 'normal', Xseedcounter);
                % distributiontype
                % 'normal' - 
                % 'laplacian' - 
                % 'powerlaw' - 
                % 'bimodalgaussian' - 
                % 'bimodalrayleigh' - 
                % 'bimodaluniform' - 
                % 'uniform' - 
                % 'bernoulli' -                 
				Xseedcounter = Xseedcounter + 1;
				y = Phi*x; % sense 
                                    
                H           = @(z) Phi*z;
                Ht          = @(z) Phi'*z;
                invHHt      = @(z) invAAt*z;
                invAAt      = inv(Phi*Phi');
           
                
                  if ~SPfail
					fprintf('SP ... \t');
					dstart = tic;
                    options.tol = tolerance;
                    Rec                         = SP(Phi, y, sparsity, options);
                    xhat = Rec.x_hat;
                    % 记录第1个算法结果
					computationtimes(1, sparsitycounter, measurementcounter,ii) = toc(dstart);
					successrecord(1, sparsitycounter, measurementcounter,ii) = success(x,xhat);
                  end
                  
                  
                  if ~DOREfail
					fprintf('DORE ... \t');
					dstart = tic;
                    
                    H           = @(z) Phi*z;
                    Ht          = @(z) Phi'*z;
                    
                    xhat = DORE(H, Ht, y, sparsity, 'Thresh',tolerance, 'visibility',0,'IsHOrthonormal',1);
                    % 记录第2个算法结果
					computationtimes(2, sparsitycounter, measurementcounter,ii) = toc(dstart);
					successrecord(2, sparsitycounter, measurementcounter,ii) = success(x,xhat);                              
                  end


                  if ~EMTPfail
					fprintf('EMTP-k ... \t');
					dstart = tic;
                    xhat                         = EMTP(Phi, y, sparsity, 1, 'Tolerance',tolerance);
                    % 记录第3个算法结果
					computationtimes(3, sparsitycounter, measurementcounter,ii) = toc(dstart);
					successrecord(3, sparsitycounter, measurementcounter,ii) = success(x,xhat);
                  end
                  
                  
                  
                  if ~EMTPbetafail
					fprintf('EMTP-beta ... \t');
					dstart = tic;
                    xhat                         = EMTPbeta(Phi, y, 0.8, 'Tolerance',tolerance);
                    % 记录第4个算法结果
					computationtimes(4, sparsitycounter, measurementcounter,ii) = toc(dstart);
					successrecord(4, sparsitycounter, measurementcounter,ii) = success(x,xhat);                                    
                  end
                  
                 
                 fprintf('\n');                  
            end          
            
			disp(sprintf('dimensions=%d: sparsity = %2.4f, ind. = %2.4f, %2.2f seconds ...', ...
			dimensions, sparsity/nummeasurements, ...
			nummeasurements/dimensions, toc(tstart)));    
        
        
            if ~sum(successrecord(1, sparsitycounter,measurementcounter,:))
				SPfail = 1;
            end


            if ~sum(successrecord(2, sparsitycounter,measurementcounter,:))
				DOREfail = 1;            
            end
            
            
            if ~sum(successrecord(3, sparsitycounter,measurementcounter,:))
				EMTPfail = 1;      
            end
            
      
            if ~sum(successrecord(4, sparsitycounter,measurementcounter,:))
				EMTPbetafail = 1;
            end
            
            
      
			sparsitycounter = sparsitycounter + 1;            
        end        
		measurementcounter = measurementcounter + 1;   
end

% end

tt = datevec(now);
str = num2str(tt(6));
filename = strcat('PhaseTrans', '_', str, '.mat'); 
save(filename)
% save(filename, 'successrecord', 'computationtimes')


%% 画图

 load PhaseTrans

%%
%  统计结果 
  stats1 = zeros(numgradations,numgradations);
  compstats = zeros(numgradations,numgradations);
  for indeterminacycounter = 1:numgradations			% indeterminacy
    for sparsitycounter = 1:numgradations	% sparsity
      for kk=1:numtrials
        stats1(sparsitycounter,indeterminacycounter) = ...
          stats1(sparsitycounter,indeterminacycounter) + ...
          sum(successrecord(1, sparsitycounter,indeterminacycounter,kk));
      end
    end
  end

  stats2 = zeros(numgradations,numgradations);
  compstats = zeros(numgradations,numgradations);
  for indeterminacycounter = 1:numgradations			% indeterminacy
    for sparsitycounter = 1:numgradations	% sparsity
      for kk=1:numtrials
        stats2(sparsitycounter,indeterminacycounter) = ...
          stats2(sparsitycounter,indeterminacycounter) + ...
          sum(successrecord(2, sparsitycounter,indeterminacycounter,kk));
      end
    end
  end

  stats3 = zeros(numgradations,numgradations);
  compstats = zeros(numgradations,numgradations);
  for indeterminacycounter = 1:numgradations			% indeterminacy
    for sparsitycounter = 1:numgradations	% sparsity
      for kk=1:numtrials
        stats3(sparsitycounter,indeterminacycounter) = ...
          stats3(sparsitycounter,indeterminacycounter) + ...
          sum(successrecord(3, sparsitycounter,indeterminacycounter,kk));
      end
    end
  end
  
 
  stats4 = zeros(numgradations,numgradations);
  compstats = zeros(numgradations,numgradations);
  for indeterminacycounter = 1:numgradations			% indeterminacy
    for sparsitycounter = 1:numgradations	% sparsity
      for kk=1:numtrials
        stats4(sparsitycounter,indeterminacycounter) = ...
          stats4(sparsitycounter,indeterminacycounter) + ...
          sum(successrecord(4, sparsitycounter,indeterminacycounter,kk));
      end
    end
  end
  

% %  画图 xlabel('Sparsity k/m');
figure;
% for ii=1:3:16
for ii=1:5:16
plot(rho,stats2(:,ii)./numtrials,'r:', 'LineWidth', 8*(20-ii)/20); 
hold on 
plot(rho,stats1(:,ii)./numtrials,'b--', 'LineWidth', 8*(20-ii)/20); 
% hold on 
% plot(rho,stats3(:,ii)./numtrials,'r-', 'LineWidth', 8*(20-ii)/20); 
hold on 
plot(rho,stats4(:,ii)./numtrials,'k-', 'LineWidth', 8*(20-ii)/20); 

% plot(rho,stats1(:,ii)./numtrials,'Color',sqrt((18-ii)/21)*ones(3,1),'LineWidth',4*(20-ii)/16,'LineStyle','--');

% plot(rho,stats2(:,ii)./numtrials,'Color',sqrt((18-ii)/21)*ones(3,1),'LineWidth',4*(20-ii)/16,'LineStyle','-.');

% plot(rho,stats3(:,ii)./numtrials,'Color',sqrt((18-ii)/21)*ones(3,1),'LineWidth',4*(20-ii)/16,'LineStyle','--');

% plot(rho,stats4(:,ii)./numtrials,'Color',sqrt((18-ii)/21)*ones(3,1), ...
%     'LineWidth',4*(20-ii)/16,'LineStyle','-');

% 	texth = gtext(sprintf('delta=%2.2f',delta(ii)));
	texth = gtext(sprintf('%2.2f',delta(ii)));
	set(texth,'FontSize',14, 'BackgroundColor','w');
    
% axis([0.05 0.5 -0.01 1.01]);
% % xlabel(['Sparsity k/m, ' sprintf('delta=%2.2f',delta(ii))]);
% xlabel(['Sparsity k/m']);
% ylabel('Probability of Exact Recovery');
% grid on;
% legend('IHT', 'IDE', 'GDE');
% str = num2str(delta(ii));
% fn = strcat('Fig_Comparison_IHT_IDE_GDE_', str, '.fig'); 
% saveas(gcf, fn) 
end
axis([0.05 0.5 -0.01 1.01]);
% xlabel(['Sparsity k/m, ' sprintf('delta=%2.2f',delta(ii))]);
xlabel(['Sparsity k/m']);
ylabel('Probability of Exact Recovery');
grid on;
% legend('SP', 'DORE', 'EMTP-k', 'EMTP-\beta=0.8');
legend('DORE', 'SP', 'EMTP-\beta=0.8',3);
% fn = strcat('Fig_PhaseTrans_DORE',  '.fig'); 
% saveas(gcf, fn) 
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
print(gcf,'-depsc2',['PhaseTrans_recovery' '.eps']);



%%
% %画相变

phase = NaN*ones(numalgorithms, numgradations);


for alg =[1 2 4]%1:numalgorithms
	stats = zeros(numgradations,numgradations);
	for indeterminacycounter = 1:numgradations % indeterminacy
        for sparsitycounter = 1:numgradations % sparsity
            for kk=1:numtrials
                stats(sparsitycounter,indeterminacycounter) = ...
				stats(sparsitycounter,indeterminacycounter) + ...
				sum(successrecord(alg, sparsitycounter,indeterminacycounter,kk)); % alg 
			end
		end
    end
    
	stats = stats./numtrials;
	% for each indeterminacy where there are non-zero values,
	% find approximate sparsity at which 50% recovery occurs
	for indeterminacycounter = 1:max(find(stats(1,:)>0))
		I = [];
		upsamplefactor = 10;
		while isempty(I)
			sparseline = interp(stats(:,indeterminacycounter),upsamplefactor);
			I = find(sparseline > 0.49 & sparseline < 0.51);
			upsamplefactor = upsamplefactor + 1;
		end
			rhoupsampled = interp(rho,upsamplefactor-1);
			[temp,J] = min(abs(sparseline(I)-0.5));
			phase(alg, indeterminacycounter) = rhoupsampled(I(J)); % alg 
	end

end

figure;
plot(delta(1:end),phase(2, 1:end),'r:', 'LineWidth',4); % alg 
hold on
plot(delta(1:end),phase(1, 1:end),'b--', 'LineWidth',4); % alg 
% hold on
% plot(delta(1:end),phase(3, 1:end),'r-', 'LineWidth',4); % alg 
hold on
plot(delta(1:end),phase(4, 1:end),'k-', 'LineWidth',4); % alg 

axis([0.05 0.55 0 0.6]);
ylabel('Sparsity k/m');
xlabel('Indeterminacy m/n');
grid on;
axis([0.05 0.5 0.05 0.6]);
% legend('SP', 'DORE', 'EMTP-k', 'EMTP-\beta=0.8');
legend('DORE', 'SP', 'EMTP-\beta=0.8');
% fn = strcat('Fig_PhaseTrans_DORE', '.fig'); 
% saveas(gcf, fn) 
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
print(gcf,'-depsc2',['PhaseTrans' '.eps']);


disp(sprintf('start at [%s]',mytime));
disp(sprintf('over at  [%s]',datestr(now)));