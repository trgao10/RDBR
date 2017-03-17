%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clearvars;%clc;
path(pathdef);
addpath(path, genpath('./utils'));
addpath(path, './data/');
scrsz = get(groot,'ScreenSize');

num_group = 2;
NoiseLevel = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% choose a regression method
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% BSFK (Knot-Free B-Spline)
% RegressionType = 'BSFK';
% RegressionOpt = struct('nknots',20,'knotremoval_factor',1.01,'order',3);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% SVM (Support Vector Machine)
%%%%% KernelScale=0.015 works even better. In any case, for unseen data one
%%%%% should cross-validate to pick an "optimal" KernelScale, probably once
%%%%% so for each iteration.
% RegressionType = 'SVM';
% RegressionOpt = struct('KernelFunction','gaussian','KernelScale',0.02,'Standardize',true);
% RegressionType = 'SVM';
% RegressionOpt = struct('KernelFunction','gaussian','KernelScale',0.015,'Standardize',true);
RegressionType = 'SVM';
RegressionOpt = struct('KernelFunction','gaussian','KernelScale',0.02,'Standardize',true);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% GPR (Gaussian Process Regression)
% RegressionType = 'GPR';
% RegressionOpt = struct('KernelFunction','squaredexponential');
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% MLS (Moving Least Squares)
% RegressionType = 'MLS';
% RegressionOpt = struct('tol', 1e-8, 'order', 3, 'Bandwidth',0.014);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% NW (Nadaraya-Watson Local Linear Regression)
% RegressionType = 'NW';
% RegressionOpt = struct('KNN',23);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% PB (Partition-Based Regression, or Nadaraya-Watson with square kernel)
% RegressionType = 'PB';
% RegressionOpt = struct('Bandwidth',0.01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% generate shape functions with Gaussian noise %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example 1: Fig11.m
% N = 2^12;
% ns = NoiseLevel*randn(1,N);
% 
% x = (0:N-1)/N;
% t = x;
% amp = 0.006;%0.01
% F1 = 150;
% F2 = 220;
% 
% sh1 = @(x) gen_shape2(x,3);
% sh2 = @(x) gen_shape2(x,2);
% 
% ins_freq = zeros(num_group,N);
% ins_amplt = zeros(num_group,N);
% ins_pre_phase = zeros(num_group,N);
% 
% xx = x + amp*sin(2*pi*x);
% ins_freq(1,:) = (1+amp*2*pi*cos(2*pi*x))*F1;
% ins_amplt(1,:) = 1+0.05*sin(2*pi*x);
% f1 = ins_amplt(1,:).*sh1(F1*xx); %%%%%ECG gen_shape(FF*xx,2)
% 
% yy = x + amp*cos(2*pi*x);
% ins_freq(2,:) = (1-amp*2*pi*sin(2*pi*x))*F2;
% ins_amplt(2,:) = 1+0.05*cos(2*pi*x);
% f2 = ins_amplt(2,:).*sh2(F2*yy);
% 
% ins_pre_phase(1,:) = (xx)*F1;
% ins_pre_phase(2,:) = (yy)*F2;
% 
% %%% only for validation purposes, not used in the RDBR algorithm
% fTrue = {f1,f2};
% shapeTrue = cell(1,2);
% shapeTrue{1} = @(x) sh1(x);
% shapeTrue{2} = @(x) sh2(x);
% 
% % %%% test permuting f1 and f2
% % ins_pre_phase = ins_pre_phase(2:-1:1,:);
% % ins_amplt = ins_amplt(end:-1:1,:);
% % ins_freq = ins_freq(end:-1:1);
% % shapeTrue = shapeTrue(end:-1:1);
% % tf = f1;
% % f1 = f2;
% % f2 = tf;
% 
% fff = f1 + f2  + ns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example 2: Fig10.m (Low Number of Periods)
N = 2^12;
ns = NoiseLevel*randn(1,N);

x = (0:N-1)/N;
t = x;
amp = 0.006;%0.01
F1 = 6;
F2 = 6;

sh1 = @(x) gen_shape(x,5);
sh2 = @(x) gen_shape(x,2);

ins_freq = zeros(num_group,N);
ins_amplt = zeros(num_group,N);
ins_pre_phase = zeros(num_group,N);

xx = x + amp*sin(2*pi*x);
ins_freq(1,:) = (1+amp*2*pi*cos(2*pi*x))*F1;
ins_amplt(1,:) = 1+0.05*sin(2*pi*x);
f1 = ins_amplt(1,:).*sh1(F1*xx); %%%%%ECG gen_shape(FF*xx,2)

yy = x + amp*cos(2*pi*x);
ins_freq(2,:) = (1-amp*2*pi*sin(2*pi*x))*F2;
ins_amplt(2,:) = 1+0.05*cos(2*pi*x);
f2 = ins_amplt(2,:).*sh2(F2*yy); %%%%%%

ins_pre_phase(1,:) = xx*F1;
ins_pre_phase(2,:) = yy*F2;

%%% only for validation purposes, not used in the RDBR algorithm
fTrue = {f1,f2};
shapeTrue = cell(1,2);
shapeTrue{1} = @(x) sh1(x);
shapeTrue{2} = @(x) sh2(x);

fff = f1 + f2  + ns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example 3: Fig6.m
% N = 2^16;
% ns = NoiseLevel*randn(1,N);
% 
% x = (0:N-1)/N;
% t = x;
% amp = 0.006;
% F1 = 6;
% F2 = 6;
% 
% sh1 = @(x) gen_shape(x,5);
% sh2 = @(x) gen_shape(x,2);
% 
% ins_freq = zeros(num_group,N);
% ins_amplt = zeros(num_group,N);
% ins_pre_phase = zeros(num_group,N);
% 
% xx = x + amp*sin(2*pi*x);
% ins_freq(1,:) = (1+amp*2*pi*cos(2*pi*x))*F1;
% ins_amplt(1,:) = 1+0.05*sin(2*pi*x);
% f1 = ins_amplt(1,:).*sh1(F1*xx); %%%%%ECG gen_shape(FF*xx,2)
% 
% yy = x + amp*cos(2*pi*x);
% ins_freq(2,:) = (1-amp*2*pi*sin(2*pi*x))*F2;
% ins_amplt(2,:) = 1+0.05*cos(2*pi*x);
% f2 = ins_amplt(2,:).*sh2(F2*yy);
% 
% ins_pre_phase(1,:) = (xx)*F1;
% ins_pre_phase(2,:) = (yy)*F2;

% %%% only for validation purposes, not used in the RDBR algorithm
% fTrue = {f1,f2};
% shapeTrue = cell(1,2);
% shapeTrue{1} = @(x) sh1(x);
% shapeTrue{2} = @(x) sh2(x);

% fff = f1 + f2  + ns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% apply RDBR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt = struct('maxIter', 10000,...
             'eps_error', 1e-6,...
             'eps_diff', 1e-6,...
             'RegressionType', RegressionType,...
             'RegressionOpt', RegressionOpt,...
             'visFlag', true);
% opt = struct('maxIter', 40,...
%              'eps_error', 1e-13,...
%              'eps_diff', 1e-13,...
%              'RegressionType', RegressionType,...
%              'RegressionOpt', RegressionOpt);

[shape,compos,rms,numIter] = RDBR(fff,ins_amplt,ins_pre_phase,opt,shapeTrue);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% post processing to get better components %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shapepp = cell(1,num_group);
for cntGroup=1:num_group
    shapepp{cntGroup} = spline(linspace(0,1,length(shape{cntGroup}))',shape{cntGroup});
    compos{cntGroup} = ppval(shapepp{cntGroup},mod(ins_pre_phase(cntGroup,:),1));
    cnt = cntGroup; cnty = cnt - sign(cnt-1.5);
    
    if cnt ==1
        fff = f2;
    else
        fff = f1;
    end
    
    compos{cntGroup} = compos{cntGroup}.* ins_amplt(cntGroup,:)...
        *(fff/(compos{cntGroup}.* ins_amplt(cntGroup,:)));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% plotting results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig_opts=struct('Position',[scrsz(4)/2,scrsz(4)-(scrsz(4)/2),scrsz(4)/2*2.5,scrsz(4)/3]);
% ind = 1:800;
% fig = figure;plot(x,fff);
% set(fig,fig_opts);axis tight
% xlabel('Time (Second)');title('A superposition of general modes');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',scrsz);
for j=1:num_group
    [trans_est_shape,min_error] = shape_phase_trans(shape{j}.',shapeTrue{j}(linspace(0,1,1000)));
    L = length(trans_est_shape);
    gd = 0:1/L:(1-1/L);
    
    subplot(1,num_group,j);
    plot(gd,trans_est_shape);
    hold on;
    plot(gd,shapeTrue{j}(linspace(0,1,1000)));
    hold off;
    legend('Est','True');
    title(sprintf('Shape %d, Residue = %.4f', j, min_error),'Interpreter','latex');
    axis square;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);
    set(b.YLabel, 'FontSize', 16);
    set(b.ZLabel, 'FontSize', 16);
    set(b.Title, 'FontSize', 16);
end

[~,h] = suplabel(RegressionType,'t');
set(h,'FontSize',30,'Interpreter','Latex');

savefig([RegressionType, '_N=' num2str(F1) '_shapes.fig']);
% close all;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;

subplot(1,2,1);
plot(rms);

subplot(1,2,2);
% delta = rms(1:end-1)-rms(2:end);
% pos = find(abs(delta)<1e-4);
% if ~isempty(pos)
%     pos = pos(1);
% else
%     pos = length(rms)-2;
% end
% pos = min(pos, length(eta))';
% plot(eta(1:pos),'bo-');
mu = log(abs(rms(1:end-1)-rms(2:end)));
eta = mu(1:end-1)-mu(2:end);
plot(eta,'bo-');
xlabel('Number of Iterations $j$','Interpreter','Latex');
ylabel('$\eta_j$','Interpreter','Latex');
title([RegressionType, ', $N=' num2str(F1) '$'],'Interpreter','Latex');

savefig([RegressionType, '_N=' num2str(F1) '_rms.fig']);
% close all;
