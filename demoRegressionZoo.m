%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clearvars;%clc;
path(pathdef);
addpath(path, genpath('./utils'));
addpath(path, '../data/');
scrsz = get(groot,'ScreenSize');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% generate shape functions with Gaussian noise %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 8192/4;
x = (0:N-1)/N;

%%% clean shape function over the unit interval
Y = gen_shape2(x,2);

%%% overlaying two shape functions ---- many more periods
amp = 0.01;
F1 = 60;
F2 = 90;
xx = x + amp*sin(2*pi*x);
insPhase = F1*xx;
am = 1+0.05*sin(4*pi*x);
f1 = am.*gen_shape2(insPhase,2);
yy = x + amp*cos(2*pi*x);
bm = 1+0.1*sin(2*pi*x);
f2 = bm.*gen_shape2(F2*yy,3);
NM = 0;
ns = NM*randn(1,N);
fff = f1 + f2 + ns;

%%% overlayed two shape functions, folded over to the unit interval
%%% YWrap is shown in subplot(2,4,1)
insAmplitude = am;
phi2 = mod(insPhase, 1);
X = phi2(:);
YWrap = fff./insAmplitude;
YWrap = YWrap(:);

%%% re-sampled points for regression
x_regression = 0:1/1000:1;
shapeTrue = gen_shape2(x_regression,2);

%%% BSFK (B-Spline with Free Knots) regression from (X,YWrap)
%%% shapeEstBSFK is shown in subplot(2,4,2)
options = struct('animation', 0, 'knotremoval_factor',1.0001);
[pp] = BSFK(X',YWrap,3,20,[],options);
shapeEstBSFK = ppval(pp,x_regression);

%%% SVM (Supprting Vector Machine) regression from (X,YWrap)
%%% shapeEstSVM is shown in subplot(2,4,3)
SVMModel = fitrsvm(X,YWrap','KernelFunction','gaussian','KernelScale',0.02,'Standardize',true);
shapeEstSVM = predict(SVMModel, x_regression')';

%%% GP (Gaussian Process) regression from (X,YWrap)
%%% shapeEstGP is shown in subplot(2,4,4)
sigma0 = Chi2Estimation(YWrap);
GPRModel = fitrgp(X,YWrap','KernelFunction','squaredexponential','Sigma',sigma0);
shapeEstGPR = predict(GPRModel, x_regression')';

%%% Random Forest regression from (X,YWrap)
%%% shapeEstRF is shwo in subplot(2,4,5)
% RFModel = fitrensemble(X,YWrap','OptimizeHyperparameters','auto',...
%     'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
%     'expected-improvement-plus'));
% RFModel = fitrensemble(X,YWrap','Method','LSBoost',...
%     'NumLearningCycles',498,'LearnRate',0.011042);
RFModel = fitrensemble(X,YWrap','Method','LSBoost',...
    'NumLearningCycles',300,'Learners',templateTree('MaxNumSplits',15),'LearnRate',0.01);
shapeEstRF = predict(RFModel, x_regression')';

% err_opt = inf;
% for lr=0.009:0.0001:0.01
%     for mns = 10:1:25
%         RFModel = fitrensemble(X,YWrap','Method','LSBoost',...
%             'NumLearningCycles',500,'Learners',templateTree('MaxNumSplits',mns),'LearnRate',lr);
%         shapeEstRF = predict(RFModel, x_regression')';
%         err = norm(shapeEstRF-shapeTrue)/length(x_regression);
%         if err < err_opt
%             err_opt = err;
%             lr_best = lr;
%             mns_best = mns;
%         end
%     end
%     fprintf('current err_opt: %f\n',err_opt);
% end

%%% MLS (Moving Least Squares) regression from (X,YWrap)
%%% shapeEstMLS is shown in subplot(2,4,6)
bandwidth = 0.014;
W = exp(-getEuclideanDistance(x_regression,X').^2/bandwidth^2);
ord = 3;
tol = 1e-8;
coeff = MLS(X,x_regression',W,ord,tol);
shapeEstMLS = (coeff*YWrap)';

%%% Local Linear Kernel (Nadaraya-Watson) regression from (X,YWrap)
%%% https://www.mathworks.com/matlabcentral/fileexchange/19564-local-linear-kernel-regression
%%% https://www.mathworks.com/matlabcentral/fileexchange/35316-kernel-regression-with-variable-window-width
%%% shapeEstNW is shown in subplot(2,4,7)
% bandwidth = 0.01;
% r = ksrlin(X,YWrap,bandwidth,length(x_regression));
% r = ksrlin_vw(X,YWrap,23,length(x_regression));
r = ksrlin_vw(X,YWrap,23,x_regression);
shapeEstNW = r.f;

%%% Paritition-Based Regression (a silly instance of kernel regression)
%%% shapeEstPBR is shown in subplot(2,4,8)
% err_opt = inf;
% for bw=0.005:0.001:0.025
%     r = ksr_characteristic(X,YWrap,bw,length(x_regression));
%     shapeEstPBR = r.f;
%     shapeEstPBR(isnan(shapeEstPBR)) = 0;
%     
%     err = norm(shapeEstPBR-shapeTrue)/length(x_regression);
%     if err < err_opt
%         err_opt = err;
%         bw_best = bw;
%     end
%     fprintf('current err_opt: %f\n',err_opt);
% end
% r = ksr_characteristic(X,YWrap,0.01,length(x_regression));
r = ksr_characteristic(X,YWrap,0.01,x_regression);
shapeEstPBR = r.f;
shapeEstPBR(isnan(shapeEstPBR)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);

% subplot(2,3,1);
% plot(x,Y,'.','LineWidth',4);
% set(gca, 'FontSize', 16);
% axis square
% title('ground truth shape');
% b = get(gca);
% set(b.XLabel, 'FontSize', 16);
% set(b.YLabel, 'FontSize', 16);
% set(b.ZLabel, 'FontSize', 16);
% set(b.Title, 'FontSize', 16);

% subplot(2,3,2);
subplot(2,4,1);
plot(X,YWrap,'.','LineWidth',4);
set(gca, 'FontSize', 16);
axis square
title('folded $f(t)=f_1(t)+f_2(t)$','Interpreter','latex');
b = get(gca);
set(b.XLabel, 'FontSize', 16);
set(b.YLabel, 'FontSize', 16);
set(b.ZLabel, 'FontSize', 16);
set(b.Title, 'FontSize', 16);

subplot(2,4,2);
plot(x_regression,shapeTrue,'r','LineWidth',2);
hold on;
plot(x_regression,shapeEstBSFK,'b','LineWidth',2);
scatter(X,YWrap,5,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
legend({'true','est'},'Location','NorthEast');
axis square;
axis([0,1,-1,6]);
set(gca, 'FontSize', 16);
title(sprintf('BSFK, rms=%f',norm(shapeEstBSFK-shapeTrue)/length(x_regression)), 'Interpreter','latex');
b = get(gca);
set(b.XLabel, 'FontSize', 16);
set(b.YLabel, 'FontSize', 16);
set(b.ZLabel, 'FontSize', 16);
set(b.Title, 'FontSize', 16);

subplot(2,4,3);
plot(x_regression,shapeTrue,'r','LineWidth',2);
hold on;
plot(x_regression,shapeEstSVM,'b','LineWidth',2);
scatter(X,YWrap,5,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
legend({'true','est'},'Location','NorthEast');
axis square;
axis([0,1,-1,6]);
set(gca, 'FontSize', 16);
title(sprintf('SVM, rms=%f',norm(shapeEstSVM-shapeTrue)/length(x_regression)), 'Interpreter','latex');
b = get(gca);
set(b.XLabel, 'FontSize', 16);
set(b.YLabel, 'FontSize', 16);
set(b.ZLabel, 'FontSize', 16);
set(b.Title, 'FontSize', 16);

subplot(2,4,4);
plot(x_regression,shapeTrue,'r','LineWidth',2);
hold on;
plot(x_regression,shapeEstGPR,'b','LineWidth',2);
scatter(X,YWrap,5,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
legend({'true','est'},'Location','NorthEast');
axis square;
axis([0,1,-1,6]);
set(gca, 'FontSize', 16);
title(sprintf('GPR, rms=%f',norm(shapeEstGPR-shapeTrue)/length(x_regression)), 'Interpreter','latex');
b = get(gca);
set(b.XLabel, 'FontSize', 16);
set(b.YLabel, 'FontSize', 16);
set(b.ZLabel, 'FontSize', 16);
set(b.Title, 'FontSize', 16);

subplot(2,4,5);
plot(x_regression,shapeTrue,'r','LineWidth',2);
hold on;
plot(x_regression,shapeEstRF,'b','LineWidth',2);
scatter(X,YWrap,5,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
legend({'true','est'},'Location','NorthEast');
axis square;
axis([0,1,-1,6]);
set(gca, 'FontSize', 16);
title(sprintf('RF, rms=%f',norm(shapeEstRF-shapeTrue)/length(x_regression)), 'Interpreter','latex');
b = get(gca);
set(b.XLabel, 'FontSize', 16);
set(b.YLabel, 'FontSize', 16);
set(b.ZLabel, 'FontSize', 16);
set(b.Title, 'FontSize', 16);

subplot(2,4,6);
plot(x_regression,shapeTrue,'r','LineWidth',2);
hold on;
plot(x_regression,shapeEstMLS,'b','LineWidth',2);
scatter(X,YWrap,5,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
legend({'true','est'},'Location','NorthEast');
axis square;
axis([0,1,-1,6]);
set(gca, 'FontSize', 16);
title(sprintf('MLS, rms=%f',norm(shapeEstMLS-shapeTrue)/length(x_regression)), 'Interpreter','latex');
b = get(gca);
set(b.XLabel, 'FontSize', 16);
set(b.YLabel, 'FontSize', 16);
set(b.ZLabel, 'FontSize', 16);
set(b.Title, 'FontSize', 16);

subplot(2,4,7);
plot(x_regression,shapeTrue,'r','LineWidth',2);
hold on;
plot(x_regression,shapeEstNW,'b','LineWidth',2);
scatter(X,YWrap,5,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
legend({'true','est'},'Location','NorthEast');
axis square;
axis([0,1,-1,6]);
set(gca, 'FontSize', 16);
title(sprintf('Nadaraya-Watson, rms=%f',norm(shapeEstNW-shapeTrue)/length(x_regression)), 'Interpreter','latex');
b = get(gca);
set(b.XLabel, 'FontSize', 16);
set(b.YLabel, 'FontSize', 16);
set(b.ZLabel, 'FontSize', 16);
set(b.Title, 'FontSize', 16);

subplot(2,4,8);
plot(x_regression,shapeTrue,'r','LineWidth',2);
hold on;
plot(x_regression,shapeEstPBR,'b','LineWidth',2);
scatter(X,YWrap,5,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
legend({'true','est'},'Location','NorthEast');
axis square;
axis([0,1,-1,6]);
set(gca, 'FontSize', 16);
title(sprintf('Partition-Based, rms=%f',norm(shapeEstPBR-shapeTrue)/length(x_regression)), 'Interpreter','latex');
b = get(gca);
set(b.XLabel, 'FontSize', 16);
set(b.YLabel, 'FontSize', 16);
set(b.ZLabel, 'FontSize', 16);
set(b.Title, 'FontSize', 16);

