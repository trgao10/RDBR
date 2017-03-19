function [shapes,compos,rms,numIter] = RDBR(fff,insAmplitude,insPhase,opt,shapeTrue)
%RDBR Recursive Diffeomorphism-Based Regression for Empirical Mode
%Decomposition (EMD)
%   This code adds various types of regression methods to the original
%   implementation of RDBR published with the following preprint:
%   
%   [1] Recursive Diffeomorphism-Based Regression for Shape Functions,
%       https://arxiv.org/abs/1610.03819
%
%   Tingran Gao (trgao10@math.duke.edu)
%   Last Modified: March 14, 2017
%   

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% parse inputs
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if nargin < 5
    shapeTrue = {};
end
if min(size(insAmplitude))~=min(size(insPhase))
    error('Inconsistent Number of Components!');
end

numGroup = min(size(insAmplitude));
maxIter = getoptions(opt, 'maxIter', 10000);
eps_error = getoptions(opt, 'eps_error', 1e-6);
eps_diff = getoptions(opt, 'eps_diff', 1e-6);
visFlag = getoptions(opt, 'visFlag', false);
xRegr = linspace(0,1,1000);

shapes = cell(1,numGroup);
compos = cell(1,numGroup);

RegressionParams = ParseRegressionOpt(opt.RegressionType,opt.RegressionOpt);

currTimeStamp = strrep(datestr(datetime('now')),' ', '-');
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% main iteration
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if visFlag
    F1Est = round(insPhase(1,end));
    F2Est = round(insPhase(2,end));
    scrsz = get(groot,'ScreenSize');
    figure('Position', scrsz);
    videoObj = VideoWriter(['videos/' opt.RegressionType, '_N1=' num2str(F1Est) '_N2=' num2str(F2Est) '_' currTimeStamp '.avi']);
    videoObj.FrameRate = 2;
    open(videoObj);
end

numIter = 0;
r_n = fff;
compoMag = 0;
errReg = Inf;
rms = [];
while numIter<maxIter
    numIter = numIter + 1;
    errRegPrev = errReg;
    
    fn = r_n;
    for cntGroup = 1:numGroup
        foldedPhase = mod(insPhase(cntGroup,:), 1);
        Y = fn./insAmplitude(cntGroup,:);
        
        [sortFoldedPhase,ia] = sort(foldedPhase);
        [shapeInLoop,YEst] = RegressionWrapper(sortFoldedPhase,Y(ia),xRegr,RegressionParams);
        [~,ic] = sort(ia);
        YEst = YEst(ic);
        
%         [shapeInLoop,YEst] = RegressionWrapper(foldedPhase,Y,xRegr,RegressionParams);
%         RegressionParams.KernelScale = RegressionParams.KernelScale * 0.995;
%         RegressionParams.KernelScale = RegressionParams.KernelScale * 1.005;
%         RegressionParams.KernelScale = RegressionParams.KernelScale * 1.01;
        
        if visFlag
            subplot(2,numGroup+1,cntGroup);
            [sortFoldedPhase,ia] = sort(foldedPhase);
            plot(sortFoldedPhase,YEst(ia),'rx','MarkerSize',5);
            hold on
            plot(xRegr,shapeInLoop,'b-','LineWidth',1);
            legend('YEst','shapeInLoop');
            scatter(foldedPhase,Y,5,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
            hold off
            title(sprintf('Shape %d, Iter %d', cntGroup, numIter),'Interpreter','Latex');
            drawnow
        end

        shapeOffset = mean(shapeInLoop);
        shapeInLoop = shapeInLoop - shapeOffset;
        %%% Can we use a different shapeOffset for component estimation?
        compoInLoop = (YEst-mean(YEst)).*insAmplitude(cntGroup,:);
%         compoInLoop = (YEst-shapeOffset).*insAmplitude(cntGroup,:);
                
        if numIter == 1
            shapes{cntGroup} = shapeInLoop;
            compos{cntGroup} = compoInLoop;
        else
            shapes{cntGroup} = shapes{cntGroup} + shapeInLoop;
            compos{cntGroup} = compos{cntGroup} + compoInLoop;
        end
        
        if visFlag
            subplot(2,numGroup+1,cntGroup+numGroup+1);
            if ~isempty(shapeTrue)
                plot(xRegr,shapeTrue{cntGroup}(xRegr),'ro');
                hold on
            end
            plot(xRegr,shapes{cntGroup},'b-','LineWidth',2);
            if ~isempty(shapeTrue)
                legend('True','Est');
                hold off
            end
            title(sprintf('Shape %d, Iter %d', cntGroup, numIter),'Interpreter','Latex');
            drawnow
        end
        
        compoMag = max(compoMag,norm(compoInLoop)/sqrt(length(compoInLoop)));
        r_n  = r_n - compoInLoop;
    end
    
    errReg = norm(r_n)/sqrt(length(r_n));
    fprintf('Iteration %03d, errReg = %f\n', numIter, errReg);
    rms = [rms errReg];
    
    if visFlag
        subplot(2,numGroup+1,numGroup+1);
        plot(rms,'b-o');
        xlabel('Number of Iterations $j$','Interpreter','Latex');
        ylabel('$\|r_j\|_2$','Interpreter','Latex');
        title([opt.RegressionType, ', $N_1=' num2str(F1Est) ', N_2=' num2str(F2Est) '$'],'Interpreter','Latex');
        
        subplot(2,numGroup+1,2*numGroup+2);
        mu = log(abs(rms(1:end-1)-rms(2:end)));
        eta = mu(1:end-1)-mu(2:end);
        plot(eta,'bo-');
        xlabel('Number of Iterations $j$','Interpreter','Latex');
        ylabel('$\eta_j$','Interpreter','Latex');
        title([opt.RegressionType, ', $N_1=' num2str(F1Est) ', N_2=' num2str(F2Est) '$'],'Interpreter','Latex');        
        drawnow
        
        writeVideo(videoObj, getframe(gcf));
    end

    if (errReg<eps_error) || (compoMag<eps_error) || (abs(errRegPrev-errReg)<eps_diff)
%     if (errReg<eps_error) || (compoMag<eps_error) || (abs(errRegPrev-errReg)<eps_diff) || (errReg>errRegPrev)
        break;
    end
    
%     if strcmpi(RegressionParams.Type,'SVM') && (errReg>errRegPrev)
%         RegressionParams.Standardize = false;
%     end
end

savefig(['results/' opt.RegressionType, '_N1=' num2str(F1Est) '_N2=' num2str(F2Est) '_' currTimeStamp '.fig']);

if visFlag
    close(videoObj);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RegressionParams = ParseRegressionOpt(Type,Opt)

RegressionParams.Type = Type;

if strcmpi(Type, 'BSFK')
    RegressionParams.nknots = Opt.nknots;
    RegressionParams.order = Opt.order;
    RegressionParams.BSFKOptions = struct('animation', 0, 'knotremoval_factor', Opt.knotremoval_factor);
elseif strcmpi(Type, 'SVM')
    RegressionParams.KernelFunction = Opt.KernelFunction;
    RegressionParams.Standardize = Opt.Standardize;
    RegressionParams.KernelScale = Opt.KernelScale;
elseif strcmpi(Type, 'GPR')
    RegressionParams.KernelFunction = Opt.KernelFunction;
elseif strcmpi(Type, 'RF')
elseif strcmpi(Type, 'MLS')
    RegressionParams.tol = Opt.tol;
    RegressionParams.order = Opt.order;
    RegressionParams.Bandwidth = Opt.Bandwidth;
elseif strcmpi(Type, 'NW')
    RegressionParams.KNN = Opt.KNN;
elseif strcmpi(Type, 'PB')
    RegressionParams.Bandwidth = Opt.Bandwidth;
else
    error('RegressionType %s not supported!\n', Type);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shapeRegr,YEst,res] = RegressionWrapper(X,Y,XGrid,RegressionParams)

Type = RegressionParams.Type;

if strcmpi(Type, 'BSFK')    
    pp = BSFK(X,Y,RegressionParams.order,RegressionParams.nknots,[],RegressionParams.BSFKOptions);
    shapeRegr = ppval(pp, XGrid);
    YEst = ppval(pp, X);
    res = norm(Y-YEst)/sqrt(length(Y));
elseif strcmpi(Type, 'SVM')
    if size(X,1) < size(X,2)
        X = X';
    end
    if isnumeric(RegressionParams.KernelScale)
        SVMModel = fitrsvm(X,Y,'KernelFunction',RegressionParams.KernelFunction,...
            'KernelScale',RegressionParams.KernelScale,...
            'Standardize',RegressionParams.Standardize);
    else
        SVMModel = fitrsvm(X,Y,'KernelFunction',RegressionParams.KernelFunction,...
            'KernelScale','auto',...
            'Standardize',RegressionParams.Standardize);
        autoKernelScale = SVMModel.KernelParameters.Scale;
        CVSVMModel = crossval(SVMModel);
        kfLoss = kfoldLoss(CVSVMModel);
        for m=[0.01,0.1,1,10]
            mSVMModel = fitrsvm(X,Y,'KernelFunction',RegressionParams.KernelFunction,...
                'KernelScale',m*autoKernelScale,...
                'Standardize',RegressionParams.Standardize);
            CVSVMModel = crossval(mSVMModel);
            mkfLoss = kfoldLoss(CVSVMModel);
            if mkfLoss < kfLoss
                SVMModel = mSVMModel;
                kfLoss = mkfLoss;
            end
        end
    end
    shapeRegr = predict(SVMModel, XGrid')';
    YEst = resubPredict(SVMModel)';
    res = norm(Y-YEst)/sqrt(length(Y));
elseif strcmpi(Type, 'GPR')
    if size(X,1) < size(X,2)
        X = X';
    end
    sigma0 = Chi2Estimation(Y);
    GPRModel = fitrgp(X,Y,'KernelFunction',RegressionParams.KernelFunction,'Sigma',sigma0);
    shapeRegr = predict(GPRModel, XGrid')';
    YEst = resubPredict(GPRModel)';
    res = norm(Y-YEst)/sqrt(length(Y));
elseif strcmpi(Type, 'RF')
elseif strcmpi(Type, 'MLS')
    W = exp(-getEuclideanDistance(XGrid,X).^2/RegressionParams.Bandwidth^2);
    coeff = MLS(X',XGrid',W,RegressionParams.order,RegressionParams.tol);
    shapeRegr = (coeff*Y')';
    W1 = exp(-getEuclideanDistance(X,X).^2/RegressionParams.Bandwidth^2);
    coeff1 = MLS(X',X',W1,RegressionParams.order,RegressionParams.tol);
    YEst = (coeff1*Y')';
    res = norm(Y-YEst)/sqrt(length(Y));
elseif strcmpi(Type, 'NW')
    r = ksrlin_vw(X,Y,RegressionParams.KNN,XGrid);
    shapeRegr = r.f;
    shapeRegr(isnan(shapeRegr)) = 0;
    rEst = ksrlin_vw(X,Y,RegressionParams.KNN,X);
    YEst = rEst.f;
    YEst(isnan(YEst)) = 0;
    res = norm(Y-YEst)/sqrt(length(Y));
elseif strcmpi(Type, 'PB')
    r = ksr_characteristic(X,Y,RegressionParams.Bandwidth,XGrid);
    shapeRegr = r.f;
    shapeRegr(isnan(shapeRegr)) = 0;
    rEst = ksr_characteristic(X,Y,RegressionParams.Bandwidth,X);
    YEst = rEst.f;
    YEst(isnan(YEst)) = 0;
    res = norm(Y-YEst)/sqrt(length(Y));
else
    error('RegressionType %s not supported!\n', Type);
end

end

