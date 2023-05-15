function [physcmp, w, h, roibloodtac, roitissuetac, nmfbloodtac, nmftissuetac,bldi,myoi] = cnmf_last_steps(V, w, h, t, par)
%%cnmf_last_steps Computes and displays the CNMF TACs from W and H and 
%   Let the user to decide if the decomposition is satisfactory or not.
%   If Liver TAC is given then it is also displayed and the user can
%   clic two points away from each other on the last 15 min of the TAC.
%   The program will fit a biexponential to the liver tail and extrapolates
%   the fit to the last blood sample if given.
%   @Author Otman Sarrhini, Sherbrooke University
% **This code comes with no guarantee or warranty of any kind**

%% Matrices size
[n,r] = size(w);
[~,m] = size(h);

%% Normalized data
Vn = normalize_NMF_input(V,par.NormalizedData);

%% Contribution of each component to the voxel -----------------------------
Contr=100*sum(w)/sum(w(:));

%% matrix of physiological components
physcmp=zeros(n,m,r);
for i=1:r
    physcmp(:,:,i)=w(:,i)*h(i,:);
end

%% nmf tacs over the whole heart -------------------------------------------
tacs=zeros(m,r);
for i=1:r
    tacs(:,i)=(mean(squeeze(physcmp(:,:,i))))';
end

%% Plot the factors as the mean of each one over the heart ----------------
leg=num2str((1:r)');
bstr='FAC ';
str{1}='Data';
for i=1:r;str{i+1}=[bstr leg(i,:) ' : ' num2str(Contr(i)) ' %'];end;
leg=char(str);

f1=findobj(0,'Tag',par.strTag);
if ~isempty(f1);    
    figure(f1);
    clf(f1);
else    
    f1=figure;
    set(f1,'Tag',par.strTag,'DoubleBuffer','on','Name',...
        'Non-Negative Matrix Factorization (by Otman SARRHINI)',...
        'NumberTitle','off','Position',[508 291 698 525],...
        'Units','Normalized');
end
% Detect if t is in seconds or minutes
tm=t;
if t(end) > 200
    tm = t/60;
end
Vtac=(mean(Vn))';
plot(tm,[Vtac tacs],'-','linewidth',2);
set(gca,'Color',[.8 .8 .8],'FontSize',14,'FontName','Helvetica','FontWeight','bold');
leg1=legend(leg);
set(leg1,'Color',[.8 .8 .8],'EdgeColor',[.8 .8 .8],...
    'Position',[0.2159    0.7123    0.3052    0.1892]);

xlabel('Time (min)','FontSize',14,'FontName','Helvetica','FontWeight','bold');
ylabel('Normalized activity','FontSize',14,'FontName','Helvetica','FontWeight','bold');

% Root mean squared error (See Briggs and Levine 1997)
RMSE = sqrt((1/(n*m))*sum(sum((Vn-w*h).^2)));
title(['CNMF. Used algo: ''' par.algo ''' RMSE = ' num2str(RMSE)]); 

% Add axes for sum of w to see if this is equal to 1 for each voxel
sumWa = axes('Parent',f1,...
    'Position',[0.667621776504298 0.735238095238095 0.260744985673352 0.188571428571429],...
    'FontWeight','bold',...
    'Color',[0.8 0.8 0.8]);
box(sumWa,'on');
hold(sumWa,'all');
% plot(1:n,sum(w,2),'b-',1:n,sum(wa,2),'r-.','Parent',sumWa);
plot(1:n,sum(w,2),'b-','Parent',sumWa);
axis([0 n -inf inf]);

% pause;      % To zoom and identify blood components
% ------------------------------------------------------------------------
%% Specify blood component
cpts=1:r;

% New method to specify blood component
bldi = [];
while isempty(bldi)
    prompt = {'Enter blood component indice'};
    dlgTitle = 'Input';
    dims = [1 numel(prompt{1})];
    
    answer = inputdlg(prompt,dlgTitle,dims);
    if isempty(answer)
        bldi=[];
    elseif isempty(answer{1})
        bldi = [];
    else
        bldi0 = str2num(answer{1});
        
        if ~isempty(bldi0)
            for ii=1:numel(bldi0)
                cbldi = bldi0(ii);
                if ~isempty(find(cpts==cbldi,1,'first'))
                     bldi = [bldi;cbldi];
                end
            end
        end
    end
end

%% ........................................................................
% Calculate the blood tac
nmfbloodtac = [];
roibloodtac = [];
if ~isempty(bldi) && ~isempty(par.indsOfLVCInHeart)
    % ROI blood TAC
    roibloodtac = mean(V(par.indsOfLVCInHeart,:))';
    nmfbloodtac = get_cnmf_tac(physcmp, V, w, par.indsOfLVCInHeart, bldi, par.NormalizedData);
end

%% .........................................................................
% Calculate the Tissue TAC
nmftissuetac = [];
roitissuetac = [];
myoi=[];
if ~isempty(bldi) && ~isempty(par.indsOfMYOInHeart)
    % ROI tissue tac
    roitissuetac = mean(V(par.indsOfMYOInHeart,:))';
    
    % NMF blood and tissue TACs
    cpts(bldi)=[];myoi=cpts;
	nmftissuetac = get_cnmf_tac(physcmp, V, w, par.indsOfMYOInHeart, myoi, par.NormalizedData);
end

%% ------------------------------------------------------------------------
% Plot TACs
f0=findobj(0,'Tag',par.strTag0);
if ~isempty(f0);
    figure(f0);
    clf(f0);
else
    f0=figure;
    set(f0,'Tag',par.strTag0,'DoubleBuffer','on','Name',...
        'TACs',...
        'NumberTitle','off','Position',[872 381 872 606],...
        'Units','Normalized');
end
warning('off','MATLAB:Axes:NegativeDataInLogAxis');

%% Measured blood tac
current_leg={};
if ~isempty(roibloodtac)
    figure(f0);semilogx(...
        tm,roibloodtac,'bp-',...
        'linewidth',2);
    current_leg{end+1} = 'MEAS LV';
end

%% NMF blood tac
if ~isempty(nmfbloodtac)
    figure(f0);hold on;semilogx(...
        tm,nmfbloodtac,'rs-',...
        'linewidth',2);
    current_leg{end+1} = 'NMF LV';
end

%% Measured tissue tac
if ~isempty(roitissuetac)
    figure(f0);hold on;semilogx(...
        tm,roitissuetac,'ko-',...
        'linewidth',2);
    current_leg{end+1} = 'MEAS MYO';
end

%% NMF tissue tac
if ~isempty(nmftissuetac)
    figure(f0);hold on;semilogx(...
        tm,nmftissuetac,'m+-',...
        'linewidth',2);
    current_leg{end+1} = 'NMF MYO';
end

%% Liver
if ~isempty(par.liver_tac)
	figure(f0);hold on;semilogx(tm,par.liver_tac,'g*-','linewidth',2);
    current_leg{end+1}='MEAS LIV';
    
    msg_txt = "In the next step, you will be invited to click on\n" + ...
        "two distinct points of the liver curve (liver_tac) to delimit\n" + ...
        "a part of it which will be fitted by a biexponential.\n" + ...
        "This will allow you to compare nmf_blood_tac to the blood sample\n" +...
        "and/or to liver_tac tail.";
    uiwait(msgbox(msg_txt,"Information","help","modal"));

    % Fit a bi-exponential function to the tail of this curve
    gg=ginput(2);gx=gg(:,1);
    i1=find(abs(tm-gx(1)) == min(abs(tm-gx(1))),1,'first');
    i2=find(abs(tm-gx(2)) == min(abs(tm-gx(2))),1,'first');
    xdata=tm(i1:i2);ydata=par.liver_tac(i1:i2);
    [BiEx_param,BiExp_fit]=bi_exponential_fit(xdata,ydata);

    
    % Extrapolation of the fit to pass through the blood sample
    if ~isempty(par.bst)
        xdata = sort([xdata;par.bst/60]);
        BiExp_fit=bi_exponential_func(BiEx_param,xdata);
    end
    
    % Scaling the fit to have the same AUC as nmf bloodtac
    if ~isempty(nmfbloodtac)
        BiExp_fit = BiExp_fit*sum(nmfbloodtac(i1:i2))/sum(par.liver_tac(i1:i2));
        figure(f0);hold on;semilogx(xdata,BiExp_fit,'c:','linewidth',2);
    end
end

%% Blood sample
if ~isempty(par.bst)
    figure(f0);hold on;semilogx(par.bst/60,par.bsa,'cd','MarkerSize',8,'MarkerFaceColor','c');
    current_leg{end+1} = 'BS';
end

if ~isempty(get(f0,'children'))
    set(gca,'Color',[.8 .8 .8],'FontSize',14,'FontName','Helvetica','FontWeight','bold');
    xlabel('Time(min)');ylabel(par.dataUnits);
    legend('String',current_leg,'Color',[0.8 0.8 0.8], 'edgecolor',[.8 .8 .8]);
else
    close(f0);
end