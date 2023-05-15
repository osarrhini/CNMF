function [Ki,yfit,HalfLives,resnorm]=bi_exponential_fit(t,ydata,P)
%%biexponentialfit3
% bi-exponential fitting of the data ydata(t)
%   @Author: Otman Sarrhini
% **This code comes with no guarantee or warranty of any kind.**


if nargin < 2
    error('Must supply at least time (t(min)) and ydata');
end

a1=max(ydata)/2;a2=a1;
b1=-log(2)*(t(end)-t(1))/(ydata(end)-ydata(1));              % 60 to convert to min
b2=b1/3;
[~,it0]=max(ydata);t0=t(it0);
k0=[a1; b1; a2; b2; t0];
lb0 = eps*ones(numel(k0),1);
ub0 = Inf*ones(numel(k0),1);


if nargin < 3
    P.IniG = k0;
    P.FixedP = 0*k0;
    P.lb=lb0;
    P.ub=ub0;
end

if ~isfield(P,'IniG')
    P.IniG = k0;
end
if ~isfield(P,'FixedP')
    P.FixedP = 0*k0;
end
if ~isfield(P,'lb')
    P.lb = lb0;
end
if ~isfield(P,'ub')
    P.ub = ub0;
end

% Work on the interpolated data
imethod='linear';
y=ydata;
ti=t(1):(min(diff(t))/5):t(end);ti(end)=t(end);
yi=interp1(t,y,ti,imethod,'extrap');

options=optimset('DiffMinChange',1e-8,'DiffMaxChange',1e-6,'Display','off',...
    'TolFun',1e-8,'TolX',1e-8,...
    'MaxIter',Inf,...
    'FunValCheck','on','MaxFunEvals',Inf,'TypicalX',P.IniG,...
    'ScaleProblem','Jacobian');

[Ki,resnorm]=lsqcurvefit(@bi_exponential_func,P.IniG,ti,yi,P.lb,P.ub,options,P);
yfit=bi_exponential_func(Ki,t,P);
HalfLives=[Ki(2);Ki(4)];

exp1=Ki(1)*exp(-log(2)*(t-Ki(5))/Ki(2));
exp2=Ki(3)*exp(-log(2)*(t-Ki(5))/Ki(4));

strTag0 = 'biexpfit_fig';
f0=findobj(0,'Tag',strTag0);
if ~isempty(f0);
    figure(f0);
    clf(f0);
else
    f0=figure;
    set(f0,'Tag',strTag0,'DoubleBuffer','on','Name',...
        'BiExponentialFit',...
        'NumberTitle','off','Position',[680 558 560 420],...
        'Units','Normalized');
end
figure(f0);plot(t,ydata,'ko',t,yfit,'b-',...
    t,exp1,'r:',t,exp2,'g--','linewidth',3,'MarkerSize',10);
set(gca,'FontWeight','bold','fontsize',14,'color',[.8 .8 .8]);
leg=legend('Data','Bi-Exp Fit','Exp1','Exp2');
set(leg,'color',[.8 .8 .8],'EdgeColor',[.8 .8 .8]);
xlabel('Time (min)','FontWeight','bold','fontsize',14);
ylabel('Data','FontWeight','bold','fontsize',14);

end
