function y=bi_exponential_func(k,t,P)
% Bi exponential funtion
% y = a1*exp(-(ln(2)/T1)*(t-t0))+a2*exp(-(ln(2)/T2)*(t-t0))
% @Author Otman Sarrhini
% **This code comes with no guarantee or warranty of any kind.**
if exist('P','var')
    if isfield(P,'FixedP')
        for i=1:numel(P.FixedP)
            if P.FixedP(i)
                k(i)=P.IniG(i);
            end
        end
    end
end
a1 = k(1);
T1 = k(2);
a2 = k(3);
T2 = k(4);
t0 = k(5);

% y=k(1)*exp(-log(2)*(t-k(5))/k(2))+k(3)*exp(-log(2)*(t-k(5))/k(4));
y=a1*exp(-log(2)*(t-t0)/T1)+a2*exp(-log(2)*(t-t0)/T2);

end