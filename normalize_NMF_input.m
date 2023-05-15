function Y = normalize_NMF_input(heart_mat, NormalizedData)
%%normalize_NMF_input LV PET data normalization
% Normalization of heart_mat so that each row sums to 1 and returns results
% in Y. Negative data are replaced by 0. Rows which sum to 0 are untouched
% @Author: Otman Sarrhini; Sherbrooke University
% -------------------------------------------------------------------------

heart_mat(heart_mat(:) < 0) = 0;
if ~NormalizedData
    Y=heart_mat;
else
    s2=sum(heart_mat,2);
    s2(s2 == 0) = 1;
    Y = heart_mat./(s2*ones(1,size(heart_mat,2)));
end


