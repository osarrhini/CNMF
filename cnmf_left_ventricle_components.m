function cnmf_left_ventricle_components(...
    heart_mat, heart_ijk, indsOfLVCInHeart, indsOfMYOInHeart, ...
    nfactors, fmt, NormalizedData, alpha, beta, dataUnits, pn, ...
    liver_tac, MaxIter, bs)
%%cnmf_left_ventricle_components Enhanced Extraction of Blood and Tissue 
% Time-activity Curves in Cardiac Mouse FDG PET Imaging by Means of 
% Constrained Nonnegative Matrix Factorization
% Inputs:
%   heart_mat: PET data from the whole left ventricle (LV) ROI (N voxels x M frames)
%   heart_ijk: LV ROI voxel coordinates (N voxels x 3)
%   indsOfLVCInHeart: Indices of left ventricle cavity (LVC) in the whole LV ROI
%   indsOfMYOInHeart: Indices of myocardium (MYO) in the whole LV ROI
%   nfactors: Number of foctors to be extracted (mainly 2)
%   fmt: Mean frame time in seconds
%   liver_tac: Liver Time Activity curve (TAC)
%   NormalizedData: If true, CNMF will use normalized data (Each row of heart_mat will be summed to 1)
%   alpha and beta: regularization parameters
%   dataUnits: Data Activity Units [Default '[prop.to cnts/sec]']
%   pn: Folder name where to save session data
%   MaxIter: Maximum number of iterations [Default 5000]
%   bs: Whould be a vector of two values [bst, bsIDg] where
%   bst is the time in seconds and bsIDg is the activity value in %ID/g
%   of a late blood sample if any
%
% Note: All session data will be saved to a .mat file
%
% @Author: Otman Sarrhini (otman.sarrhini@usherbrooke.ca)
% Copyright: 2023-2025 by Otman Sarrhini,
%       Molecular Imaging Center of Sherbrooke,
%       Université de Sherbrooke, QC, Canada

%% Check if mandatory arguments are provided
nbr_ma = 11;
max_nbr_ia = 14;
if nargin < nbr_ma
    error('Please provide the first %d mondatory arguments', nbr_ma)
end
switch nargin
    case nbr_ma
        liver_tac = [];
        MaxIter = 5e3;
        bst = [];
        bsIDg = [];
    case nbr_ma + 1
        MaxIter = 5e3;
        bst = [];
        bsIDg = [];
    case nbr_ma + 2
        bst = [];
        bsIDg = [];
    case nbr_ma + 3
        bst = bs(1);
        bsIDg = bs(2);
    otherwise
        error('Number of input arguments should be at least %d and not exceed %d', nbr_ma, max_nbr_ia)
end

%% Initialization of w and h by ROIs data
[w0, h0, D, M, bldi, myoi] = ...
    init_wh(heart_mat, heart_ijk,indsOfLVCInHeart, indsOfMYOInHeart, nfactors, NormalizedData);

%% ------------------------------------------------------------------------
%  - Provide alpha = beta = 0 for Unconstrained NMF
%  - ----------------------------------------------------------------------
%  - For Constrained NMF, try differents values of alpha and beta.
%  -   I suggest varying alpha first (beta fixed to 0), then try to 
%  -   change beta with fixed alpha. The tail of the liver tac or a later 
%  -   blood sample should guide you during the tests.
[phys_cnmf,W_cnmf,H_cnmf, roi_btac, roi_myotac, cnmf_btac,cnmf_myotac] = ...
        call_cnmf_hp(heart_mat,nfactors,fmt,'alpha', alpha, 'beta',beta, ...
        'W_INIT',w0, 'H_INIT',h0, 'Blood_inds',indsOfLVCInHeart, ...
        'Tissue_inds', indsOfMYOInHeart, 'MAX_ITER', MaxIter, ...
        'LIV_TAC',liver_tac, 'BST', bst, 'BSA', bsIDg, ...
        'NormalizedData', NormalizedData, 'Data_units',dataUnits);

%% ------------------------------------------------------------------------
%  - Saving TACs to text files
%  ------------------------------------------------------------------------
time_col_header = 'Time (sec)';                         % Header for the time column
roi_plasma_col_header = ['ROI Plasma ' dataUnits];      % Header for roi plasma
roi_tissue_col_header = ['ROI MYO ' dataUnits];         % Header for roi tissue
cnmf_plasma_col_header = ['CNMF Plasma ' dataUnits];    % Header for cnmf plasma component
cnmf_tissue_col_header = ['CNMF MYO ' dataUnits];       % Header for cnmf tissue component

cnmf_ptac = plasma_from_blood_activity(cnmf_btac, fmt);
roi_ptac = plasma_from_blood_activity(roi_btac, fmt);
TACs_path = fullfile(pn, 'TACs.tac');
T_data = array2table([fmt(:) roi_ptac(:) roi_myotac(:) cnmf_ptac(:) cnmf_myotac(:)]);
T_data.Properties.VariableNames={time_col_header, roi_plasma_col_header, roi_tissue_col_header, cnmf_plasma_col_header, cnmf_tissue_col_header};
writetable(T_data, TACs_path, 'FileType', 'text', 'Delimiter', '\t');

%% ------------------------------------------------------------------------
%  - Calculation of Recovery coefficients and spill over fractions
%  ------------------------------------------------------------------------
W_TT = mean(squeeze(W_cnmf(indsOfMYOInHeart, myoi)));   % Estimation of Recovery coef in MYOCARDIUM ROI
W_TB = mean(squeeze(W_cnmf(indsOfLVCInHeart, myoi)));   % Spill over from MYO to LVC
W_BB = mean(squeeze(W_cnmf(indsOfLVCInHeart, bldi)));   % Recovery coef in LVC ROI
W_BT = mean(squeeze(W_cnmf(indsOfMYOInHeart, bldi)));   % Spill over from MYO to LVC
W_data = [W_TT W_BB W_TB W_BT];
T_W_data = array2table(W_data);
T_W_data.Properties.VariableNames = {'W_TT', 'W_BB', 'W_TB', 'W_BT'};
disp(T_W_data);

%% ------------------------------------------------------------------------
%  - Saving the session variables to .mat file
%  ------------------------------------------------------------------------
fncurrent = ['CNMF_analysis_' datestr(now(), 'yyyy-mm-dd_HHMMSS') '.mat'];             % Output file where session analysis data are stored
save(fullfile(pn, fncurrent));

end