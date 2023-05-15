function [physcmp, w, h, roibloodtac, roitissuetac, nmfbloodtac, nmftissuetac, bldi, myoi]=call_cnmf_hp(V,r,t,varargin)
%%cnmf_hp Decomposition of the matrix V (n x m) into 2 matrices W(n x r) and 
% H(r x m) so that V = W * H using the non-negative matrix factorization
% V: Matrix from an ROI intensity values (n voxels and m frames)
% r: Number of the physiological factors in the ROI
% t: time (s) (for plot purpose)
%   Optional parameter must be by pairs 'Parameter_Name',Parameter_value
%             'ALGO',               Used algorithme ('regularized')
%             'ALPHA',              alpha
%             'BETA',               beta
%             'MAX_ITER',           MaxIter
%             'W_INIT',             Initial guess of W
%             'H_INIT',             Initial guess of H
%             'BLOOD_INDS',         Inds of blood ROI voxels in Heart ROI to estimate nmf blood tac
%             'TISSUE_INDS',        Inds of tissue ROI voxels in Heart ROI to estimate nmf tissue tac
%             'LIV_TAC',            Liver TAC if any
%             'BST',                Blood sample time (s)
%             'BSA',                Blood sample activity in data_units
%             'NORMALIZEDDATA',     true (recommended to use normalized data) or false
%             'DATA_UNITS',         data units
% To call:
% [physcmp, w, h, roibloodtac, roitissuetac, nmfbloodtac, nmftissuetac, bldi, myoi]=cnmf_hp(V,r,t,varargin)
% @Author: Otman Sarrhini
% Date of creation : 2014/06/11
% Date of last update:  2019/09/11
% **This code comes with no guarantee or warranty of any kind.**

% Required parameters: V,r,t
if nargin < 3
    error('Must supply roi data, number of components and time in seconds')
end

%% Check the input matrix
[n,m]=size(V);
if any([n,m] == 1)
    error(' Only matices are allowed ');
elseif ~ismatrix(V)
    error(' Only 2D matrices are allowed ');
end
if numel(t) ~= m
    error('Size of the input matrix V along the second dimesion must be the same as time (t) ');
end

%% Negative values get replaced by eps ------------------------------------
if min(V(:)) < 0
    warning('WarnNMF:negVal','Negative values were found in input matrix.\nThey will be replaced by matlab eps');
    V(V(:) < 0) = eps;
end

%% rank r criterion must be respected --------------------------------------
if (n+m)*r >= n*m
    error('r should be choosen so that (n+m)*r << n*m');
end

%% Min activity should not be less than eps -------------------------------
V(abs(V) < eps)=eps;

%% Read optional parameters
if (rem(numel(varargin),2) == 1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(numel(varargin)-1)
        switch upper(varargin{i})
            case 'ALGO',                par.algo = varargin{i+1};
            case 'ALPHA',               par.alpha = varargin{i+1};
            case 'BETA',                par.beta = varargin{i+1};
            case 'MAX_ITER',            par.MaxIter = varargin{i+1};
            case 'W_INIT',              par.w0 = varargin{i+1};
            case 'H_INIT',              par.h0 = varargin{i+1};
            case 'BLOOD_INDS',          par.indsOfLVCInHeart = varargin{i+1};
            case 'TISSUE_INDS',         par.indsOfMYOInHeart = varargin{i+1};
            case 'LIV_TAC',             par.liver_tac = varargin{i+1};
            case 'BST',                 par.bst = varargin{i+1};
            case 'BSA',                 par.bsa = varargin{i+1};
            case 'NORMALIZEDDATA',      par.NormalizedData = varargin{i+1};
            case 'DATA_UNITS',          par.dataUnits = varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end

%% Apply default configuration for missed parameters
if ~exist('par','var')
    error('Optional parameters are required');
else
    % Default values
    if ~isfield(par,'algo')
        par.algo='regularized';
    end
    if ~isfield(par,'MaxIter')
        par.MaxIter = 5e3;
    end
    if ~isfield(par,'liver_tac')
        par.liver_tac = [];
    end
    if ~isfield(par,'bst')
        par.bst = [];
    end
    if ~isfield(par,'bsa')
        par.bsa = [];
    end
    
    % Requiered values
    if ~isfield(par,'indsOfLVCInHeart')
        error('indsOfLVCInHeart is required');
    end
    if ~isfield(par,'indsOfMYOInHeart')
        error('indsOfMYOInHeart is required');
    end
    if ~isfield(par,'alpha')
        error('alpha is required');
    end
    if ~isfield(par,'beta')
        error('beta is required');
    end
    if ~isfield(par,'w0')
        error('Initial guess w0 is required');
    end
    if ~isfield(par,'h0')
        error('Initial guess h0 is required');
    end
    if ~isfield(par,'NormalizedData')
        error('Whether or not to use normalized data must be specified');
    end
    if ~isfield(par,'dataUnits')
        error('Data units is required');
    end
end

%% The nmf curves figure tag
strTag='nnmfroifigHP';
par.strTag = strTag;
par.strTag0 = [par.strTag '_' mfilename];

%% Decomposition of V in w and h ------------------------------------------
tol = 1e-6;
Vn = normalize_NMF_input(V, par.NormalizedData);
[w,h] = cnmf_hp(Vn,r,'W_INIT',par.w0,'H_INIT',par.h0,'TYPE',par.algo, ....
    'ALPHA',par.alpha,'BETA',par.beta,'MAX_ITER',par.MaxIter,'MIN_ITER',par.MaxIter,'TOL',tol);

%% Scale constraint (Bodvarsson_2007.pdf; WeiMu_NMF_2013.pdf) so that sum(w,2)=1 for each row
[w,h] = bodvarsson_wh_scaling( w, h );

%% Remaining CNMF steps are completed in the following funciton
[physcmp, ~, ~, roibloodtac, roitissuetac, nmfbloodtac, nmftissuetac, bldi, myoi] = cnmf_last_steps(V, w, h, t, par);
