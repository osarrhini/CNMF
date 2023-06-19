%% ########################################################################
% Load your data here to use CNMF algorithm
% ------------------------------------------------------------------------
% You need to a .mat file containing the following variables:
%   heart_mat: a N by M matrix where N is the number of the voxels and M is 
%   the number of the frame time. This matrix contains in each row the Time
%   Activity Curve (TAC) of each image voxel from the Left ventricle.
%   heart_ijk: a N by 3 matrix containing the (i, j, k) coordinates of the 
%   left ventricle voxels.
%   indsOfLVCInHeart: a vector containing the indices of the Left ventricle
%   Cavity (blood) voxels in the left ventricle ROI
%   indsOfMYOInHeart: a vector containing the indices of the myocardium in
%   the left ventricle ROI
%   nfactors: number of factors to be extracted (only 2 are supported)
%   fmt: Frame mean time in seconds
%   NormalizeData: Must be true to use normalized voxel's TAC
%   alpha and beta are the regularization parameters
%   dataUnits: Data units (example %ID/g, kBq/uL, ...)
%   pn: Folder name where the session data will be stored
%   Optional arguments:
%   liver_tac: Liver TAC (optional argument but strongly recommanded)
%   MaxIter: Maximum number of iterations (default 5000)
%   bs: vector containing two values: blood sample activity and time in
%   seconds (optional argument)
%
%   Here is my suggestion for using alpha and beta
%       Start with alpha = beta = 0, then try to modify alpha first 
%       from 1e-6 for example and see how the CNMF blood TAC is relative to
%       Liver TAC or to a late blood sample. 
%       In most case beta doesn't need to be changed from 0.
%       If needed, set beta to non zero value (1e-2 for example) and see the result.
%       You can try many values of alpha and beta in the range published in
%       the manuscript.
% -------------------------------------------------------------------------       

% load(<path to a file containing data as described above>)
alpha = 0;
beta = 0;

cnmf_left_ventricle_components(...
    heart_mat, heart_ijk, indsOfLVCInHeart, indsOfMYOInHeart, ...
    nfactors, fmt, NormalizedData, alpha, beta, dataUnits, pn, ...
    liver_tac, MaxIter, bs);

% If Unconstrained NMF did not provide satisfactory results then :
% change alpha and beta and call cnmf_left_ventricle_components until
% you get satisfactory results (See the above suggestion).
%  ########################################################################

