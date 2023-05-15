function nmftac = get_cnmf_tac(phys, V, w, roi_indices_in_heart, component, bNormalizedData)
%%get_cnmf_tac: calculate nmf tac from physiological components 
%   (phys(m voxels x n frames x nbr of components)
%   V is the ROI data (m voxels x n frames)
%   w is the matrix of weight
%   roi_indices_in_heart, component index
%   Since the PET images are noisy and some adjacent tissues can contribute
%   to the decomposed ROI then we use this function to calculate the nmf
%   tac given the indices of the roi in the whole heart ROI
%   @Author: Otman Sarrhini, Sherbrooke University
% **This code comes with no guarantee or warranty of any kind.**

if ~exist('phys','var') || ...
        ~exist('V','var') || ...
        ~exist('w','var') || ...
        ~exist('roi_indices_in_heart','var') || ...
        ~exist('component','var')
    error('All input arguments are required');
else
    % W*H of the component at the roi indices
    % Since the decomposition was done on normalized data then we multiply
    % the result wh_roi by the sum of the corresponding rows
    s2=sum(V(roi_indices_in_heart,:),2);    % as a column vector
    s2=s2*ones(1,size(phys,2));             % Repeat the column to get a matrix
    if numel(component) == 1
        W_roi = squeeze(w(roi_indices_in_heart,component));
        wh_roi=squeeze(phys(roi_indices_in_heart,:,component));
    else
        W_roi = squeeze(sum(w(roi_indices_in_heart,component),2));
        wh_roi=squeeze(sum(phys(roi_indices_in_heart,:,component),3));
    end
    if bNormalizedData
        wh_roi = wh_roi.*s2;
    end
    % Fix the problem with zeros in W_roi
    W_roi_zeros = find(W_roi == 0);
    W_roi(W_roi_zeros) = [];
    wh_roi(W_roi_zeros,:) = [];
    nmftac = mean(wh_roi./(W_roi*ones(1,size(phys,2))))';
end