function [w0, h0, D, M, bldi, myoi] = init_wh(heart_mat, heart_ijk, ...
    indsOfLVCInHeart, indsOfMYOInHeart, nfactors, normalizedData)
%%init_wh: Initialisation of w and h using curves obtained with crude ROIs
% Inputs:
%   heart_mat: Whole left ventricle (LV) data (n_voxels x nbr_frames)
%   heart_ijk: LV voxels's coordinates as [i j k] values
%   indsOfLVCInHeart: Indices of LVC voxels in LV
%   indsOfMYOInHeart: Indices of MYO voxels in LV
%   nfactors: number of physiological factors (2 in most case)
%   normalizedData: If true (1) then data will be normalized before.
% Outputs:
%   w0: Guess of W
%   h0: Guess of H
%   D: Matrix of Euclidian distance of each voxel to each ROI
%   M: Matrix of masks (M(i, j) = 1 if and only if i=j else 0
%   bldi: blood component number
%   myoi: Tissue component number
% **This code comes with no guarantee or warranty of any kind.**

%% Normalize data
Y = normalize_NMF_input(heart_mat, normalizedData);

%% Heart ROI data dims
[nbPix, nbIm]=size(Y);

%% Mask of ROIs col_1 is tissue mask, col_2 is blood mask
myoi=1;
bldi=2;
M = zeros(nbPix,nfactors);
M(indsOfMYOInHeart, myoi) = 1;
M(indsOfLVCInHeart, bldi) = 1;

%% Matrix of euclidian distance to each ROI
D = compute_D(M, heart_ijk);

%% Set w0 for each voxel
w0 = 1./(1+D);
w0 = w0./(sum(w0,2)*ones(1,nfactors));

%% h0 is initialized with crude rois curves
h0=zeros(nfactors, nbIm);
h0(myoi,:) = mean(Y(indsOfMYOInHeart, :));
h0(bldi,:) = mean(Y(indsOfLVCInHeart, :));

end
