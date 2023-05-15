function priorA = compute_D(M, XYZ)
%%compute_D: Euclidian distance between each voxel to each ROI
% @input :
%   M : Matrix of masks (ROIs)
%   XYZ: Matrix of coordinates
% @return : 
%   priorA : Matrix of euclidean distance of each pixel to each mask
%   (vectorized)
% Adapted from Filippi et al:
%   M. Filippi, M. Desvignes, and E. Moisan, 
%   “Robust Unmixing of Dynamic Sequences Using Regions of Interest,” 
%   IEEE Trans. Med. Imaging, vol. 37, no. 1, pp. 306–315, 2018, 
%   doi: 10.1109/TMI.2017.2759661.
% **This code comes with no guarantee or warranty of any kind.**

[nbPix, K]=size(M);
priorA = 0*M;
for i = 1:nbPix
    xyzi = XYZ(i, :);                   % Current voxel coordiantes
    for k=1:K
        Mk = M(:,k);                    % current ROI
        inds_k = find(Mk == 1);
    
        xyzi_matrix = ones(numel(inds_k),1) * xyzi;
        priorA(i, k) = min(sqrt(sum((xyzi_matrix - XYZ(inds_k,:)).^2,2)));
    end
end

end