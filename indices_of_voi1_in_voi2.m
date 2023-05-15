function indsOfVOI1InVOI2 = indices_of_voi1_in_voi2(voi1_ijk, voi2_ijk)
%indices_of_voi1_in_voi2 Find indices VOI1 in VOI2 and returns them
%   This function is part of NMF toolbox
%   @Author Otman Sarrhini
%   To call: indsOfVOI1InVOI2 = indices_of_voi1_in_voi2(voi1_ijk, voi2_ijk);
% Inputs:
%   voi1_ijk: a (N1, 3) matrix of voxels's position of VOI1 in the PET image
%   voi2_ijk: a (N2, 3) matrix of voxels's position of VOI2 in the PET image
% Outputs:
%   indsOfVOI1InVOI2: Indices of VOI1 voxels in VOI2 matrix
% **This code comes with no guarantee or warranty of any kind.**

EPS1 = 1e-5;
indsOfVOI1InVOI2=zeros(size(voi1_ijk, 1), 1);
for i=1:size(voi1_ijk,1)
    xyzi_mat = ones(size(voi2_ijk,1),1) * voi1_ijk(i,:);
    dd=sum((voi2_ijk-xyzi_mat).^2,2);
    ii=find(dd <= EPS1,1,'first');
    if ~isempty(ii)
        indsOfVOI1InVOI2(i) = ii;
    end
end

end