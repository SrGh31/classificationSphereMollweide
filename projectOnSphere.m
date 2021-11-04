function [onSphere] = projectOnSphere(original)
% This reshapes the original input space into data projected on the surface
% of an unit sphere (radius =1)
    dimension=size(original,2); radius=1;
    onSphere=original.*repmat(radius*(ones(size(original,1),1).^(1/dimension))./sqrt(sum(original.^2,2)),1,dimension);
end