function indx=centerprop(d,D,ind)
% where, 
% d is the step/jump propose
% D is the matrix with distances between each pair of voxel in the source region. 
% ind is the current CF center position in the source region

    D=D(:,ind);
    absD=abs(D-d);
    indx=find(absD==min(absD));
    indx=indx(randperm(size(indx,1),1));
    
end