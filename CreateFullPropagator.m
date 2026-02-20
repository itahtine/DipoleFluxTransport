% Creates full propagator matrices from latitude slices

%%
load('VectorMapSliceDiff350Mer11CR')

VectorMap = cat(1,repmat(VectorMapSlice,1,360,1,1),flip(repmat(VectorMapSlice,1,360,1,1),1));
VectorMap(91:180,:,2,:) = -VectorMap(91:180,:,2,:);

for i = 1:360
    VectorMap(:,i,3,:) = VectorMap(:,i,3,:)+(i-1);
end

%save('VectorMapFullDiff350Mer11CR','VectorMap')
%%
load('VectorMapSliceDiff350Mer11Daily')

VectorMap = cat(1,repmat(VectorMapSlice,1,360,1,1),flip(repmat(VectorMapSlice,1,360,1,1),1));
VectorMap(91:180,:,2,:) = -VectorMap(91:180,:,2,:);

for i = 1:360
    VectorMap(:,i,3,:) = VectorMap(:,i,3,:)+(i-1);
end

%save(VectorMapFullDiff350Mer11Daily','VectorMap')

%%
