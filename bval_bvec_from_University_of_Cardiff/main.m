%% 
%  Loading bvec and bval for both longer scan (60 directions: in 1_DWI directory) and shorter scan (30 directions: in 3_DWI directory)
%    Here we only process data with b = 1200
%    If you are interested in the imaging data, please contact the authors of
%           https://doi.org/10.1016/j.neuroimage.2019.01.077
%  

subjectID = "./A_prisma"

bvecLong = dlmread(fullfile(subjectID, "scans/1_DWI/TXT/dwi.bvec"));
bvalLong = dlmread(fullfile(subjectID, "scans/1_DWI/TXT/dwi.bval"));
bvecShort= dlmread(fullfile(subjectID, "scans/3_DWI/TXT/dwi.bvec"));
bvalShort= dlmread(fullfile(subjectID, "scans/3_DWI/TXT/dwi.bval"));
find1200Long = find(bvalLong == 1200);
bvec1200Long = bvecLong(:, find1200Long);
find1200Short = find(bvalShort == 1200);
bvec1200Short = bvecShort(:, find1200Short);

%% Under-sampling through matching directions

bvecLongMatchWithShort1200 = zeros(3,length(find1200Short));
bvecLongMatchWithShort1200Index = zeros(1,length(find1200Short));
for cnt = 1:length(find1200Short)
    tmp = bvec1200Short(:,cnt);
    tmp = repmat(tmp,[1 length(find1200Long)]);
    tmp2 = dot(tmp,bvec1200Long);
    tmp2 = abs(tmp2); % Note that we take absolute value here (e.g., [1,0,0] and [-1,0,0] have the same effect on diffusion signals)
    tmpL = find(tmp2 == max(tmp2(:)));
    bvecLongMatchWithShort1200Index(cnt)= find1200Long(tmpL(1));
    bvecLongMatchWithShort1200(:,cnt) = bvec1200Long(:,tmpL(1));
end
%% Visualization

tmpx =[bvec1200Short(1,:)' ;-bvec1200Short(1,:)'];
tmpy =[bvec1200Short(2,:)' ;-bvec1200Short(2,:)'];
tmpz =[bvec1200Short(3,:)' ;-bvec1200Short(3,:)'];
DTShort = delaunayTriangulation(tmpx,tmpy,tmpz);

tmpx =[bvecLongMatchWithShort1200(1,:)' ;-bvecLongMatchWithShort1200(1,:)'];
tmpy =[bvecLongMatchWithShort1200(2,:)' ;-bvecLongMatchWithShort1200(2,:)'];
tmpz =[bvecLongMatchWithShort1200(3,:)' ;-bvecLongMatchWithShort1200(3,:)'];
DTLongMatchWithShort = delaunayTriangulation(tmpx,tmpy,tmpz);

[KShort,~] = convexHull(DTShort);
figure(1); 
trisurf(KShort,DTShort.Points(:,1),DTShort.Points(:,2),DTShort.Points(:,3))
axis equal; view([-27.00 1.72])
title("Delaunay triangles from 30 directions")

[KLongMatchWithShort,~] = convexHull(DTLongMatchWithShort);
figure(2);
trisurf(KLongMatchWithShort,DTLongMatchWithShort.Points(:,1),DTLongMatchWithShort.Points(:,2),DTLongMatchWithShort.Points(:,3))
axis equal; view([-27.00 1.72]); 
title("Delaunay triangles from 30 directions down-sampled from 60 directions")

%% Calculating areas of triangles

areaShort = zeros(size(KShort,1),1);
areaLongMatchWithShort = zeros(size(KLongMatchWithShort,1),1);

for cnt = 1:size(KShort,1)
    tmp = KShort(cnt,:);
    tmpA = DTShort.Points(tmp(1),:);
    tmpB = DTShort.Points(tmp(2),:);
    tmpC = DTShort.Points(tmp(3),:);
    tmpArea = 0.5 * norm(cross(tmpB-tmpA,tmpC-tmpA));
    areaShort(cnt,1)=tmpArea;
end
for cnt = 1:size(KLongMatchWithShort,1)
    tmp = KLongMatchWithShort(cnt,:);
    tmpA = DTLongMatchWithShort.Points(tmp(1),:);
    tmpB = DTLongMatchWithShort.Points(tmp(2),:);
    tmpC = DTLongMatchWithShort.Points(tmp(3),:);
    tmpArea = 0.5 * norm(cross(tmpB-tmpA,tmpC-tmpA));
    areaLongMatchWithShort(cnt,1)=tmpArea;
end

areaRef = mean(areaShort);
areaShort = areaShort/areaRef;
areaLongMatchWithShort = areaLongMatchWithShort/areaRef;

figure(3);hold off;
plot(1:size(KShort,1),sort(areaShort),'bo-');
figure(3);hold on;
plot(size(KShort,1)+1:size(KShort,1)+size(KShort,1),sort(areaLongMatchWithShort),'ro-');
display(std(areaLongMatchWithShort)/std(areaShort))
hold off;

%% Randomly choosing a subset of sampling vectors

randomList = randperm(length(find1200Long));
randomList = randomList(1:length(find1200Short)); 
vecLongMatchWithShortRandom = bvec1200Long(:,randomList);

%% Visualization

tmpx =[vecLongMatchWithShortRandom(1,:)' ;-vecLongMatchWithShortRandom(1,:)'];
tmpy =[vecLongMatchWithShortRandom(2,:)' ;-vecLongMatchWithShortRandom(2,:)'];
tmpz =[vecLongMatchWithShortRandom(3,:)' ;-vecLongMatchWithShortRandom(3,:)'];
DTLongMatchWithShortRandom = delaunayTriangulation(tmpx,tmpy,tmpz);

figure(4);
[KLongMatchWithShortRandom,~] = convexHull(DTLongMatchWithShortRandom);
trisurf(KLongMatchWithShortRandom,DTLongMatchWithShortRandom.Points(:,1),DTLongMatchWithShortRandom.Points(:,2),DTLongMatchWithShortRandom.Points(:,3))
axis equal; view([-27.00 1.72]); 
title("Delaunay triangles from 30 directions randomly down-sampled from 60 directions")

%% Calculating areas of triangles

areaLongMatchWithShortRandom = zeros(size(KLongMatchWithShortRandom,1),1);

for cnt = 1:size(KLongMatchWithShortRandom,1)
    tmp = KLongMatchWithShortRandom(cnt,:);
    tmpA = DTLongMatchWithShortRandom.Points(tmp(1),:);
    tmpB = DTLongMatchWithShortRandom.Points(tmp(2),:);
    tmpC = DTLongMatchWithShortRandom.Points(tmp(3),:);
    tmpArea = 0.5 * norm(cross(tmpB-tmpA,tmpC-tmpA));
    areaLongMatchWithShortRandom(cnt,1)=tmpArea;
end
areaLongMatchWithShortRandom = areaLongMatchWithShortRandom/areaRef;

figure(5);clf; hold off;
plot(1:size(KShort,1),sort(areaShort),'bo-');
figure(5);hold on;
plot(size(KShort,1)+1:size(KShort,1)+size(KShort,1),sort(areaLongMatchWithShortRandom),'ro-');
display(std(areaLongMatchWithShortRandom)/std(areaShort))

%% Under-sampling with random process (1000 times)
%
% stdSaveMe = zeros(1000,1);
% randomListSaveMe = zeros(1000,length(find1200Short));
% for cnt2 = 1:1000
%     randomList = randperm(length(find1200Long));
%     randomList = randomList(1:length(find1200Short));
%     vecLongMatchWithShortRandom = bvec1200Long(:,randomList);
%     tmpx =[vecLongMatchWithShortRandom(1,:)' ;-vecLongMatchWithShortRandom(1,:)'];
%     tmpy =[vecLongMatchWithShortRandom(2,:)' ;-vecLongMatchWithShortRandom(2,:)'];
%     tmpz =[vecLongMatchWithShortRandom(3,:)' ;-vecLongMatchWithShortRandom(3,:)'];
%     DTLongMatchWithShortRandom = delaunayTriangulation(tmpx,tmpy,tmpz);
%     [KLongMatchWithShortRandom,~] = convexHull(DTLongMatchWithShortRandom);
%     areaLongMatchWithShortRandom = zeros(size(KLongMatchWithShortRandom,1),1);
%     for cnt = 1:size(KLongMatchWithShortRandom,1)
%         tmp = KLongMatchWithShortRandom(cnt,:);
%         tmpA = DTLongMatchWithShortRandom.Points(tmp(1),:);
%         tmpB = DTLongMatchWithShortRandom.Points(tmp(2),:);
%         tmpC = DTLongMatchWithShortRandom.Points(tmp(3),:);
%         tmpArea = 0.5 * norm(cross(tmpB-tmpA,tmpC-tmpA));
%         areaLongMatchWithShortRandom(cnt,1)=tmpArea;
%     end
%     areaLongMatchWithShortRandom = areaLongMatchWithShortRandom/areaRef;
%     tmp10 = std(areaLongMatchWithShortRandom)/std(areaShort);
%     stdSaveMe(cnt2,1) = tmp10;
%     randomListSaveMe(cnt2,:) = randomList;
% end
% figure;
% plot(stdSaveMe);
