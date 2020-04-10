%% A software procedure for guiding and evaluating diffusion MRI down-sampling schemes
% The codes were implemented and tested with Matlab R2019
% Nan-kuei Chen


%% Loadng diffusion encoding directions from two protocols that need to be harmonized in post-processing
% Protocol P64: A protocol with 64 directions at b = 800 s/mm2
% Protocol P30: A protocol with 30 directions at b = 800 s/mm2

fileLong = 'P64.txt';   % Lines #28 to #91 describes 64 directions;
                        % Lines #1 to #27 could be ignored (b<800)

fileShort = 'P30.txt';  % Lines #2 to #31 describes 30 direction
                        % Line #1 could be ignored (b<800)

fileID = fopen(fileLong);
infoLong = textscan(fileID, '%f %f %f');
fclose(fileID);
lenLong = length(infoLong{1});

fileID = fopen(fileShort);
infoShort = textscan(fileID, '%f %f %f');
fclose(fileID);
lenShort = length(infoShort{1});


%% Down-sampling P64 to 30 directions, to match P30 as much as possible
%
% Specifically, for each of the 30 diffusion-encoding gradient directions 
% of scan protocol P30, we use the absolute value of inner product 
% to identify the most similar direction among the 64 diffusion-encoding 
% directions of protocol P64. 
%
% Note that two diffusion encoding directions are considered the same, 
% when the inner product is either n or -n (ignoring eddy current effect)
% 
% The codes were implemented to ignore the vectors with length<1 
%
% The output of this section include "vecLongMatchWithShortIndex" 
%                                and "vecLongMatchWithShort"

vecShort = zeros(3,lenShort);
vecLong = zeros(3,lenLong);
for cnt = 1: 3
    vecShort(cnt,:)=infoShort{cnt};
    vecLong(cnt,:)=infoLong{cnt};
end

ampVectorLong = zeros(lenLong,1);
for cnt = 1:lenLong
    tmp1 = vecLong(:,cnt);
    tmp2 = norm(tmp1);
    ampVectorLong(cnt,1) = tmp2;
end
L_Long = find(ampVectorLong>0.99 & ampVectorLong<1.01);
vecLong = vecLong(:,L_Long);

ampVectorShort = zeros(lenShort,1);
for cnt = 1:lenShort
    tmp1 = vecShort(:,cnt);
    tmp2 = norm(tmp1);
    ampVectorShort(cnt,1) = tmp2;
end
L_Short = find(ampVectorShort>0.99 & ampVectorShort<1.01);
vecShort = vecShort(:,L_Short);

vecLongMatchWithShort = zeros(3,length(L_Short));
vecLongMatchWithShortIndex = zeros(1,length(L_Short));
for cnt = 1:length(L_Short)
    tmp = vecShort(:,cnt);
    tmp = repmat(tmp,[1 length(L_Long)]);
    tmp2 = dot(tmp,vecLong);
    tmp2 = abs(tmp2); % Note that we take absolute value here (e.g., assuming [1,0,0] and [-1,0,0] have the same effect on diffusion signals)
    tmpL = find(tmp2 == max(tmp2(:)));
    vecLongMatchWithShortIndex(cnt)= L_Long(tmpL(1));
    vecLongMatchWithShort(:,cnt) = vecLong(:,tmpL(1));
end


%% Visualization of Delaunay triangulation representation 

tmpx =[vecShort(1,:)' ;-vecShort(1,:)'];
tmpy =[vecShort(2,:)' ;-vecShort(2,:)'];
tmpz =[vecShort(3,:)' ;-vecShort(3,:)'];
DTShort = delaunayTriangulation(tmpx,tmpy,tmpz);

tmpx =[vecLongMatchWithShort(1,:)' ;-vecLongMatchWithShort(1,:)'];
tmpy =[vecLongMatchWithShort(2,:)' ;-vecLongMatchWithShort(2,:)'];
tmpz =[vecLongMatchWithShort(3,:)' ;-vecLongMatchWithShort(3,:)'];
DTLongMatchWithShort = delaunayTriangulation(tmpx,tmpy,tmpz);

figure(1);
[KShort,~] = convexHull(DTShort);
trisurf(KShort,DTShort.Points(:,1),DTShort.Points(:,2),DTShort.Points(:,3))
axis equal; view([-27.00 1.72])
title('Delaunay triangles for P30')

figure(2);
[KLongMatchWithShort,~] = convexHull(DTLongMatchWithShort);
trisurf(KLongMatchWithShort,DTLongMatchWithShort.Points(:,1),DTLongMatchWithShort.Points(:,2),DTLongMatchWithShort.Points(:,3))
axis equal; view([-27.00 1.72]); 
title('Delaunay triangles for 30 directions down-sampled from P64')


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
hold off;
title('Triangle areas for P30 (blue) and down-sampled version of P64 (red)');


%% Calculating the spatial uniformity index: SUindex
SUindex = std(areaLongMatchWithShort)/std(areaShort);
display(SUindex);

