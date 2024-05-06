 filename = 'C:\Users\zhaoy\Desktop\jstmaltab\99.stl';    

theta =115; % 角度输入

[Y, X, Z, C] = stlread(filename);



    % 将X, Y, Z矩阵转换为一个单一的M-by-3矩阵
    points = [X(:) Y(:) Z(:)];

% 提取满足条件的点
finalPoints = points;

% 创建新的点云对象
finalPointCloud = pointCloud(finalPoints);

% 点云数据清洗和优化
% 1. 去噪
finalPointCloud = pcdenoise(finalPointCloud);

% 2. 降采样
finalPointCloud = pcdownsample(finalPointCloud, 'gridAverage', 0.01);

% 3. 移除重复点
finalPointCloud = pointCloud(unique(finalPointCloud.Location, 'rows', 'stable'));

basicPoints =finalPointCloud ;



alpha =160;  % alpha 值的选择取决于点云的密度和分布
shp = alphaShape(finalPointCloud.Location, alpha);


% 获取表面三角形
[triangles, vertices] = boundaryFacets(shp);

% 初始化质心和总面积
totalCentroid = [0 0 0];
totalArea = 0;

% 遍历每个三角形，计算质心和面积
for i = 1:size(triangles, 1)
    % 获取三角形的顶点
    triPoints = vertices(triangles(i, :), :);
    
    % 计算三角形的面积
    edge1 = triPoints(2, :) - triPoints(1, :);
    edge2 = triPoints(3, :) - triPoints(1, :);
    crossProduct = cross(edge1, edge2);
    area = 0.5 * norm(crossProduct);
    
    % 计算三角形的质心
    centroid = mean(triPoints, 1);
    
    % 累加计算加权质心
    totalCentroid = totalCentroid + centroid * area;
    
    % 累加总面积
    totalArea = totalArea + area;
end

% 计算整体质心
overallCentroid = totalCentroid / totalArea;

overallCentroid =[-317.025 , -0.297 , -103.952];

% 显示结果
disp('Overall Centroid:');
disp(overallCentroid);









if theta <90

% 使用 fzero 函数计算 c 的值
c = fzero(@(c) volumeDifference1(finalPointCloud.Location, theta, c,  4000000), 0);

% 根据 c 过滤点云
final = finalPointCloud.Location(:,3) < finalPointCloud.Location(:,2) * tand(theta) + c;

elseif theta >90

c = fzero(@(c) volumeDifference2(finalPointCloud.Location, theta, c,  4000000), 0);

% 根据 c 过滤点云
final = finalPointCloud.Location(:,3) > finalPointCloud.Location(:,2) * tand(theta) + c;

else

end



finalPointCloud = finalPointCloud.Location(final, :);




% 创建 alphaShape 对象
alpha = 1000000;  
shp2 = alphaShape(finalPointCloud, alpha);

% 获取三角化的顶点索引
[elements, vertices] = alphaTriangulation(shp2);

% 初始化质心和总体积
totalCentroid2 = [0 0 0];
totalVolume2 = 0;

% 遍历每个四面体，计算质心和体积
for i = 1:size(elements, 1)
    % 获取四面体的顶点
    tetraPoints = vertices(elements(i, :), :);
    
    % 计算四面体的体积
    % 体积公式是 |a·(b×c)| / 6，其中 a, b, c 是四面体的三条边向量
    volume2 = abs(det([tetraPoints(2, :) - tetraPoints(1, :); ...
                       tetraPoints(3, :) - tetraPoints(1, :); ...
                       tetraPoints(4, :) - tetraPoints(1, :)]) / 6);
    
    % 计算四面体的质心
    centroid = mean(tetraPoints, 1);
    
    % 累加计算加权质心
    totalCentroid2 = totalCentroid2 + centroid * volume2;
    
    % 累加总体积
    totalVolume2 = totalVolume2 + volume2;
end

% 计算整体质心
overallCentroid2 = totalCentroid2 / totalVolume2;

% 输出质心
disp('质心坐标:');
disp(overallCentroid2);

% 可视化 alphaShape
figure;
pcshow(finalPointCloud);
hold on;
plot3(overallCentroid(1), overallCentroid(2), overallCentroid(3), 'r*', 'MarkerSize', 10);
plot3(overallCentroid2(1), overallCentroid2(2), overallCentroid2(3), 'b*', 'MarkerSize', 10);
title('Alpha Shape and its Centroid');
xlabel('X');
ylabel('Y');
zlabel('Z');

figure;
plot(shp2);
hold on;
plot3(overallCentroid2(1), overallCentroid2(2), overallCentroid2(3), 'b*', 'MarkerSize', 10);
title('Alpha Shape and its Centroid');
xlabel('X');
ylabel('Y');
zlabel('Z');



figure;
pcshow(basicPoints.Location);
hold on;

plot3(overallCentroid(1), overallCentroid(2), overallCentroid(3), 'r*', 'MarkerSize', 10);
plot3(overallCentroid2(1), overallCentroid2(2), overallCentroid2(3), 'b*', 'MarkerSize', 10);
title('切割后的点云和质心');
xlabel('X轴');
ylabel('Y轴');
zlabel('Z轴');

Mass = 150;    %船体自重



finalans = liju(overallCentroid,overallCentroid2,Mass,theta);



disp(finalans);

function diff = volumeDifference1(points, alpha, c, targetVolume)
    % 根据给定的 c 过滤点云
    filteredPoints = points(:,3) <  points(:,2) * tand(alpha) + c;
    
    if sum(filteredPoints) < 4
        diff = -targetVolume;  % 如果点不够，假设体积为0，差异为负目标体积
        return;
    end

    % 创建点云对象
    pc = pointCloud(points(filteredPoints, :));
    volume = volumeOfConvexHull(pc.Location);
    diff = volume - targetVolume;
end

function diff = volumeDifference2(points, alpha, c, targetVolume)
    % 根据给定的 c 过滤点云
    filteredPoints = points(:,3) >  points(:,2) * tand(alpha) + c;
    
    if sum(filteredPoints) < 4
        diff = -targetVolume;  % 如果点不够，假设体积为0，差异为负目标体积
        return;
    end

    % 创建点云对象
    pc = pointCloud(points(filteredPoints, :));
    volume = volumeOfConvexHull(pc.Location);
    diff = volume - targetVolume;
end






function V = volumeOfConvexHull(points)
% 创建 alphaShape 对象
alpha = 10000;  % alpha 值的选择取决于点云的密度和分布
sshp = alphaShape(points, alpha);

[elements, vertices] = alphaTriangulation(sshp);
V = 0;

% 遍历每个四面体，计算质心和体积
for i = 1:size(elements, 1)
    % 获取四面体的顶点
    tetraPoints = vertices(elements(i, :), :);
    
    % 计算四面体的体积
    % 体积公式是 |a·(b×c)| / 6，其中 a, b, c 是四面体的三条边向量
    volume = abs(det([tetraPoints(2, :) - tetraPoints(1, :); ...
                       tetraPoints(3, :) - tetraPoints(1, :); ...
                       tetraPoints(4, :) - tetraPoints(1, :)]) / 6);
    V = V + volume;
end   
end

function a = liju(A,B,mass,theta)
   D=[0,B(1,2)-A(1,2),B(1,3)-A(1,3)];
   theta=theta+90;
   theta=deg2rad(theta);
   fuli=[0,mass*cos(theta),mass*sin(theta)];
   a=cross(D,fuli);
end






