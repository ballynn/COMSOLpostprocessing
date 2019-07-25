clear
clc
%% Data Injestion
%Ingest Nodes
Data_ingest = readtable('test_mesh.csv');
Data_ingest = table2cell(Data_ingest);
Data = cell2mat(Data_ingest);
x = Data(:,2);
y = Data(:,3);
z = Data(:,4);

%Ingest Results
Data_ingest = readtable('R1.6T_Results.csv');
Data_ingest = table2cell(Data_ingest);
Data2 = cell2mat(Data_ingest);
x1 = Data2(:,2);
y1 = Data2(:,3);
z1 = Data2(:,4);
color1 = Data2(:,5);

%% User Defined Parameters
Threshold = 5; %Particles per unit meter

%% Create Mesh
DT = delaunayTriangulation(y,z);
Center = incenter(DT);

%% Plot Mesh and Centers
figure(1)
plot(y1,z1,'.')
xlabel('Meters (m)')
ylabel('Meters (m)')

figure(2)
plot(y1,z1,'.')
hold on
plot(y,z,'o')
xlabel('Meters (m)')
ylabel('Meters (m)')

figure(3)
tri = delaunay(y,z);
trimesh(tri,y,z)
xlabel('Meters (m)')
ylabel('Meters (m)')

figure(4)
plot(y1,z1,'.')
hold on
plot(Center(:,1),Center(:,2),'x')
tri = delaunay(y,z);
trimesh(tri,y,z)
xlabel('Meters (m)')
ylabel('Meters (m)')
hold off

%Things for later
y1z1 = [y1,z1];
n = length(Center);
Center_Color = zeros([n,1]);

%% KNN Algorithm
%Define training and truth sets
trainingset = y1z1;
truth_train = color1;
testingset = Center;
truth_test = Center_Color;

%Implement K = 1 algorithm
K = 1;
[Final_Centroid_Value,Error] = knnclassify(testingset, trainingset, truth_train, truth_test, K);

%% Area of Triangles
%Looping variable initialization
i = 1;
Area = zeros([n,1]);

%Area Calculation Loop
while i < n
    AX = DT.Points(DT.ConnectivityList(i,1),1);
    AY = DT.Points(DT.ConnectivityList(i,1),2);
    BX = DT.Points(DT.ConnectivityList(i,2),1);
    BY = DT.Points(DT.ConnectivityList(i,2),2);
    CX = DT.Points(DT.ConnectivityList(i,3),1);
    CY = DT.Points(DT.ConnectivityList(i,3),2);
    Area(i) = abs((AX*(BY-CY)+BX*(CY-AY)+CX*(AY-BY))/2);
    i = i+1;
end

%% Particles Per Unit Area
Particles_Area = Final_Centroid_Value./Area;
Final = [Area, Final_Centroid_Value, Particles_Area];
i = 1;
z2 = z;
y2 = y;
tri2 = tri;
Area2 = Area;

while i <= n
    if Particles_Area(i) < Threshold
        Final_Centroid_Value(i) = 0;
        Area(i) = 0;
        tri2(i,:) = [0 0 0];
        y2(i) = 0;
        z2(i) = 0;
        i = i+1;
    else
        i = i+1;
    end
end
    
a = Area2 - Area;
TotalArea = pi*(15^2);
Percent = sum(a)/TotalArea*100;
P1 = 'The total Astronaut Area is ';
P2 =  'm^2.';
P3 = 'Or ';
P4 = '% of the total area studied.';
disp([P1 ,num2str(sum(sum(a))),P2]) 
disp([P3,num2str(Percent),P4])
P5 = 'With ';
P6 = ' particle interactions with the total area studied';
disp([P5,num2str(sum(Final(:,2))),P6])

%% Mesh Area Display

figure(5)
trimesh(tri,y,z,x,'facecolor', 'r')
hold on
isspecial = Final_Centroid_Value > 0;
trimesh(tri(isspecial,:),y,z,x)
title('Max 1.6T, 90Deg, Elements with Particles/Unit Area below Threshold of 5 Particles/m')
xlabel('Meters (m)')
ylabel('Meters (m)')
zlabel('Meters (m)')
hold off

%% Mesh Area Display  - TOTAL GARBAGE      
%Remove Rows with zero values


% tri2( ~any(tri2,2), :) = [];
% y3 = y2;
% z3 = z2;
% y2( ~any(y3-z3,2),:) = [];
% z2( ~any(y3-z3,2),:) = [];
% i = 1;
% q = zeros([n,1]);
% f = zeros([n,1]);
% r = [q, (y3 - z3), f];
% k = 1;
% 
% while i <= n
%     r(i,1) = i;
%     if r(i,2) ~=0
%         r(i,3) = k;
%         k = k + 1;
%         i = i+1;
%     else
%         i = i+1;
%     end
% end
% 
% n = length(tri2);
% i = 1;
% k = 1;
% 
% while i < n
%     while k < 4
%         [tf, rowElement] = ismember(tri2(i,k),r(:,1));
%         output = r(rowElement,3);
%         tri2(i,k) = output;
%         k = k+1;
%     end
%     k = 1;
%     i = i+1;
% end

% figure(5)
% triplot(tri,y,z)
% hold on
% triplot(tri2,y3,z3,'r')
