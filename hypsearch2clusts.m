%1/10/2022
%Let's do a greedy algorithm to to the hypersearch
%i.e. iterate over all points, find the biggest improvement, remove
%repeat until...(no improvment?)
%
%28/9/2022
%Modified to use k-means instead of CMA-ES. Also has the multiple removal
%of points mod

% 6/9/2022
% Modifying the cec19n5.m to do different datasets, etc.
%------

% Used to generate results in my CEC'19 paper
% Clustering on the Ruspini dataset, using CMA-ES to optimize
% A simple local/greedy meta optimizer to try removing single data points
% to find problems that CMA-ES takes longer to solve.
%---
%ptsrem = 100;
%%load s1_data;
%ruspinirangemin = 4;
%ruspinirange = 156-4;
%%s1n = normalize(s1,'range');
%%datarange=1;
%%datamin=0;
%Try sorting the data to help the hypsearch
%s1n = sortrows(s1n);
%%data = s1n;
%


%-------------
% %Generate data from 2 circles overlapping
% dim=2;
% q=200;
% data = zeros(q,dim);
% r = rand(q,1);
% theta = 2*pi*rand(q,1);
% data(:,1) = sqrt(r).*cos(theta);
% data(:,2) = sqrt(r).*sin(theta);
%-------------


dim=1;
q = 100;
sobolGen = sobolset(dim); % 1 dimension
%interval 1
points1 = rescale(net(sobolGen,q),0,2)
%interval 2
points2 = rescale(net(sobolGen,q),3,5)
%combined intervals
data= [points1; points2];

%-------------



datarange=1;
datamin=0;


k = 2;
dim = k*size(data,2);
n = size(data,1);
numreps = 5;
%Record index of points as they are removed
%recind = zeros(size(data,1),1);
%Record the dataset every so often - how often?
recfreq = 30;
recdata = zeros(n,size(data,2),round(n/recfreq));

%rusind = ones(n,1);
%Randomly shuffle data?
permdata = data(randperm(size(data,1)),:);
wdata = permdata;

clusterror = 0;
%for z=1:numreps
    %[iterations, clusterrortrial, centers, runflag] = mykmeans(wdata,k);
    %clusterror = clusterror + clusterrortrial;
%end

for z=1:numreps
    [XMIN,FMIN,ceval] = cmaes('fitnessclustd', (datarange/2)*ones(dim,1)+datamin, (datarange/3), [], wdata);
    clusterror = clusterror + FMIN;
end
clusterror = clusterror/numreps;
%Change SSE to MSE
clusterror = clusterror/n;
%for i=1:3
%    [XMIN,FMIN,ceval] = cmaes('fitnessclustd', (datarange/2)*ones(dim,1)+datamin, (datarange/3), [], wdata);
%    cit(i) = ceval;
%end
%fminprev = FMIN
%cevalprev = mean(cit);
%cevalprev = iterations;
cevalprev = clusterror;
figure(1);
hold on;

j = n;
t = 0;
while j>1
    figure(1);
    plot(t,cevalprev,'*');
    drawnow;
    %figure(2)
    %colormap(winter(size(wdata,1)));
    %hold on;
    %plot(wdata(:,1),wdata(:,2),'o');
    %scatter(wdata(:,1),wdata(:,2),25,j*ones(size(wdata,1),1),"filled");
    %drawnow;
    %Remove point(s) from the dataset
    % Note that we don't want to run over the end of the dataset
    % So our picking point only goes up to n-ptsrem
    %rusind(i)=0;
    %i = randi(size(wdata,1));
    %i = randi(n-ptsrem);
    %Remove points i,i+1,... from dataset
    %permdata = wdata(randperm(size(wdata,1)),:);
    %wdata = permdata;
    %wdatanew = [wdata(1:i-1,:); wdata(i+ptsrem:n,:)];
    %--
    % Loop over the dataset, removing each point and seeing what happens
    errtobeat = cevalprev;
    rmptidx = 0; %Index of point to remove
    for i=1:j
        wdatanew = [wdata(1:i-1,:); wdata(i+1:j,:)];
        clusterror = 0;
        %for z=1:numreps
            %Now we need to shuffle the data for the kmeans initialization
            %wdatanew = wdatanew(randperm(size(wdatanew,1)),:);
            %[iterations, clusterrortrial, centers, runflag] = mykmeans(wdatanew,k);
            %clusterrortrial
            %clusterror = clusterror + clusterrortrial;
        %end
        
        for z=1:numreps
           [XMIN,FMIN,ceval] = cmaes('fitnessclustd', (datarange/2)*ones(dim,1)+datamin, (datarange/3), [], wdatanew);
           clusterror = clusterror + FMIN;
        end
        
        clusterror = clusterror/numreps;
        %Change SSE to MSE
        clusterror = clusterror/(j-1);
        figure(1);
        plot(t,clusterror,'.');
        drawnow;
    %for i=1:3
    %    [XMIN,FMIN,ceval] = cmaes('fitnessclustd', (datarange/2)*ones(dim,1)+datamin, (datarange/3), [], wdata);
    %    cit(i) = ceval;
    %end
    %if FMIN > fminprev
    %if mean(cit) > cevalprev
    %if iterations > cevalprev
        %i
        %clusterror
        if clusterror > errtobeat
            %****************
            %Point is the new contender to be removed
            %wdata = wdatanew;
            %n = n-ptsrem;
        %fminprev = FMIN
        %cevalprev = mean(cit)
        %cevalprev = iterations
            errtobeat = clusterror
            rmptidx = i
        end
        %Else we move on and try the next point
    end
    %If we didn't find any improving point, remove a random one!
    if rmptidx == 0
        rmptidx = randi(size(wdata,1));
        disp('Removed a random point');
    end
    %Now we have the point to remove
    figure(2)
    colormap(copper(size(wdata,1)));
    hold on;
    %plot(wdata(:,1),wdata(:,2),'o');
    scatter(wdata(rmptidx,1),wdata(rmptidx,2),50,j,"filled");
    drawnow;

    %if (mod((t+1),recfreq) == 0)
        %Record the current dataset
        recdata(1:size(wdata,1),:,floor(t/recfreq)+1) = wdata;
    %end

    wdata = [wdata(1:rmptidx-1,:); wdata(rmptidx+1:j,:)];
    cevalprev = errtobeat;
    j = j-1
    t = t+1;
    %recidx(t) = rmptidx;
    %n = n-ptsrem;
end

csvwrite('recdata.csv', recdata)
