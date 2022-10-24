clear;
clc;
load('rho_TDOA_experiment1.mat');
%% 
%samples per tag
samples{1} = 1:800;
samples{2}= 1:2008;
samples{3} = 1:2006;
samples{4} = 1:2014;

%% 
close all

for j = 1:4
    for i = 1:5 %measurements
        signal = rho{j}(i,:);
    
        signal =inpaint_nans(signal,5);
    
        [B tf]= rmoutliers(signal,'movmedian',50,'SamplePoints',samples{j});
        signal(tf) = NaN;
        signal =inpaint_nans(signal);
        if j == 1 && i == 1
            figure;
            plot(signal);
            hold on;
        end
            

        [B tf]= rmoutliers(signal,'movmedian',10,'SamplePoints',samples{j});
        signal(tf) = NaN;
        signal =inpaint_nans(signal);
       if j == 1 && i == 1
            
            plot(signal);
        end
    
        order = 3;
        framelen = 5; 
    
        sgf = sgolayfilt(signal,order,framelen);
    
        [B tf]= rmoutliers(sgf,'movmedian',5,'SamplePoints',samples{j});
        sgf(tf) = NaN;
        sgf = inpaint_nans(sgf);
        if j == 1 && i == 1
            
            plot(sgf);
        end
    
        filtered{j}(i,:) = sgf;%smoothdata(sgf,'sgolay');

        if j == 1 && i == 1
            
            plot(filtered{1}(1,:));
        end

    end
end
%% 


% for i = 1:5
%    
%     for j = 400:2014
%         smoothData{4}(i,j) = smoothData{4}(i,399);
%     end
% end


%% NLS and EKF
load('AP.mat')
load('smoothData.mat');
NiterMax = 100;
sigma_driving = 0.01:0.01:0.1;

for j = 1:4
    
    position_ekf{j} = zeros(length(samples{j}), 4);
    position_nls{j} = zeros(length(samples{j}), 2);

    
    signal_dummy = filtered{1,j};
   
    R = buildCovarianceMatrix( signal_dummy);

    for i = 1:length(samples{j})
        
        [ uHat{j} , numberOfPerformedIterations{j}] = iterativeNLS(NiterMax, AP(:,1:2),filtered{j}(:,i),i ,R); %this function should return a NiterMax x 2 vector of the estimates provided by the NLS
        iterations{j}(i) = numberOfPerformedIterations{j};
        position_nls{j}(i,:) = uHat{j}(numberOfPerformedIterations{j}, :);
       
    end

    position_ekf{j} = KFTrack(AP(:,1:2),filtered{j},length(samples{j}),position_nls{j}(1,:),sigma_driving(4));%sigma_driving(4)
    position_groundroot{j} = KFTrack(AP(:,1:2),smoothData{j},length(samples{j}),position_nls{j}(1,:),sigma_driving(4));%sigma_driving(4)
    
    
    plot(position_nls{j}(:,1),position_nls{j}(:,2));
    hold on
    plot(position_ekf{j}(:,1),position_ekf{j}(:,2),'--');
    grid on
    
%drawnow;
    
end
%% 
for j = 1:4
for i = 1:length(samples{j})
error{j}(i) = norm(position_groundroot{j}(i,1:2) - position_ekf{j}(i,1:2));

end
end

plot(error{1});

%% 
pos_ekf1(:,1) = position_ekf{1}(:,1);
pos_ekf1(:,2) = position_ekf{1}(:,2);
pos_ekf2(:,1) = position_ekf{2}(:,1);
pos_ekf2(:,2) = position_ekf{2}(:,2);
pos_ekf3(:,1) = position_ekf{3}(:,1);
pos_ekf3(:,2) = position_ekf{3}(:,2);
pos_ekf4(:,1) = position_ekf{4}(:,1);
pos_ekf4(:,2) = position_ekf{4}(:,2);

xi1 = linspace(1,numel(samples{1}),numel(samples{4}));
xi2 = linspace(1,numel(samples{2}),numel(samples{4}));
xi3 = linspace(1,numel(samples{3}),numel(samples{4}));

pos_ekf1p(:,1) = interp1(samples{1},pos_ekf1(:,1),xi1);
pos_ekf1p(:,2) = interp1(samples{1},pos_ekf1(:,2),xi1);
pos_ekf2p(:,1) = interp1(samples{2},pos_ekf2(:,1),xi2);
pos_ekf2p(:,2) = interp1(samples{2},pos_ekf2(:,2),xi2);
pos_ekf3p(:,1) = interp1(samples{3},pos_ekf3(:,1),xi3);
pos_ekf3p(:,2) = interp1(samples{3},pos_ekf3(:,2),xi3);

averageEKF = sum(pos_ekf1p + pos_ekf2p + pos_ekf3p + pos_ekf4,3)/4;


pos_nls1(:,1) = position_nls{1}(:,1);
pos_nls1(:,2) = position_nls{1}(:,2);
pos_nls2(:,1) = position_nls{2}(:,1);
pos_nls2(:,2) = position_nls{2}(:,2);
pos_nls3(:,1) = position_nls{3}(:,1);
pos_nls3(:,2) = position_nls{3}(:,2);
pos_nls4(:,1) = position_nls{4}(:,1);
pos_nls4(:,2) = position_nls{4}(:,2);



pos_nls1p(:,1) = interp1(samples{1},pos_nls1(:,1),xi1);
pos_nls1p(:,2) = interp1(samples{1},pos_nls1(:,2),xi1);
pos_nls2p(:,1) = interp1(samples{2},pos_nls2(:,1),xi2);
pos_nls2p(:,2) = interp1(samples{2},pos_nls2(:,2),xi2);
pos_nls3p(:,1) = interp1(samples{3},pos_nls3(:,1),xi3);
pos_nls3p(:,2) = interp1(samples{3},pos_nls3(:,2),xi3);

averageNLS = sum(pos_nls1p + pos_nls2p + pos_nls3p + pos_nls4,3)/4;


plot(averageNLS(:,1),averageNLS(:,2));
hold on
plot(averageEKF(:,1),averageEKF(:,2));


%%
h = animatedline;
xlim([6 13]);
ylim([16.5 18.5]);
grid on;
%axis square;


samplenum = numel(samples{4});

for i = 1:samplenum
addpoints(h,averageEKF(i,1),averageEKF(i,2));
pp =text(0.01,0.01, strcat(num2str(i*0.1),' sec'), ...
     'Units', 'normalized', ...   % Not depending on the data
     'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'bottom');

drawnow
    
   pause(0.02);

   if i~= 2014
        delete(pp);
    end
end

%% 
%k -sigma driving
%j -tags
%i - 
for k = 1:9
for j = 1:4
    for i = 1:2

    [cdf , bin] = histcounts(error{k,j}(:,i),'Normalization','cdf');
    area(k,j,i) = trapz(linspace(0,1,numel(cdf)),cdf);

    end
end
end

test  = sum(sum(area,3),2);
