function [x_hat] = KFTrack(AP,rho,sampleNumber,init,sigma_driving)


R = buildCovarianceMatrix( rho );
numberOfAP = 6;
samplingTime = 0.1;
UE_init = [init(1), init(2), 0, 0];
UE_init_COV = diag([0.04^2, 0.04^2, 0.04^2, 0.04^2]);
%figure;
%hold on
%plotCovariance( UE_init_COV(1:2,1:2) , UE_init(1,1) , UE_init(1,2)  , 3 ,'Initialization');

x_hat = NaN( sampleNumber , 4);
P_hat = NaN( 4 , 4 , sampleNumber );



L = [0.5*samplingTime^2, 0; 
                  0,samplingTime^2;
                  samplingTime,0; 
                  0,samplingTime];
              
% Q = sigma_driving^2 *(L * L');
% 
% F = [1, 0, samplingTime, 0;
%      0, 1, 0, samplingTime;
%      0, 0, 1, 0;
%      0, 0, 0, 1];

 %update over time
for time = 1 : sampleNumber

    if time > 51
        sampleVar = mean(var(rho(:,time-50:time)'));
        disp(sampleVar);
    else
        sampleVar = 1000;
    end

 if sampleVar < 10e-3

     F = [1, 0, 0, 0;
          0, 1, 0, 0;
          0, 0, 0, 0;
          0, 0, 0, 0];

     Q = 0;

 else

     F = [1, 0, samplingTime, 0;
     0, 1, 0, samplingTime;
     0, 0, 1, 0;
     0, 0, 0, 1];
     
     Q = sigma_driving^2 *(L * L');

 end

    %prediction
    if time == 1
        x_pred =  UE_init';
        P_pred = UE_init_COV;
    else
        x_pred = F * x_hat(time-1,:)';
        P_pred = F * P_hat(:,:,time-1) *F' + Q;
    end

    H = buildJacobianMatrixH(x_pred(1:2)' , AP);
    H = [H zeros(numberOfAP-1, 2)];%no direct measurements of the velocity

    %update
    G = P_pred * H' * inv( H*P_pred*H' + R);
    x_hat(time,:) = x_pred + G * (rho( : , time )'  - measurementModel(x_pred(1:2)',AP) )';
    P_hat(:,:,time) = P_pred - G * H * P_pred;

    
    
   
end
 
end

 %plotCovariance( P_pred(1:2,1:2)  , x_pred(1) , x_pred(2)  , 3 , 'Prior');
    %plotCovariance( P_hat(1:2,1:2,time)  , x_hat(time,1) , x_hat(time,2)  , 3 , 'Update');

%     %plot evolution
%     fig = figure(11);
%     fig.WindowState = 'maximized';
%     plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255), hold on
%     plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.30,0.75,0.93] ), 
%     plot( UE(time,1) , UE(time,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0.30,0.75,0.93],'MarkerFaceColor',[0.50,0,0] ),
%     plotCovariance( P_pred(1:2,1:2)  , x_pred(1) , x_pred(2)  , 3 , 'Prior');
% %     axis equal
%     xlim([UE(time,1)-100 UE(time,1)+100])
%     ylim([UE(time,2)-50 UE(time,2)+50])
%     xlabel('[m]'), ylabel('[m]');
%     legend('AP','True UE (all)','True UE (current)','Cov. pred.')
%     title('Time: ' , num2str(time) )
%     pause(0.3)
%     plot( x_hat(:,1) , x_hat(:,2) , '-g','Marker','s')
%     plotCovariance( P_hat(1:2,1:2,time)  , x_hat(time,1) , x_hat(time,2)  , 3 , 'Update');
%     legend('AP','True UE (all)','True UE (current)','Cov. pred.','KF est.','Cov. upd.')
% %     axis equal
%     xlim([UE(time,1)-100 UE(time,1)+100])
%     ylim([UE(time,2)-50 UE(time,2)+50])
%     
%     pause(0.3)
%     hold off
