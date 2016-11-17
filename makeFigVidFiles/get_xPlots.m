function get_xPlots(data, mainTitle, Legend)
%% load univariate data from estimation
% cd('/Users/AlliMatzelle/Manuscripts/Matzelleetal_LO/DEBmodel/DEBMODEL_Mcalifornianus_033116')
% run_data = load('data.mat'); run_data = run_data.data;
% cd('/Users/AlliMatzelle/Manuscripts/Matzelleetal_LO/AFTER_REESTIMATION_033116')

%% plot DEB quantities
% data: structure containing n number of tables created in analysis.m to be compared

dataFields = fieldnames(data);
datetime = daten2datet(data(1).t);
dtnum = datenum(data(1).t);
data = data(1);

% figure
% for g = 1:size(data,2)
subplot(3,5,1)
	plot(dtnum, data.Tb_out, 'LineWidth', 1)
    hold on
%     ylim([0 30])
    xlabel('Date')
    ylabel('T_b, C')
    title('Body Temperature')
    
subplot(3,5,2)
    plot(dtnum, data.X_out, 'LineWidth', 1)
    hold on
%     ylim([0 80])
    xlabel('Date')
    ylabel('Chl-a (ug/L)')
    title('food density')
    
subplot(3,5,3)
    plot(dtnum, data.f, 'LineWidth', 1)
    hold on
    ylim([0 1])
    xlabel('Date')
    ylabel('-')
    title('functional response')

subplot(3,5,4)
    plot(dtnum, data.E, 'LineWidth', 1)
    hold on
%     ylim([0 8e4])
    xlabel('Date')
    ylabel('E (J)')
    title('Reserve Energy')

subplot(3,5,5)
    plot(dtnum, data.V, 'LineWidth', 1)
    hold on
%     ylim([0 150])
    xlabel('Date')
    ylabel('V (cm^3)')
    title('Structural Volume')    
    
subplot(3,5,6)
    plot(dtnum, data.E_H, 'LineWidth', 1)
    hold on
%     ylim([0 45])
    xlabel('Date')
    ylabel('E_H (J)')
    title('Maturity')     
    
subplot(3,5,7)
    plot(dtnum, data.E_R, 'LineWidth', 1)
    hold on
%     ylim([0 1600])
    xlabel('Date')
    ylabel('E_R (J)')
    title('Reproduction')     
    
subplot(3,5,8)
    plot(dtnum, data.e, 'LineWidth', 1)
    hold on
%     ylim([0 1.5])
    xlabel('Date')
    ylabel('e (-)')
    title('Scaled reserve density')    
% end
% subplot(3,5,5)
legend(Legend, 'Location', 'Best')
% hold on
sp1 = mtit(mainTitle, 'fontsize', 18);
%% Plot observable quantities
% figure
% for j = 1:length(dataFields)
% subplot(3,5,1)
%     plot(dtnum, data.Tb_out, 'LineWidth', 1)
%     hold on
% %     ylim([0 30])
%     xlabel('Date')
%     ylabel('T_b, C')
%     title('Body Temperature')
%     
%     dspan = xlim;
%     
% subplot(3,5,2)
%     plot(dtnum, data.X_out, 'LineWidth', 1)
%     hold on
% %     ylim([0 80])
%     xlabel('Date')
%     ylabel('Chl-a (ug/L)')
%     title('food density')
    
subplot(3,5,9)
    plot(dtnum, data.Lw, 'LineWidth', 1)
    hold on
%     ylim([0 18])
    xlabel('Date')
    ylabel('Lw (cm)')
    title('Physical (shell) length') 
    
subplot(3,5,10)
    plot(dtnum, data.W_w, 'LineWidth', 1)
    hold on
%     ylim([0 150])
    xlabel('Date')
    ylabel('Ww (g)')
    title('Wet weight') 
    
subplot(3,5,11)
    plot(dtnum, data.W_d, 'LineWidth', 1)
    hold on
%     ylim([0 15])
    xlabel('Date')
    ylabel('Wd (g)')
    title('Dry weight')     
    
subplot(3,5,12)
    plot(data.Lw, data.W_w, 'LineWidth', 1)
    hold on
%     ylim([0 100])
    xlabel('Lw (cm)')
    ylabel('Ww (g)')
    title('Length-Wet weight')   
    
subplot(3,5,13)
    plot(data.Lw, data.W_d, 'LineWidth', 1)
    hold on
%     ylim([0 15])
    xlabel('Lw (cm)')
    ylabel('Wd (g)')
    title('Length-Dry weight') 
    
subplot(3,5,14)
E_Hp = 40.78; kap_R = 0.95;   % FROM SET_PARS_NEW.M
% Fraw = zeros(size(data.rawE_R,1),size(data.rawE_R,2));
E_0 = data.E(1);
E_R = data.E_R;
E_H = data.E_H;
i_sp = find(and((E_R == 0),(E_H>=E_Hp)));
i_sp = i_sp - 1; % look at preceeding line for E_R value before spawning
F = kap_R .* E_R ./ E_0;  % fecundity = egg number

    plot(dtnum(i_sp), F(i_sp), '*')
    hold on
%     xlim(dspan)
%     ylim([0 7e6])
    xlabel('Date')
    ylabel('F (#)')
    title('Spawning events') 
    
    subplot(3,5,15)
    ylim([0 5])
    text(0,x,strcat('+',int2str(muPlus(x)),'; TRO = ',int2str(data.TRO),'; #repyrs = ',int2str(data.nrepyrs),'; time 2 maturity = ',int2str(data.time2mat)))
    set(gca, 'visible', 'off')
    hold on
% end  

% subplot(3,5,3)
%     t_tL = run_data.tL(:,1); 
%     t_tL = t_tL + datenum(datevec(dtnum(1)));
%     t_tL = datetime(datevec(t_tL));  % add to start date of current datetime
%     plot(t_tL, run_data.tL(:,2), '.r', 'MarkerSize', 20) 
%     
% subplot(3,5,6)
%     plot(run_data.LWw(:,1), run_data.LWw(:,2), '.k', 'MarkerSize', 20)    
%     
% subplot(3,5,7) 
%     plot(run_data.LWd(:,1), run_data.LWd(:,2), '.k', 'MarkerSize', 20)
%     hold on
%     plot(run_data.LWd2(:,1), run_data.LWd2(:,2), '.k', 'MarkerSize', 20)
%     
% subplot(3,5,8)
%     legend(Legend, 'Location', 'Best')
%     
% sp2 = mtit(mainTitle, 'fontsize', 18);

    
