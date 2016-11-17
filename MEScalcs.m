function MEScalcs(seasonal_amp, muPlus, foodDensity)
% p = path;
% path(p, '/Users/AlliMatzelle/Manuscripts/Matzelleetal_LO/050716/DEB_out_3sdNoise')
% figure
% hold on

W_w = []; T_b = []; s_a = []; f_d = [];
for l = 1:length(foodDensity)
    fd = repmat(foodDensity(l), 10, 1);
    for k = 1:length(seasonal_amp)
        sa = repmat(seasonal_amp(k), 10, 1);
        for j = 1:length(muPlus)
            Tb = repmat(muPlus(j), 10, 1);
            fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f.mat', ...
                muPlus(j), seasonal_amp(k), foodDensity(l));
            load(fname)
            
            % pull out Ww data from each replicate model run
            for i = 1:10
                data = DEB_out(i); mainTitle = '';
                Ww(i) = real(max(data.W_w)); % store in matrix
            end
            
            W_w = [W_w; Ww'];
            T_b = [T_b; Tb];
            s_a = [s_a; sa];
            f_d = [f_d; fd];
            
%             plot(sa,Ww,'o')
%             hold on
        end
    end
end
A = [s_a f_d T_b W_w];
B = sortrows(A,[1 2 3]);

s_a = B(:,1);
f_d = B(:,2);
T_b = B(:,3);
W_w = B(:,4);  
%% hedgesg (wrong)
[h,p,ci,stats] = ttest2(W_w(s_a == 0 & f_d == 2.62) , ...
    W_w(s_a == 10& f_d == 2.62))
%
MES = mes(W_w(s_a == 0 & f_d == 2.62) , ...
    W_w(s_a == 10 & f_d == 2.62), 'hedgesg') %, 'doPlot','true');
% hedgesg = MES.hedgesg
% hedgesgCi = MES.hedgesgCi


%% 2-way MES

% main comparison: compare levels of amplitude
cW = [1; 0-1/3; 0-1/3; 0-1/3];
[stats, resTable] = mes2way(W_w(T_b == 8), [s_a(T_b == 8) f_d(T_b == 8)], {}, 'fName', {'amplitude','food'}, ...
    'doDataPlot',1);

[stats, resTable] = mes2way(W_w, [s_a f_d], 'eta2', 'fName', {'amplitude','food'}, ...
    'doDataPlot',1);
% 'cWeight', cW,


% path(p)



end

