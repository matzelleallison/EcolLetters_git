%% Movie Test.

%%


%% Set up the movie.


% load parameters
load('metaData.mat')
load('data.mat')
estimData.metaData = metaData;  estimData.data = data;
% set par structure
[par, txtPar] = set_pars(smpl_interval, estimData);
cPar = parscomp_st(par);
vars_pull(cPar);

% set par structure
[par, txtPar] = set_pars(smpl_interval, estimData);

for z = 1:length(foodDensity)
    for y = 1:length(seasonal_amp)
        vfname = sprintf('pX-all_envplus-env_samp_%1.3f-env_X_%1.3f.avi', seasonal_amp(y), foodDensity(z));
        writerObj = VideoWriter(vfname); % Name it.
        writerObj.FrameRate = 60; % How many frames per second.
        open(writerObj);
        for x = 1:length(muPlus)
            fname = sprintf('DEB_out-env_plus_%1.3f-env_samp_%1.3f-env_X_%1.3f.mat', muPlus(x), seasonal_amp(y), foodDensity(z));
            
            load(fname)
            data = DEB_out(1); mainTitle = '';
            f = data.f;
            E = data.E;
            V = data.V;
            T = data.Tb_out;
            X = data.X_out;
            Lw = data.Lw;
            
            initpars = get_initpars(X(1), par);
            M = shapecorr(Lw, initpars.L_b, initpars.L_j); % shape correction factor
            
            
            % TP parameters
            T_1 = par.T1;
            Tpars = [par.T_A par.T_L par.T_H par.T_AL par.T_AH];
            c_T = tempcorr(T+273.15,T_1,Tpars);
            
            % correct for temperature and shape
            p_AmT = c_T .* cPar.p_Am .* M; % specific maximum assimilation rate
            p_XmT = p_AmT ./ par.kap_X;  % max specific ingestion rate
            p_MT = c_T .* par.p_M;
            v_T = c_T .* par.v .* M;          % energy conductance
            
            % get rates
            pX = f .* p_XmT .* V.^(2/3);  pX = real(pX); % ingestion rate
            pX = interp1(data.t, pX, 1:nt);
            p_X(:,x) = pX;
            Tb = interp1(data.t, data.Tb_out, 1:nt);
            T_b(:,x) = Tb;
            pA = par.kap_X .* pX;    % assmilation rate
            pM = p_MT .* V;             % somatic maintenance flux
            pS = pM;	% volume linked somatic maintenance, pS = pM + pT pT = 0 for bvalves)
            %                 pT = p_TT .* V.^(2/3);   % surface area-linked somatic maintenance; maintenance of osmotic regulatory boundaries
            pC = (E./V) .* (par.E_G .* v_T .* V.^(2/3) + pS ) ./ (par.kap .* E./V + par.E_G ); %eq. 2.12 p.37 Kooijman 2010
            pG = (par.kap .* pC - pS)/ par.E_G;        % growth (in structural volume)
            G = pG.^(1/3); % conert to structural length
            Gw = G ./ par.del_M;
            
            for i = 1:length(data.W_w)-1
                GR(i) = (real(data.W_w(i+1)) - real(data.W_w(i)))./real(data.W_w(i));
            end
        end
        
        %% plots
        n_t = 200;
        nnt = 1:n_t;
        
        makepXvid
        makeTbVid
        makeTPCVid
        
        
    end
end