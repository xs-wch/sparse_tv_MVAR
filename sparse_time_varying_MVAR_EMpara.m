classdef sparse_time_varying_MVAR_EMpara < handle
    %   This version estimate the model using SAEM algorithm
    %   This version for parallel pomputing on BME cluster
    %   2015-07-01 this version do not include sparse prior
    %   2015-09-18 measurement error covariance matrix varys in this version
    
    properties
        chan %number of channels
        len % data length of each trial 
        trial % number of trials

        alpha2
        beta2

        lambda2_new
        poster_precisex

        poster_mux
        m_order
        filename_full
        EEGdata
        predict_EEG
        inovation_process
        Sigma % covriance matrix, sparse matrix
        Sigma0 % prior in Wishart distribution. Initialized with MVAR or adaptive MVAR
        X % vecterized MVAR coefficients, sparse matrix 
        
        R
        Q
        D

        varx2
        W
        
        threshold_X
        threshold_all
        iteration_threshold

        
        synthetic_A1
        synthetic_A2
        synthetic_Atv
        synthetic_SNR
        changeamp
        sin_frq
        
        %%%%%%%%%
        MDDM
        err_matrix
        WVARXW
        EX
        PRESX
        HX
        
        ELambda
        ElnLambda
        HLambda
        nu0
        Lambda0
        
        Elambda
        Elnlambda 
        Hlambda 
        alpha0
        beta0
        beta_for_cal
        lbconst
    end
    
    methods
        function obj = sparse_time_varying_MVAR_EMpara(working_mode,SNR,amp,sin_frq, data_to_analyze,order)
            % working_mode 1: real data
            %              2: synthetic data.It may do not work when order > 2
            %              3: programes for my thesis
        
           
            obj.alpha2 = 10^-3;

            obj.beta2 = 10;
           
            obj.threshold_X = 10^-3;

            obj.threshold_all = 10^-7;
            obj.iteration_threshold =20;

            obj.synthetic_SNR = SNR;
            obj.changeamp = amp;
            if working_mode == 1
                obj.m_order = order;

                obj.EEGdata = data_to_analyze;
                [obj.chan, obj.len, obj.trial] = size(obj.EEGdata);

                obj.synthetic_A1 = [];obj.synthetic_A2 = [];
            end
            
            if working_mode == 2
                obj.m_order = order;
                obj.chan = input('number of channels:');
                obj.len = input('data length:');
                obj.trial = input('number of trials:');
                obj.filename_full = [];
                
                obj.EEGdata = obj.synthtic_generate(200);
                
                
            end
            
            if working_mode == 3
                % simulation using 1 order model, it is difficult to find a higher order stable system

                %obj.EEGdata = obj.synthtic_generate(200);
                %obj.m_order = order;
                if order == 1;
                    obj.chan = 4;
                    obj.len = 500;
                    obj.trial = 30;
                    obj.sin_frq = sin_frq;
                
                    obj.filename_full = [];
                    obj.m_order = 1;
                    [obj.synthetic_Atv, obj.EEGdata] = obj.AR_generation_thesis(200,obj.m_order);
                else if order == 2 % Only one 2-order model used, since higher order may be unstable  
                    obj.chan = 2;
                    obj.len = 500;
                    obj.trial = 30;
                    obj.m_order = order;
                    [obj.synthetic_Atv, obj.EEGdata] = obj.AR_generation_thesis(200,obj.m_order);   
                        
                    end
                end
            end
                
                obj.matrix_D();

                


        end
        
        function EEGdata = synthtic_generate(obj, MVAR_break_point)
            % only support one-order model in this verison
            EEGdata = randn(obj.chan, obj.len, obj.trial);
            
            if (MVAR_break_point <0) || (MVAR_break_point > obj.len) 
                error('break point shoul betwween 1 and length');
            end
            
            [A1, A2, Md1, Md2] = obj.AR_generation_oneorder();
            obj.synthetic_A1 = A1;
            obj.synthetic_A2 = A2;
            
            EEG0 = randn(obj.m_order, obj.chan, obj.trial);
            noise0 = randn(obj.m_order, obj.chan, obj.trial); 
            for itrial = 1:obj.trial
                           
                [EEGdata1, noise1] = vgxsim(Md1, obj.len, [], EEG0(:,:,itrial), noise0(:,:,itrial));
         
                EEGdata(:,1: MVAR_break_point,itrial) = EEGdata1(end-MVAR_break_point+1:end, :)';
                
                [EEGdata2, ~] = vgxsim(Md2, obj.len, [], EEGdata1(end-obj.len+1:end,:),  noise1(end-obj.len+1:end,:));
                
                EEGdata(:,1+MVAR_break_point:end,itrial) = EEGdata2(1:(obj.len - MVAR_break_point),:)';
                
                %%%%%% add white noise
                
                pow_EEGdata = sum(var(squeeze(EEGdata(:,:,itrial)),0,2));
                
                pow_noise = pow_EEGdata/(10^(obj.synthetic_SNR/10));
                
                noise_factor = sqrt(pow_noise/obj.chan);
                
                EEGdata(:,:,itrial) = EEGdata(:,:,itrial)+noise_factor.*randn(obj.chan, obj.len,1);
                
                
                
                
            end

        end
        
        
        
        function [AR_coef_tv,EEGdata] = AR_generation_thesis(obj,MVAR_break_point,p,sigman)
            
            % high order model is not stable - to be solved 
            % p: synthetic model order
                Isstable = false;
                while ~Isstable
                    if nargin <4
                        sigman = eye(obj.chan);
                        
                    end
                    sigman_sqrt =  sigman^-0.5;
                    if p == 1
                        AR_coef1 = 2*rand(obj.chan)-1;
                        [U,S,V] = svd(AR_coef1);
                        sprd = max(abs(diag(S)));
                        if max(abs(diag(S)))>=1
                            S = S./sprd*0.9;
                        end
                        AR_coef1 = U*S*V';
                        AR_coef_temp = repmat(AR_coef1, [1,1,obj.len]);
                        chan_tv1 = randi(obj.chan,1,2);
                        AR_coef_temp(chan_tv1(1),chan_tv1(2),:) = 0.5*sin(2*pi*obj.sin_frq*[1:obj.len]);
                        issamechan = true;
                        while issamechan
                            chan_tv2 = randi(obj.chan,1,2);
                            if ~isequal(chan_tv1,chan_tv2)
                                issamechan = false;
                            end
                        end
                        [varingcoefamp,Ivaringcoef] = max([squeeze(AR_coef_temp(chan_tv2(1),chan_tv2(2),1))+1,...
                            1-squeeze(AR_coef_temp(chan_tv2(1),chan_tv2(2),1))]);  
                        if Ivaringcoef == 1
                            varingcoef = AR_coef_temp(chan_tv2(1),chan_tv2(2),1) - ((varingcoefamp-obj.changeamp)*rand(1)+obj.changeamp);
                        else
                            varingcoef = AR_coef_temp(chan_tv2(1),chan_tv2(2),1) +((varingcoefamp-obj.changeamp)*rand(1)+obj.changeamp);
                        end
                        AR_coef_temp(chan_tv2(1),chan_tv2(2),MVAR_break_point:end) = varingcoef;
                    else if p == 2 
                            AR_coef_temp = [0.4,0.3, 0.35, -0.4;1.2,0.7,-0.3,-0.5 ];
                            AR_coef_temp = repmat(AR_coef_temp, [1,1,obj.len]);
                            AR_coef_temp (2,1,:) = 0.8*sin(2*pi*0.003*[1:obj.len]);


                        else
                            error('higher order simulation in development');
                        end
                    end
                        AR_coef_tv = AR_coef_temp;

                    EEGdata = zeros(obj.chan, obj.len+p, obj.trial);

                    for it = 1:obj.trial
                        for i = p+1:size(EEGdata,2)
                            tempEEGdata = zeros(obj.chan,1);
                            for ip = 1:p
                                tempEEGdata = tempEEGdata+squeeze(AR_coef_tv(:,(ip-1)*obj.chan+1:ip*obj.chan,i-p))*squeeze(EEGdata(:,i-ip,it));
                            end

                            EEGdata(:,i,it) = tempEEGdata +sigman_sqrt*randn(obj.chan,1);
                        end
                    end
                    EEGdata = EEGdata(:,p+1:end,:);
                    if max(abs(EEGdata(1,:,1))) <= 20
                        Isstable = true;
                    end

                    pow_EEGdata = sum(var(EEGdata(:,:),0,2));

                    pow_noise = pow_EEGdata/(10^(obj.synthetic_SNR/10));

                    noise_factor = sqrt(pow_noise/obj.chan);
                    %noise_tensor = noise_factor.*randn(obj.chan, obj.len,);
                    EEGdata = EEGdata+noise_factor.*randn(obj.chan, obj.len,obj.trial);
                end

        end
        
            
        
        
        function [AR_coef1, AR_coef2, Md1, Md2] = AR_generation_oneorder(obj)
            %AR_coef1 = 2*sprand(obj.chan, obj.chan*obj.m_order, 0.2)-1;% may be constrained
            %%%% most nonzero coefficients on the diagnoal
            isStable =true;
            num_try = 0;
            while(isStable)
                num_try = num_try+1;
                AR_coef1 = 2*rand(obj.chan)-1;
                [U,S,V] = svd(AR_coef1);
                sprd = max(abs(diag(S)));
                if max(abs(diag(S)))>=1
                    S = S./sprd*0.9;
                end
                AR_coef1 = U*S*V';
                    
                AR_coef2 = 2*rand(obj.chan)-1;
                [U,S,V] = svd(AR_coef2);
                sprd = max(abs(diag(S)));
                if max(abs(diag(S)))>=1
                    S = S./sprd*0.9;
                end
                AR_coef2 = U*S*V';
                
                clear U S V sprd


                A1 = cell(1,1);
               
                A2 = cell(1,1);

                A1{1} = AR_coef1;
                A2{1} = AR_coef2;
                
                
                epislon = randn(obj.chan);
                epislon = epislon*epislon';
                Md1 = vgxset('n',obj.chan,'nAR',1,'AR',A1,'Q',5*epislon);
                Md2 = vgxset('n',obj.chan,'nAR',1,'AR',A2,'Q',5*epislon);

                [isStable1,~,eigAR1,~] = vgxqual(Md1);
                [isStable2,~,eigAR2,~] = vgxqual(Md2);

                isStable = ~(isStable1&&isStable2);
                disp(['the ',num2str(num_try),' trial to generate synthetic data']);
            end

        end
        
        

        
        
        
        function get_init_new(obj, method)
            % method 1: randomly generated matrix
            % method 2: using all data to estimate AR
            addpath /home/wch/Documents/tv_AR/arfit

            p = obj.m_order;
            c = obj.chan;
            N = obj.len;
            if method == 1
                obj.lambda2_new = 1*ones((N-p-1),1);

                EEG_perm = permute(obj.EEGdata,[2 1 3]);                
                [~, ~, obj.Sigma0, ~, ~, ~] = arfit(EEG_perm, ceil(p/2), p, 'sbc', 'zero');
                
                obj.Sigma = repmat(obj.Sigma0,[1,1,(N-p)]);
            end
            
            if method == 2 
             
                EEG_perm = permute(obj.EEGdata,[2 1 3]);
                
                [~, ~, obj.Sigma0, ~, ~, ~] = arfit(EEG_perm, ceil(p/2), 15, 'sbc', 'zero');
                obj.Sigma = repmat(obj.Sigma0,[1,1,(N-p)]);
                %%%% initialize the EM algorithm
                stv_mvar_initmodel = sparse_time_varying_MVAR(1,obj.EEGdata,p);
 %               stv_mvar_initmodel.Sigma = obj.Sigma;
                init_X = stv_mvar_initmodel.estimate_model();
                
                figure_flag = false;
                if figure_flag
                    A = reshape(full(init_X),[c,c*p,N-p]);
                    
                    figure
                    for i = 1:4
                        for j = 1:4
                            subplot(4,4,(i-1)*4+j)
                            chanindex = [i,j];
                            plot(1:497, squeeze(A(chanindex(1),chanindex(2)+8,:)))
                            hold on
                            plot([1:500], [repmat(obj.synthetic_A1(chanindex(1),chanindex(2)),1,200), repmat(obj.synthetic_A2(chanindex(1),chanindex(2)),1,300)],'r-')
                            ylim([-1, 1])
                            xlim([0, 500])
                            xlabel('time(ms)')
                            
                        end
                    end
                end
                
                
                
  %              obj.lambda1_new  =  ((c^2*p*(N-p)+2*obj.alpha1-2)/(init_X'*init_X+2*obj.beta1)).*ones((N-p),1);
                
                tempD = obj.D'*obj.D;
                tempxDx = zeros(1);
                for i = 1:N-p-1
                    
                    samplexn = init_X((i-1)*c^2*p+1:(i+1)*c^2*p);
                    %%%%%% modified at 9.6
                    %obj.lambda2_new(i) = (c^2*p+2*obj.alpha2-2)/(samplexn'*tempD*samplexn+2*obj.beta2);
                    tempxDx = tempxDx + samplexn'*tempD*samplexn;
                end
                obj.lambda2_new = ((c^2*p*(N-p-1)+2*obj.alpha2-2)/(tempxDx+2*obj.beta2))*ones(N-p-1,1);
                %%%% this version use lambda2 for all data points
                
%                disp(['initial lambda1: ',num2str(obj.lambda1_new(1))]);
                disp(['initial lambda2: ',num2str(obj.lambda2_new(1))]);

            end
            
        end
        

        
        function matrix_RQ(obj)
            p = obj.m_order;
            c = obj.chan;
            obj.R = sparse(c^2*p*(obj.len-p),c^2*p*(obj.len-p));
            obj.Q = sparse(1, c^2*p*(obj.len-p));

            for i = p+1:obj.len
          
                WW = sparse(c^2*p,c^2*p);
                sW = sparse(1,c^2*p);
                I = false(c^2*p*(obj.len-p),1);
                I((i-p-1)*c^2*p+1:(i-p)*c^2*p) = true;
                Sigman = squeeze(obj.Sigma(:,:,i-p));
                for k = 1:obj.trial
                    W_nk = obj.W{i-p,k};
                    Sigma_W = Sigman\W_nk;
                    WW = WW + W_nk'*Sigma_W;
                    sW = sW + sparse(squeeze(obj.EEGdata(:,i,k)))'*Sigma_W;
                end

                obj.R(I,I) = obj.R(I,I) + WW;
                obj.Q(I) = obj.Q(I) + sW;
            end
           
        end
        
        function matrix_W(obj)
            p = obj.m_order;
            c = obj.chan;
            N = obj.len;
            
            obj.W = cell(N-p,obj.trial);
            
            tempW_nk = zeros(c, c^2*p);
            for iwin = p+1:obj.len
                EEG_in_win = obj.EEGdata(:,iwin-1:-1:iwin-p,:);
                for k = 1:obj.trial
                    for i = 1:p
                        tempW_nk(:,(i-1)*c^2+1:i*c^2) = kron(squeeze(EEG_in_win(:,i,k))',speye(c));
                    end
                    obj.W{iwin-p,k} = sparse(tempW_nk);
                end
            end
       
        end
        
        function matrix_D(obj)
            p = obj.m_order;
            c = obj.chan;
            obj.D = [sparse(-eye(c^2*p)), sparse(eye(c^2*p))];
            
        end
        
        function M_n = matrix_Mn(obj,i)
            p = obj.m_order;
            c = obj.chan;
            
            M_n = sparse(c^2*p,c^2*p*(obj.len-p));
            M_n(:,(i-p-1)*c^2*p+1:(i-p)*c^2*p) = sparse(eye(c^2*p));
        end
        

        function sigma = matrix_Sigma(obj,mode,poster_varx)
            
            p = obj.m_order;
            c = obj.chan;
            np = obj.len-p;

            tempmux = obj.poster_mux;
          %  tempvarx = poster_varx;
            
            obj.predict_EEG = zeros(c,np,obj.trial);
            obj.inovation_process = zeros(c,np,obj.trial);
            
%             sigma1 = zeros(c);
%             sigma2 = zeros(c);
            for i = p+1:obj.len

                EEG_target = squeeze(obj.EEGdata(:,i,:));

                X_n = sparse(tempmux((i-p-1)*c^2*p+1:(i-p)*c^2*p));

                epsilon_n = zeros(c, obj.trial);
                poster_varxn = zeros(c);
              
                temp = sparse(squeeze(poster_varx(:,:,i-p)));
               
                for k = 1:obj.trial
                    W_nk = obj.W{i-p,k};
                    obj.predict_EEG(:,i-p,k) =  W_nk*X_n;
                    epsilon_n(:,k) = EEG_target(:,k) -obj.predict_EEG(:,i-p,k);
                    obj.inovation_process(:,i-p,k) = epsilon_n(:,k);
                    poster_varxn = poster_varxn + W_nk*temp*W_nk';
                    
                end
             %   sigma1 = sigma1 + epsilon_n*epsilon_n';
             %   sigma2 = sigma2 + poster_varxn;
            
                sigma =epsilon_n*epsilon_n'+poster_varxn;
            
                if mode == 1 % the bayesian version with conjugate prior
                    if obj.trial == 1
                        
                        obj.Sigma(:,:,i-p) = (sigma+obj.Sigma0);
                    else
                        obj.Sigma(:,:,i-p) = (sigma+obj.Sigma0)/(obj.trial-1);
                    end
                end
                
                if mode == 2 % the ML version
                    obj.Sigma(:,:,i-p) = (sigma)/(obj.trial);
                end
                
                if mode == 3 % the bayesian version with conjugate prior
                    if obj.trial == 1
                        if i == p+1
                            obj.Sigma(:,:,i-p) = (sigma+obj.Sigma0);
                        else
                            obj.Sigma(:,:,i-p) = (sigma+obj.Sigma(:,:,i-p-1));
                        end
                    else
                        if i == p+1
                            obj.Sigma(:,:,i-p) = (sigma+obj.Sigma0)/(obj.trial - 1);
                        else
                            obj.Sigma(:,:,i-p) = (sigma+obj.Sigma(:,:,i-p-1))/(obj.trial - 1);
                        end
                    end
                end
            
            end
            
        end
        

        
        function matrix_varx2_new(obj) % \sum_n lambda_2n*Mn,n+1'*D'*D*Mn,n+1
            p = obj.m_order;
            c = obj.chan;
            temp_varx = sparse(c^2*p*(obj.len-p),c^2*p*(obj.len-p));
            tempD = obj.D'*obj.D;
            
            
            for i = p+1:obj.len-1
                temp_varx((i-p-1)*c^2*p+1:(i-p+1)*c^2*p,(i-p-1)*c^2*p+1:(i-p+1)*c^2*p)...
                = temp_varx((i-p-1)*c^2*p+1:(i-p+1)*c^2*p,(i-p-1)*c^2*p+1:(i-p+1)*c^2*p)+obj.lambda2_new(i-p)*tempD;
                
            end
            obj.varx2 = temp_varx;
        end
        

        function E_step_new(obj)
            
            obj.matrix_RQ();
 
            obj.matrix_varx2_new();
            
            obj.poster_precisex = obj.varx2+obj.R;
            obj.poster_mux = obj.poster_precisex\obj.Q';
            
            
            
        end
        

        
        function M_step_SAEM(obj,niter)
            
            %%%% the parameter lambda1 is identical across the time period
            %%%% using SAEM to update the parameter
            p = obj.m_order;
            c = obj.chan;
            N = obj.len;
        
               
            r = obj.generate_r(niter);
            r = 0;
            if r~=0
                [L, prank]= chol(obj.poster_precisex);
                sample = L\randn(c^2*p*(N-p),1) + obj.poster_mux;
            end
            

            tempD = obj.D'*obj.D;
            tempxDx_SEM = zeros(N-p-1,1);
            tempxDx_EM = zeros(N-p-1,1);
            poster_varx = zeros(c^2*p,c^2*p,N-p);
            vxn = cell(N-p-1,1);

            parfor i = 1:N-p-1

                xn = obj.poster_mux((i-1)*c^2*p+1:(i+1)*c^2*p);
                tempM = sparse(1:c^2*p*2,(i-1)*c^2*p+1:(i+1)*c^2*p,ones(c^2*p*2,1), c^2*p*2, c^2*p*(N-p));
                vxn{i}= tempM*(obj.poster_precisex\tempM');
   
                tempxDx_EM(i) = xn'*tempD*xn+trace(vxn{i}*tempD);

                
            end

            
            for i = 1:N-p-1
                
                if i ==1
                    poster_varx(:,:,1) = vxn{1}(1:c^2*p,1:c^2*p);
                    poster_varx(:,:,2) = vxn{1}(c^2*p+1:2*c^2*p,c^2*p+1:2*c^2*p);
                else
                    poster_varx(:,:,i+1) = vxn{i}(c^2*p+1:2*c^2*p,c^2*p+1:2*c^2*p);
                end
            end
            
            tempxDx_SEM = sum(tempxDx_SEM);
            tempxDx_EM  = sum(tempxDx_EM);
            
            if r~=0
            lambda2_SEM = (c^2*p*(N-p-1)+2*obj.alpha2-2)/( tempxDx_SEM +2*obj.beta2);
            end
            lambda2_EM = (c^2*p*(N-p-1)+2*obj.alpha2-2)/(tempxDx_EM +2*obj.beta2);
            
            
            if r~=0
            obj.lambda2_new = ((1-r)*lambda2_EM+r*lambda2_SEM)*ones(N-p-1,1);
            else
              obj.lambda2_new = (lambda2_EM)*ones(N-p-1,1);  
            end
            
            disp(['lambda',num2str(lambda2_EM)]);
            obj.matrix_Sigma(1,poster_varx);
            
        end
        
        
        function r = generate_r(obj, niter)
            if niter < 5
                %r = 1/(0.1*niter + 0.9);
                r = 1;
            else if niter < 18
                    r = 1/(niter-4);
                else
                    r = 0;
                end
                
            end
         
          
        end
        
        
        function A = estimate_model_EM(obj)
            
            obj.get_init(1);
%            obj.get_initX(1);
            obj.matrix_RQ();

            disp('parameters:')
            disp(['alpha1: ',num2str(obj.alpha1)])
            disp(['alpha2: ',num2str(obj.alpha2)])
            disp(['beta1: ',num2str(obj.beta1)])
            disp(['beta2: ',num2str(obj.beta2)])
            disp(['model order: ', num2str(obj.m_order)])
            disp(['EEGdata size: ',num2str(obj.chan),'*',num2str(obj.len),'*',num2str(obj.trial)])

            p = obj.m_order;
            c = obj.chan;
            cp = c*p;
            np = obj.len-p;
            deltaX = 1;
            X_old = sparse(ones(c^2*p*np,1));
            niter = 0;
            while (deltaX > obj.threshold_all) && (niter < obj.iteration_threshold) 
                niter = niter+1;
                tic;
                obj.E_step();
                t1 = toc;
                tic;
                obj.M_step();
                t2 = toc;
                tic;
                obj.matrix_RQ();
                t3 = toc;
                
                deltaX = norm(obj.poster_mux-X_old)/norm(X_old);
                X_old = obj.poster_mux;
                
                disp(['iteration: ', num2str(niter), ' deltaX: ',num2str(deltaX)]);
                disp(['E-step takes ',num2str(t1),'s']);
                disp(['M-step takes ',num2str(t2),'s']);
                disp(['Update Q and R takes ',num2str(t3),'s']);
                
            end
            
            A = obj.poster_mux;
            
            A(abs(A)<obj.threshold_X) = 0;
            
            A = reshape(full(A),[c,cp,np]);
 
            
        end
        
        
        
         function A = estimate_model_EM_new(obj)
            
            obj.get_init_new(2);

            obj.matrix_RQ();
      

            
            disp('parameters:')
            disp(['alpha1: ',num2str(obj.alpha1)])
            disp(['alpha2: ',num2str(obj.alpha2)])
            disp(['beta1: ',num2str(obj.beta1)])
            disp(['beta2: ',num2str(obj.beta2)])
            disp(['model order: ', num2str(obj.m_order)])
            disp(['EEGdata size: ',num2str(obj.chan),'*',num2str(obj.len),'*',num2str(obj.trial)])

            p = obj.m_order;
            c = obj.chan;
            cp = c*p;
            np = obj.len-p;
            deltaX = 1;
            X_old = sparse(ones(c^2*p*np,1));
            niter = 0;
            while (deltaX > obj.threshold_all) && (niter < obj.iteration_threshold)
                niter = niter+1;
                tic;
                obj.E_step_new();
                t1 = toc;
                tic;
                obj.M_step_constlambda1();
                t2 = toc;
                tic;
                obj.matrix_RQ();
                t3 = toc;
                
                deltaX = norm(obj.poster_mux-X_old)/norm(X_old);
                X_old = obj.poster_mux;
                
                disp(['iteration: ', num2str(niter), ' deltaX: ',num2str(deltaX)]);
                disp(['E-step takes ',num2str(t1),'s']);
                disp(['M-step takes ',num2str(t2),'s']);
                disp(['Update Q and R takes ',num2str(t3),'s']);
                
            end
            
            A = obj.poster_mux;
            
            A(abs(A)<obj.threshold_X) = 0;
            
            A = reshape(full(A),[c,cp,np]);
 
            
         end
        
         
         function A = estimate_model_SAEM(obj)
             
             obj.get_init_new(1);
             obj.matrix_W();
        
%             obj.matrix_RQ();

             disp('parameters:')
%             disp(['alpha1: ',num2str(obj.alpha1)])
             disp(['alpha2: ',num2str(obj.alpha2)])
%             disp(['beta1: ',num2str(obj.beta1)])
             disp(['beta2: ',num2str(obj.beta2)])
             disp(['model order: ', num2str(obj.m_order)])
             disp(['EEGdata size: ',num2str(obj.chan),'*',num2str(obj.len),'*',num2str(obj.trial)])
             
             p = obj.m_order;
             c = obj.chan;
             cp = c*p;
             np = obj.len-p;
             deltaX = 1;
             X_old = sparse(ones(c^2*p*np,1));
             niter = 0;
             fig_flag = false;
             matlabpool open 3
             while (deltaX > obj.threshold_all) && (niter < obj.iteration_threshold)
                 niter = niter+1;
                 tic;
                 obj.E_step_new();
                 t1 = toc;
                 tic;
                 obj.M_step_SAEM(niter);
                 t2 = toc;
%                  tic;
%                  obj.matrix_RQ();
%                  t3 = toc;
                 
                 deltaX = norm(obj.poster_mux-X_old)/norm(X_old);
                 X_old = obj.poster_mux;
                 
                 disp(['iteration: ', num2str(niter), ' deltaX: ',num2str(deltaX)]);
                 disp(['E-step takes ',num2str(t1),'s']);
                 disp(['M-step takes ',num2str(t2),'s']);
                 
                 if fig_flag
                    figure
                    for ichan = 1:obj.chan
                    subplot(obj.chan,2,2*ichan-1);
                    plot(squeeze(obj.predict_EEG(ichan,:,1)),'r-','LineWidth',2);
                    hold on
                    plot(squeeze(obj.EEGdata(ichan,p+1:end,1)),'b-')
                    subplot(obj.chan,2,2*ichan);
                    plot(squeeze(obj.inovation_process(ichan,:,1)),'b-')
                    
                    end
                     
                     
                 end
        %         disp(['Update Q and R takes ',num2str(t3),'s']);
                 
             end
             matlabpool close
             A = obj.poster_mux;
             
             A(abs(A)<obj.threshold_X) = 0;
             
             A = reshape(full(A),[c,cp,np]);
             
             
         end
         
         function [A, lowerbound_new]= estimate_model_VB(obj)
             
             p = obj.m_order;
             c = obj.chan;
             np = obj.len-p;
             obj.init_VB(1,p);
             
             dl = 1;
             lowerbound_old = 1;
             update_count = 0;
             minupdate_flag = true;
             %matlabpool open
             while ((dl > 5*10^-3)&& (update_count < 50)) || minupdate_flag
             update_count =  update_count+1;
             if update_count > 15 %%%% maximum iteration times 15
                 minupdate_flag = false;
             end
             tic    
             obj.updateX();
             obj.updatelambda();
             obj.updateLambda();
             lowerbound_new = obj.lowerbound_VB();
             dl = abs((lowerbound_new-lowerbound_old)/ lowerbound_old );
             disp(['the ',num2str(update_count),' update: \Delta LB: ', num2str(lowerbound_new-lowerbound_old),'  LB: ',num2str(lowerbound_new)]);
             lowerbound_old =  lowerbound_new ;
             toc
             end
             %matlabpool close
             A = obj.EX;
             
             A(abs(A)<obj.threshold_X) = 0;
             
             A = reshape(full(A),[c,c*p,np]);
             
         end
         
         function init_VB(obj,pmin,pmax) % initialize the VB algorithm
             obj.alpha0 = 0;
             obj.beta0 = 0;
                     
             addpath /home/wch/Documents/tv_AR/arfit
             addpath /home/wch/Documents/SSVEP/WOSSPA_Mathworks_v2
             EEG_perm = permute(obj.EEGdata,[2 1 3]);
             [~, ~, obj.Sigma0, sbc, ~, ~] = arfit(EEG_perm, pmin, pmax);
             obj.ELambda = inv(obj.Sigma0);
             obj.nu0 = obj.chan;
             obj.Lambda0 = obj.ELambda;
             %%%%%% init matrix D/matrix Mn/matrix Mn'D'*D*Mn
             obj.matrix_D();
             obj.matrix_W();
             obj.matrix_MDDM();
             
            UC = 0.01;% Update mode of the process noise covariance matrix Q
            %%%%%% initlize time varying MVAR model with AAR
            [A_TV,~,~,~] = mvaar(squeeze(EEG_perm(:,:,1)),obj.m_order,UC,2); % Multivariate Adaptive AR estimation using Kalman filter --> BioSig toolbox
            CH = 2;
            A_TV_reshape = reshape(A_TV', CH*obj.m_order, CH, size(A_TV,1));
            A_TV_reshape =  A_TV_reshape(:,:,obj.m_order+1:end);
            for i = 1 : size(A_TV_reshape,3)
                A_TV2(:,:,i) = A_TV_reshape(:,:,i)';
            end
            A = A_TV2(:);
            
            obj.Elambda = 0.5*obj.chan^2*obj.m_order*(obj.len-obj.m_order-1)/(A'*obj.MDDM*A );
            
             p = obj.m_order;
             c = obj.chan;
             N = obj.len;
             K = obj.trial;
              %%%% calculate the constant in Lower bound
            obj.lbconst = -K*(N-p)*c/2*log(2*pi)-gammaln(10^-5)-0.5*obj.nu0*obj.logdet(obj.Lambda0,'chol')...
                 -obj.nu0*c/2*log(2) - sum(gammaln((obj.nu0+1-[1:c])/2));
         end
         
         function matrix_MDDM(obj) % \sum_n lambda_2n*Mn,n+1'*D'*D*Mn,n+1
            p = obj.m_order;
            c = obj.chan;
            temp_varx = sparse(c^2*p*(obj.len-p),c^2*p*(obj.len-p));
            tempD = obj.D'*obj.D;
            
            
            for i = p+1:obj.len-1
                temp_varx((i-p-1)*c^2*p+1:(i-p+1)*c^2*p,(i-p-1)*c^2*p+1:(i-p+1)*c^2*p)...
                = temp_varx((i-p-1)*c^2*p+1:(i-p+1)*c^2*p,(i-p-1)*c^2*p+1:(i-p+1)*c^2*p)+tempD;
                
            end
            obj.MDDM = temp_varx;
         end
         
         function updateX(obj)
             
            p = obj.m_order;
            c = obj.chan;
            obj.R = sparse(c^2*p*(obj.len-p),c^2*p*(obj.len-p));
            obj.Q = sparse(1, c^2*p*(obj.len-p));
     

            for i = p+1:obj.len
            
                WW = sparse(c^2*p,c^2*p);
                sW = sparse(1,c^2*p);
                I = false(c^2*p*(obj.len-p),1);
                I((i-p-1)*c^2*p+1:(i-p)*c^2*p) = true;
              
                for k = 1:obj.trial
                    W_nk = obj.W{i-p,k};
                    tempLambdaW_nk = sparse(obj.ELambda*W_nk);
                    WW = WW + W_nk'*tempLambdaW_nk;
                    sW = sW + sparse(squeeze(obj.EEGdata(:,i,k)))'*tempLambdaW_nk;
                end


                obj.R(I,I) = obj.R(I,I) + 0.5*(WW+WW');
          
                obj.Q(I) = obj.Q(I) + sW;
            end

            
            obj.PRESX = obj.Elambda*obj.MDDM + obj.R;  
            obj.EX =  obj.PRESX \ obj.Q';
            obj.HX = -0.5*obj.logdet(obj.PRESX,'chol')+c^2*p*(obj.len-p)/2*(1+log(2*pi));

            if size(obj.PRESX)< 10^4
                obj.beta_for_cal = obj.EX'*obj.MDDM*obj.EX + trace(obj.PRESX \ obj.MDDM);
            else    
                obj.beta_for_cal = obj.EX'*obj.MDDM*obj.EX + obj.trace_invA_B();
            end
            

            [obj.err_matrix, obj.WVARXW] = obj.matrix_err_VB();
         end
         
         function updatelambda(obj)
             p = obj.m_order;
             c = obj.chan;
             N = obj.len;
           
             alpha_new = obj.alpha0 + 0.5*c^2*p*(N-p-1);
             
             beta_new = 0.5*(obj.beta_for_cal)+obj.beta0;
          
             obj.Elambda = alpha_new/beta_new;
             obj.Elnlambda = psi(alpha_new)-log(beta_new);
             obj.Hlambda = gammaln(alpha_new)-(alpha_new-1)*psi(alpha_new)-log(beta_new)+alpha_new;
         end
         
         
         function updateLambda(obj)
             p = obj.m_order;
             c = obj.chan;
             N = obj.len;
             K = obj.trial;
             nu_new = K*(N-p)+obj.nu0;
             
            % [err_matrix, WVARXW] = obj.matrix_err_VB();
             W_new = inv(inv(obj.Lambda0)+obj.err_matrix+obj.WVARXW);
             obj.ELambda = nu_new*W_new;
             obj.ElnLambda = obj.logdet(W_new,'chol') + c*log(2) + sum(psi((nu_new+1-[1:c])/2));
             obj.HLambda = nu_new/2*obj.logdet(W_new,'chol') +  nu_new*c/2*log(2) + sum(gammaln((nu_new+1-[1:c])/2))...
                 -(nu_new-c-1)/2* obj.ElnLambda+ nu_new*c/2;
             
         end
         
         function lowerbound = lowerbound_VB(obj)
             p = obj.m_order;
             c = obj.chan;
             N = obj.len;
             K = obj.trial;
             
       
             term1 = K*(N-p)/2 * obj.ElnLambda;
           
             term2 = -0.5*(trace((obj.err_matrix + obj.WVARXW)*obj.ELambda));
             term3 = -c^2*p*(N-p-1)/2*log(2*pi)+c^2*p*(N-p-1)/2*obj.Elnlambda-obj.Elambda/2*(obj.beta_for_cal);
             term4 = (obj.alpha0-1)*obj.Elnlambda - obj.beta0*obj.Elambda+...
                 (obj.nu0-c-1)/2*obj.ElnLambda-0.5*trace(obj.Lambda0\obj.ELambda);
             lowerbound = term1+term2+term3+term4+obj.HX+obj.Hlambda+obj.HLambda+obj.lbconst;
         end
         

         function [err_matrix, WVARXW] = matrix_err_VB(obj)
            p = obj.m_order;
            c = obj.chan;
            np = obj.len-p;

            obj.predict_EEG = zeros(c,np,obj.trial);
            obj.inovation_process = zeros(c,np,obj.trial);
            err_matrix = zeros(c);
            WVARXW = zeros(c);
            
            for i = p+1:obj.len 

                EEG_target = squeeze(obj.EEGdata(:,i,:));

                X_n = sparse(obj.EX((i-p-1)*c^2*p+1:(i-p)*c^2*p));

                epsilon_n = zeros(c, obj.trial);
                temp_varxn = zeros(c);
                M_n = sparse(obj.matrix_Mn(i));
                temp = M_n*(obj.PRESX \ M_n');
               
                for k = 1:obj.trial
                    W_nk = obj.W{i-p,k};
                    obj.predict_EEG(:,i-p,k) =  W_nk*X_n;
                    epsilon_n(:,k) = EEG_target(:,k) -obj.predict_EEG(:,i-p,k);
                    obj.inovation_process(:,i-p,k) = epsilon_n(:,k);
                    temp_varxn = temp_varxn + W_nk*temp*W_nk';
                    
                end

                sigma_n = epsilon_n*epsilon_n';
                err_matrix = err_matrix + sigma_n;   
                WVARXW = WVARXW + temp_varxn;
                
       
            
            end
         end
         
         function t= trace_invA_B(obj) % trace(A\B).When B is very large, A\B have to be calculated seperately
             
             p = obj.m_order;
             c = obj.chan;
             np = obj.len-p;
             ncol= size(obj.MDDM,2);
             tempt = zeros(ncol,1);
             for i = 1:c^2*p
                 temp = obj.PRESX\obj.MDDM(:,(i-1)*np+1:i*np);
                 tempt((i-1)*np+1:i*np) = diag(temp,-(i-1)*np);
             end
             t = sum(tempt);
         end
         
         function t= trace_invA_B_para(obj,A,B)
             
             p = obj.m_order;
             c = obj.chan;
             np = obj.len-p;
             tempt = zeros(c^2*p,1);
             parfor i = 1:c^2*p
                 temp = A\B(:,(i-1)*np+1:i*np);
                 tempt(i) = sum(diag(temp,-(i-1)*np));
             end
             t = sum(tempt);
         end 
         
        function v = logdet(obj,A, op)
        %LOGDET Computation of logarithm of determinant of a matrix
        %
        %   v = logdet(A);
        %       computes the logarithm of determinant of A. 
        %
        %       Here, A should be a square matrix of double or single class.
        %       If A is singular, it will returns -inf.
        %
        %       Theoretically, this function should be functionally 
        %       equivalent to log(det(A)). However, it avoids the 
        %       overflow/underflow problems that are likely to 
        %       happen when applying det to large matrices.
        %
        %       The key idea is based on the mathematical fact that
        %       the determinant of a triangular matrix equals the
        %       product of its diagonal elements. Hence, the matrix's
        %       log-determinant is equal to the sum of their logarithm
        %       values. By keeping all computations in log-scale, the
        %       problem of underflow/overflow caused by product of 
        %       many numbers can be effectively circumvented.
        %
        %       The implementation is based on LU factorization.
        %
        %   v = logdet(A, 'chol');
        %       If A is positive definite, you can tell the function 
        %       to use Cholesky factorization to accomplish the task 
        %       using this syntax, which is substantially more efficient
        %       for positive definite matrix. 
        %
        %   Remarks
        %   -------
        %       logarithm of determinant of a matrix widely occurs in the 
        %       context of multivariate statistics. The log-pdf, entropy, 
        %       and divergence of Gaussian distribution typically comprises 
        %       a term in form of log-determinant. This function might be 
        %       useful there, especially in a high-dimensional space.       
        %
        %       Theoretially, LU, QR can both do the job. However, LU 
        %       factorization is substantially faster. So, for generic
        %       matrix, LU factorization is adopted. 
        %
        %       For positive definite matrices, such as covariance matrices,
        %       Cholesky factorization is typically more efficient. And it
        %       is STRONGLY RECOMMENDED that you use the chol (2nd syntax above) 
        %       when you are sure that you are dealing with a positive definite
        %       matrix.
        %
        %   Examples
        %   --------
        %       % compute the log-determinant of a generic matrix
        %       A = rand(1000);
        %       v = logdet(A);
        %
        %       % compute the log-determinant of a positive-definite matrix
        %       A = rand(1000);
        %       C = A * A';     % this makes C positive definite
        %       v = logdet(C, 'chol');
        %

        %   Copyright 2008, Dahua Lin, MIT
        %   Email: dhlin@mit.edu
        %
        %   This file can be freely modified or distributed for any kind of 
        %   purposes.
        %

        % argument checking

        assert(isfloat(A) && ndims(A) == 2 && size(A,1) == size(A,2), ...
        'logdet:invalidarg', ...
        'A should be a square matrix of double or single class.');

        if nargin <= 2
        use_chol = 0;
        else
        assert(strcmpi(op, 'chol'), ...
        'logdet:invalidarg', ...
        'The second argument can only be a string ''chol'' if it is specified.');
        use_chol = 1;
        end

        % computation

        if use_chol
        v = 2 * sum(log(diag(chol(A))));
        else
        [~, U, P] = lu(A);
        du = diag(U);
        c = det(P) * prod(sign(du));
        v = log(c) + sum(log(abs(du)));
        end
 
        end
    end
    
end

