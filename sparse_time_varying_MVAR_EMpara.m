classdef sparse_time_varying_MVAR_EMpara < handle
    %   This version estimate the model using SAEM algorithm
    %   This version for parallel pomputing on BME cluster
    %   2015-07-01 this version do not include sparse prior
    %   2015-09-18 measurement error covariance matrix varys in this version
    
    properties
        chan
        len
        trial
%        alpha1
        alpha2
%        beta1
        beta2
%        lambda1
%        lambda2
%        lambda1_new
        lambda2_new
        poster_precisex
%        poster_varx
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
%        varx1
        varx2
        W
        
        threshold_X
        threshold_all
        iteration_threshold
%        epsilon
        
        synthetic_A1
        synthetic_A2
        synthetic_SNR
        changeamp
        
    end
    
    methods
        function obj = sparse_time_varying_MVAR_EMpara(working_mode,SNR,amp,data_to_analyze,order)
            % working_mode 1: real data
            %              2: synthetic data
            %              3: random generated data, just for test
        
            
 %           obj.alpha1 = 10^-3;
            obj.alpha2 = 10^-3;
 %           obj.beta1 = 10^15;
            obj.beta2 = 10;
           
            obj.threshold_X = 10^-3;

            obj.threshold_all = 10^-7;
            obj.iteration_threshold =50;

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
                obj.m_order = 1; % simulation using 1 order model, it is difficult to find a higher order stable system
                obj.chan = 4;
                obj.len = 500;
                obj.trial = 30;
                
                obj.filename_full = [];
                
                obj.EEGdata = obj.synthtic_generate(200);
                %obj.m_order = order;
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
        
        
        
        function [AR_coef1, AR_coef2, Md1, Md2] = AR_generation_highorder(obj)
            
            % high order model is not stable - to be solved 
            
            AR_coef1 = cell(obj.m_order,1);
            AR_coef2 = cell(obj.m_order,1);
  %          A1 = cell(obj.m_order,1);
  %          A2 = cell(obj.m_order,1);
            
            [AR_coef1{1}, AR_coef2{1}, ~, ~] = obj.AR_generation_oneorder();
            for i = 2:obj.m_order
                [temp_AR1, temp_AR2, ~, ~] = obj.AR_generation_oneorder();
                
                AR_coef1{i} = AR_coef1{i-1}*temp_AR1;
                AR_coef2{i} = AR_coef2{i-1}*temp_AR1;
            end
            

            
            epislon = eye(obj.chan);
            Md1 = vgxset('n',obj.chan,'nAR',obj.m_order,'AR',AR_coef1 ,'Q',epislon);
            Md2 = vgxset('n',obj.chan,'nAR',obj.m_order,'AR',AR_coef2 ,'Q',epislon);
            
            [isStable1,~,eigAR1,~] = vgxqual(Md1);
            [isStable2,~,eigAR2,~] = vgxqual(Md2);

                
            if ~(isStable1 && isStable2)
                  disp(['eigAR1: ',num2str(eigAR1)]);
                  disp(['eigAR2: ',num2str(eigAR2)]);
                  error('model is not stable');
                
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
%                 AR_coef1 = diag(2*rand(obj.chan,1)-1);
%                 index_zero = randi(obj.chan,[round(obj.chan*0.25),1]);
%          
%                 for i = 1:length(index_zero)
%                     tempAR_coef = AR_coef1(:,index_zero(i));
%                     tempAR_coef = circshift(tempAR_coef, randi(obj.chan-1,1));
%                     AR_coef1(:,index_zero(i)) = tempAR_coef;
%                 end
% 
%                 %%%% some nonzero coefficients off the diagnoal
%                 AR_coef1 = AR_coef1+(2*rand(obj.chan)-1)...
%                     .*(sprand(obj.chan, obj.chan, 0.3/obj.chan) ~= 0);
%                 %%%% keep all the elements in the matrix lower than 0.9
%                 %%%% this help improve stability
%                 I = find(abs(AR_coef1) > 0.90);
%                 if ~isempty(I)
%                     AR_coef1(I) = sign(randn(length(I),1)).*(0.90*rand(length(I),1));
%                 end
%                 
%                 %%%% some elements changes after 200ms
%                 %number_nonzero = sum(sum(AR_coef1 ~= 0));
%                 number_change = 1;
% 
% 
%                 AR_coefI = (2*rand(obj.chan)-1).*(AR_coef1~=0);
%                 temp_threshold = sort(AR_coefI(:));
%                 change_threshold = temp_threshold(number_change);
%                 change_I = find(AR_coefI <= change_threshold);
%            
%                 if ~isempty(change_I)
%                     tempARcoef = AR_coef1(:);
%                     for i = 1:length(change_I)
%                         if (tempARcoef(change_I(i)) +obj.changeamp) >=1
%                             tempARcoef(change_I(i)) = tempARcoef(change_I(i)) - obj.changeamp;
%                         end
%                         if (tempARcoef(change_I(i)) -obj.changeamp) <= -1
%                             tempARcoef(change_I(i)) = tempARcoef(change_I(i)) + obj.changeamp;
%                         end
%                         if (tempARcoef(change_I(i)) -obj.changeamp >= -1) && (tempARcoef(change_I(i)) +obj.changeamp<= 1)
%                             temp = shuffle([tempARcoef(change_I(i))+obj.changeamp,tempARcoef(change_I(i))-obj.changeamp]);
%                             tempARcoef(change_I(i)) = temp(1);
%                         end
%                     end
%                     AR_coef2 = reshape(tempARcoef,[obj.chan,obj.chan]);
%                 else
%                     pause
%                     AR_coef2 = AR_coef1;
%                             
%                 end

                A1 = cell(1,1);
               
                A2 = cell(1,1);

                A1{1} = AR_coef1;
                A2{1} = AR_coef2;
                
                
                epislon = randn(obj.chan);
                epislon = epislon*epislon';
                Md1 = vgxset('n',obj.chan,'nAR',1,'AR',A1,'Q',5*epislon);
                Md2 = vgxset('n',obj.chan,'nAR',1,'AR',A2,'Q',10*epislon);

                [isStable1,~,eigAR1,~] = vgxqual(Md1);
                [isStable2,~,eigAR2,~] = vgxqual(Md2);

                isStable = ~(isStable1&&isStable2);
                disp(['the ',num2str(num_try),' trial to generate synthetic data']);
            end

        end
        
%         function [A1, A2] = coef_refactor(obj, AR_coef1, AR_coef2)
%             sign_A1 = sign(AR_coef1);
%             sign_A2 = sign(AR_coef2);
%       
%             abs_A1 = abs(AR_coef1);
%             abs_A2 = abs(AR_coef2);
%             
%             max_A = max(max(abs_A1(AR_coef1 ~=0), abs_A2(AR_coef2 ~=0)));
%             min_A = min(min(abs_A1(AR_coef1 ~=0), abs_A2(AR_coef2 ~=0)));
%             
%             abs_A1 = (abs_A1-min_A)/(max_A-min_A)*0.8+0.15;
%             abs_A2 = (abs_A2-min_A)/(max_A-min_A)*0.8+0.15;
%             
%             A1 = sign_A1.*abs_A1;
%             A2 = sign_A2.*abs_A2;
% 
%             
%             
%         end
        
% 
%         function get_init(obj, method)
%             % method 1: randomly generated matrix
%             % method 2: using all data to estimate AR
%             addpath /home/wch/Documents/tv_AR/arfit
%             p = obj.m_order;
%             c = obj.chan;
%             N = obj.len;
%             if method == 1
%                 obj.X = sparse(randn(c^2*p*(N  - p),1));
%             end
%             
%             if method == 2 
%              
%                 EEG_perm = permute(obj.EEGdata,[2 1 3]);
%                 
%                 [~, A, obj.Sigma, ~, ~, ~] = arfit(EEG_perm, ceil(p/2), p, 'sbc', 'zero');
%                 [A_d1, A_d2] = size(A);
%                 if A_d1 ~= c
%                     error('dimension do not match');
%                 end
%                 est_order = floor(A_d2/A_d1);
%                 if est_order <p
%                     initX = zeros(c, c*p);
%                     initX(:,1: c*est_order) = A;
%                     
%                 else
%                     initX = A(:, 1:c*p);                    
%                     
%                 end
%                 %%%% initialize with totally non-smooth image
%                
%                 obj.poster_mux = sparse(repmat(initX(:), N - p,1));
%                 obj.poster_precisex = sparse(eye(c^2*p*(N -p)));
%                 obj.poster_varx = obj.poster_precisex;
%                 
%                 temp1 = trace(obj.varx1*obj.poster_varx) + obj.poster_mux'*obj.varx1*obj.poster_mux+2*obj.beta1;
%                 obj.lambda1 = (c^2*p*(N-p)+2*obj.alpha1-2)/temp1;
%             
%                 temp2 = trace(obj.varx2*obj.poster_varx) + obj.poster_mux'*obj.varx2*obj.poster_mux+2*obj.beta2;
%                 obj.lambda2 = (c^2*p*(N-p-1)+2*obj.alpha2-2)/temp2;
%                 
%                 obj.matrix_varx1();
%                 obj.matrix_varx2();
%                 
%                 disp(['initial lambda1: ',num2str(obj.lambda1)]);
%                 disp(['initial lambda2: ',num2str(obj.lambda2)]);
% 
%             end
%             
%         end
        
        
        
        function get_init_new(obj, method)
            % method 1: randomly generated matrix
            % method 2: using all data to estimate AR
            addpath /home/wch/Documents/tv_AR/arfit
      %      obj.matrix_D();
            p = obj.m_order;
            c = obj.chan;
            N = obj.len;
            if method == 1

               % obj.lambda1_new = 30*rand((N-p),1);
                obj.lambda2_new = 100*rand((N-p-1),1);

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
            
          %  Sigma_inv = sparse(inv(obj.Sigma));
            for i = p+1:obj.len
                %EEG_in_win = obj.EEGdata(:,i-1:-1:i-p,:);
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
            tempvarx = poster_varx;
            
            obj.predict_EEG = zeros(c,np,obj.trial);
            obj.inovation_process = zeros(c,np,obj.trial);
            
%             sigma1 = zeros(c);
%             sigma2 = zeros(c);
            for i = p+1:obj.len

                EEG_target = squeeze(obj.EEGdata(:,i,:));

                X_n = sparse(tempmux((i-p-1)*c^2*p+1:(i-p)*c^2*p));

                epsilon_n = zeros(c, obj.trial);
                poster_varxn = zeros(c);
              
                temp = sparse(tempvarx((i-p-1)*c^2*p+1:(i-p)*c^2*p,(i-p-1)*c^2*p+1:(i-p)*c^2*p));
               
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
        
%         function matrix_varx1(obj) % \sum_n Mn'*Mn
%             p = obj.m_order;
%             c = obj.chan;
% %             temp_varx = sparse(zeros(c^2*p*(obj.len-p)));
% %             
% %             for i = p+1:obj.len
% %                 M_n = obj.matrix_Mn(i);
% %                 temp_varx = temp_varx+M_n'*M_n;
% %                 
% %             end
% %             obj.varx1 = temp_varx;
%              obj.varx1 = sparse(eye(c^2*p*(obj.len-p)));
%             
%         end
        
%         function matrix_varx1_new(obj) % \sum_n lambda_1n*Mn'*Mn
%             p = obj.m_order;
%             c = obj.chan;
%             temp = repmat(obj.lambda1_new', c^2*p,1);
%             obj.varx1 = spdiags(temp(:),0,c^2*p*(obj.len-p),c^2*p*(obj.len-p));
%             
%         end
        
%         function matrix_varx2(obj) % \sum_n Mn,n+1'*D'*D*Mn,n+1
%             p = obj.m_order;
%             c = obj.chan;
%             temp_varx = sparse(c^2*p*(obj.len-p),c^2*p*(obj.len-p));
%             
%             for i = p+1:obj.len-1
%                 M_n_n1 = [obj.matrix_Mn(i);obj.matrix_Mn(i+1)];
%                 temp_varx = temp_varx+M_n_n1'*obj.D'*obj.D*M_n_n1;
%                 
%             end
%             obj.varx2 = temp_varx;
%         end
        
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
        
%         function E_step(obj)
%         
%             obj.poster_precisex = obj.lambda1*obj.varx1+obj.lambda2*obj.varx2+obj.R;
%             obj.poster_varx = inv(obj.poster_precisex);
%             obj.poster_mux = obj.poster_precisex\obj.Q';
%             
%         end
        
        function E_step_new(obj)
            
            obj.matrix_RQ();
        %    obj.matrix_varx1_new();
            obj.matrix_varx2_new();
            
            obj.poster_precisex = obj.varx2+obj.R;
            obj.poster_mux = obj.poster_precisex\obj.Q';
            
            
            
        end
        
%         function M_step(obj)
%             p = obj.m_order;
%             c = obj.chan;
%             N = obj.len;
%             
% %             obj.matrix_varx1();
% %             obj.matrix_varx2();
%             
%             temp1 = trace(obj.varx1*obj.poster_varx) + obj.poster_mux'*obj.varx1*obj.poster_mux+2*obj.beta1;
%             obj.lambda1 = (c^2*p*(N-p)+2*obj.alpha1-2)/temp1;
%             
%             temp2 = trace(obj.varx2*obj.poster_varx) + obj.poster_mux'*obj.varx2*obj.poster_mux+2*obj.beta2;
%             obj.lambda2 = (c^2*p*(N-p-1)+2*obj.alpha2-2)/temp2;
%             
%             obj.matrix_Sigma(2);
%             
%         end
        
%         function M_step_new(obj)
%             p = obj.m_order;
%             c = obj.chan;
%             N = obj.len;
%             %%%%% to be continue
%             for i = 1:N-p
%                xn = obj.poster_mux((i-1)*c^2*p+1:i*c^2*p);
%                vxn =  obj.poster_varx((i-1)*c^2*p+1:i*c^2*p,(i-1)*c^2*p+1:i*c^2*p);
%                obj.lambda1_new(i) = (c^2*p+2*obj.alpha1-2)/(xn'*xn+trace(vxn)+2*obj.beta1);
%             end
%             
%             tempD = obj.D'*obj.D;
%             for i = 1:N-p-1
%                xn = obj.poster_mux((i-1)*c^2*p+1:(i+1)*c^2*p);
%                vxn =  obj.poster_varx((i-1)*c^2*p+1:(i+1)*c^2*p,(i-1)*c^2*p+1:(i+1)*c^2*p);
%                obj.lambda2_new(i) = (c^2*p+2*obj.alpha2-2)/(xn'*tempD*xn+trace(vxn*tempD)+2*obj.beta2);
%             end
%             
% 
%             obj.matrix_Sigma(2);
%             
%         end
%         
%         
%         
%         function M_step_constlambda1(obj)
%             
%             %%%% the parameter lambda1 is identical across the time period
%             p = obj.m_order;
%             c = obj.chan;
%             N = obj.len;
%             %%%%% to be continue
% %             for i = 1:N-p
% %                xn = obj.poster_mux((i-1)*c^2*p+1:i*c^2*p);
% %                vxn =  obj.poster_varx((i-1)*c^2*p+1:i*c^2*p,(i-1)*c^2*p+1:i*c^2*p);
% %                obj.lambda1_new(i) = (c^2*p+2*obj.alpha1-2)/(xn'*xn+trace(vxn)+2*obj.beta1);
% %             end
%             
%             temp1 = trace(obj.poster_varx) + obj.poster_mux'*obj.poster_mux+2*obj.beta1;
%             obj.lambda1_new = (c^2*p*(N-p)+2*obj.alpha1-2)/temp1*ones((N-p),1);
% 
% 
% 
% 
%             tempD = obj.D'*obj.D;
%             for i = 1:N-p-1
%                xn = obj.poster_mux((i-1)*c^2*p+1:(i+1)*c^2*p);
%                vxn =  obj.poster_varx((i-1)*c^2*p+1:(i+1)*c^2*p,(i-1)*c^2*p+1:(i+1)*c^2*p);
%                obj.lambda2_new(i) = (c^2*p+2*obj.alpha2-2)/(xn'*tempD*xn+trace(vxn*tempD)+2*obj.beta2);
%             end
%             
% 
%             obj.matrix_Sigma(2);
%             
%         end
        
        
        function M_step_SAEM(obj,niter)
            
            %%%% the parameter lambda1 is identical across the time period
            %%%% using SAEM to update the parameter
            p = obj.m_order;
            c = obj.chan;
            N = obj.len;
            
            if (c^2*p*(N-p))^2 < 10*(1024^3*8) 
                poster_varx = obj.poster_precisex\eye(c^2*p*(N-p));
            else
                % to construct a sparse covariance matrix using matlab
                % sparse matrix format
                tempi1 = repmat((1:2*c^2*p)',1,c^2*p);
                tempi2 = zeros(3*c^2*p,c^2*p,N-p-2);
                tempi = repmat((1:3*c^2*p)',1,c^2*p);
                tempi3 = tempi1 + c^2*p*(N-p-2);
                
                %index_i = zeros((c^2*p)^2*(3*(N-p)-2));
                %index_j = zeros((c^2*p)^2*(3*(N-p)-2));
                tempij2 = zeros(3*c^2*p,c^2*p,N-p-2);
                
                tempj1 = repmat(1:c^2*p,2*c^2*p,1);
                tempj2 = repmat(c^2*p+1:c^2*p*(N-p-1),3*c^2*p,1);
                tempj3 = repmat(c^2*p*(N-p-1)+1:c^2*p*(N-p),2*c^2*p,1);
                index_j = [tempj1(:);tempj2(:);tempj3(:)];
                clear tempj1 tempj2 tempj3
                
                tic
                for i = 1:N-p
                    tempM = sparse((i-1)*c^2*p+1:i*c^2*p,1:c^2*p,ones(c^2*p,1), c^2*p*(N-p), c^2*p);
                    %tempM(:,(i-1)*c^2*p+1:i*c^2*p) = eye(c^2*p,c^2*p);
                    tempposter_varx = obj.poster_precisex\tempM;
 
                    
                    if i ==1
                        %tempi = repmat((1:2*c^2*p)',1,c^2*p);
                        tempij1 = tempposter_varx(1:2*c^2*p,:);
                        %index_i(1:(c^2*p)^2*2) = tempi(:);
                        %value_ij(1:(c^2*p)^2*2) = tempv(:);
                        
                    end
                    
                    if i == N-p
                        tempij3 = tempposter_varx(c^2*p*(N-p-2)+1:c^2*p*(N-p),:);
                        %index_i(1:(c^2*p)^2*2) = tempi(:);
                        %value_ij(end-(c^2*p)^2*2+1:end) = tempv(:);
                    end
                    
                    if (i < N-p) && (i > 1)
                       tempi2(:,:,i-1) = tempi+c^2*p*(i-2);
                       tempij2(:,:,i-1)= tempposter_varx(c^2*p*(i-2)+1:c^2*p*(i+1),:);
                       
                       
                    end
                    
                    
%                     if i ==1
%                         %poster_varx = sparse(c^2*p*(N-p),c^2*p);
%                         poster_varx = [tempposter_varx(1:2*c^2*p,:);sparse(c^2*p*(N-p-2),c^2*p)];
%                     end
%                     
%                     if i == N-p
%                         tempposter_varx = [sparse(c^2*p*(N-p-2),c^2*p);tempposter_varx(c^2*p*(N-p-2)+1:c^2*p*(N-p),:)];
%                         poster_varx = [poster_varx,tempposter_varx];
%                     end
%                     
%                     if (i < N-p) && (i > 1)
%                         tempposter_varx = [sparse(c^2*p*(i-2),c^2*p);tempposter_varx(c^2*p*(i-2)+1:c^2*p*(i+1),:);sparse(c^2*p*(N-p-i-1),c^2*p)];
%                         poster_varx = [poster_varx,tempposter_varx];
%                     end
                    
                end
                toc
                
                index_i = [tempi1(:);tempi2(:);tempi3(:)];
                clear tempi1 tempi2 tempi3 tempi
                value_ij = [tempij1(:);tempij2(:);tempij3(:)];
                clear tempij1 tempij2 tempij3
                
                poster_varx = sparse(index_i,index_j,value_ij,c^2*p*(N-p),c^2*p*(N-p));
                clear index_i index_j value_ij
                
                
            end
               
            r = obj.generate_r(niter);
            r = 0;
            if r~=0
                [L, prank]= chol(obj.poster_precisex);
                sample = L\randn(c^2*p*(N-p),1) + obj.poster_mux;
            end
            
%            lambda1_SEM =  (c^2*p*(N-p)+2*obj.alpha1-2)/(sample'*sample+2*obj.beta1);
            
%             if r == 1
%                 obj.lambda1_new  =  lambda1_SEM*ones((N-p),1);
%             else
%                 
%                 temp1 = trace(poster_varx) + obj.poster_mux'*obj.poster_mux+2*obj.beta1;
%                 lambda1_EM = (c^2*p*(N-p)+2*obj.alpha1-2)/temp1*ones((N-p),1);
%                 
%                 obj.lambda1_new = (1-r)*lambda1_EM+r*lambda1_SEM;
%             end
            
%%%%%%%%%%%%%%%%%%%% modified at 9.6            
%             tempD = obj.D'*obj.D;
%             for i = 1:N-p-1
%                 
%                 samplexn = sample((i-1)*c^2*p+1:(i+1)*c^2*p);
%                 
%                 lambda2_SEM = (c^2*p+2*obj.alpha2-2)/(samplexn'*tempD*samplexn+2*obj.beta2);
%                 if r == 1
%                     obj.lambda2_new(i) = lambda2_SEM;
%                 else
%                     
%                     xn = obj.poster_mux((i-1)*c^2*p+1:(i+1)*c^2*p);
%                     vxn =  poster_varx((i-1)*c^2*p+1:(i+1)*c^2*p,(i-1)*c^2*p+1:(i+1)*c^2*p);
%                     lambda2_EM = (c^2*p+2*obj.alpha2-2)/(xn'*tempD*xn+trace(vxn*tempD)+2*obj.beta2);
%                     obj.lambda2_new(i) = (1-r)*lambda2_EM+r*lambda2_SEM;
%                 end
%                 
%                 
%             end

            tempD = obj.D'*obj.D;
            tempxDx_SEM = 0;
            tempxDx_EM = 0;
            for i = 1:N-p-1
                if r~=0
                samplexn = sample((i-1)*c^2*p+1:(i+1)*c^2*p);
                
                tempxDx_SEM = tempxDx_SEM+samplexn'*tempD*samplexn;
                end
                    
                xn = obj.poster_mux((i-1)*c^2*p+1:(i+1)*c^2*p);
                vxn =  poster_varx((i-1)*c^2*p+1:(i+1)*c^2*p,(i-1)*c^2*p+1:(i+1)*c^2*p);
                tempxDx_EM = tempxDx_EM+xn'*tempD*xn+trace(vxn*tempD);
             
              
                
                
            end
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
             
             A = obj.poster_mux;
             
             A(abs(A)<obj.threshold_X) = 0;
             
             A = reshape(full(A),[c,cp,np]);
             
             
         end
        
 
    end
    
end

