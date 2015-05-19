classdef sparse_time_varying_MVAR_EM < handle
    %   This version estimate the model using EM algorithm
    %   Detailed explanation goes here
    
    properties
        chan
        len
        trial
        alpha1
        alpha2
        beta1
        beta2
        lambda1
        lambda2
        lambda1_new
        lambda2_new
        poster_precisex
        poster_varx
        poster_mux
        m_order
        filename_full
        EEGdata
        Sigma % covriance matrix, sparse matrix
        X % vecterized MVAR coefficients, sparse matrix 
        
        R
        Q
        D
        varx1
        varx2
        
        threshold_X
        threshold_all
        iteration_threshold
%        epsilon
        
        synthetic_A1
        synthetic_A2
        synthetic_SNR
        
    end
    
    methods
        function obj = sparse_time_varying_MVAR_EM(working_mode,data_to_analyze)
            % working_mode 1: real data
            %              2: synthetic data
            %              3: random generated data, just for test
        
            
            obj.alpha1 = 10^-10;
            obj.alpha2 = 10^-10;
            obj.beta1 = 10^-10;
            obj.beta2 = 10^-10;
           
            obj.threshold_X = 10^-3;
          %  obj.threshold_MM2 = 5*10^-4;
          %  obj.threshold_MM1 = 5*10^-4;
            obj.threshold_all = 5*10^-6;
            obj.iteration_threshold =50;
%            obj.epsilon = 10^-10;
            obj.synthetic_SNR = 15;
            if working_mode == 1
                 obj.m_order = 1;
%                 obj.filename_full = input('Please type in the full-path filename', 's');
%                 if isempty(obj.filename_full)
%                     obj.EEGdata = load(obj.filename_full);
%                     [obj.chan, obj.len, obj.trial] = size(obj.EEGdata);
%                     
%                 end
%                 
%                 obj.synthetic_A1 = [];obj.synthetic_A2 = [];
                obj.EEGdata = data_to_analyze;
                [obj.chan, obj.len, obj.trial] = size(obj.EEGdata);

                obj.synthetic_A1 = [];obj.synthetic_A2 = [];
            end
            
            if working_mode == 2
                obj.m_order = 1;
                obj.chan = input('number of channels:');
                obj.len = input('data length:');
                obj.trial = input('number of trials:');
                obj.filename_full = [];
                
                obj.EEGdata = obj.synthtic_generate(200);
                
                
            end
            
            if working_mode == 3
                obj.m_order = 1; % simulation using 1 order model, it is difficult to find a higher order stable system
                obj.chan = 6;
                obj.len = 500;
                obj.trial = 60;
                
                obj.filename_full = [];
                
                obj.EEGdata = obj.synthtic_generate(200);
            end
%                 obj.poster_precisex = sparse(zeros(obj.chan^2*obj.m_order*(obj.len-obj.m_order)));
%                 obj.poster_varx = sparse(zeros(obj.chan^2*obj.m_order*(obj.len-obj.m_order)));
%                 obj.poster_mux = sparse(zeros(obj.chan^2*obj.m_order*(obj.len-obj.m_order),1));   
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
                AR_coef1 = diag(2*rand(obj.chan,1)-1);
                index_zero = randi(obj.chan,[round(obj.chan*0.25),1]);
                %AR_coef1(:,index_zero) = 0;
                for i = 1:length(index_zero)
                    tempAR_coef = AR_coef1(:,index_zero(i));
                    tempAR_coef = circshift(tempAR_coef, randi(obj.chan-1,1));
                    AR_coef1(:,index_zero(i)) = tempAR_coef;
                end

                %%%% some nonzero coefficients off the diagnoal
                AR_coef1 = AR_coef1+(2*rand(obj.chan)-1)...
                    .*(sprand(obj.chan, obj.chan, 0.3/obj.chan) ~= 0);


                number_nonzero = sum(sum(AR_coef1 ~= 0));
                number_change = ceil(number_nonzero*0.1);
                number_setzero = ceil(number_nonzero*0.05);
                number_setnonzero = ceil(number_nonzero*0.05);

                AR_coefI = (2*rand(obj.chan)-1).*(AR_coef1~=0);
                temp_threshold = sort(AR_coefI(:));
                change_threshold = temp_threshold(number_change);
                setzero_threshold = temp_threshold(end - number_setzero);
                AR_coef2 = (2*rand(obj.chan)-1).*(AR_coefI <= change_threshold) + AR_coef1;
                AR_coef2(AR_coefI > setzero_threshold) = 0;

                AR_coefI = (2*rand(obj.chan)-1).*(AR_coef1==0);
                temp_threshold = sort(AR_coefI(:));
                setnonzero_threshold = temp_threshold(number_setnonzero);
                AR_coef3 = (2*rand(obj.chan)-1).*(AR_coefI <=  setnonzero_threshold);
                AR_coef2 = AR_coef2+AR_coef3;

                A1 = cell(1,1);
               
                A2 = cell(1,1);
                
                [A1{1},A2{1}] = obj.coef_refactor(AR_coef1, AR_coef2);
                AR_coef1 = A1{1};
                AR_coef2 = A2{1};

                epislon = eye(obj.chan);
                Md1 = vgxset('n',obj.chan,'nAR',1,'AR',A1,'Q',epislon);
                Md2 = vgxset('n',obj.chan,'nAR',1,'AR',A2,'Q',epislon);

                [isStable1,~,eigAR1,~] = vgxqual(Md1);
                [isStable2,~,eigAR2,~] = vgxqual(Md2);

                isStable = ~(isStable1&&isStable2);
                disp(['the ',num2str(num_try),' trial to generate synthetic data']);
            end

        end
        
        function [A1, A2] = coef_refactor(obj, AR_coef1, AR_coef2)
            sign_A1 = sign(AR_coef1);
            sign_A2 = sign(AR_coef2);
      
            abs_A1 = abs(AR_coef1);
            abs_A2 = abs(AR_coef2);
            
            max_A = max(max(abs_A1(AR_coef1 ~=0), abs_A2(AR_coef2 ~=0)));
            min_A = min(min(abs_A1(AR_coef1 ~=0), abs_A2(AR_coef2 ~=0)));
            
            abs_A1 = (abs_A1-min_A)/(max_A-min_A)*0.8+0.15;
            abs_A2 = (abs_A2-min_A)/(max_A-min_A)*0.8+0.15;
            
            A1 = sign_A1.*abs_A1;
            A2 = sign_A2.*abs_A2;

            
            
        end
        

        function get_init(obj, method)
            % method 1: randomly generated matrix
            % method 2: using all data to estimate AR
            addpath /home/wch/Documents/tv_AR/arfit
            p = obj.m_order;
            c = obj.chan;
            N = obj.len;
            if method == 1
                obj.X = sparse(randn(c^2*p*(N  - p),1));
            end
            
            if method == 2 
             
                EEG_perm = permute(obj.EEGdata,[2 1 3]);
                
                [~, A, obj.Sigma, ~, ~, ~] = arfit(EEG_perm, ceil(p/2), p, 'sbc', 'zero');
                [A_d1, A_d2] = size(A);
                if A_d1 ~= c
                    error('dimension do not match');
                end
                est_order = floor(A_d2/A_d1);
                if est_order <p
                    initX = zeros(c, c*p);
                    initX(:,1: c*est_order) = A;
                    
                else
                    initX = A(:, 1:c*p);                    
                    
                end
                %%%% initialize with totally non-smooth image
               
                obj.poster_mux = sparse(repmat(initX(:), N - p,1));
                obj.poster_precisex = sparse(eye(c^2*p*(N -p)));
                obj.poster_varx = obj.poster_precisex;
                
                temp1 = trace(obj.varx1*obj.poster_varx) + obj.poster_mux'*obj.varx1*obj.poster_mux+2*obj.beta1;
                obj.lambda1 = (c^2*p*(N-p)+2*obj.alpha1-2)/temp1;
            
                temp2 = trace(obj.varx2*obj.poster_varx) + obj.poster_mux'*obj.varx2*obj.poster_mux+2*obj.beta2;
                obj.lambda2 = (c^2*p*(N-p-1)+2*obj.alpha2-2)/temp2;
                
                obj.matrix_varx1();
                obj.matrix_varx2();
                
                disp(['initial lambda1: ',num2str(obj.lambda1)]);
                disp(['initial lambda2: ',num2str(obj.lambda2)]);

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
                %obj.X = sparse(randn(c^2*p*(N  - p),1
                obj.lambda1_new = ones((N-p),1);
                obj.lambda2_new = ones((N-p-1),1);
                
                obj.matrix_varx1_new();
                obj.matrix_varx2_new();
                
                EEG_perm = permute(obj.EEGdata,[2 1 3]);                
                [~, ~, obj.Sigma, ~, ~, ~] = arfit(EEG_perm, ceil(p/2), p, 'sbc', 'zero');
            end
            
            if method == 2 
             
                EEG_perm = permute(obj.EEGdata,[2 1 3]);
                
                [~, A, obj.Sigma, ~, ~, ~] = arfit(EEG_perm, ceil(p/2), p, 'sbc', 'zero');
                [A_d1, A_d2] = size(A);
                if A_d1 ~= c
                    error('dimension do not match');
                end
                est_order = floor(A_d2/A_d1);
                if est_order <p
                    initX = zeros(c, c*p);
                    initX(:,1: c*est_order) = A;
                    
                else
                    initX = A(:, 1:c*p);                    
                    
                end
                %%%% initialize with totally non-smooth image
               
                obj.poster_mux = sparse(repmat(initX(:), N - p,1));
                obj.poster_precisex = sparse(eye(c^2*p*(N -p)));
                obj.poster_varx = obj.poster_precisex;
                
                temp1 = trace(obj.varx1*obj.poster_varx) + obj.poster_mux'*obj.varx1*obj.poster_mux+2*obj.beta1;
                obj.lambda1 = (c^2*p*(N-p)+2*obj.alpha1-2)/temp1*ones(c^2*p*(N-p),1);
            
                temp2 = trace(obj.varx2*obj.poster_varx) + obj.poster_mux'*obj.varx2*obj.poster_mux+2*obj.beta2;
                obj.lambda2 = (c^2*p*(N-p-1)+2*obj.alpha2-2)/temp2*ones(c^2*p*(N-p-1),1);
                
%                 obj.matrix_varx1_new();
%                 obj.matrix_varx2_new();
                
                disp(['initial lambda1: ',num2str(obj.lambda1(1))]);
                disp(['initial lambda2: ',num2str(obj.lambda2(1))]);

            end
            
        end
        
        
        
        
        
%         function get_initSigma(obj)
%             tempEEG = obj.EEGdata(:,:);
%             obj.Sigma = (tempEEG*tempEEG')/(obj.trial*obj.len).*eye(obj.chan);
%             
%         end 
        
        function matrix_RQ(obj)
            p = obj.m_order;
            c = obj.chan;
            obj.R = sparse(c^2*p*(obj.len-p),c^2*p*(obj.len-p));
            obj.Q = sparse(1, c^2*p*(obj.len-p));
            
            Sigma_sqrt = obj.Sigma^(-0.5);
            for i = p+1:obj.len
                M_n = obj.matrix_Mn(i);

                EEG_in_win = obj.EEGdata(:,i-1:-1:i-p,:);
                %EEGn = obj.EEGdata(:,i,:);
                WW = sparse(c^2*p,c^2*p);
                sW = sparse(1,c^2*p);
                
                for k = 1:obj.trial
                    W_nk = obj.matrix_W(EEG_in_win,k);
                    tilde_W_nk = sparse(Sigma_sqrt)*W_nk;  
                    WW = WW + tilde_W_nk'*tilde_W_nk;
                    % import change on May 11
                    sW = sW + sparse(Sigma_sqrt*squeeze(obj.EEGdata(:,i,k)))'*tilde_W_nk;
                end

                obj.R = obj.R + M_n'*WW*M_n;
                obj.Q = obj.Q + sW*M_n;
            end
           
            %obj.R = obj.R/2;
           
        end
        
        function W_nk = matrix_W(obj,EEG_in_win,k)
            p = obj.m_order;
            W_nk = zeros(obj.chan, obj.chan^2*p);
            for i = 1:p
                W_nk(:,(i-1)*obj.chan^2+1:i*obj.chan^2) = kron(squeeze(EEG_in_win(:,i,k))',eye(obj.chan));
            end
            W_nk = sparse(W_nk);
        end
        
        function matrix_D(obj)
            p = obj.m_order;
            c = obj.chan;
            %D1 = [sparse(-eye(c^2*p*(obj.len - p -1))), sparse(zeros(c^2*p*(obj.len - p -1),c^2*p ))];
            %D2 = [sparse(zeros(c^2*p*(obj.len - p -1),c^2*p )), sparse(eye(c^2*p*(obj.len - p -1)))];
            obj.D = [sparse(-eye(c^2*p)), sparse(eye(c^2*p))];
            
        end
        
        function M_n = matrix_Mn(obj,i)
            p = obj.m_order;
            c = obj.chan;
            
            M_n = sparse(c^2*p,c^2*p*(obj.len-p));
            M_n(:,(i-p-1)*c^2*p+1:(i-p)*c^2*p) = sparse(eye(c^2*p));
        end
        

        function sigma = matrix_Sigma(obj,mode)
            
            p = obj.m_order;
            c = obj.chan;
 %           cp = c*p;
            np = obj.len-p;
           % A = reshape(full(obj.poster_mux),[c,cp,np]);
            
%             [c,cp,np] = size(A);
%             if (c ~= obj.chan) || (cp/c ~= p) || (np ~= obj.len-p)
%                 error('dimension do not match');
%             end
            tempmux = obj.poster_mux;
            tempvarx = obj.poster_varx;
            sigma1 = zeros(c);
            sigma2 = zeros(c);
            for i = p+1:obj.len
   
                EEG_in_win = obj.EEGdata(:,i-1:-1:i-p,:);
                EEG_target = squeeze(obj.EEGdata(:,i,:));

                X_n = sparse(tempmux((i-p-1)*c^2*p+1:(i-p)*c^2*p));
               
                %X_n = A(:,:,i-p);
                %X_n = sparse(X_n(:));
                %EEGn = obj.EEGdata(:,i,:);
%                 WW = sparse(zeros(c^2*p));
%                 sW = sparse(zeros(1,c^2*p));
                epsilon_n = zeros(c, obj.trial);
                poster_varxn = zeros(c);
              
                temp = sparse(tempvarx((i-p-1)*c^2*p+1:(i-p)*c^2*p,(i-p-1)*c^2*p+1:(i-p)*c^2*p));
               
                for k = 1:obj.trial
                    W_nk = obj.matrix_W(EEG_in_win,k);
                    epsilon_n(:,k) = EEG_target(:,k) - W_nk*X_n;
%                     WW = WW + tilde_W_nk'*tilde_W_nk;
%                     sW = sW + sparse(squeeze(obj.EEGdata(:,i,k)))'*tilde_W_nk;
                    poster_varxn = poster_varxn + W_nk*temp*W_nk';
                    
                end
                sigma1 = sigma1 + epsilon_n*epsilon_n';
                sigma2 = sigma2 + poster_varxn;
            end
                sigma = sigma1+sigma2;
            
            if mode == 1 % the bayesian version with conjugate prior
               obj.Sigma = (sigma+eye(c))/(np*obj.trial -1);
            end
            if mode == 2 % the bayesian version with uniform prior
               obj.Sigma = (sigma)/(np*obj.trial);
            end
            
        end
        
        function matrix_varx1(obj) % \sum_n Mn'*Mn
            p = obj.m_order;
            c = obj.chan;
%             temp_varx = sparse(zeros(c^2*p*(obj.len-p)));
%             
%             for i = p+1:obj.len
%                 M_n = obj.matrix_Mn(i);
%                 temp_varx = temp_varx+M_n'*M_n;
%                 
%             end
%             obj.varx1 = temp_varx;
             obj.varx1 = sparse(eye(c^2*p*(obj.len-p)));
            
        end
        
        function matrix_varx1_new(obj) % \sum_n lambda_1n*Mn'*Mn
            p = obj.m_order;
            c = obj.chan;
            temp = repmat(obj.lambda1_new', c^2*p,1);
            obj.varx1 = sparse(diag(temp(:)));
            
        end
        
        function matrix_varx2(obj) % \sum_n Mn,n+1'*D'*D*Mn,n+1
            p = obj.m_order;
            c = obj.chan;
            temp_varx = sparse(c^2*p*(obj.len-p),c^2*p*(obj.len-p));
            
            for i = p+1:obj.len-1
                M_n_n1 = [obj.matrix_Mn(i);obj.matrix_Mn(i+1)];
                temp_varx = temp_varx+M_n_n1'*obj.D'*obj.D*M_n_n1;
                
            end
            obj.varx2 = temp_varx;
        end
        
        function matrix_varx2_new(obj) % \sum_n lambda_2n*Mn,n+1'*D'*D*Mn,n+1
            p = obj.m_order;
            c = obj.chan;
            temp_varx = sparse(c^2*p*(obj.len-p),c^2*p*(obj.len-p));
            tempD = obj.D'*obj.D;
            
            
            for i = p+1:obj.len-1
                % M_n_n1 = [obj.matrix_Mn(i);obj.matrix_Mn(i+1)];
                tempMD = sparse(c^2*p*(obj.len-p),c^2*p*(obj.len-p));
                tempMD((i-p-1)*c^2*p+1:(i-p+1)*c^2*p,(i-p-1)*c^2*p+1:(i-p+1)*c^2*p) = obj.lambda2_new(i-p)*tempD;
                temp_varx = temp_varx+tempMD;
                
            end
            obj.varx2 = temp_varx;
        end
        
        function E_step(obj)
        
            obj.poster_precisex = obj.lambda1*obj.varx1+obj.lambda2*obj.varx2+obj.R;
            obj.poster_varx = inv(obj.poster_precisex);
            obj.poster_mux = obj.poster_precisex\obj.Q';
            
        end
        
        function E_step_new(obj)
        
            obj.poster_precisex = obj.varx1+obj.varx2+obj.R;
            obj.poster_varx = inv(obj.poster_precisex);
            obj.poster_mux = obj.poster_precisex\obj.Q';
            
        end
        
        function M_step(obj)
            p = obj.m_order;
            c = obj.chan;
            N = obj.len;
            
            obj.matrix_varx1();
            obj.matrix_varx2();
            
            temp1 = trace(obj.varx1*obj.poster_varx) + obj.poster_mux'*obj.varx1*obj.poster_mux+2*obj.beta1;
            obj.lambda1 = (c^2*p*(N-p)+2*obj.alpha1-2)/temp1;
            
            temp2 = trace(obj.varx2*obj.poster_varx) + obj.poster_mux'*obj.varx2*obj.poster_mux+2*obj.beta2;
            obj.lambda2 = (c^2*p*(N-p-1)+2*obj.alpha2-2)/temp2;
            
            obj.matrix_Sigma(2);
            
        end
        
        function M_step_new(obj)
            p = obj.m_order;
            c = obj.chan;
            N = obj.len;
            %%%%% to be continue
            for i = 1:N-p
               xn = obj.poster_mux((i-1)*c^2*p+1:i*c^2*p);
               vxn =  obj.poster_varx((i-1)*c^2*p+1:i*c^2*p,(i-1)*c^2*p+1:i*c^2*p);
               obj.lambda1_new(i) = (c^2*p+2*obj.alpha1-2)/(xn'*xn+trace(vxn)+2*obj.beta1);
            end
            
            tempD = obj.D'*obj.D;
            for i = 1:N-p-1
               xn = obj.poster_mux((i-1)*c^2*p+1:(i+1)*c^2*p);
               vxn =  obj.poster_varx((i-1)*c^2*p+1:(i+1)*c^2*p,(i-1)*c^2*p+1:(i+1)*c^2*p);
               obj.lambda2_new(i) = (c^2*p+2*obj.alpha2-2)/(xn'*tempD*xn+trace(vxn*tempD)+2*obj.beta2);
            end
            
            obj.matrix_Sigma(2);
            
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
            
            obj.get_init_new(1);
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
                obj.M_step_new();
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
        
 
    end
    
end

