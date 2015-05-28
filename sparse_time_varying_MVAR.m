classdef sparse_time_varying_MVAR < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        chan
        len
        trial
        alpha1
        alpha2
        beta1
        beta2
        theta1
        theta2
        m_order
        filename_full
        EEGdata
        Sigma % covriance matrix, sparse matrix
        X % vecterized MVAR coefficients, sparse matrix 
        
        R
        Q
        D
        H
        
        threshold_X
        threshold_MM2
        threshold_MM1
        threshold_all
        iteration_threshold
        epsilon
        
        synthetic_A1
        synthetic_A2
        synthetic_SNR
        
    end
    
    methods
        function obj = sparse_time_varying_MVAR(working_mode,data_to_analyze)
            % working_mode 1: real data
            %              2: synthetic data
            %              3: random generated data, just for test
            obj.alpha1 = 10^-4;
            obj.alpha2 = 10^-4;
            obj.beta1 = 1;
            obj.beta2 = 1;
            obj.m_order = 3;
            obj.threshold_X = 10^-3;
            obj.threshold_MM2 = 5*10^-4;
            obj.threshold_MM1 = 5*10^-4;
            obj.threshold_all = 5*10^-4;
            obj.iteration_threshold =30;
            obj.epsilon = 10^-10;
            obj.synthetic_SNR = 15;
            if working_mode == 1
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
           2 
            abs_A1 = abs(AR_coef1);
            abs_A2 = abs(AR_coef2);
            
            max_A = max(max(abs_A1(AR_coef1 ~=0), abs_A2(AR_coef2 ~=0)));
            min_A = min(min(abs_A1(AR_coef1 ~=0), abs_A2(AR_coef2 ~=0)));
            
            abs_A1 = (abs_A1-min_A)/(max_A-min_A)*0.8+0.15;
            abs_A2 = (abs_A2-min_A)/(max_A-min_A)*0.8+0.15;
            
            A1 = sign_A1.*abs_A1;
            A2 = sign_A2.*abs_A2;

            
            
        end
        

        function get_initX(obj, method)
            % method 1: randomly generated matrix
            % method 2: using all data to estimate AR
            addpath /home/wch/Documents/tv_AR/arfit
            p = obj.m_order;
            c = obj.chan;
            if method == 1
                obj.X = sparse(randn(c^2*p*(obj.len - p),1));
            end
            
            if method == 2 
             
                EEG_perm = permute(obj.EEGdata,[2 1 3]);
                
                [~, A, ~, ~, ~, ~] = arfit(EEG_perm, ceil(p/2), p, 'sbc', 'zero');
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
                noise_factor = 0.1*std(initX(:));
                obj.X = sparse(repmat(initX(:), obj.len - p,1) + noise_factor*randn(c^2*p*(obj.len - p),1));

            end
            
        end
        
        
        function get_initSigma(obj)
            tempEEG = obj.EEGdata(:,:);
            obj.Sigma = (tempEEG*tempEEG')/(obj.trial*obj.len).*eye(obj.chan);
            
        end 
        
        function matrix_RQ(obj)
            p = obj.m_order;
            c = obj.chan;
            obj.R = sparse(zeros(c^2*p*(obj.len-p)));
            obj.Q = sparse(zeros(1, c^2*p*(obj.len-p)));
            
            Sigma_sqrt = obj.Sigma^(-0.5);
            for i = p+1:obj.len
                M_n = sparse(zeros(c^2*p,c^2*p*(obj.len-p)));
                M_n(:,(i-p-1)*c^2*p+1:(i-p)*c^2*p) = sparse(eye(c^2*p));
                EEG_in_win = obj.EEGdata(:,i-1:-1:i-p,:);
                %EEGn = obj.EEGdata(:,i,:);
                WW = sparse(zeros(c^2*p));
                sW = sparse(zeros(1,c^2*p));
                
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
           
            obj.R = obj.R/2;
           
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
            D1 = [sparse(-eye(c^2*p*(obj.len - p -1))), sparse(zeros(c^2*p*(obj.len - p -1),c^2*p ))];
            D2 = [sparse(zeros(c^2*p*(obj.len - p -1),c^2*p )), sparse(eye(c^2*p*(obj.len - p -1)))];
            obj.D = D1+D2;
            
        end
        

        function sigma = matrix_Sigma(obj,A,mode)
            
            p = obj.m_order;
            c = obj.chan;
            cp = c*p;
            np = obj.len-p;
            A = reshape(full(A),[c,cp,np]);
            
%             [c,cp,np] = size(A);
%             if (c ~= obj.chan) || (cp/c ~= p) || (np ~= obj.len-p)
%                 error('dimension do not match');
%             end
%             
            sigma = zeros(c);
            for i = p+1:obj.len

                EEG_in_win = obj.EEGdata(:,i-1:-1:i-p,:);
                EEG_target = squeeze(obj.EEGdata(:,i,:));
                X_n = A(:,:,i-p);
                X_n = sparse(X_n(:));
                %EEGn = obj.EEGdata(:,i,:);
%                 WW = sparse(zeros(c^2*p));
%                 sW = sparse(zeros(1,c^2*p));
                epsilon_n = zeros(c, obj.trial);
                for k = 1:obj.trial
                    W_nk = obj.matrix_W(EEG_in_win,k);
                    epsilon_n(:,k) = EEG_target(:,k) - W_nk*X_n;
%                     WW = WW + tilde_W_nk'*tilde_W_nk;
%                     sW = sW + sparse(squeeze(obj.EEGdata(:,i,k)))'*tilde_W_nk;
                    
                    
                end
                sigma = sigma + epsilon_n*epsilon_n';

            end
            
            if mode == 1 % the bayesian version with conjugate prior
               obj.Sigma = (sigma+eye(c))/(np*obj.trial -1);
            end
            if mode == 2 % the bayesian version with uniform prior
               obj.Sigma = sigma/(np*obj.trial);
            end
            
        end
        
        
        function X_update = get_X_MM2(obj, rho1,  rho2, X_old)
            s = warning('error', 'MATLAB:nearlySingularMatrix');
      %      X_old = obj.X;
            delta_t = 1;
            target_old = 10^-10; 
            iteration_time = 0;
            while (delta_t > obj.threshold_MM2) &&  (iteration_time < obj.iteration_threshold)
                iteration_time = iteration_time+1;
                Lambda1 = sparse(diag((2*abs(X_old)+obj.epsilon).^-1));
                Lambda2 = sparse(diag((2*abs(obj.D*X_old)+obj.epsilon).^-1));
                obj.H = obj.R+rho1*Lambda1+rho2*obj.D'*Lambda2*obj.D;
                % when dignal element goes to infinity, how to simplify the
                % calculation
                % some elements can be set to zeros directly
                % both the dimension of H and Q can be reduced before
                % calculation
                
               
%                 if condest(H) > 10^18
%                     I_nonzeros = (X_old ~= 0);
%                     H_reduce = H(:,I_nonzeros);
%                     
%                     X_reduce = H_reduce\(obj.Q./2)';
%                     X_update = sparse(zeros(size(X_old)));
%                     X_update(I_nonzeros) = X_reduce;
%                 else
%                     X_update = H\(obj.Q./2)';
%                 end
                    
               % X_update(~I_nonzeros) = 0;
               
                %condest(H)
                try
                    X_update = obj.H\(obj.Q./2)';
                catch 
                    X_update = X_old;
                    disp('matrix close to sigular');
                    return
                    
                end
                
                %X_update = cgs(H, (obj.Q/2)', 10^-6, 15);
                %
                
                X_update = sparse(X_update);
                
                target_update = X_update'*obj.R*X_update-obj.Q*X_update+rho1*sum(abs(X_update))+rho2*sum(abs(obj.D*X_update));
                
                delta_t_sign = (target_update - target_old)/abs(target_old);
                delta_t = -(delta_t_sign);
                abs(norm(X_update-X_old)/norm(X_old));
                %(target_update - target_old)/abs(target_old)
                
                target_old = target_update;
                X_old = X_update; 
                
                disp(['    MM2: iteration ',  num2str(iteration_time), ' delta t: ', num2str(delta_t_sign)]);

            end
            
            %X_update(abs(X_update) < obj.threshold_X) = 0;
            
            warning(s);
            if iteration_time >= obj.iteration_threshold
                disp('iteration time over threshold, the result may be not accrute');
            end
        
        end
        
        function X_update = get_X_MM1(obj, mu1, mu2, X_old)
         %   X_old = obj.X;
            delta_t = 1;
            target_old = 10^-10; 
            iteration_time = 0;
            while (delta_t > obj.threshold_MM1) &&  (iteration_time < obj.iteration_threshold)
                iteration_time = iteration_time+1;
                rho1 = mu1/(sum(abs(X_old))+obj.beta1);
                rho2 = mu2/(sum(abs(obj.D*X_old))+obj.beta2);

                X_update = obj.get_X_MM2(rho1,  rho2, X_old);

                %target_update = X_update'*obj.R*X_update - obj.Q*X_update + rho1*sum(abs(X_update)) + rho2*sum(abs(obj.D*X_update));

                target_update = X_update'*obj.R*X_update - obj.Q*X_update + mu1*log(sum(abs(X_update))+obj.beta1) + mu2*log(sum(abs(obj.D*X_update))+obj.beta2);
    
                delta_t_sign = (target_update - target_old)/abs(target_old);
                
                delta_t = abs(delta_t_sign);
                
                target_old = target_update;
                X_old = X_update;
                
                 disp(['  MM1: iteration ',  num2str(iteration_time), ' delta t: ', num2str(delta_t_sign)]);
            
            end
            
            if iteration_time >= obj.iteration_threshold
                disp('iteration time over threshold, the result may be not accrute');
            end
            
            
        end
        
        
        function update_theta(obj, X)
            %%%% update theta1 and theta2 using stochastic EM
            c = obj.chan;
            p = obj.m_order;
            N= obj.len;
            %%%%E_step: M-H algorithm
            %e1 = log(sum(abs(X))+obj.beta1);
            %e2 = log(sum(abs(obj.D*X))+obj.beta2);
            MCsize = 1000;
            [e1,e2] = M_H_algorithm(obj, X, MCsize);
            
            %%%%M_step: digamma(theta1) = e
            psi1_lowerbound = psi(obj.alpha1);
            psi1_upperbound = psi(c^2*p*(N-p)+obj.alpha1);
            if e1 > psi1_upperbound
                obj.theta1 = 0;
                disp('e1>psi1')
            end
            if e1 < psi1_lowerbound 
                obj.theta1 = 1;
                disp('e1<psi1')
            end
            if (e1 < psi1_upperbound) && (e1 > psi1_lowerbound )
               obj.theta1 = (invpsi(e1)-obj.alpha1)/(c^2*p*(N-p));               
                
            end
            
            psi2_lowerbound = psi(obj.alpha2);
            psi2_upperbound = psi(c^2*p*(N-p-1)+obj.alpha2);
            
            if e2 > psi2_upperbound
                obj.theta2 = 0;
                disp('e2>psi2')
            end
            if e2 < psi2_lowerbound 
                obj.theta2 = 1;
                disp('e2<psi2')
            end
            if (e2 < psi2_upperbound) && (e2 > psi2_lowerbound )
               obj.theta2 = (invpsi(e2)-obj.alpha2)/(c^2*p*(N-p-1));               
                
            end

            
        end
        
        function [e1,e2] = M_H_algorithm(obj, X, MCsize)
            
            if isempty(MCsize)
                MCsize = 1000;
            end
            c = obj.chan;
            p = obj.m_order;
            N= obj.len;
            %%%% proposal distributuion : laplace approximate
            sign1 = sign(X);
            DX = obj.D*X;
            sign2 = obj.D'*sign(DX);
            
            mu1 = obj.theta1*c^2*p*(N-p)+obj.alpha1;
            mu2 = obj.theta2*c^2*p*(N-p-1)+obj.alpha2;
            
            A1 = sparse((mu1/((sum(abs(X))+obj.beta1)^2))*(sign1*sign1'));
            A2 = sparse((mu2/((sum(abs(DX))+obj.beta2)^2))*(sign2*sign2'));
            
            A = 2*obj.R-A1-A2;
            
            %%%% Monte Carlo
            [L,err] = cholcov(A);
            if err ~= 0
                error(message('stats:mvnrnd:BadCovariance2DSymPos'));
            end
            
            randvec = randn(c^2*p*(N-p),MCsize);
            MCsample = L\randvec+repmat(X,1,MCsize);
            
            %%%% M-H
            p_tilde = @(x) exp(-x'*obj.R*x+obj.Q*x-mu1*log(sum(abs(x))+obj.beta1)-mu2*log(sum(abs(obj.D*x))+obj.beta2));
            naccept = 0;
            for i = 1:MCsize
                ak = min(1,p_tilde(MCsample(:,i))/p_tilde(X));
                if rand(1) > ak
                    MCsample(:,i) = X;
                else
                    naccept = naccept + 1;
                end
             
            end
            
            e1 = mean(log(sum(abs(MCsample),1)+obj.beta1));
            e2 = mean(log(sum(abs(obj.D*MCsample),1)+obj.beta2));
            
            
            
            %q_tilde = @(z)
            
            
        end
        
        function Y=invpsi(X)
        % Y = INVPSI(X)
        %
        % Inverse digamma (psi) function.  The digamma function is the
        % derivative of the log gamma function.  This calculates the value
        % Y > 0 for a value X such that digamma(Y) = X.
        %
        % This algorithm is from Paul Fackler:
        % http://www4.ncsu.edu/~pfackler/
        %
          L = 1;
          Y = exp(X);
          while L > 10e-8
            Y = Y + L*sign(X-psi(Y));
            L = L / 2;
          end
        end
        
        function A = estimate_model(obj,sp_mode)
            
            obj.get_initSigma();
            obj.get_initX(1);
            obj.matrix_RQ();
            obj.matrix_D();
            obj.theta1 = 0.3;
            obj.theta2 = 0.3;
            disp(['parameters: \n','alpha1: ',num2str(obj.alpha1)])
            disp(['alpha2: ',num2str(obj.alpha2)])
            disp(['beta1: ',num2str(obj.beta1)])
            disp(['beta2: ',num2str(obj.beta2)])
            disp(['model order: ', num2str(obj.m_order)])
            disp(['EEGdata size: ',num2str(obj.chan),'*',num2str(obj.len),'*',num2str(obj.trial)])
            %data = obj.EEGdata;
            c = obj.chan;
            p = obj.m_order;
            N= obj.len;
            for it = 1:5 
            
                if sp_mode == 1 % with sparse
                    mu1 = obj.theta1*c^2*p*(N-p)+obj.alpha1;

                end
                if sp_mode == 2 % without sparse
                    mu1 = 0;
                end

                mu2 = obj.theta2*c^2*p*(N-p-1)+obj.alpha2;

                disp(['mu1: ',num2str(mu1),' mu2', num2str(mu2)]);
                i = 0;
                if it == 1
                    X_old = obj.X;
                else
                    X_old = X_update;
                end
                L_old = 1;
                deltaL = 1;
                while (deltaL > obj.threshold_all) && (i<20)
                    i = i+1;

                    X_update = obj.get_X_MM1(mu1, mu2, X_old);
                    %%%%%%%% some problem to be sloved 
                    %%%%%%%% update the algorithm, sigma is updated in each
                    %%%%%%%% step
                    matrix_Sigma(obj,X_update,2);

                    obj.matrix_RQ();

                    L_update = (obj.len-p)*obj.trial/2*log(det(obj.Sigma)) + mu1*log(sum(abs(X_update))+obj.beta1) + mu2*log(sum(abs(obj.D*X_update))+obj.beta2);


                    X_old = X_update;

                    dL = (L_update - L_old)/abs( L_old);
                    L_old = L_update;
                    disp(['Iteration ', num2str(i), 'delta t:', num2str(dL)]);

                    deltaL = abs(dL);

                end

                update_theta(obj, X_update);
            end
            
            X_update(abs(X_update) < obj.threshold_X) = 0;
            A = reshape(full(X_update),[c,c*p,N-p]);
            
            
            
        end
 
    end
    
end

