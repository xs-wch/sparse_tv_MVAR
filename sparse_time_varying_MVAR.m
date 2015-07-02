classdef sparse_time_varying_MVAR < handle
    % this version for the initialization of sparse_time_varying_MVAR_EM
    % 2014-07-01, this version ignore sparse prior, use TV regularization
    
    properties
        chan
        len
        trial
        %alpha1
        alpha2
        %beta1
        beta2
        %theta1
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
        W
        
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
        function obj = sparse_time_varying_MVAR(working_mode,data_to_analyze,order)
            % working_mode 1: real data
            %              2: synthetic data
            %              3: random generated data, just for test
            %obj.alpha1 = 10^-3;
            obj.alpha2 = 10^-3;
            %obj.beta1 = 10^-3;
            obj.beta2 = 10^-3;
            obj.threshold_X = 10^-3;
            obj.threshold_MM2 = 5*10^-3;
            obj.threshold_MM1 = 5*10^-2;
            obj.threshold_all = 5*10^-2;
            obj.iteration_threshold =5;
            obj.epsilon = 10^-10;
            obj.synthetic_SNR = 15;
            if working_mode == 1

                obj.EEGdata = data_to_analyze;
                [obj.chan, obj.len, obj.trial] = size(obj.EEGdata);
                obj.m_order = order;
                obj.synthetic_A1 = [];
                obj.synthetic_A2 = [];
            end

            
        end
        


        function get_initX(obj, method)
            % method 1: randomly generated matrix
            % method 2: using all data to estimate AR
            addpath /home/wch/Documents/tv_AR/arfit
            p = obj.m_order;
            c = obj.chan;
            if method == 1
                obj.X = sparse(randn(c^2*p*(obj.len - p),1));
                EEG_perm = permute(obj.EEGdata,[2 1 3]);                
                [~, ~, obj.Sigma, ~, ~, ~] = arfit(EEG_perm, ceil(p/2), p, 'sbc', 'zero');
            end
%             
%             if method == 2 
% 
%                 EEG_perm = permute(obj.EEGdata,[2 1 3]);
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
%                 noise_factor = 0.1*std(initX(:));
%                 obj.X = sparse(repmat(initX(:), obj.len - p,1) + noise_factor*randn(c^2*p*(obj.len - p),1));
% 
%             end
            
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
            
            Sigma_inv = sparse(inv(obj.Sigma));
            for i = p+1:obj.len
                %EEG_in_win = obj.EEGdata(:,i-1:-1:i-p,:);
                WW = sparse(c^2*p,c^2*p);
                sW = sparse(1,c^2*p);
                I = false(c^2*p*(obj.len-p),1);
                I((i-p-1)*c^2*p+1:(i-p)*c^2*p) = true;
                for k = 1:obj.trial
                    W_nk = obj.W{i-p,k};
                    WW = WW + W_nk'*Sigma_inv*W_nk;
                    sW = sW + sparse(squeeze(obj.EEGdata(:,i,k)))'*Sigma_inv*W_nk;
                end

                obj.R(I,I) = obj.R(I,I) + WW;
                obj.Q(I) = obj.Q(I) + sW;
            end
                obj.R = obj.R/2;
           
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
                        tempW_nk(:,(i-1)*c^2+1:i*c^2) = kron(squeeze(EEG_in_win(:,i,k))',eye(c));
                    end
                    obj.W{iwin-p,k} = sparse(tempW_nk);
                end
            end
       
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
       
            sigma = zeros(c);
            for i = p+1:obj.len

                EEG_target = squeeze(obj.EEGdata(:,i,:));
                X_n = A(:,:,i-p);
                X_n = sparse(X_n(:));

                epsilon_n = zeros(c, obj.trial);
                for k = 1:obj.trial
                    W_nk =  obj.W{i-p,k};
                    epsilon_n(:,k) = EEG_target(:,k) - W_nk*X_n;

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
        
        
        function X_update = get_X_MM2(obj, rho2, X_old)
            s = warning('error', 'MATLAB:nearlySingularMatrix');
      %      X_old = obj.X;
            delta_t = 1;
            target_old = 10^-10; 
            iteration_time = 0;
            while (delta_t > obj.threshold_MM2) &&  (iteration_time < obj.iteration_threshold)
                iteration_time = iteration_time+1;
               % Lambda1 = sparse(diag((2*abs(X_old)+obj.epsilon).^-1));
                Lambda2 = sparse(diag((2*abs(obj.D*X_old)+obj.epsilon).^-1));
                obj.H = obj.R+rho2*obj.D'*Lambda2*obj.D;
                % when dignal element goes to infinity, how to simplify the
                % calculation
                % some elements can be set to zeros directly
                % both the dimension of H and Q can be reduced before
                % calculation

                try
                    X_update = obj.H\(obj.Q./2)';
                catch 
                    X_update = X_old;
                    disp('matrix close to sigular');
                    return
                    
                end

                X_update = sparse(X_update);
                
                target_update = X_update'*obj.R*X_update-obj.Q*X_update+rho2*sum(abs(obj.D*X_update));
                
                delta_t_sign = (target_update - target_old)/abs(target_old);
                delta_t = -(delta_t_sign);
                abs(norm(X_update-X_old)/norm(X_old));
                
                target_old = target_update;
                X_old = X_update; 
                
                disp(['    MM2: iteration ',  num2str(iteration_time), ' delta t: ', num2str(delta_t_sign)]);

            end
            
           
            warning(s);
            if iteration_time >= obj.iteration_threshold
                disp('iteration time over threshold, the result may be not accrute');
            end
        
        end
        
        function X_update = get_X_MM1(obj, mu2, X_old)

            delta_t = 1;
            target_old = 10^-10; 
            iteration_time = 0;
            while (delta_t > obj.threshold_MM1) &&  (iteration_time < obj.iteration_threshold)
                iteration_time = iteration_time+1;
               
         %       rho1 = 0;
          
                rho2 = mu2/(sum(abs(obj.D*X_old))+obj.beta2);

                X_update = obj.get_X_MM2(rho2, X_old);

                %target_update = X_update'*obj.R*X_update - obj.Q*X_update + rho1*sum(abs(X_update)) + rho2*sum(abs(obj.D*X_update));

                target_update = X_update'*obj.R*X_update - obj.Q*X_update + mu2*log(sum(abs(obj.D*X_update))+obj.beta2);
    
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
%            e1 = log(sum(abs(X))+obj.beta1);
            e2 = log(sum(abs(obj.D*X))+obj.beta2);
% 
%             %%%%M_step: digamma(theta1) = e
%             psi1_lowerbound = psi(obj.alpha1);
%             psi1_upperbound = psi(c^2*p*(N-p)+obj.alpha1);
%             if e1 > psi1_upperbound
%                 obj.theta1 = 0;
%                 disp('e1>psi1')
%             end
%             if e1 < psi1_lowerbound 
%                 obj.theta1 = 1;
%                 disp('e1<psi1')
%             end
%             if (e1 < psi1_upperbound) && (e1 > psi1_lowerbound )
%                obj.theta1 = (obj.invpsi(e1)-obj.alpha1)/(c^2*p*(N-p));               
%                 
%             end
            
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
               obj.theta2 = (obj.invpsi(e2)-obj.alpha2)/(c^2*p*(N-p-1));               
                
            end

            
        end
        
        function Y=invpsi(obj,X)
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
        

        
        function X_update = estimate_model(obj)
            
%            obj.get_initSigma();
            obj.get_initX(1);
            obj.matrix_W();
            obj.matrix_RQ();            
            obj.matrix_D();
%            obj.matrix_RnQn();
%            obj.theta1 = 0.1;
            obj.theta2 = 0.1;
            disp(['parameters: \n','alpha1: ',num2str(obj.alpha2)])
%            disp(['alpha2: ',num2str(obj.alpha2)])
%            disp(['beta1: ',num2str(obj.beta1)])
            disp(['beta2: ',num2str(obj.beta2)])
            disp(['model order: ', num2str(obj.m_order)])
            disp(['EEGdata size: ',num2str(obj.chan),'*',num2str(obj.len),'*',num2str(obj.trial)])
            %data = obj.EEGdata;
            c = obj.chan;
            p = obj.m_order;
            N= obj.len;
            
            delta_theta = 1;
            niterouter = 0;
            while (delta_theta > 10^-4) && (niterouter < 3)
                niterouter = niterouter + 1;

                mu1 = 0;

                mu2 = obj.theta2*c^2*p*(N-p-1)+obj.alpha2;

                disp(['mu1: ',num2str(mu1),' mu2: ', num2str(mu2)]);
                i = 0;
                if niterouter == 1
                    X_old = obj.X;
                else
                    X_old = X_update;
                end
                L_old = 1;
                deltaL = 1;
                while (deltaL > obj.threshold_all) && (i < 5)
                    i = i+1;

                    X_update = obj.get_X_MM1(mu2, X_old);
                    %%%%%%%% some problem to be sloved 
                    %%%%%%%% update the algorithm, sigma is updated in each
                    %%%%%%%% step
                    matrix_Sigma(obj,X_update,2);

                    obj.matrix_RQ();

                    L_update = (obj.len-p)*obj.trial/2*log(det(obj.Sigma)) + mu2*log(sum(abs(obj.D*X_update))+obj.beta2);


                    X_old = X_update;

                    dL = (L_update - L_old)/abs( L_old);
                    L_old = L_update;
                    disp(['Iteration ', num2str(i), 'delta t:', num2str(dL)]);

                    deltaL = abs(dL);

                end
           %     temptheta1 = obj.theta1;
                temptheta2 = obj.theta2;
                obj.update_theta(X_update);
                
                delta_theta = abs((obj.theta2-temptheta2)/temptheta2);
                
                disp(['outer iteration: ',num2str(niterouter),' delta theta: ',num2str(delta_theta)]);
            end
            
           
        end
 
    end
    
end

