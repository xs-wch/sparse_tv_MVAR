 classdef sparse_TVMVAR_realdata < handle
    % This class is to apply sparse time varying MVAR model on real EEG
    % data analysis
    % 
    
    properties
        rawdata
        tfiltdata
        EEGdata
        chan
        len
        trial
        srate
        event
                
        channel_used
        chanlocs
        
        sp_pattern
        sp_filter
        
        fit_mvar
        
        
        
    end
    
    methods
        function obj = sparse_TVMVAR_realdata(EEG)
            obj.rawdata = EEG.data;
            obj.srate = EEG.srate;
            obj.event = EEG.event;
            obj.chanlocs = EEG.chanlocs;
            obj.chan = length(obj.chanlocs);
            
        end
        
        function time_filter(obj, f_low,f_up,figure_flag)
            
            %%%% detrend & remove spike
            obj.tfiltdata = detrend(obj.rawdata')';
            
            figure;plot(obj.tfiltdata(1,:));
            num_bp = input('Number of breaking points: '); 
            diffdata = cat(2,zeros(obj.chan,1),diff(obj.tfiltdata,1,2));
            
            for ii = 1:num_bp
                [~,I] = max(abs(diffdata(1,:)));
                obj.tfiltdata(:,I-1:I+1) = repmat(obj.tfiltdata(:,I+2),[1,3]);
                diffdata(:,I-1:I+1) = 0;
            end
            
            obj.tfiltdata = detrend(obj.tfiltdata')';
            %%%% filter
            obj.tfiltdata = eegfilt(obj.tfiltdata, obj.srate, f_low, 0,0,2*fix(obj.srate/f_low));
            if figure_flag
                figure
                plot([0:10^4-1]./10^4.*obj.srate, abs(fft(obj.tfiltdata(1,1:10^4))));
            end
            
            obj.tfiltdata = eegfilt(obj.tfiltdata, obj.srate, 0, f_up);
            
            if figure_flag
                figure
                plot([0:10^4-1]./10^4.*obj.srate, abs(fft(obj.tfiltdata(1,1:10^4))));
            end
        end
        
        function time_section(obj, label, pts_before, pts_after)
            
            num_event = length(obj.event);
            I = false(num_event,1);
            for ie = 1:num_event
                if strcmp(obj.event(ie).type,label)
                    I(ie) = true;
                end
            end
            obj.trial = sum(I);
            obj.len = pts_before+pts_after+1;
            obj.chan = size(obj.tfiltdata,1);
            
            obj.EEGdata = zeros(obj.chan, obj.len, obj.trial);
            event_use = obj.event(I);
            for it = 1:obj.trial
                templatency = event_use(it).latency;
                obj.EEGdata(:,:,it) = obj.tfiltdata(:,[templatency-pts_before,templatency+pts_after]);                               
            end
            
            
        end
        
        function time_section_specific(obj)
            %%%% this function to extract trials on Yijun's data
            len_EC = floor(1*obj.srate);% effector cue
            len_MP = floor(0.7*obj.srate);% movement planning
            len_prepare = floor(0.5*obj.srate);% prepare time
            
            obj.EEGdata.hand_left_EC = [];
            obj.EEGdata.hand_right_EC = [];
            obj.EEGdata.hand_center_EC = [];
            obj.EEGdata.eye_left_EC = [];
            obj.EEGdata.eye_right_EC = [];
            obj.EEGdata.eye_center_EC = [];
            obj.EEGdata.both_left_EC = [];
            obj.EEGdata.both_right_EC = [];
            obj.EEGdata.both_center_EC = [];
            
            obj.EEGdata.hand_left_MP = [];
            obj.EEGdata.hand_right_MP = [];
            obj.EEGdata.hand_center_MP = [];
            obj.EEGdata.eye_left_MP = [];
            obj.EEGdata.eye_right_MP = [];
            obj.EEGdata.eye_center_MP = [];
            obj.EEGdata.both_left_MP = [];
            obj.EEGdata.both_right_MP = [];
            obj.EEGdata.both_center_MP = [];
            
            obj.EEGdata.hand_left_pre = [];
            obj.EEGdata.hand_right_pre = [];
            obj.EEGdata.hand_center_pre = [];
            obj.EEGdata.eye_left_pre = [];
            obj.EEGdata.eye_right_pre = [];
            obj.EEGdata.eye_center_pre = [];
            obj.EEGdata.both_left_pre = [];
            obj.EEGdata.both_right_pre = [];
            obj.EEGdata.both_center_pre = [];

            obj.EEGdata.hand_left_all = [];
            obj.EEGdata.hand_right_all = [];
            obj.EEGdata.hand_center_all = [];
            obj.EEGdata.eye_left_all = [];
            obj.EEGdata.eye_right_all= [];
            obj.EEGdata.eye_center_all = [];
            obj.EEGdata.both_left_all= [];
            obj.EEGdata.both_right_all = [];
            obj.EEGdata.both_center_all = [];
            
            num_event = length(obj.event);
            ie = 1;
            while ie < num_event
                switch obj.event(ie).type
                    case '22'
                        templatency = int32(obj.event(ie).latency);
                        ie = ie+1;
                        switch obj.event(ie).type
                            case '10'
                                obj.EEGdata.eye_left_EC = cat(3,...
                                    obj.EEGdata.eye_left_EC,obj.tfiltdata(:,templatency:templatency+len_EC));
                                obj.EEGdata.eye_left_MP = cat(3,...
                                    obj.EEGdata.eye_left_MP,obj.tfiltdata(:,int32(obj.event(ie).latency):int32(obj.event(ie).latency)+len_MP));
                                obj.EEGdata.eye_left_pre = cat(3,...
                                    obj.EEGdata.eye_left_pre,obj.tfiltdata(:,templatency-len_prepare:templatency));
                                obj.EEGdata.eye_left_all = cat(3,...
                                    obj.EEGdata.eye_left_all,obj.tfiltdata(:,templatency-len_prepare:templatency+len_EC+len_MP));
                                ie = ie+1;
                            case '11'
                                obj.EEGdata.eye_center_EC = cat(3,...
                                    obj.EEGdata.eye_center_EC,obj.tfiltdata(:,templatency:templatency+len_EC));
                                obj.EEGdata.eye_center_MP = cat(3,...
                                    obj.EEGdata.eye_center_MP,obj.tfiltdata(:,int32(obj.event(ie).latency):int32(obj.event(ie).latency)+len_MP));
                                obj.EEGdata.eye_center_pre = cat(3,...
                                    obj.EEGdata.eye_center_pre,obj.tfiltdata(:,templatency-len_prepare:templatency));
                                obj.EEGdata.eye_center_all = cat(3,...
                                    obj.EEGdata.eye_center_all,obj.tfiltdata(:,templatency-len_prepare:templatency+len_EC+len_MP));

                                ie = ie+1;

                            case '12'
                                obj.EEGdata.eye_right_EC = cat(3,...
                                    obj.EEGdata.eye_right_EC,obj.tfiltdata(:,templatency:templatency+len_EC));
                                obj.EEGdata.eye_right_MP = cat(3,...
                                    obj.EEGdata.eye_right_MP,obj.tfiltdata(:,int32(obj.event(ie).latency):int32(obj.event(ie).latency)+len_MP));
                                obj.EEGdata.eye_right_pre = cat(3,...
                                    obj.EEGdata.eye_right_pre,obj.tfiltdata(:,templatency-len_prepare:templatency));
                                obj.EEGdata.eye_right_all = cat(3,...
                                    obj.EEGdata.eye_right_all,obj.tfiltdata(:,templatency-len_prepare:templatency+len_EC+len_MP));

                                ie = ie+1;
                            otherwise
                                ie = ie+1;
                        end
                        
                    case '21'
                        templatency = int32(obj.event(ie).latency);
                        ie = ie+1;
                        switch obj.event(ie).type
                            case '10'
                                obj.EEGdata.hand_left_EC = cat(3,...
                                    obj.EEGdata.hand_left_EC,obj.tfiltdata(:,templatency:templatency+len_EC));
                                obj.EEGdata.hand_left_MP = cat(3,...
                                    obj.EEGdata.hand_left_MP,obj.tfiltdata(:,int32(obj.event(ie).latency):int32(obj.event(ie).latency)+len_MP));
                                obj.EEGdata.hand_left_pre = cat(3,...
                                    obj.EEGdata.hand_left_pre,obj.tfiltdata(:,templatency-len_prepare:templatency));
                                obj.EEGdata.hand_left_all = cat(3,...
                                    obj.EEGdata.hand_left_all,obj.tfiltdata(:,templatency-len_prepare:templatency+len_EC+len_MP));

                                ie = ie+1;
                            case '11'
                                obj.EEGdata.hand_center_EC = cat(3,...
                                    obj.EEGdata.hand_center_EC,obj.tfiltdata(:,templatency:templatency+len_EC));
                                obj.EEGdata.hand_center_MP = cat(3,...
                                    obj.EEGdata.hand_center_MP,obj.tfiltdata(:,int32(obj.event(ie).latency):int32(obj.event(ie).latency)+len_MP));
                                obj.EEGdata.hand_center_pre = cat(3,...
                                    obj.EEGdata.hand_center_pre,obj.tfiltdata(:,templatency-len_prepare:templatency));
                                obj.EEGdata.hand_center_all = cat(3,...
                                    obj.EEGdata.hand_center_all,obj.tfiltdata(:,templatency-len_prepare:templatency+len_EC+len_MP));

                                ie = ie+1;

                            case '12'
                                obj.EEGdata.hand_right_EC = cat(3,...
                                    obj.EEGdata.hand_right_EC,obj.tfiltdata(:,templatency:templatency+len_EC));
                                obj.EEGdata.hand_right_MP = cat(3,...
                                    obj.EEGdata.hand_right_MP,obj.tfiltdata(:,int32(obj.event(ie).latency):int32(obj.event(ie).latency)+len_MP));
                                obj.EEGdata.hand_right_pre = cat(3,...
                                    obj.EEGdata.hand_right_pre,obj.tfiltdata(:,templatency-len_prepare:templatency));
                                obj.EEGdata.hand_right_all = cat(3,...
                                    obj.EEGdata.hand_right_all,obj.tfiltdata(:,templatency-len_prepare:templatency+len_EC+len_MP));

                                ie = ie+1;
                            otherwise
                                ie = ie+1;
                        end


                    case '20'
                        templatency = int32(obj.event(ie).latency);
                        ie = ie+1;
                        switch obj.event(ie).type
                            case '10'
                                obj.EEGdata.both_left_EC = cat(3,...
                                    obj.EEGdata.both_left_EC,obj.tfiltdata(:,templatency:templatency+len_EC));
                                obj.EEGdata.both_left_MP = cat(3,...
                                    obj.EEGdata.both_left_MP,obj.tfiltdata(:,int32(obj.event(ie).latency):int32(obj.event(ie).latency)+len_MP));
                                obj.EEGdata.both_left_pre = cat(3,...
                                    obj.EEGdata.both_left_pre,obj.tfiltdata(:,templatency-len_prepare:templatency));
                                obj.EEGdata.both_left_all = cat(3,...
                                    obj.EEGdata.both_left_all,obj.tfiltdata(:,templatency-len_prepare:templatency+len_EC+len_MP));

                                ie = ie+1;
                            case '11'
                                obj.EEGdata.both_center_EC = cat(3,...
                                    obj.EEGdata.both_center_EC,obj.tfiltdata(:,templatency:templatency+len_EC));
                                obj.EEGdata.both_center_MP = cat(3,...
                                    obj.EEGdata.both_center_MP,obj.tfiltdata(:,int32(obj.event(ie).latency):int32(obj.event(ie).latency)+len_MP));
                                obj.EEGdata.both_center_pre = cat(3,...
                                    obj.EEGdata.both_center_pre,obj.tfiltdata(:,templatency-len_prepare:templatency));
                                obj.EEGdata.both_center_all = cat(3,...
                                    obj.EEGdata.both_center_all,obj.tfiltdata(:,templatency-len_prepare:templatency+len_EC+len_MP));

                                ie = ie+1;

                            case '12'
                                obj.EEGdata.both_right_EC = cat(3,...
                                    obj.EEGdata.both_right_EC,obj.tfiltdata(:,templatency:templatency+len_EC));
                                obj.EEGdata.both_right_MP = cat(3,...
                                    obj.EEGdata.both_right_MP,obj.tfiltdata(:,int32(obj.event(ie).latency):int32(obj.event(ie).latency)+len_MP));
                                obj.EEGdata.both_right_pre = cat(3,...
                                    obj.EEGdata.both_right_pre,obj.tfiltdata(:,templatency-len_prepare:templatency));
                                obj.EEGdata.both_right_all = cat(3,...
                                    obj.EEGdata.both_right_all,obj.tfiltdata(:,templatency-len_prepare:templatency+len_EC+len_MP));

                                ie = ie+1;
                            otherwise
                                ie = ie+1;
                        end


                    otherwise
                        ie = ie+1;
                        
                end
            end
            
            

            
            
        end
        
        function data = EEGdownsample(obj,data_to_downsample,n)
            
            [tempchan,~,temptrial] = size(data_to_downsample);
            templen = size(downsample(squeeze(data_to_downsample(:,:,1))',n)',2);
            tempEEGdata = zeros(tempchan,templen,temptrial);
            for it = 1:temptrial
                temp = downsample(squeeze(data_to_downsample(:,:,it))',n)';
                %%%% demean, keep zeros mean
                tempEEGdata(:,:,it) = detrend(temp', 'constant')';
            end
            
            data = tempEEGdata;
           % obj.srate = obj.srate/n;
            
        end
        
        function EEGICA(obj,percent_kept)
            tempdata = obj.EEGdata(:,:);
            
            e = eig(tempdata*tempdata');
            
            I_kept = cumsum(e) >= ((1-percent_kept)*sum(e)) ;
            obj.chan = sum(I_kept)+1;
            [weight, sphere] = runica(tempdata,'pca',obj.chan);
            obj.sp_filter = weight*sphere;
            obj.sp_pattern = pinv(obj.sp_filter);
            obj.EEGdata = obj.sp_filter*tempdata;
            obj.EEGdata = reshape(obj.EEGdata,[obj.chan,obj.chan,obj.trial]);
                       
        end
        
        function data = channel_select(obj,data)
            %%%% channel_used is a cell,each cell contains the name of
            %%%% each channel
            if ~isempty(obj.channel_used)
                
                I_chan = false(obj.chan,1);
                for ic_use = 1:length(obj.channel_used)
                    tempchan = obj.channel_used{ic_use};
                    for ic = 1:obj.chan
                        if strcmp(obj.chanlocs(ic).labels, tempchan)
                            I_chan(ic) = true;
                        end
                        
                    end
                    
                end
                obj.chanlocs = obj.chanlocs(I_chan);
                obj.chan = sum(I_chan);
                data = data(I_chan,:,:);
            end
            
            
        end
        
        
        function estimate_MVAR(obj,data,preprocess_method,order)
            
            [obj.chan,obj.len,obj.trial] = size(data);
            ERP = mean(data,3);
            ERPtensor  = reshape(repmat(ERP,1,obj.trial),[obj.chan,obj.len,obj.trial]);
            data = data-ERPtensor;
            
            %%%% PCA to make the two components uncorrelated
            if preprocess_method == 1
                tempdata = data(:,:);
                [V,D] = eig(tempdata*tempdata');
                tempdata = (D^-0.5*V')* tempdata;
                data = reshape(tempdata,[obj.chan,obj.len,obj.trial]);
                stv_MVAR = sparse_time_varying_MVAR_EMpara(1,[],[],data,order);
                %  stv_MVAR.m_order = 3;
                obj.fit_mvar = stv_MVAR.estimate_model_SAEM();
                [c,cp,l] = size(obj.fit_mvar);
                temp_mvar = obj.fit_mvar(:,:);
                obj.fit_mvar = reshape((V*D^0.5)*temp_mvar,[c,cp,l]);
            end
            %%%% 
            if preprocess_method == 2
                tempdata = data(:,:);
                amp_factor = std(tempdata,0,2);
                tempdata = diag(amp_factor.^-1)* tempdata;
                data = reshape(tempdata,[obj.chan,obj.len,obj.trial]);
                stv_MVAR = sparse_time_varying_MVAR_EMpara(1,[],[],data,order);
                %  stv_MVAR.m_order = 3;
                obj.fit_mvar = stv_MVAR.estimate_model_SAEM();
               % [c,cp,l] = size(obj.fit_mvar);
                %temp_mvar = obj.fit_mvar(:,:);
                %obj.fit_mvar = reshape((V*D^0.5)*temp_mvar,[c,cp,l]);
            end
            
            if preprocess_method == 3

                stv_MVAR = sparse_time_varying_MVAR_EMpara(1,[],[],data,order);
                %  stv_MVAR.m_order = 3;
                obj.fit_mvar = stv_MVAR.estimate_model_SAEM();
                % [c,cp,l] = size(obj.fit_mvar);
                %temp_mvar = obj.fit_mvar(:,:);
                %obj.fit_mvar = reshape((V*D^0.5)*temp_mvar,[c,cp,l]);
            end
            
        end
        
        function [PDC, DTF] = PDC_DTF(obj)
            
            m_order = size(obj.fit_mvar,2)/obj.chan; 
            p = obj.len-m_order;
  
            Fmax = 40;
            Nf = 40;
        
            [PDC, DTF] = obj.PDC_DTF_matrix(obj.fit_mvar, m_order , obj.srate , Fmax , Nf );
      
        end
            
        
        
        
        
    end
    
end

