function showAlignment(src,event)
aligned_data = src.Parent.Parent.UserData.aligned_data;
data_sources = src.Parent.Parent.UserData.alignment_source;

%% 5 Hz stim pulse
if isfield(aligned_data,'left_LFP_table') || isfield(aligned_data,'right_LFP_table')
    fig_hand = findobj('Name','Alignment Results - 5 Hz stim pulse');
    if ~isempty(fig_hand)
        close(fig_hand,'reset');
        figure('Name','Alignment Results - 5 Hz stim pulse');
    else
        figure('Name','Alignment Results - 5 Hz stim pulse');
    end
    if isfield(aligned_data,'left_LFP_table')
        ax(1) = subplot(4,1,1);
        if ~isempty(data_sources{1})
            plot(ax(1),aligned_data.left_taxis,aligned_data.left_LFP_table.(data_sources{1}),'-k');
            title(['Left RCS ',data_sources{1}]);
        else
            plot(ax(1),aligned_data.left_taxis,aligned_data.left_LFP_table.key0,'-k');
            title('Left RCS key0');
        end
    end
    
    if isfield(aligned_data,'right_LFP_table')
        ax(2) = subplot(4,1,2);
        if ~isempty(data_sources{2})
            plot(ax(2),aligned_data.right_taxis,aligned_data.right_LFP_table.(data_sources{2}),'-k');
            title(['Right RCS ',data_sources{2}]);
        else
            plot(ax(2),aligned_data.right_taxis,aligned_data.right_LFP_table.key0,'-k');
            title('Right RCS key0');
        end
    end
    
%     if isfield(aligned_data,'Delsys')
%         ax(3) = subplot(4,1,3);
%         left_name_ind = find(cellfun(@(x) contains(x,'Left','IgnoreCase',true) && contains(x,'EMG','IgnoreCase',true),aligned_data.Delsys.Chan_names));
%         if ~isempty(left_name_ind)
%             plot(ax(3),aligned_data.Delsys.Time.(aligned_data.Delsys.Chan_names{left_name_ind}),aligned_data.Delsys.Data.(aligned_data.Delsys.Chan_names{left_name_ind}),'-k');
%             title('Delsys Left EMG');
%         end
%         
%         ax(4) = subplot(4,1,4);
%         right_name_ind = find(cellfun(@(x) contains(x,'Right','IgnoreCase',true) && contains(x,'EMG','IgnoreCase',true),aligned_data.Delsys.Chan_names));
%         if ~isempty(right_name_ind)
%             plot(ax(4),aligned_data.Delsys.Time.(aligned_data.Delsys.Chan_names{right_name_ind}),aligned_data.Delsys.Data.(aligned_data.Delsys.Chan_names{right_name_ind}),'-k');
%             title('Delsys Right EMG');
%         end
%     end
    linkaxes(ax,'x');
end

%% Acceleration
if isfield(aligned_data,'left_Accel_table') || isfield(aligned_data,'right_Accel_table')
    fig_hand = findobj('Name','Alignment Results - RCS vs Delsys acceleration movement');
    if ~isempty(fig_hand)
        close(fig_hand);
        figure('Name','Alignment Results - RCS vs Delsys acceleration movement');
    else
        figure('Name','Alignment Results - RCS vs Delsys acceleration movement');
    end
    if isfield(aligned_data,'left_Accel_table')
        ax(1) = subplot(4,1,1);
        if ~isempty(data_sources{3})
            plot(ax(1),aligned_data.left_accel_taxis,aligned_data.left_Accel_table.(data_sources{3}),'-k');
            title(['Left RCS Accel ',data_sources{3}]);
        else
            plot(ax(1),aligned_data.left_accel_taxis,aligned_data.left_Accel_table.XSamples,'-k');
            title('Left RCS Accel XSamples');
        end
    end
    
    if isfield(aligned_data,'right_Accel_table')
        ax(2) = subplot(4,1,2);
        if ~isempty(data_sources{4})
            plot(ax(2),aligned_data.right_accel_taxis,aligned_data.right_Accel_table.(data_sources{4}),'-k');
            title(['Right RCS Accel ',data_sources{4}]);
        else
            plot(ax(2),aligned_data.right_accel_taxis,aligned_data.right_Accel_table.XSamples,'-k');
            title('Right RCS Accel XSamples');
        end
    end
    
    if isfield(aligned_data,'Delsys')
        ax(3) = subplot(4,1,3);
        if contains(data_sources{5},'Acc','IgnoreCase',true)
            if contains(data_sources{5},'left','IgnoreCase',true) || contains(data_sources{5},'L','IgnoreCase',true) 
                left_name_ind = find(strcmp(data_sources{5},aligned_data.Delsys.Chan_names));
            else
                if contains(data_sources{5},'X','IgnoreCase',true)
                    axis_name = 'X';
                elseif contains(data_sources{5},'Y','IgnoreCase',true)
                    axis_name = 'Y';
                elseif contains(data_sources{5},'Z','IgnoreCase',true)
                    axis_name = 'Z';
                end

                if contains(data_sources{5},'IM')
                    axis_name = [axis_name,'_IM'];
                else
                    axis_name = [axis_name,'_G'];
                end

                left_name_ind = find(cellfun(@(x) contains(x,'Left','IgnoreCase',true) && contains(x,'DBS','IgnoreCase',true) && contains(x,'Acc','IgnoreCase',true) && contains(x,axis_name,'IgnoreCase',true),aligned_data.Delsys.Chan_names));
            end
        else
            left_name_ind = find(cellfun(@(x) contains(x,'Left','IgnoreCase',true) && contains(x,'Acc','IgnoreCase',true) && contains(x,'Y','IgnoreCase',true),aligned_data.Delsys.Chan_names));
        end
        if ~isempty(left_name_ind)
            plot(ax(3),aligned_data.Delsys.Time.(aligned_data.Delsys.Chan_names{left_name_ind}),aligned_data.Delsys.Data.(aligned_data.Delsys.Chan_names{left_name_ind}),'-k');
            if exist('axis_name','var')
                title(['Delsys Left IPG ',axis_name,' Acceleration']);
            else
                title('Delsys Left IPG Acceleration');
            end
        end
        
        
        ax(4) = subplot(4,1,4);
        if contains(data_sources{5},'Acc','IgnoreCase',true)
            if contains(data_sources{5},'right','IgnoreCase',true) || contains(data_sources{5},'R','IgnoreCase',true)
                right_name_ind = find(strcmp(data_sources{5},aligned_data.Delsys.Chan_names));
            else
                if contains(data_sources{5},'X','IgnoreCase',true)
                    axis_name = 'X';
                elseif contains(data_sources{5},'Y','IgnoreCase',true)
                    axis_name = 'Y';
                elseif contains(data_sources{5},'Z','IgnoreCase',true)
                    axis_name = 'Z';
                end

                if contains(data_sources{5},'IM')
                    axis_name = [axis_name,'_IM'];
                else
                    axis_name = [axis_name,'_G'];
                end

                right_name_ind = find(cellfun(@(x) contains(x,'Right','IgnoreCase',true) && contains(x,'DBS','IgnoreCase',true) && contains(x,'Acc','IgnoreCase',true) && contains(x,axis_name,'IgnoreCase',true),aligned_data.Delsys.Chan_names));
            end
        else
            right_name_ind = find(cellfun(@(x) contains(x,'Right','IgnoreCase',true) && contains(x,'Acc','IgnoreCase',true) && contains(x,'Y','IgnoreCase',true),aligned_data.Delsys.Chan_names));
        end
        if ~isempty(right_name_ind)
            plot(ax(4),aligned_data.Delsys.Time.(aligned_data.Delsys.Chan_names{right_name_ind}),aligned_data.Delsys.Data.(aligned_data.Delsys.Chan_names{right_name_ind}),'-k');
            if exist('axis_name','var')
                title(['Delsys Right IPG ',axis_name,' Acceleration']);
            else
                title('Delsys Right IPG Acceleration');
            end
        end
    end
    linkaxes(ax,'x');
end

%% Delsys vs Xsens
if isfield(aligned_data,'Delsys') && isfield(aligned_data,'Xsens')
    fig_hand = findobj('Name','Alignment Results - Delsys vs Xsens accleration movement');
    if ~isempty(fig_hand)
        close(fig_hand);
        figure('Name','Alignment Results - Delsys vs Xsens accleration movement');
    else
        figure('Name','Alignment Results - Delsys vs Xsens accleration movement');
    end
    
    if isfield(aligned_data,'Delsys')
        ax(1) = subplot(4,1,1);
        if contains(data_sources{5},'Acc','IgnoreCase',true)
            if contains(data_sources{5},'left','IgnoreCase',true) || contains(data_sources{5},'L','IgnoreCase',true)
                 left_name_ind = find(strcmp(data_sources{5},aligned_data.Delsys.Chan_names));
            else
                if contains(data_sources{5},'X','IgnoreCase',true)
                    axis_name = 'x';
                elseif contains(data_sources{5},'Y','IgnoreCase',true)
                    axis_name = 'y';
                elseif contains(data_sources{5},'Z','IgnoreCase',true)
                    axis_name = 'z';
                end
                left_name_ind = find(cellfun(@(x) contains(x,'Left','IgnoreCase',true) && contains(x,'DBS','IgnoreCase',true) && contains(x,'Acc','IgnoreCase',true) && contains(x,axis_name,'IgnoreCase',true),aligned_data.Delsys.Chan_names));
            end
        else
            left_name_ind = find(cellfun(@(x) contains(x,'Left','IgnoreCase',true) && contains(x,'Acc','IgnoreCase',true) && contains(x,'Y','IgnoreCase',true),aligned_data.Delsys.Chan_names));
        end
        if ~isempty(left_name_ind)
            plot(ax(1),aligned_data.Delsys.Time.(aligned_data.Delsys.Chan_names{left_name_ind}),aligned_data.Delsys.Data.(aligned_data.Delsys.Chan_names{left_name_ind}),'-k');
            if exist('axis_name','var')
                title(['Delsys Left IPG ',axis_name,' Acceleration']);
            else
                title('Delsys Left IPG Y Acceleration');
            end
        end
        
        
        ax(2) = subplot(4,1,2);
        if contains(data_sources{5},'Acc','IgnoreCase',true)
            if contains(data_sources{5},'right','IgnoreCase',true) || contains(data_sources{5},'R','IgnoreCase',true)
                right_name_ind = find(strcmp(data_sources{5},aligned_data.Delsys.Chan_names));
            else
                if contains(data_sources{5},'X','IgnoreCase',true)
                    axis_name = 'x';
                elseif contains(data_sources{5},'Y','IgnoreCase',true)
                    axis_name = 'y';
                elseif contains(data_sources{5},'Z','IgnoreCase',true)
                    axis_name = 'z';
                end
                right_name_ind = find(cellfun(@(x) contains(x,'Right','IgnoreCase',true) && contains(x,'DBS','IgnoreCase',true) && contains(x,'Acc','IgnoreCase',true) && contains(x,axis_name,'IgnoreCase',true),aligned_data.Delsys.Chan_names));
            end
        else
            right_name_ind = find(cellfun(@(x) contains(x,'Right','IgnoreCase',true) && contains(x,'Acc','IgnoreCase',true) && contains(x,'Y','IgnoreCase',true),aligned_data.Delsys.Chan_names));
        end
        if ~isempty(right_name_ind)
            plot(ax(2),aligned_data.Delsys.Time.(aligned_data.Delsys.Chan_names{right_name_ind}),aligned_data.Delsys.Data.(aligned_data.Delsys.Chan_names{right_name_ind}),'-k');
            if exist('axis_name','var')
                title(['Delsys Right IPG ',axis_name,' Acceleration']);
            else
                title('Delsys Right IPG Y Acceleration');
            end
        end
    end
    
    if isfield(aligned_data,'Xsens')
        ax(3) = subplot(4,1,3);
        if ~isempty(data_sources{6})
            var_name = replace(data_sources{6},'Right','Left');
            plot(ax(3),aligned_data.Xsens.Time,aligned_data.Xsens.(var_name),'-k');
            title(['Xsens ',var_name]);
        else
            plot(ax(3),aligned_data.Xsens.Time,aligned_data.Xsens.LeftShoulder_AccZ,'-k');
            title('Xsens Left Shoulder Z Acceleration');
        end
        
        ax(4) = subplot(4,1,4);
        if ~isempty(data_sources{6})
            var_name = replace(data_sources{6},'Left','Right');
            plot(ax(4),aligned_data.Xsens.Time,aligned_data.Xsens.(var_name),'-k');
            title(['Xsens ',var_name]);
        else
            plot(ax(4),aligned_data.Xsens.Time,aligned_data.Xsens.RightShoulder_AccZ,'-k');
            title('Xsens Right Shoulder Z Acceleration');
        end
    end
end
linkaxes(ax,'x');

%% Delsys TTL vs Force plate force values
if isfield(aligned_data,'FP') && isfield(aligned_data,'Delsys')
    fig_hand = findobj('Name','Alignment Results - Force plate vs Delsys TTL');
    if ~isempty(fig_hand)
        close(fig_hand);
        figure('Name','Alignment Results - Force plate vs Delsys TTL');
    else
        figure('Name','Alignment Results - Force plate vs Delsys TTL');
    end
    
    TTL_ind = find(cellfun(@(x) contains(x,'TTL','IgnoreCase',true),aligned_data.Delsys.Chan_names));
    
%     plot(aligned_data.FP.Time,aligned_data.FP.Fx_N);
%     hold on;
%     plot(aligned_data.FP.Time,aligned_data.FP.Fy_N);
%     plot(aligned_data.FP.Time,aligned_data.FP.Fz_N-aligned_data.FP.Fz_N(1));
%     ylim([min(min(aligned_data.FP.Fx_N),min(aligned_data.FP.Fy_N))*1.05,max(max(aligned_data.FP.Fx_N),max(aligned_data.FP.Fy_N))*1.05]);
%     yyaxis right;
%     plot(aligned_data.Delsys.Time.(aligned_data.Delsys.Chan_names{TTL_ind}),aligned_data.Delsys.Data.(aligned_data.Delsys.Chan_names{TTL_ind}),'-k');
    
    
end
end