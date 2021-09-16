function calcRCS_AvgGaitCycleCoherence(aligned_data,varargin)

for i = 1:2:nargin-2
    switch varargin{i}
        case 'cycle_start_event'
            cycle_start_event = varargin{i+1};
        case 'n_percent_bins'
            n_percent_bins = varargin{i+1};
        case 'baseline_data'
            baseline_data = varargin{i+1};
        case 'baseline_normalization'
            baseline_normalization = varargin{i+1};
        case 'subjectID'
            subjectID = varargin{i+1};
        case 'save_flag'
            save_flag = varargin{i+1};
        case 'save_dir'
            save_dir = varargin{i+1};
    end
end

elements_to_compare = [];
if isfield(aligned_data,'left_LFP_table')
    elements_to_compare = [elements_to_compare,1:4];
end

if isfield(aligned_data,'right_LFP_table')
    elements_to_compare = [elements_to_compare,5:8];
end

if ~exist('cycle_start_event','var')
    cycle_start_event = 'LHS';
end

if ~exist('n_percent_bins','var')
    n_percent_bins = 100;
end

if ~exist('subjectID','var')
    subjectID = 'RCSXX';
end

if ~exist('save_flag','var')
    save_flag = 0;
end

if ~exist('save_dir','var') && save_flag == 1
    save_dir = uigetdir();
end

gait_events_sorted = sortGaitEvents(aligned_data.gait_events,cycle_start_event);

n_percent_bins = 100;

fig_vec = [];
pairs = nchoosek(elements_to_compare,2);
gait_cycle_avg = cell(1,size(pairs,1));

for i = 1:size(pairs,1)
    
    % Determine the indexing
    switch pairs(i,1)
        case {1,5}
            contact1 = 'key0';
        case {2,6}
            contact1 = 'key1';
        case {3,7}
            contact1 = 'key2';
        case {4,8}
            contact1 = 'key3';
    end
    
    switch pairs(i,2)
        case {1,5}
            contact2 = 'key0';
        case {2,6}
            contact2 = 'key1';
        case {3,7}
            contact2 = 'key2';
        case {4,8}
            contact2 = 'key3';
    end
    
    if pairs(i,1) <= 4
        side1 = 'left_LFP_table';
    else
        side1 = 'right_LFP_table';
    end
    
    if pairs(i,2) <= 4
        side2 = 'left_LFP_table';
    else
        side2 = 'right_LFP_table';
    end
    
    % Determine index limits
    n_val1 = height(aligned_data.(side1));
    n_val2 = height(aligned_data.(side2));
    [n_val,min_ind] = min([n_val1,n_val2]);
    
    if min_ind == 1
        if contains(side1,'left')
            time_vec = aligned_data.left_taxis;
        else
            time_vec = aligned_data.right_taxis;
        end
    else
        if contains(side2,'left')
            time_vec = aligned_data.left_taxis;
        else
            time_vec = aligned_data.right_taxis;
        end
    end
    
    % Get sample rate
    sr = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);
    
    % Calculate wavelet coherence and extract gait cycle values
    [x,~,y] = wcoherence(aligned_data.(side1).(contact1)(1:n_val),aligned_data.(side2).(contact2)(1:n_val),sr,'VoicesPerOctave',10);
    
    gait_cycle_mat = zeros(length(y),n_percent_bins,1);
    count = 1;
    for j = 1:height(gait_events_sorted)-1
        if ~isnan(gait_events_sorted.(cycle_start_event)(j)) && ~isnan(gait_events_sorted.(cycle_start_event)(j+1)) && (diff(gait_events_sorted.(cycle_start_event)(j:j+1)) < 2)
            [~,start_ind] = min(abs(time_vec-gait_events_sorted.(cycle_start_event)(j)));
            [~,end_ind] = min(abs(time_vec-gait_events_sorted.(cycle_start_event)(j+1)));
            data_snip = x(:,start_ind:end_ind);
            
            if sum(isinf(data_snip),'all') == 0
                percent_inds = round(linspace(1,size(data_snip,2),n_percent_bins+1));
                for k = 1:length(percent_inds)-1
                    if k == 1
                        gait_cycle_mat(:,k,count) = mean(data_snip(:,percent_inds(k):percent_inds(k+1)),2);
                    else
                        gait_cycle_mat(:,k,count) = mean(data_snip(:,percent_inds(k)+1:percent_inds(k+1)),2);
                    end
                end
                count = count + 1;
            end
        end
    end
    gait_cycle_avg{i} = mean(gait_cycle_mat,3,'omitnan');
    
    fig_vec(end+1) = figure();
    h = pcolor(1:100,log2(y),gait_cycle_avg{i});
    h.EdgeColor = 'none';
    yticks(log2([2.5,2.^(2:6)]));
    curr_y_ticks = yticks;
    yticklabels(cellstr(num2str(2.^curr_y_ticks')));
    ylim(log2([2.5,64]));
    shading interp;
    colormap jet;
    title({subjectID;createFigureTitle(side1,side2,contact1,contact2)});
    ylabel('Frequency (Hz)');
    xlabel('% Gait Cycle');
end

% Format and save figures
if save_flag
    figure_format(6,6,12,[],'painters');
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'Coherence'))
        mkdir(fullfile(save_dir,'Coherence'));
    end
    
    if ~isfolder(fullfile(save_dir,'Coherence',[aligned_data.stim_condition,'_STIM']))
        mkdir(fullfile(save_dir,'Coherence',[aligned_data.stim_condition,'_STIM']))
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if ~isfolder(fullfile(save_dir,'Coherence',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{n}))
            mkdir(fullfile(save_dir,'Coherence',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{n}));
        end
    end
    
    for j = 1:length(fig_vec)
        curr_axes = gca(fig_vec(j));
        save_name = [];
        for k = 1:length(curr_axes.Title.String)
            if isempty(save_name)
                save_name = curr_axes.Title.String{k};
            else
                save_name = [save_name,' ', curr_axes.Title.String{k}];
            end
        end
        
        savefig(fig_vec(i),fullfile(save_dir,'Coherence',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{1},strrep(save_name,' ','_')));
        
        for m = 2:length(folders_to_check)
            print(fig_vec(j),[fullfile(save_dir,'Coherence',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{m},strrep(save_name,' ','_')),extension{m}],'-r300',['-d',extension{m}(2:end)]);
        end
    end
end
end

function out_str = createFigureTitle(side1,side2,contact1,contact2)
out_str = [];
if contains(side1,'left')
    out_str = [out_str,'Left '];
    switch contact1
        case 'key0'
            out_str = [out_str, '+2-0'];
        case 'key1'
            out_str = [out_str, '+3-1'];
        case 'key2'
            out_str = [out_str, '+9-8'];
        case 'key3'
            out_str = [out_str, '+11-10'];
    end
else
    out_str = [out_str,'Right '];
    switch contact1
        case 'key0'
            out_str = [out_str, '+2-0'];
        case 'key1'
            out_str = [out_str, '+3-1'];
        case 'key2'
            out_str = [out_str, '+9-8'];
        case 'key3'
            out_str = [out_str, '+11-10'];
    end
end

out_str = [out_str,' to '];

if contains(side2,'left')
    out_str = [out_str,'Left '];
    switch contact2
        case 'key0'
            out_str = [out_str, '+2-0'];
        case 'key1'
            out_str = [out_str, '+3-1'];
        case 'key2'
            out_str = [out_str, '+9-8'];
        case 'key3'
            out_str = [out_str, '+11-10'];
    end
else
    out_str = [out_str,'Right '];
    switch contact2
        case 'key0'
            out_str = [out_str, '+2-0'];
        case 'key1'
            out_str = [out_str, '+3-1'];
        case 'key2'
            out_str = [out_str, '+9-8'];
        case 'key3'
            out_str = [out_str, '+11-10'];
    end
end

end