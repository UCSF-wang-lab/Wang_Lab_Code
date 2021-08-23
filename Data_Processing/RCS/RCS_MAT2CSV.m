function RCS_MAT2CSV(filename)
%% Create output variables
out_table = [];
subject_id = [];
time = [];
side = [];
pair1 = [];
pair2 = [];
pair3 = [];
pair4 = [];

%% Load data and check for left and right data
load(filename);
left_contact_names = [];
right_contact_names = [];

if isfield(aligned_data,'left_LFP_table')
    time = [time;aligned_data.left_taxis];
    side = [side;repelem('L',length(aligned_data.left_taxis),1)];
    pair1 = [pair1;aligned_data.left_LFP_table.key0];
    pair2 = [pair2;aligned_data.left_LFP_table.key1];
    pair3 = [pair3;aligned_data.left_LFP_table.key2];
    pair4 = [pair4;aligned_data.left_LFP_table.key3];
    
    left_contact_names = {aligned_data.DeviceSettings.Left.timeDomainSettings(end,:).TDsettings{1}.chanOut};
end

if isfield(aligned_data,'right_LFP_table')
    time = [time;aligned_data.right_taxis];
    side = [side;repelem('R',length(aligned_data.right_taxis),1)];
    pair1 = [pair1;aligned_data.right_LFP_table.key0];
    pair2 = [pair2;aligned_data.right_LFP_table.key1];
    pair3 = [pair3;aligned_data.right_LFP_table.key2];
    pair4 = [pair4;aligned_data.right_LFP_table.key3];
    
    right_contact_names = {aligned_data.DeviceSettings.Right.timeDomainSettings(end,:).TDsettings{1}.chanOut};
end

% Check if the contact names match, if they do great. If they don't...poop
if(~isempty(left_contact_names) && ~isempty(right_contact_names)) && sum(strcmp(left_contact_names,right_contact_names)) == length(left_contact_names)
    table_contacts = left_contact_names;
elseif ~isempty(left_contact_names) && isempty(right_contact_names)
    table_contacts = left_contact_names;
elseif isempty(left_contact_names) && ~isempty(right_contact_names)
    table_contacts = right_contact_names;
else
    error('Recording contact pairs do not match across hemispheres.');
end

% Create subject id cell vector
name_idx = strfind(filename,'RCS');
if ~isempty(name_idx)
    subject_id = repelem({filename(name_idx(1):name_idx(1)+4)},length(time),1);
else
    error('Filename does not contain the subject name');
end

%% Create output table
table_contacts = strrep(strrep(table_contacts,'+','p'),'-','m');
out_table = table(subject_id,side,time,pair1,pair2,pair3,pair4,'VariableNames',[{'SubjectID','Side','Time'},table_contacts]);
writetable(out_table,strrep(filename,'.mat','.csv'));
end