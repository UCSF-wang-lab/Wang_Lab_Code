function changeTrialNum(src,obj)

switch src.Tag
    case 'trial_num_edit'
        src.Parent.Parent.UserData.trial_num = str2num(src.String);
end
end