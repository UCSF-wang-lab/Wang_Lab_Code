function changeStimMedCondition(src,obj)

switch src.Tag
    case 'Stim_condition_drop_down'
        src.Parent.Parent.UserData.stim_condition = src.String{src.Value};
    case 'Med_condition_drop_down'
        src.Parent.Parent.UserData.med_condition = src.String{src.Value};
end
end