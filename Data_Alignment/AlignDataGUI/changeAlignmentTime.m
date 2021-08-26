function changeAlignmentTime(src,event)
if contains(src.Tag,'pre')
    src.Parent.Parent.UserData.pre_alignment_time = str2num(src.String);
elseif contains(src.Tag,'post')
    src.Parent.Parent.UserData.post_alignment_time = str2num(src.String);
end

end