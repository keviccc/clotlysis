function IR = infusionRate_func_general(time, tPAdose_array, timePoints, t_ramping)
% IR, infusion rate in [g/s]
% time, current time in [s]
% tPAdose_array, tPA dose in each stage [g]
% timePoints: t0,tb,td,tc in time
% t_ramping: ramping duration when switching infusion modes

t_inteval = timePoints(2:end)-timePoints(1:end-1);
I_array = tPAdose_array./t_inteval; % rate of infusion [g/s]
I_array_new = I_array./(1-0.5*t_ramping./t_inteval);

tempVar = timePoints-time;
yy = 1; currentStage = [];
while (yy<length(tempVar))
    if tempVar(yy)*tempVar(yy+1)<=0
        currentStage = yy;
        break
    else
        yy = yy+1;
    end
end

if isempty(currentStage)
    IR = 0;
else
    t_lb = timePoints(currentStage);
    t_ub = timePoints(currentStage+1);
    
    slp = I_array_new(currentStage)/(0.5*t_ramping);
    if time<=(t_lb+0.5*t_ramping)
        IR = slp*time - slp*t_lb;
    elseif time>=(t_ub-0.5*t_ramping)
        IR = -slp*time + slp*t_ub;
    else
        IR = I_array_new(currentStage);
    end
end

