function err_C = calculateError_btw_simNexp(expData,simData)
global fileID
% expData: [time min, conc]
% simData: [time sec, conc]

timeSec = round(60*expData(:,1)); % [s]
err_C = [];
for kk=1:length(timeSec)
    fprintf(fileID, 'Data compare: %d out of %d\n',kk,length(timeSec));
    curTime = timeSec(kk);
%     indTime = find(simData(:,1)==curTime);
    indTime = curTime+1;
    C_sim = simData(indTime,2);
    % ---- Error between simulation and experiment ------------
    if curTime==0
        err_C_temp = [];
    else
        err_C_temp = (C_sim - expData(kk,2))/expData(kk,2);
    end
    err_C = [err_C, err_C_temp];
    % ---------------------------------------------------------
end