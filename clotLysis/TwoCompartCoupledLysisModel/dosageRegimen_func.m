function [tPAdose_array, timePoints, t_ramping] = dosageRegimen_func(tPAdose,N_stage,t_infusion,t_delay,perc_infusion)
global Mw_tPA V_plasma Vc

%% Drug and patient specs
% --- Patient spec
% =========================================================================
weight = 80;                    % weight of patient [kg]
% =========================================================================
V_plasma = Vc*weight;           % volume of patient [L]3.9
Mw_tPA = 59042.3;               % molecular weight of tPA [g/mol]

% --- tPA specs
% =========================================================================
% tPAdose = 0.75;                     % total tPA dose in [mg]
% =========================================================================
tPAdose_tot = 1e-3*tPAdose;	% tPA dose [g]

%% Treatment protocol
% =========================================================================
% N_stage = 1;
t0 = 0;                     % start of infusion [s]
% t_infusion = [75]*60;    % infusion times [s]
% t_delay = []*60;           % time delay between infusion modes [s]
t_ramping = 10;             % ramping period when switching infusion mode [s]
% perc_infusion = [1];    % percentage of each infusion of total amount [-]
% =========================================================================
tPAdose_each = tPAdose_tot*perc_infusion;   % tPA dose in each dosing stage [g]
% input data check
if length(t_infusion)~=N_stage; error('No. of input for "infusion times" not matching!!!'); end
if length(t_delay)~=(N_stage-1); error('No. of input for "time delay" not matching!!!'); end
if sum((t_infusion-t_ramping)>0)==0; error('Adjust the ramping span, t_ramp < t_infusion.'); end
if length(perc_infusion)~=N_stage; error('No. of input for "ratio infusion" not matching!!!'); end
if sum(perc_infusion)~=1; error('Sum of infusion percentage should be 1!!!!'); end

% create doses for each stage
N_dose = N_stage + sum(t_delay>0);
tPAdose_array = zeros(1,N_dose);
indDel = 0;
for ii=1:N_stage
    if isempty(t_delay) || ii==1
        tPAdose_array(ii) = tPAdose_each(ii);
        if ~isempty(t_delay) && t_delay(ii)>0
            indDel = indDel+1;
        end
    else
        indd = ii + indDel;
        tPAdose_array(indd) = tPAdose_each(ii);
        if t_delay(ii-1)>0
            indDel = indDel+1;
        end
    end
end

% create time points
N_timePoints = (N_stage+1) + sum(t_delay>0);
timePoints = zeros(1,N_timePoints);
timePoints(1) = t0; pp = 1; tt = 1;
while (pp<N_timePoints)
    pp = pp+1; tt = tt+1;
    timePoints(pp) = timePoints(pp-1) + t_infusion(tt-1);
    if length(t_delay)>=(tt-1) && t_delay(tt-1)>0
        pp = pp+1;
        timePoints(pp) = timePoints(pp-1) + t_delay(tt-1);
    end
end

