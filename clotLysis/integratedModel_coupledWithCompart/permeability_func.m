function k_ret = permeability_func(Rf,phi,H,RG)
% Rf: fibre radius in nm
% epsilon: porosity [-]
% phi: fibre volume fraction
% H: haematocrit [-]
% alp: rho_fibrinogen in plasma/rop_fibrinogen in fibrin fibre
% RG = retraction (0 for no retraction, 0.5 for 50% volume reduction)

if nargin==2
    % fibrin clot
    k_ret = (1e-18*Rf*Rf*4)./(70*(phi.^1.5).*( 1+52*phi.^1.5 ));
else 
    % Alpha value, source from Diamond's review paper, "Engineering...."
    if Rf>=300 % [nm]
        alp = 0.0171;
    elseif Rf<=100 % [nm]
        alp = 0.0143;
    else
        alp = 0.5*(0.0143+0.0171);
    end
    
    % compaction of firin within the volume surrounding the RBC
    C = ( (RG) - (H+(1-H)*alp))/((1-H)*(1-alp));
    epsilon = (C*(1-alp))/(alp+C*(1-alp));
    phi = 1-epsilon;
    
    %% Permeability
    % fibrin clot
    % k_fib = (1e-18*Rf*Rf*4)./(70*(phi.^1.5).*( 1+52*phi.^1.5 ));
    k_fib = (Rf*Rf*4)./(70*(phi.^1.5).*( 1+52*phi.^1.5 ));
    
    % retracted whole blood clot
    k_ret = k_fib*(( (1-H)*alp + C*(1-H)*(1-alp) )/( H + (1-H)*alp + C*(1-H)*(1-alp) ));
end

