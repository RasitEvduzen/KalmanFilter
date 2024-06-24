function [Phi] = state_transition(F,Ts,model_order)
% Input Parameters
% F  -> System State Space Model
% Ts -> Sampling Periode
% model_order -> Taylor Series Truncation Order
% Output Parameters
% Phi -> State Transition Matrix

dim = size(F,1);

switch model_order
    case 1
        Phi = (eye(dim,dim)) + (F*Ts);        
    case 2
        Phi = (eye(dim,dim)) + (F*Ts) + (0.5*F^2*Ts^2);    
    case 3
        Phi = (eye(dim,dim)) + (F*Ts) + (0.5*F^2*Ts^2) + (1/6*F^3*Ts^3);
    % case 4
    % Phi = eye(dim,dim) + F*Ts + 0.5*F^2*Ts^2 + 1/6*F^3*Ts^3 + 1/24*F^4*Ts^4;
    % case 5
    % Phi = eye(dim,dim) + F*Ts + 0.5*F^2*Ts^2 + 1/6*F^3*Ts^3 + 1/24*F^4*Ts^4 + 1/120*F^5*Ts^5;
end

end