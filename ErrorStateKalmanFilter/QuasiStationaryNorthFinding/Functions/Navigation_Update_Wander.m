function [C_b_w_plus,V_eb_w_plus,R_eb_w_plus,W_ie_w,W_ang_plus] = Navigation_Update_Wander(Ts,L_b,Lambda_b,H_b,C_bw_minus,V_eb_w_minus,R_eb_w_minus,F_ib_b,W_ib_b,W_ang)
% -------------------------------------------------
% ------- Sub Scripts & Input Parameters ----------
% -------------------------------------------------
% (i) -> (ECI), (e) -> (ECEF), (n) -> (NED), (b) -> (Body), (w) -> (Wander)
% -------------------------------------------------
% Ts       -> Sampling Periode [Sn]
% L_b      -> Latitude [Rad]
% Lambda_b -> Longitude [Rad]
% H_b      -> Height [M]
% C_bw     -> C_bw Coordinate Transform Matrix b -> w
% V_eb_n   -> V_eb^n  Velocity vector [Vnorth Veast Vdown] b to n [M/Sn]
% F_ib_b   -> Accelerometer Data [ax ay az] [M/Sn^2]
% w_ib_b   -> Gyroscope Data  [wx wy wz]  [Rad/Sn]
% W_ang    -> Wander Angle[Rad]
% -------------------------------------------------
% -------------- Output Parameters ----------------
% -------------------------------------------------
% V_eb_w_plus   -> 
% C_b_w_plus    -> 
% R_eb_w_plus   -> 
% W_ie_w        -> Earth Angular Velocity Sensing From Wander Frame
% -------------------------------------------------
% -------------------------------------------------
global W_ie

% Earth Angular vel transform Eq [15.23]  % TRUE!
W_ie_w = W_ie * [cos(W_ang)*cos(L_b) -sin(W_ang)*cos(L_b) -sin(L_b)]';  

% Attitude Update  eq [15.18]   % TRUE!
C_b_w_plus = ((C_bw_minus) * (eye(3) + skew_symmetric(W_ib_b) * Ts)) - ((W_ie*Ts)*([0 sin(L_b) -sin(W_ang)*cos(L_b); -sin(L_b) 0 -cos(W_ang)*cos(L_b); sin(W_ang)*cos(L_b) cos(W_ang)*cos(L_b) 0]*C_bw_minus)); 
C_b_w_plus = ortho_normalisation(C_b_w_plus); % Orthonormalization decomposition

% Force Transform 5.28 Update
F_ib_w = .5 * (C_b_w_plus + C_bw_minus)* F_ib_b;

% Velocity Update eq[15.19] 
grav_vec = Gravity_NED(L_b,H_b);
V_eb_w_plus = (V_eb_w_minus) + (F_ib_w * Ts) + ([0 0 grav_vec(3,1)]' * Ts);

% Pose Update eq[15.20]
R_eb_w_plus = (R_eb_w_minus) + ((Ts/2) * (V_eb_w_minus + V_eb_w_plus));

W_ang_plus = W_ang; % Wander Angle Update
end