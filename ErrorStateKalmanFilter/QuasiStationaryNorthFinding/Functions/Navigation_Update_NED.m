function [L_b_plus,lambda_b_plus,H_b_plus,C_b_n_plus,V_eb_n_plus] = Navigation_Update_NED(Ts,L_b_minus,Lambda_b_minus,H_b_minus,C_bn_minus,V_eb_n_minus,f_ib_b,w_ib_b,wander_angle)
% -------------------------------------------------
% ------- Sub Scripts & Input Parameters ----------
% -------------------------------------------------
% i -> ECI, e -> ECEF, n -> NED, b -> Body, w -> Wander
% Ts -> Sampling Periode [Sn]
% L_b -> Latitude [Rad]
% Lambda_b -> Longitude [Rad]
% H_b -> Height [M]
% C_bn -> C_b^n Coordinate Transform Matrix b to n
% V_eb_n -> V_eb^n  Velocity vector [Vnorth Veast Vdown] b to n [M/Sn]
% f_ib_b -> Accelerometer Data [ax ay az] [M/Sn^2]
% w_ib_b -> Gyroscope Data  [wx wy wz]  [Rad/Sn]
% wander_ange -> [Rad]
% W_ie -> Earth Angular Rate ECEF to ECI [Rad/Sn]

% -------------------------------------------------
% -------------- Output Parameters ----------------
% -------------------------------------------------
% L_b_plus      -> 
% lambda_b_plus ->
% h_b_plus      ->
% v_eb_n_plus   ->
% C_b_n_plus    ->
% -------------------------------------------------
% -------------------------------------------------

global W_ie Eq_rad Pol_rad Flat Eccen mu

alpha_ib_b = w_ib_b * Ts;
mag_alpha = sqrt(alpha_ib_b' * alpha_ib_b);
Alpha_ib_b = skew_symmetric(alpha_ib_b); 

% angular rate of the ECEF to ECI frame, resolved Wander frame [eq 15.23]
% omega_ie_w = W_ie * [cos(wander_angle)*cos(L_b_minus); sin(wander_angle)*cos(L_b_minus); sin(L_b_minus)];

% angular rate of the ECEF to ECI frame, resolved NED frame [eq 5.41]
omega_ie_n = skew_symmetric(W_ie * [cos(L_b_minus); 0; -sin(L_b_minus)]); 

% Eq(5.44) w.r.t the ECEF frame, resolved about NED
[R_N_minus,R_E_minus] = Rad_Curv(L_b_minus);
omega_en_n_minus = [ V_eb_n_minus(2)/(R_E_minus+H_b_minus);
                    -V_eb_n_minus(1)/(R_N_minus+H_b_minus);
                    -V_eb_n_minus(2)*tan(L_b_minus)/(R_E_minus+H_b_minus)];

% Force Frame Transform
% avarage Cbn [eq 5.84 & 5.86]  !Maybe remove *Ts therm in below equations!
if mag_alpha>1.e-8
    Cbb = eye(3) + (((1-cos(mag_alpha))/mag_alpha^2)*Alpha_ib_b) + (((1-(sin(mag_alpha)/mag_alpha))/mag_alpha^2)*Alpha_ib_b*Alpha_ib_b); 
    C_b_n_avarage = C_bn_minus*Cbb-0.5*(skew_symmetric(omega_en_n_minus)+omega_ie_n)*C_bn_minus*Ts;
else
     C_b_n_avarage = C_bn_minus-0.5*(skew_symmetric(omega_en_n_minus)+omega_ie_n)*C_bn_minus*Ts;
end     
f_ib_n = C_b_n_avarage*f_ib_b';

% Velocity Update [eq 5.54]
V_eb_n_plus = V_eb_n_minus + Ts*(f_ib_n+Gravity_NED(L_b_minus,H_b_minus)-(skew_symmetric(omega_en_n_minus)+2*omega_ie_n)*V_eb_n_minus);

% Curvilinear Pose Update [eq 5.56]
H_b_plus = H_b_minus+0.5*Ts*(V_eb_n_minus(3) + V_eb_n_plus(3));

L_b_plus = L_b_minus+0.5*Ts*(V_eb_n_minus(1)/(R_N_minus+H_b_minus) + V_eb_n_plus(1)/(R_N_minus+H_b_plus));

[R_N_plus,R_E_plus]= Rad_Curv(L_b_plus);
lambda_b_plus = Lambda_b_minus+0.5*Ts*(V_eb_n_minus(2)/((R_E_minus+H_b_minus)*cos(L_b_minus))+V_eb_n_plus(2)/((R_E_plus+H_b_plus)*cos(L_b_plus))); 

% Attitude Update [eq 5.44]
omega_en_n_plus = [ V_eb_n_plus(2)/(R_E_plus+H_b_plus);
                   -V_eb_n_plus(1)/(R_N_plus+H_b_plus);
                   -V_eb_n_plus(2)*tan(L_b_plus)/(R_E_plus+H_b_plus)];
%  [eq 5.73]
if mag_alpha>1.e-8
    C_Bplus_Bminus = [eye(3)] + [sin(mag_alpha)/mag_alpha*Alpha_ib_b]+...
        [(1-cos(mag_alpha))/mag_alpha^2*Alpha_ib_b*Alpha_ib_b];
else
    C_Bplus_Bminus =eye(3)+Alpha_ib_b;
end     
    
% Update attitude  [eq 5.77]
C_b_n_plus = (eye(3)-(omega_ie_n+0.5*skew_symmetric(omega_en_n_minus)+0.5*skew_symmetric(omega_en_n_plus))*Ts)*C_bn_minus*C_Bplus_Bminus;
C_b_n_plus = ortho_normalisation(C_b_n_plus);

end