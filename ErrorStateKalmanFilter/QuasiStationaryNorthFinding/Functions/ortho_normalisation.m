function [C_plus] = ortho_normalisation(C_minus)

[C_plus,R] = qr(C_minus);
% [C_plus] = orth(C_minus);
% norm(C_plus(:, 1))  


% C_1_minus = C_minus(:,1); 
% C_2_minus = C_minus(:,2);
% C_3_minus = C_minus(:,3);
% 
% Delta_12 = C_1_minus'*C_2_minus;
% Delta_13 = C_1_minus'*C_3_minus;
% Delta_23 = C_2_minus'*C_3_minus;
% 
% C_1_plus = C_1_minus - .5 * Delta_12 * C_2_minus - .5 * Delta_13 * C_3_minus;
% C_2_plus = C_2_minus - .5 * Delta_12 * C_1_minus - .5 * Delta_23 * C_3_minus;
% C_3_plus = C_3_minus - .5 * Delta_13 * C_1_minus - .5 * Delta_23 * C_2_minus;
% 
% C_1_plus_nrm = C_1_plus / sqrt(C_1_plus'*C_1_plus);
% C_2_plus_nrm = C_2_plus / sqrt(C_2_plus'*C_2_plus);
% C_3_plus_nrm = C_3_plus / sqrt(C_3_plus'*C_3_plus);
% 
% C_plus = [C_1_plus_nrm C_2_plus_nrm C_3_plus_nrm];

end