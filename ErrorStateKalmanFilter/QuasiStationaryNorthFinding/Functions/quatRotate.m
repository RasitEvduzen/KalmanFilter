function v = quatRotate(v, q)


[row col] = size(v);
qRot = quatMultiplication(quatMultiplication(q, [zeros(row, 1) v]), quatConjugate(q));
v = qRot(:, 2:4);


end

