skew = @skew;
q2R  = @q2R;

%robot toolbox
quat2rotm([1 0 0 0])
rotm2quat(eye(3))
eul2rotm([1 2 3],"XYZ")
eul2rotm(eye(3),"XYZ")

%[V,D] = eig(A) %V: 3x3 (eigenvalues on diagonal)