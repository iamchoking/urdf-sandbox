function R = elemR(th,axis)
    if axis == 'x'
        R = [1 0 0;0 cos(th) -sin(th);0 sin(th) cos(th)];
    elseif axis == 'y'
        R = [cos(th) 0 sin(th);0 1 0;-sin(th) 0 cos(th)];
    elseif axis == 'z'
        R = [cos(th) -sin(th) 0;sin(th) cos(th) 0; 0 0 1];
    else
        R = eye(3);
    end
end

