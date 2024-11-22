function R = q2R(q)
    th = acos(q(1))*2;
    if(th == 0)
        R = eye(3);
    else
        n  = q(2:4)/sin(th/2);

        R = eye(3) + sin(th) * skew(n) + (1-cos(th)) * (skew(n) * skew(n));
    end
    
%     R = [
%         1 - 2 * (q(3)*q(3) + q(4)*q(4)), 2 * (q(2)*q(3) - q(4)*q(1))  , 2 * (q(2)*q(4) + q(3)*q(1));
%         2 * (q(2)*q(3) + q(4)*q(1))    , 1 - 2 * (q(2)*q(2) + q(4)*q(4)), 2 * (q(3)*q(4) - q(2)*q(1));
%         2 * (q(2)*q(4) - q(3)*q(1))    , 2 * (q(3)*q(4) + q(2)*q(1))    , 1 - 2 * (q(2)*q(2) + q(3)*q(3))
%     ];

end

