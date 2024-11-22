function q = randq()
    n = rand(3,1);
    n = n/norm(n);
    th = rand(1,1)*2*pi;
    q = [cos(th/2) ; sin(th/2)*n];
    if q(1) < 0
        q = -q;
    end
end

