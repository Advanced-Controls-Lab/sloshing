
function [solution, its] = stewart_forward(a,l,num)
    tolf = 10^-3;
    tola = 10^-3;
    platform_parameters=[130,25,90,170]; %p r, p theta, b theta, b r

    x=a(1);
    y=a(2);
    z=a(3);
    alpha=a(4);
    beta=a(5);
    gamma=a(6);

    sumf=0; % to sum fi(a)
    n=0;m=0;

    % define func, B for i 1-6
    func = zeros(1,6);
    A = zeros(6,6);
    B = zeros(6,1);

    % compute rotation matrix
    R = [cos(gamma)*cos(beta), cos(gamma)*sin(beta)*sin(alpha)-sin(gamma)*cos(alpha), cos(gamma)*sin(beta)*cos(alpha)+sin(gamma)*sin(alpha);
        sin(gamma)*cos(beta), sin(gamma)*sin(beta)*sin(alpha)+cos(gamma)*cos(alpha), sin(gamma)*sin(beta)*cos(alpha)-cos(gamma)*sin(alpha);
        -sin(beta), cos(beta)*sin(alpha), cos(beta)*cos(alpha)];

    for i = 1:6
        b=mod(i,2);
        if b==1
            n=(60*i)-(platform_parameters(2)/2);
            m=(60*i)-(platform_parameters(3)/2); 
        else
            n=n+platform_parameters(2);
            m=m+platform_parameters(3);
        end
        
        
        % calculate xb = db - bb
        bB=[platform_parameters(4)*cos((m/180.0)*pi);
            platform_parameters(4)*sin((m/180.0)*pi);
            0];
        xB=[x-bB(1),y-bB(2),z-bB(3)];

        % calculate pb = R*pp
        pP=[platform_parameters(1)*cos((n/180.0)*pi);
            platform_parameters(1)*sin((n/180.0)*pi);
            0];
        pB = R*pP;


        % compute fi(a)
        func(i) = (xB(1)+pB(1))^2 + (xB(2)+pB(2))^2 + (xB(3)+pB(3))^2 - l(i)^2;
        sumf = sumf + abs(func(i)); % add to sum of fi(a)
        B(i) = -func(i);

        % compute A partial derivatives
        A(i,1)=2*(xB(1)+pB(1));
        A(i,2)=2*(xB(2)+pB(2));
        A(i,3)=2*(xB(3)+pB(3));
        A(i,4)=2*(-xB(1)*pB(2)+xB(2)*pB(1));
        A(i,5)=2*(pB(3)*(-xB(1)*cos(gamma)+xB(2)*sin(gamma))-xB(3)*(pP(1)*cos(beta)+pP(2)*sin(beta)*sin(alpha)));
        A(i,6)=2*pP(2)*(xB(1)*R(1,3)+xB(2)*R(2,3)+xB(3)*R(3,3));
    
    end

    if sumf < tolf
        solution = a;
        its = num;
        return
    end

    % LU decomp step
    % solve for delta a (j): A*deltaa=B
    deltaa = B\A;
    if sum(abs(deltaa)) < tola
        solution = a;
        its = num;
        return
    end

    anew = a + deltaa
    num=num+1;
    [solution, its] = stewart_forward(anew,l,num);

end