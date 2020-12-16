function [ke,fe] = weakform(xe,Psie,porosity)

  % 1 point formula - degree of precision 1
  gp =  [ 1/4, 1/4, 1/4 ];
  w =  1;

  ngp = length(w);
  
  % initialize stiffness matrix
  ke = zeros(4,4);
  
  % right hand size
  fe = zeros(4,1);
  
  % stress-strain displacement matrix
  B = zeros(1,4);
  % loop over gauss points
  for i=1:ngp
    [N,dN,jac] = shape(gp(i,:),xe);
    z = N * xe(:,3);
    if z > 30
        porosity = 0.0;
    end
    for j=1:4 % loop over local nodes
      B(j) = dN(j,3);
    end
    %fprintf("B is %f %f %f %f\n",B(1),B(2),B(3),B(4));
    % fill k
    ke = ke + N' * B * w(i) * jac;
    % fill f
    fe = fe - N' * ( ( porosity * Psie ) / ( 1 - porosity * ( 1 - Psie ) ) ) ...
        * w(i) * jac;
  end

end
