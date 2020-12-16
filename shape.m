function [N,dN,jac] = shape(gp,xe)

  % local coordinate
  r = gp(1);
  s = gp(2);
  t = gp(3);
  
  % shape functions
  N = [ r, s, t, 1-r-s-t ];
  
  N_r = [ 1, 0, 0, -1 ];
  N_s = [ 0, 1, 0, -1 ];
  N_t = [ 0, 0, 1, -1 ];
  
  x_r = N_r * xe(:,1);
  x_s = N_s * xe(:,1);
  x_t = N_t * xe(:,1);
  %
  y_r = N_r * xe(:,2);
  y_s = N_s * xe(:,2);
  y_t = N_t * xe(:,2);
  %
  z_r = N_r * xe(:,3);
  z_s = N_s * xe(:,3);
  z_t = N_t * xe(:,3);

  %fprintf("xr %f xs %f xt %f,yr %f ys %f yt %f,zr %f zs %f zt %f\n",x_r,x_s,x_t,y_r,y_s,y_t,z_r,z_s,z_t);
 
  jacobian = [x_r, x_s, x_t;  y_r, y_s, y_t; z_r, z_s, z_t];
  jac = abs(det(jacobian));

  % check jacobian
  if jac <= 0.0
    error('Negative jacobian, element too distorted!');
  end

  inv_jac = inv(jacobian);
  
  dN = zeros(4,3);
  for i=1:4
    dN(i,:) = [N_r(i), N_s(i), N_t(i)] * inv_jac; %#ok<MINV>
  end
  
  %fprintf("N is %f %f %f %f\n",N(1),N(2),N(3),N(4));
  %fprintf("dN is %f %f %f %f %f %f %f %f %f %f %f %f\n",dN(1),dN(2),dN(3)...
      %,dN(4),dN(5),dN(6),dN(7),dN(8),dN(9),dN(10),dN(11),dN(12));
  %fprintf("jacobian is %f\n",jac);
  end
