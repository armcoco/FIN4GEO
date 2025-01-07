function [u,v] = Relaxation(u,v,dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv)

[u,v] = GaussSeidel(u,v,dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
% [u,v] = Jacobi(u,v,dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
