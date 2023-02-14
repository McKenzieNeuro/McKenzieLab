function [R] = rot2d(a)
%ROT2D Return 2D rotation matrix
%   Returns matrix that rotates a column vector by angle a (in radians).
%   Example 1: v=[1;1]; w=rot2d(pi/4)*v; {w now = [0;1.414]}
%   Example 2: Rotate 3 (or more) [x;y] pairs at a time: 
%       v1=[1,0]; v2=[1;1]; v3=[0;1]; vc=[v1,v2,v3]; wc=rot2d(pi/6)*vc;
%   The examples show the correct multiplication order.
%   Note that this is not the matrix to use if you want to express a
%   vector in terms of a coordinate system rotated by angle a. 
%   Use w=inv(rot2d(a))*v to express vector v in a coordinate system
%   rotated by a.
R = [cos(a),-sin(a);sin(a),cos(a)];
end