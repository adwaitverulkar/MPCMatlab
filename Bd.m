function Bd = Bd(in1,in2)
%Bd
%    Bd = Bd(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    12-May-2023 11:55:44

delta = in1(7,:);
vx = in1(4,:);
Bd = reshape([0.0,0.0,0.0,3.333333333333333e-5,delta./6.6e+4,delta./8.91e+4,0.0,0.0,0.0,0.0,0.0,vx./4.4e+1,vx.*(5.0./2.97e+2),1.0./2.0e+1],[7,2]);
