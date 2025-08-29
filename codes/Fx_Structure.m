function [costValue] = Fx_Structure(A,B)
C1=1e-10;
costValue=(Fx_COV(A,B)+C1)/(Fx_STD(A)*Fx_STD(B)+C1);
end

