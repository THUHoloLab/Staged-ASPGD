function [STD_A] = Fx_STD(A)
% 标准差
STD_A=sqrt(mean((A-mean(A,"all")).^2,'all'));
end

