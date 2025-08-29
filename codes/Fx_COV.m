function [COV_AB] = Fx_COV(A,B)
% 协方差
%   
COV_AB=mean((A-mean(A,"all")).*(B-mean(B,"all")),"all");
end