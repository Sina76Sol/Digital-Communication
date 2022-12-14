function [Rx]=Rx_est(X,M)
% [Rx]=Rx_est(X,M)
% RX_EST  estimates the autocorrelation of the sequence of random
% variables given in X. For input M: Only Rx(0), Rx(1), ... , Rx(M) are computed.
% **Note: because index can only start at 1, output Rx(1) actually means Rx(0), and Rx(m) actually means Rx(m-1). 
N=length(X);
Rx=zeros(1,M+1);                      
for m=1:M+1,                          
  for n=1:N-m+1,
    Rx(m)=Rx(m)+X(n)*X(n+m-1);        
  end;
  Rx(m)=Rx(m)/(N-m+1);            
end;
