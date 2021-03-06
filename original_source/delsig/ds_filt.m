function y = ds_filt(H,x)
% function y = ds_filt(H,x)
% H is an LxM array of zpk transfer functions
% x is an MxN matrix consisting of M time sequences of length N
% y is an LxN matrix consisting of L time sequences of length N
L = size(H,1);
M = size(H,2);
if M ~= size(x,1)
    error('size mismatch');
end
N = size(x,2);
y = zeros(L,N);
for i = 1:L
    for j = 1:M
        b = H(i,j).k * poly(H(i,j).z{1});
        a = poly(H(i,j).p{1});
        m = length(b);
        n = length(a);
        if n>m
            b = [zeros(1,n-m) b];
        elseif n<m
%             a = [zeros(m-n,1) b];
            error( 'H must have at least as many poles as zeros.' );
        end
        y(i,:) = y(i,:) + filter(b,a,x(j,:));
    end
end
