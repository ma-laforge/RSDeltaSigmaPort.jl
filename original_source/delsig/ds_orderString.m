function s=ds_orderString(n,capitalize)
if nargin<2
    capitalize=0;
end
strings = { 'zeroth' 'first' 'second' 'third' 'fourth' 'fifth' 'sixth' 'seventh' 'eighth' 'nineth' 'tenth' 'twelfth' 'thirteenth'};
if n<0 | n+1>length(strings)
    s = sprintf('%dth',n);
else
    s = strings{n+1};
    if capitalize
        s(1) = s(1) + 'A'-'a';
    end
end

    
