function y = wspair(x)

if('A'==x || 'a'==x)
    y = 'T';
elseif('G'==x || 'g'==x)
    y = 'C';
elseif('C'==x || 'c'==x)
    y = 'G';
elseif('T'==x || 't'==x)
    y = 'A';
else
    error('Illegal base.');
end

end