function [resultat] = PmdPn(m,n)

if (n >= m+1)
    resultat = 1-(-1)^(m+n);
else
    resultat = 0;
end