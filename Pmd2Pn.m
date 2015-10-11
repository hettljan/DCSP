% Pour LambAnisotrope3DLegendre

function [resultat] = Pmd2Pn(mm,nn)

if (nn >= mm+2)
    resultat = 0.5*(1-(-1)^(mm+nn+1))*(nn*(nn+1)-mm*(mm+1));
else
    resultat = 0;
end