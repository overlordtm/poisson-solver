## Copyright (C) 2013 az
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## cg (A, b, x)
## Metoda konjugiranih gradientov za resevanje sistema Ax=b
## vhodni parametri:
##   A: matrika
##   b: vektor
##   x: zacetni priblizek
## izhod:
##   x: resitev enacbe
##   err: gibanje napake

## Author: az <az@ares>
## Created: 2013-03-31

function [ x err ] = cg (A, b, x)

    r=b-A*x;
    p=r;
    rsold=r'*r;

    seps = sqrt(eps) * 10e4; 
    err = zeros(1000, 1);

    for i=1:1000
        Ap=A*p;
        alpha=rsold/(p'*Ap);
        x=x+alpha*p;
        r=r-alpha*Ap;
        rsnew=r'*r;
        err(i) = sqrt(rsnew);
        if sqrt(rsnew) < seps
              break;
        end
        p=r+rsnew/rsold*p;
        rsold=rsnew;
    end

    err = nonzeros(err);

endfunction
