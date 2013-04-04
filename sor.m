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

## sor

## Author: az <az@ares>
## Created: 2013-04-04

function [ x err ] = sor (A, b, w, x0)
  D = diag(diag(A));
  L = -tril(A,-1);
  U = -triu(A,1);

  DL = D - w*L;
  UE = (1-w)*D + w*U;
  bE = w*b;
  
  err = zeros(1000, 1);
  seps = sqrt(eps) * 10e4;

  for it = 1:1000
    rhv = UE*x0 + bE;
    x = DL\rhv;
    e = norm(x-x0, Inf);
    err(it) = e;
    if e < seps
        return
    else
        x0 = x;
    end
  end
endfunction
