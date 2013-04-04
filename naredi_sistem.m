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

## naredi_sistem

## Author: az <az@ares>
## Created: 2013-03-31

function [ Z, b ] = naredi_sistem (A, fun, h, border)

  if nargin < 4
    border = [0, 1; 0, 1]
  endif
  # m ... vrstice
  # n ... stolpci
  [m, n] = size(A);
  m = m - 2;
  n = n - 2;
  nm = n * m;

  tmp = ones(1, nm-1);
  idx = find(rem(1:nm-1, m) == 0);
  tmp(idx) = 0;

  i = [1:nm 2:nm 1:nm-1 1:nm-m m+1:nm];
  j = [1:nm 1:nm-1 2:nm m+1:nm 1:nm-m];
  s = [ones(1, nm) .* -4 tmp tmp ones(1, nm-m) ones(1, nm-m)];
  Z = sparse(i, j, s, nm, nm);

  xx = linspace(border(1, 1), border(1, 2), n);
  yy = linspace(border(2, 1), border(2, 2), m);
  [xxx, yyy] = meshgrid(xx, yy);
  b = 2 .* h .* fun(xxx, yyy);
  b = reshape(b, nm, 1);

  b(1:m) = b(1:m) - A(2:m+1, 1); #levi rob
  b(nm-m+1:nm) = b(nm-m+1:nm) - A(2:m+1, n+2); #desni rob
  idx_up = find(rem(1:nm, m)==1);
  idx_down = find(rem(1:nm, m)==0);
  b(idx_up) = b(idx_up) - A(1, 2:n+1)';
  b(idx_down) = b(idx_down) - A(m+2, 2:n+1)';

endfunction
