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

## naloga

## Author: az <az@ares>
## Created: 2013-03-31

function [ ret ] = test ()
  sizex = 500;
  sizey = 1000;
  border = [0, 1; 0, 2];
  x = linspace(border(1, 1), border(1, 2), sizex);
  y = linspace(border(2, 1), border(2, 2), sizey);
  
  A = zeros(sizey, sizex);
  #A(1, :) = sin(x);
  #A(sizey, :) = sin(x);
  #A(:, 1) = cos(y) .- 1;
  #A(:, sizex) = cos(y) .- 1;

  hx = (border(1,2) - border(1,1)) / sizex;
  hy = (border(2,2) - border(2,1)) / sizey;

  [Z, b] = naredi_sistem(A, @function3, hx, border);

  tic;
  [B err] = cg(Z, b, zeros(size(b)));
  #[B err] = sor(Z, b, 1.9, zeros(size(b)));
  time = toc;

  printf("Porabil sem %f sekund\n", time)

  B = reshape(B, sizey-2, sizex-2);
 
  A(2:sizey-1, 2:sizex-1) = B;
  figure(1);
	semilogy(linspace(0, length(err), length(err)), err);
  figure(2);
  surf(x, y, A, 'EdgeColor','none','LineStyle','none','FaceLighting','phong');

endfunction

function [ ret ] = function1(x, y)
  ret =	6.*x.*y.*(1.-y)-2.*x.^3;
endfunction

function [ ret ] = function2(x, y)
	ret = x.*0;
endfunction

function [ ret ] = function3(x, y)
  ret =	x.*0 + 1;
endfunction
