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

## Pozene testni primer
## Za urejanje parametrov je potrebno poeditirati file

## Author: az <az@ares>
## Created: 2013-03-31

function [ ret ] = testi ()

  # nastavi velikost mreze
  sizex = 250; 
  sizey = 250;
  # nastavi meje (pazi da so v razmerju z mrezo)
  border = [0, 1; 0, 1];
  x = linspace(border(1, 1), border(1, 2), sizex);
  y = linspace(border(2, 1), border(2, 2), sizey);
  
  A = zeros(sizey, sizex);
  # nastavi robove (ce noces da so 0)
  #A(1, :) = sin(x);
  #A(sizey, :) = sin(x);
  #A(:, 1) = cos(y) .- 1;
  #A(:, sizex) = cos(y) .- 1;

  # izberi funkcijo: function[1-3]
  fun = @function4;

  # izberi metodo method=["sor"|"cg"]
  method = "sor"

  # ne diraj lava dok spava spodaj

  hx = (border(1,2) - border(1,1)) / sizex;
  hy = (border(2,2) - border(2,1)) / sizey;

  assert(abs(hx-hy)<eps, "h ni enak v obeh dimenzijah!");

  [Z, b] = naredi_sistem(A, fun, hx, border);

  tic;
  if strcmp(method, "cg")
    [B err] = cg(Z, b, zeros(size(b)));
  else
    [wy wx] = opt_relax_fact(A);
    [B err] = sor(Z, b, (wx+wy)/2, zeros(size(b)));
  end
  time = toc;

  printf("Porabil sem %f sekund\n", time)

  B = reshape(B, sizey-2, sizex-2);
 
  A(2:sizey-1, 2:sizex-1) = B;
  figure()
  semilogy(linspace(0, length(err), length(err)), err);

  figure();
  surf(x, y, A, 'EdgeColor','none','LineStyle','none','FaceLighting','phong');

endfunction

function [wy, wx] = opt_relax_fact(A)
  
  # najde optimalni relaxacijski fator za obe dimenzije matrike

  [sizey, sizex] = size(A);
  wy = 2/(1+pi/sizey)
  wx = 2/(1+pi/sizex)
end

function [ ret ] = function1(x, y)
  # testna funckija 1
  ret =	6.*x.*y.*(1.-y)-2.*x.^3;
endfunction

function [ ret ] = function2(x, y)
  # testna funkcija 2
	ret = x.*0;
endfunction

function [ ret ] = function3(x, y)
  # testna funckija 3
  ret =	x.*0 + 1;
endfunction

function [ ret ] = function4(x, y)
  # testna funckija 3
  zoki =	zeros(length(y), length(x));
  zoki(length(y)/2, length(x)/2) = 1;
  ret = zoki;
endfunction
