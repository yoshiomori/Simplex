function [ind v] = simplex(A, b, c, m, n, x)
 disp ("Simplex: Fase 2")
 nIt = 0;
 # Selecionando os indicies básicos e não básicos
 b = find (x != 0);
 nb = find (x == 0);
 invB = A(1:m, b) \ eye (m);
 # Fazendo as iterações do método simplex
 while (1)
  disp (cstrcat ("Iterando ", num2str (nIt++)));
  disp ([ b, x(b) ]);
  disp (cstrcat ("Valor da função objetivo: ", num2str (c' * x)));
  # Procurando algum custo reduzido negativo
  p = c(b)' * invB;
  cr = c(nb) - (p * A(1:m, nb))';
  disp ("Custos reduzidos");
  disp ([ nb, cr ]);
  j = nb(find (cr < 0, 1));
  if isempty (j)
   # Não foi encontrado custo reduzido negativo para um índice não básico, terminando o algoritmo
   ind = 0;
   v = x;
   disp (cstrcat ("Solução ótima foi encotrada com custo ", num2str(c' * v), ":"));
   disp (v);
   break;
  endif
  disp (cstrcat ("Entra na base: ", num2str (j)));
  # Foi encontrado custo reduzido negativo para um índice não básico, calculando o u, inverso da direção básica.
  u = invB * A(1:m, j);
  disp ("Direção");
  disp ([ b, -u ]);
  if max (u) <= 0
   # Nenhuma componete de u é positiva, terminando o algoritmo
   ind = -1;
   v = zeros (n, 1);
   v(b) = u;
   disp ("Problema é ilimitado com custo menos infinito:");
   disp (v);
   break;
  endif
  # Uma componente positiva de u foi encontrada, calculando theta mínimo e l, índice da componete de u que foi usado no calculo do theta mínimo.
  l = find(u>0);
  [ thetamin, i ] = min (x(b(l))./u(l));
  l = l(i);
  disp ("Theta*");
  disp (thetamin);
  disp (cstrcat ("Sai da base: ", num2str(b(l))));
  [ x(b(l)), b(l), nb(find(nb == j)) ] = deal (0, j, b(l));
  x(j) = thetamin;
  x(b(1:m!=l)) -= thetamin * u(1:m!=l);
  # Calculando a invB da próxima iteração
  I = eye(m);
  I(1:m,l) = u / -u(l);
  I(l,l) /= -u(l);
  invB = I * invB;
 endwhile
endfunction
