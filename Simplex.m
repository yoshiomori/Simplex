function Q = pivoting(u, m, l)
  # função que recebe um vetor coluna u, que vai ser pivoteado,
  # sua dimenção m e a l-ésima componente de u que será o pivo.
  Q = eye(m);
  Q(1:m,l) = u / -u(l);
  Q(l,l) /= -u(l);
endfunction

function [B, N, invB, v, ind] = simplex_iteration(A, c, m, n, B, N, invB, x)
 nIt = 0;
 # Fazendo as iterações do método simplex
 while (1)
  disp (cstrcat ("Iterando ", num2str (nIt++)));
  disp ([ B, x(B) ]);
  disp (cstrcat ("Valor da função objetivo: ", num2str (c' * x)));
  # Procurando algum custo reduzido negativo
  p = c(B)' * invB;
  cr = c(N) - (p * A(1:m, N))';
  disp ("Custos reduzidos");
  disp ([ N, cr ]);
  j = N(find (cr < 0, 1));
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
  disp ([ B, -u ]);
  if max (u) <= 0
   # Nenhuma componete de u é positiva, terminando o algoritmo
   ind = -1;
   v = zeros (n, 1);
   v(B) = u;
   disp ("Problema é ilimitado com custo menos infinito:");
   disp (v);
   break;
  endif
  # Uma componente positiva de u foi encontrada, calculando theta mínimo e l, índice da componete de u que foi usado no calculo do theta mínimo.
  l = find(u>0);
  [ thetamin, i ] = min (x(B(l))./u(l));
  l = l(i);
  disp ("Theta*");
  disp (thetamin);
  disp (cstrcat ("Sai da base: ", num2str(B(l))));
  [ x(B(l)), B(l), N(N == j) ] = deal (0, j, B(l));
  x(j) = thetamin;
  x(B(1:m!=l)) -= thetamin * u(1:m!=l);
  # Calculando a invB da próxima iteração
  Q = pivoting(u, m, l);
  invB = Q * invB;
 endwhile
endfunction

function [ind v] = simplex(A, b, c, m, n)
 # Inicializando problema auxiliar
 A_aux = [ A, eye(m) ];
 c_aux = [ zeros(n, 1); ones(m, 1) ];
 n_aux = m + n;
 x_aux = [ zeros(n, 1); b ];
 B = [ n + 1:n + m ]';
 N = [ 1:n ]';
 invB = eye(m);
 
 # Fazendo as iterações do método simplex da fase 1
 [B, N, invB, v, ind] = simplex_iteration(A_aux, c_aux, m, n_aux, B, N, invB, x_aux);
 
 # Verificando a viabilidade do problema
 if c_aux'*v
  # Problema é inviável, terminando o algoritmo
  ind = 1;
  break
 endif
 
 # Conduzindo a variável artificial para fora da base
 for l = find(B>n)'
  j = find(A(l, N(N<=n)'), 1);
  [ B(l), N(N == j) ] = deal(j, B(l));
  u = invB * A(1:m, j);
  Q = pivoting(u, m, l);
  invB = Q * invB;
 endfor
 
 # Voltando ao problema original
 x = v(1:n);
 N = N(N<=n);
 
 disp ("Simplex: Fase 2")
 # Fazendo as iterações do método simplex da fase 2
 [B, N, invB, v, ind] = simplex_iteration(A, c, m, n, B, N, invB, x);
endfunction

%!test # Solução ótima encotrada
%! A = [ 1, 1, 1, 1; 2, 0, 3, 4 ];
%! b = [ 2; 2 ];
%! c = [ 2; 2; 2; 2 ];
%! m = 2;
%! n = 4;
%! [ind, v] = simplex(A, b, c, m, n);
%! assert (ind, 0);
%! assert (v, [ 1; 1; 0; 0 ]);

%!test # Custo ótimo menos infinito
%! A = [ -1, 1, 1, 0; 1, -2, 0, 1 ];
%! b = [ 1; 2 ];
%! c = [ -2; -1; 0; 0 ];
%! m = 2;
%! n = 4;
%! [ind, v] = simplex(A, b, c, m, n);
%! assert (ind, -1);
%! assert (v, [ -2; 0; -1; 0 ]);
