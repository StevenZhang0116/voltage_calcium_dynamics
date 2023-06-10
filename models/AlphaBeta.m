function [An, Bn, Am, Bm, Ah, Bh] = AlphaBeta(V, T)
   global Vr

   V = V*1000;
   dV = (V - Vr);
   phi = 3^((T-6.3)/10);

   An = phi * (eps + 0.10 - 0.01 .* dV) ./ (eps + exp(1 - 0.1 .* dV) - 1);
   Am = phi * (eps + 2.5 - 0.1 .* dV) ./ (eps + exp(2.5 - 0.1 .* dV) - 1);
   Ah = phi * 0.07 .* exp(-dV ./ 20);

   Bn = phi * 0.125 .* exp(-dV ./ 80);
   Bm = phi * 4 .* exp(-dV/18);
   Bh = phi * 1 ./ (exp(3.0 - 0.1 .* dV) + 1);
end
  