using MacroModelling
import StatsPlots

# 1) Modelo (apenas equações)
@model NK_COCHRANE begin
  # IS
  y[0]  = y[1] - (1/sigma)*( i[0] - pi[1] )

  # NKPC
  pi[0] = beta * pi[1] + kappa * y[0]

  # Regra de MP
  i[0]  = phi_pi * pi[0] + phi_y * y[0] + s_i[0]

  # v (predeterminada / backward)
  rho_b * v[0] = v[-1] + (omega * q[0] - q[-1]) - pi[0] - s[0]

  # Long-bond (q, forward) — com choque eps_q
  omega * q[1] - q[0] = i[0] + psi * std_q * eps_q[x]

  # Processos AR(1)
  s[0]   = rho_b * s[-1]   + std_b * eps_b[x]
  s_i[0] = rho_i * s_i[-1] + std_i * eps_i[x]
end;

# 2) Parâmetros (bloco separado)
@parameters NK_COCHRANE begin
  beta   = 0.99
  sigma  = 1.0
  kappa  = 0.1
  phi_pi = 0.8
  phi_y  = 0.0
  rho_b  = 0.9
  rho_i  = 0.8
  omega  = 0.9
  psi    = 0.0
  # desvios-padrão dos choques
  std_b  = 1e-12   # quase zero (como no teu .mod); se der problema, usa 0.0 ou 1e-6
  std_i  = 1.0
  std_q  = 1.0
end;

# 3) Resolver e IRFs
plot_irf(NK_COCHRANE)

