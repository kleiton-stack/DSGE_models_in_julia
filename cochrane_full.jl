using MacroModelling
import StatsPlots

####### Equações do modelo ###########
@model NK_COCHRANE_FULL begin
  # IS (gap ~ y)
  y[0]  = y[1] - sigma * (i[0] - pi[1])

  # Oferta de trabalho (MRS): w = σ c + ϕ n
  w[0]  = sigma * c[0] + varphi * n[0]

  # Produção (sem capital): y = a + n
  y[0]  = a[0] + n[0]

  # Restrição de recursos: c = y
  c[0]  = y[0]

  # Custo marginal real: mc = w - a
  mc[0] = w[0] - a[0]

  # NKPC
  pi[0] = beta * pi[1] + kappa * mc[0]

  # Regra de política monetária
  i[0]  = phi_pi * pi[0] + phi_y * y[0] + s_i[0]

  # Bloco “Cochrane” — dívida/termo
  rho_b * v[0] = v[-1] + (omega * q[0] - q[-1]) - pi[0] - s[0]

  # Precificação do long-bond (q, forward)
  omega * q[1] - q[0] = i[0] + psi * std_q * eps_q[x]

  # Processos AR(1)
  s[0]   = rho_b * s[-1]   + std_b * eps_b[x]
  s_i[0] = rho_i * s_i[-1] + std_i * eps_i[x]
  a[0]   = rho_a * a[-1]   + std_a * eps_a[x]
end;


###### Parâmetros do modelo ##############

@parameters NK_COCHRANE_FULL begin
  # Preferências / tecnologia / precificação
  beta   = 0.99
  sigma  = 1.0      # 1/IES no seu comentário original; aqui entra como multiplicador na IS
  varphi = 2.0      # inversa da Frisch
  kappa  = 0.1

  # Regra de política e persistências
  phi_pi = 0.8
  phi_y  = 0.0
  rho_b  = 0.9
  rho_i  = 0.8
  omega  = 0.9
  psi    = 0.0
  rho_a  = 0.95

  # Desvios-padrão dos choques (Dynare: stderr = 1)
  std_b  = 1.0
  std_i  = 1.0
  std_q  = 1.0
  std_a  = 1.0
end;



# Resolve o modelo

solve!(NK_COCHRANE_FULL)

# Plota IRFs 

plot_irf(NK_COCHRANE_FULL)

