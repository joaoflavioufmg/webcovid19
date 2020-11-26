##############################################################################
# Webcovid - Pandemia COVID-19 (para o Brasil - estados e municipios)
# Nescon e DEP (UFMG)
# Equipe: João Flávio de Freitas Almeida <joao.flavio@dep.ufmg.br>
#         Francisco Cardoso (Chico) <cardoso@nescon.medicina.ufmg.br>
#         Luiz Ricardo Pinto <luiz@dep.ufmg.br>
#         Samuel Vieira Conceição <svieira@dep.ufmg.br>
#         Virginia Magalhães <vmagalhaes@nescon.medicina.ufmg.br>
# Fonte:
# Dados:  Wesley Cota (https://covid19br.wcota.me/)
# SEIR:   Ryan McGee(https://github.com/ryansmcgee/seirsplus)
# Leitos: Array Advisors (https://www.healthleadersmedia.com/
# welcome-ad?toURL=/covid-19/see-when-states-will-face-hospital
# -bed-capacity-shortages-during-covid-19-outbreak)
##############################################################################
##############################################################################
# Estados e Municípios
##############################################################################
########################################################################
# Simulação por Estado ou Município
########################################################################
# Define se simula estado ou municpipio
is_state = True

# Função deve definir o Nome do 'estado'
estado = 'Minas Gerais'

# Função deve definir o Nome do 'Município'
municpipio = 'Belo Horizonte'

Ndias = 1825
# Mudança de parâmetros ao longo da simulação: políticas de distanciamento
ini_dist = 20   # Data de início do distanciamento
fim_dist = 40  # Data do término do distanciamento

Taxa_sobrev_LG = 1.0
Taxa_obitos_LG = 1.0 - Taxa_sobrev_LG
Taxa_sobrev_UTI = 0.5
Taxa_obitos_UTI = 1.0 - Taxa_sobrev_UTI
Taxa_ocup_LG_exp  = 0.0 # Leitos eletivos (Gerais)

# Tempos
Tempo_LG_sobrev = 8 # Dias
Tempo_LG_obito  = 12 # Dias
Tempo_UTI_sobrev_inicia_LG  = 1 # Dias
Tempo_UTI_sobrev_durante  = 7 # Dias
Tempo_UTI_sobrev_fim_em_LG = 6 # Dias
Tempo_UTI_sobrev = Tempo_UTI_sobrev_inicia_LG + Tempo_UTI_sobrev_durante + \
                   Tempo_UTI_sobrev_fim_em_LG # Dias
Tempo_UTI_obito_inicia_LG  = 7 # Dias
Tempo_UTI_obito_fim_em_LG = 8 # Dias
Tempo_UTI_obito = Tempo_UTI_obito_inicia_LG + Tempo_UTI_obito_fim_em_LG
