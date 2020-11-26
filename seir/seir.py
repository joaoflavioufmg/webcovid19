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
# O modelo SEIRS+ implementa um modelo SEIR genérico incluindo complementos
# de estrutura de população, distanciamento social, nível de testes,
# rastreamento de contatos e casos em quarentena detectados.
# O pacote possui a implementação do modelo estocástico em redes dinâmicas

# A dinâmica do SEIRS:
# SEIR é um modelo de infecção. A população é dividida em indivíduos:
#   (S) Susceptível: Essa pessoa fica ...
#   (E) Exposta    : ...ao contato com infectados, ficando...
#   (I) Infectada  : ... e depois, retorna (ou não) para o estado de ...
#   (R) Recuperada. No modelo SEIRS, essa pessoa talvez retorne ao estado ...
#   (S) Susceptível. (Isso pode ser excluído do modelo, se necessário)
##############################################################################

##############################################################################
import sys # Importa modulos do sistema
print('Vesão do python: ', sys.version) # Versao do python em uso
import os # Local (path: caminho) do sistema
import pathlib  # Local do diretório (pasta)

from pathlib import Path
print('Running SEIR...' if __name__ == '__main__' else 'Importing SEIR...',
      Path(__file__).resolve())


########################################################################
# Simulação por Estado ou Município
########################################################################
if config.is_state:
    # Função receber arg 'Estado'
    # e retorna lista: [população, casos, mortes]
    local = config.estado
    list_local = states_pop_cases_deaths(local)

else:
    # Função receber arg 'Município'
    # e retorna lista: [população, casos, mortes]
    local = config.municpipio
    list_local = mun_pop_cases_deaths(local)

Ndias = config.Ndias
# Mudança de parâmetros ao longo da simulação: políticas de distanciamento
ini_dist = config.ini_dist   # Data de início do distanciamento
fim_dist = config.fim_dist  # Data do término do distanciamento
TResus = 0.0  # Taxa de resusceptibilidade
Testes = 0.0  # Taxa de testes de detectados


if __name__ == '__main__':

#######################################################################
Modelo SEIR (basico)
#######################################################################
modelo_SEIR(local=local, ndias=Ndias, list_local=list_local,
                imprime=True, pop_real=True)

#######################################################################
Modelo SEIRS (com re-susceptibilidade) - Ondas
#######################################################################
modelo_SEIRS(local=local, ndias=Ndias, list_local=list_local,
             imprime=True, pop_real=True)

#######################################################################
Modelo SEIRS+
#######################################################################
ref_model = modelo_SEIR(local=local, ndias=Ndias, list_local=list_local,
                        pop_real=True)
modelo_SEIRS_plus(ref_model, local=local, ndias=Ndias,
                      list_local=list_local,T_resus=TResus,
                      id=ini_dist, fd=fim_dist, tst=Testes,
                      imprime=True, pop_real=True)

#######################################################################
Modelo SEIR Network
#######################################################################
modelo_SEIR_Network_base(local=local, ndias=Ndias, list_local=list_local,
                         imprime=True)


#######################################################################
Modelo SEIRS+ Network
#######################################################################
ref_model        = modelo_SEIR_Network_base(local=local, ndias=Ndias,
                                            list_local=list_local,)
ref_model_determ = modelo_SEIR(local=local, ndias=Ndias,
                               list_local=list_local,)
modelo_SEIRS_network(ref_model, ref_model_determ, local=local,
                         ndias=Ndias, list_local=list_local,
                         id=ini_dist, fd=fim_dist,
                         imprime=True)
