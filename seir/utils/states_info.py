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

import sys # Importa modulos do sistema
import os # Local (path) do sistema
import pathlib  # Local do diretório (pasta)
import csv
from urllib import request
import json
from seir.utils.folders import cd
import urllib.error
import pandas as pd


# DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'utils')
# OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'output')

# DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '', 'utils')
# print(f'importado! Módulo: {__name__}\tPacote: {__package__}')

# Retorna o código dos estados brasileiros selecionado por sigla.
def states_codes(state):
    states_cod = {# Estado : Cód
                "AC" : 12,
                "AL" : 27,
                "AM" : 13,
                "AP" : 16,
                "BA" : 29,
                "CE" : 23,
                "ES" : 32,
                "GO" : 52,
                "MA" : 21,
                "MG" : 31,
                "MS" : 50,
                "MT" : 51,
                "PA" : 15,
                "PB" : 25,
                "PE" : 26,
                "PI" : 22,
                "PR" : 41,
                "RJ" : 33,
                "RN" : 24,
                "RO" : 11,
                "RR" : 14,
                "RS" : 43,
                "SC" : 42,
                "SE" : 28,
                "SP" : 35,
                "TO" : 17,
                "DF" : 53,
                }
    return states_cod.get(state, '** Estado inválido **')

# Retorna o número de municípios por estado selecionado por sigla.
def states_mun(state_mun = "BR"):
    states_mun = {
                "AC" : 22,
                "AL" : 102,
                "AM" : 62,
                "AP" : 16,
                "BA" : 417,
                "CE" : 184,
                "ES" : 78,
                "GO" : 247,
                "MA" : 217,
                "MG" : 853,
                "MS" : 79,
                "MT" : 141,
                "PA" : 144,
                "PB" : 223,
                "PE" : 184,
                "PI" : 224,
                "PR" : 399,
                "RJ" : 92,
                "RN" : 167,
                "RO" : 52,
                "RR" : 15,
                "RS" : 497,
                "SC" : 295,
                "SE" : 75,
                "SP" : 645,
                "TO" : 139,
                "DF" : 1,
                "BR" : 5569,
                }
    return states_mun.get(state_mun, '** Estado inválido **')

# Retorna a população (em 2018) de cada estado selecionado por sigla.
def states_pop(state_pop):
    states_pop = {# Estado : Pop    : Cód
                "AC" : 869265,      # 12
                "AL" : 3322820,     # 27
                "AM" : 4080610,     # 13
                "AP" : 829494,      # 16
                "BA" : 14812600,    # 29
                "CE" : 9075650,     # 23
                "ES" : 3972390,     # 32
                "GO" : 9895860,     # 52
                "MA" : 7035060,     # 21
                "MG" : 21040700,    # 31
                "MS" : 2748020,     # 50
                "MT" : 3442000,     # 51
                "PA" : 8513500,     # 15
                "PB" : 3996500,     # 25
                "PE" : 9493270,     # 26
                "PI" : 3264530,     # 22
                "PR" : 11348900,    # 41
                "RJ" : 17160000,    # 33
                "RN" : 3479010,     # 24
                "RO" : 1757590,     # 11
                "RR" : 576568,      # 14
                "RS" : 11329600,    # 43
                "SC" : 7075490,     # 42
                "SE" : 2278310,     # 28
                "SP" : 45538900,    # 35
                "TO" : 1555230,     # 17
                "DF" : 3015268,     # 53
                "TOTAL" : 210147125,     #
                }
    return states_pop.get(state_pop, '** Estado inválido **')




# Leitos UTI dedicados ao COVID19
def states_leitos_UTI(estado):
    st_leitos_UTI = {# Estado : Pop    : Cód
                "AC" : 74 ,
                "AL" : 269 ,
                "AM" : 242 ,
                "AP" : 93 ,
                "BA" : 1268 ,
                "CE" : 831 ,
                "ES" : 608 ,
                "GO" : 547 ,
                "MA" : 393 ,
                "MG" : 2168 ,
                "MS" : 252 ,
                "MT" : 383 ,
                "PA" : 486 ,
                "PB" : 249 ,
                "PE" : 992 ,
                "PI" : 335 ,
                "PR" : 886 ,
                "RJ" : 1946 ,
                "RN" : 443 ,
                "RO" : 158 ,
                "RR" : 25 ,
                "RS" : 917 ,
                "SC" : 822 ,
                "SE" : 169 ,
                "SP" : 5185 ,
                "TO" : 91 ,
                "DF" : 426 ,
                "TOTAL" : 20258 ,
                }
    return st_leitos_UTI.get(estado, '** Estado inválido **')


def states_novos_leitos_UTI_covid(estado):
    st_novos_leitos_UTI_covid = {# Estado : Pop    : Cód
                "AC" : 0 ,
                "AL" : 0 ,
                "AM" : 0 ,
                "AP" : 0 ,
                "BA" : 0 ,
                "CE" : 0 ,
                "ES" : 0 ,
                "GO" : 0 ,
                "MA" : 0 ,
                "MG" : 0 ,
                "MS" : 0 ,
                "MT" : 0 ,
                "PA" : 0 ,
                "PB" : 0 ,
                "PE" : 0 ,
                "PI" : 0 ,
                "PR" : 0 ,
                "RJ" : 0 ,
                "RN" : 0 ,
                "RO" : 0 ,
                "RR" : 0 ,
                "RS" : 0 ,
                "SC" : 0 ,
                "SE" : 0 ,
                "SP" : 0 ,
                "TO" : 0 ,
                "DF" : 0 ,
                "TOTAL" : 0 ,
                }
    return st_novos_leitos_UTI_covid.get(estado, '** Estado inválido **')


# Leitos clínicos por estado
def states_leitos_Gerais(estado):
    st_leitos_Gerais = {# Estado : Pop    : Cód
                "AC" : 738 ,
                "AL" : 2639 ,
                "AM" : 2518 ,
                "AP" : 657 ,
                "BA" : 11573 ,
                "CE" : 7682 ,
                "ES" : 3382 ,
                "GO" : 6528 ,
                "MA" : 6509 ,
                "MG" : 20496 ,
                "MS" : 2265 ,
                "MT" : 3132 ,
                "PA" : 6366 ,
                "PB" : 3566 ,
                "PE" : 10622 ,
                "PI" : 3219 ,
                "PR" : 11229 ,
                "RJ" : 15143 ,
                "RN" : 3111 ,
                "RO" : 2280 ,
                "RR" : 996 ,
                "RS" : 13910 ,
                "SC" : 6775 ,
                "SE" : 1466 ,
                "SP" : 36741 ,
                "TO" : 1022 ,
                "DF" : 2906 ,
                "TOTAL" : 187471 ,
                }
    return st_leitos_Gerais.get(estado, '** Estado inválido **')


def states_novos_leitos_Gerais_covid(estado):
    st_novos_leitos_Gerais_covid = {# Estado : Pop    : Cód
                "AC" : 0 ,
                "AL" : 0 ,
                "AM" : 0 ,
                "AP" : 0 ,
                "BA" : 0 ,
                "CE" : 0 ,
                "ES" : 0 ,
                "GO" : 0 ,
                "MA" : 0 ,
                "MG" : 0 ,
                "MS" : 0 ,
                "MT" : 0 ,
                "PA" : 0 ,
                "PB" : 0 ,
                "PE" : 0 ,
                "PI" : 0 ,
                "PR" : 0 ,
                "RJ" : 0 ,
                "RN" : 0 ,
                "RO" : 0 ,
                "RR" : 0 ,
                "RS" : 0 ,
                "SC" : 0 ,
                "SE" : 0 ,
                "SP" : 0 ,
                "TO" : 0 ,
                "DF" : 0 ,
                "TOTAL" : 0 ,
                }
    return st_novos_leitos_Gerais_covid.get(estado, '** Estado inválido **')


# https://www.medrxiv.org/content/10.1101/2020.04.25.20077396v1.full.pdf
def taxa_internacao_hospitalar(sigla):
    st_tx_intern_hosp = {# Estado : Pop    : Cód
                "AC"   :   0.0559,
                "AL"   :   0.0802,
                "AM"   :   0.0729,
                "AP"   :   0.0606,
                "BA"   :   0.0681,
                "CE"   :   0.0920,
                "DF"   :   0.0633,
                "ES"   :   0.0705,
                "GO"   :   0.0695,
                "MA"   :   0.0804,
                "MG"   :   0.0774,
                "MS"   :   0.0581,
                "MT"   :   0.0556,
                "PA"   :   0.0842,
                "PB"   :   0.0715,
                "PE"   :   0.1117,
                "PI"   :   0.0704,
                "PR"   :   0.0761,
                "RJ"   :   0.0970,
                "RN"   :   0.0715,
                "RO"   :   0.0604,
                "RR"   :   0.0503,
                "RS"   :   0.0740,
                "SC"   :   0.0653,
                "SE"   :   0.0554,
                "SP"   :   0.0882,
                "TO"   :   0.0545,
                }
    return st_tx_intern_hosp.get(sigla, '** Estado inválido **')


def states_name(sigla):
    states_name = {# Estado : nome
                'RO' : 'Rondônia',
                'AC' : 'Acre',
                'AM' : 'Amazonas',
                'RR' : 'Roraima',
                'PA' : 'Pará',
                'AP' : 'Amapá',
                'TO' : 'Tocantins',
                'MA' : 'Maranhão',
                'PI' : 'Piauí',
                'CE' : 'Ceará',
                'RN' : 'Rio Grande do Norte',
                'PB' : 'Paraíba',
                'PE' : 'Pernambuco',
                'AL' : 'Alagoas',
                'SE' : 'Sergipe',
                'BA' : 'Bahia',
                'MG' : 'Minas Gerais',
                'ES' : 'Espírito Santo',
                'RJ' : 'Rio de Janeiro',
                'SP' : 'São Paulo',
                'PR' : 'Paraná',
                'SC' : 'Santa Catarina',
                'RS' : 'Rio Grande do Sul',
                'MS' : 'Mato Grosso do Sul',
                'MT' : 'Mato Grosso',
                'GO' : 'Goiás',
                'DF' : 'Distrito Federal',
                # "TOTAL" : 'Brasil Escala Nacional',
                }
    return states_name.get(sigla, '** Estado inválido **')


def states_sigla(nome):
    states_sigla = {# Estado : sigla
			 'Rondônia': 'RO' ,
			 'Acre': 'AC' ,
			 'Amazonas': 'AM' ,
			 'Roraima': 'RR' ,
			 'Pará': 'PA' ,
			 'Amapá': 'AP' ,
			 'Tocantins': 'TO' ,
			 'Maranhão': 'MA' ,
			 'Piauí': 'PI' ,
			 'Ceará': 'CE' ,
			 'Rio Grande do Norte': 'RN' ,
			 'Paraíba': 'PB' ,
			 'Pernambuco': 'PE' ,
			 'Alagoas': 'AL' ,
			 'Sergipe': 'SE' ,
			 'Bahia': 'BA' ,
			 'Minas Gerais': 'MG' ,
			 'Espírito Santo': 'ES' ,
			 'Rio de Janeiro': 'RJ' ,
			 'São Paulo': 'SP' ,
			 'Paraná': 'PR' ,
			 'Santa Catarina': 'SC' ,
			 'Rio Grande do Sul': 'RS' ,
			 'Mato Grosso do Sul': 'MS' ,
			 'Mato Grosso': 'MT' ,
			 'Goiás': 'GO' ,
			 'Distrito Federal': 'DF' ,
                }
    return states_sigla.get(nome, '** Estado inválido **')
