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

# Diretório atual
import sys # Importa modulos do sistema
print('Vesão do python: ', sys.version) # Versao do python em uso
import pathlib  # Local do diretório (pasta)
# print('Local da pasta: ', pathlib.Path().absolute())
# Onde começa o script
# print('Local da pasta: ',pathlib.Path(__file__).parent.absolute())
import os # Local (path: caminho) do sistema
# print('Local do módulo: ', os.path)
# print('Locais mapeados: ', sys.path) # Local do pyhton e do arquivo

import numpy as np  # Estatísticas (np é um 'alias' para numpy)
import scipy.stats  # Estatísticas

###############################################################
# from pacote (pasta) import modulos (arquivos) # Customizados
from math import *
import urllib.error
import pandas as pd

from seir.utils.states_info import states_pop
from seir.utils.states_info import states_leitos_UTI
from seir.utils.states_info import states_leitos_Gerais
from seir.utils.states_info import states_novos_leitos_UTI_covid
from seir.utils.states_info import states_novos_leitos_Gerais_covid
from seir.utils.states_info import taxa_internacao_hospitalar

from intro import load_states
# from intro import read_est_cases
from intro import states_pop_cases_deaths
from intro import get_EST_data_by_day
from seir.utils.folders import cd
from seir import config

from intro import get_EST_data

from seir.seirsbr import modelo_SEIRS_plus
from seir.seirsbr import taxa_de_transmissao_beta

from seir.utils.states_info import states_codes
from seir.utils.states_info import states_sigla
from seir.utils.states_info import states_pop
from seir.utils.states_info import states_name


IMAG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),'', 'images/')
###############################################################
try:
    Estados = get_EST_data_by_day()
    # Municipios = get_MUN_data_by_day()
except Exception as e:
    raise

T_dias = len(Estados.columns)
# print('T_dias: ', T_dias)

# def modelo_infeccao(Ndias=config.n_dias, sigla=config.estado_sigla, Perc_pop_afetada=config.Perc_pop_afetada):
def modelo_infeccao(Ndias, sigla, Perc_pop_afetada, Est=Estados):

    # Numero de >>> casos do passado registrados <<<
    Estados = Est

    Ncasos_est_pa = Estados[(Estados.index == sigla)].to_numpy()[0]

    # Numero de novos casos do passado registrados
    Nnovca_est_pa=[]

    for dia in range(1,len(Ncasos_est_pa)+1):
        if dia == 1:
            n_casos_est_pa = 0
            Nnovca_est_pa.append(n_casos_est_pa)
        else:
            n_casos_est_pa = max(0,(Ncasos_est_pa[dia-1] - Ncasos_est_pa[dia-2]))
            Nnovca_est_pa.append(n_casos_est_pa)




    # Previsão de >>> casos do futuro planejados <<<
    Ncasos_dia_um = Estados[(Estados.index == sigla)].iloc[:,-1].to_numpy()[0]

    # Filtro por linha e coluna (última = -1)
    y7 = Estados[(Estados.index == sigla)].iloc[:,-1].to_numpy()[0]
    y1 = Estados[(Estados.index == sigla)].iloc[:,-7].to_numpy()[0]
    Taxa_cresc_dia = exp((log(y7/y1))/7)-1
    Pop_atingida = states_pop(sigla) * Perc_pop_afetada


    Tot_casos_previstos = []
    Prev_novos_casos = []

    # #########################################
    # Novo modelo de previsão SEIR
    # #########################################
    local = states_name(sigla)
    list_local = states_pop_cases_deaths(local)
    TResus = 0
    ini_dist = 30
    fim_dist = T_dias
    Testes = 0
    N_transmissao = round((1 - 0.42*0.74),2)
    slider_pop = Perc_pop_afetada

    print('T_dias: ', T_dias)
    print('Ndias: ', Ndias)

    Prev_novos_casos = modelo_SEIRS_plus(local=local, tdias=T_dias, ndias=Ndias,
                      list_local=list_local, T_resus=TResus,
                      id=ini_dist, fd=fim_dist, tst=Testes,
                      ntrans=N_transmissao,
                      sl_pop=slider_pop,
                      imprime=False, pop_real=True)

    # Cumulativo (passado e futuro)
    Nov_casos = Nnovca_est_pa + Prev_novos_casos
    Tot_casos = list(Ncasos_est_pa)
    for i in range(len(Tot_casos), len(Nov_casos)):
        Tot_casos.append(Tot_casos[i-1] + Nov_casos[i])

    with cd(IMAG_DIR):
        with open('chart_leitos_casos.csv', 'w') as saida:
            print("dia,Tot_casos,Nov_casos", file=saida)
            for i in range(len(Tot_casos)):
                print('{0:d},{1:d},{2:d}'.format(int(i),
                                                 int(Tot_casos[i]),
                                                int(Nov_casos[i])),
                      file=saida)
    return [Tot_casos, Nov_casos]


def modelo_admissao(Ndias,
                    sigla,
                    Perc_pop_afetada,
                    Taxa_internacao_UTI,
                    Aumento_taxa_Intern,
                    Mais_leitos_UTI_p_COVID19,
                    Mais_leitos_Gerais_p_COVID19,
                    Taxa_ocup_UTI,
                    Taxa_ocup_LG,
                    Taxa_ocup_LG_exp=config.Taxa_ocup_LG_exp,
                    Taxa_sobrev_LG=config.Taxa_sobrev_LG,
                    Taxa_obitos_LG=config.Taxa_obitos_LG,
                    Taxa_sobrev_UTI=config.Taxa_sobrev_UTI,
                    Taxa_obitos_UTI=config.Taxa_obitos_UTI,
                    Tempo_LG_sobrev=config.Tempo_LG_sobrev,
                    Tempo_LG_obito=config.Tempo_LG_obito,
                    Tempo_UTI_sobrev_inicia_LG=config.Tempo_UTI_sobrev_inicia_LG,
                    Tempo_UTI_sobrev_durante=config.Tempo_UTI_sobrev_durante,
                    Tempo_UTI_sobrev_fim_em_LG=config.Tempo_UTI_sobrev_fim_em_LG,
                    Tempo_UTI_sobrev=config.Tempo_UTI_sobrev,
                    Tempo_UTI_obito_inicia_LG=config.Tempo_UTI_obito_inicia_LG,
                    Tempo_UTI_obito_fim_em_LG=config.Tempo_UTI_obito_fim_em_LG,
                    Tempo_UTI_obito=config.Tempo_UTI_obito,
                    Est=Estados):

    Estados = Est
    sigla = sigla
    Taxa_internacao_UTI = Taxa_internacao_UTI
    Taxa_internacao_LG = (1 - Taxa_internacao_UTI)

    # Leitos
    Leitos_UTI = states_leitos_UTI(sigla)
    Novos_leitos_UTI_p_COVID19 = states_novos_leitos_UTI_covid(sigla)

    Leitos_LG = states_leitos_Gerais(sigla)
    Novos_leitos_Gerais_p_COVID19 = states_novos_leitos_Gerais_covid(sigla)

    Novas_admissoes = []

    # #############################################
    # Demanda por leitos Gerais
    # #############################################
    LG_novas_admissoes = []
    # Sobrevive
    LGS_internacoes = []
    LGS_dia_saida = []
    # Óbito
    LGO_internacoes = []
    LGO_dia_saida = []

    # #############################################
    # Demanda por leitos UTI
    # #############################################
    UTI_novas_admissoes = []
    # Sobrevive
    UTIS_internacoes_LGS = []
    UTIS_dia_transf_de_LG_para_UTI = []
    UTIS_dia_transf_de_UTI_para_LG = []
    UTIS_dia_saida_LGS = []
    # Óbito
    UTIO_novas_admissoes = []
    UTIO_dia_transf_de_LG_para_UTI = []
    UTIO_dia_saida = []

    # #############################################
    # USO DE leitos (Por dia)
    # #############################################
    UTI = []
    GERAL = []

    # #############################################
    # Leitos destinados a pacientes de COVI-19 (dia)
    # #############################################
    UTI_p_COVID19 = Leitos_UTI * (1 - Taxa_ocup_UTI) + \
        Novos_leitos_UTI_p_COVID19 + Mais_leitos_UTI_p_COVID19

    LG_p_COVID19 = Leitos_LG * (1 - Taxa_ocup_LG + Taxa_ocup_LG_exp) + \
        Novos_leitos_Gerais_p_COVID19 + Mais_leitos_Gerais_p_COVID19

    # #############################################
    # Leitos disponíveis (por dia)
    # #############################################
    UTI_disp = []
    LG__disp = []

    Nd = Ndias
    Sg = sigla
    Pp = Perc_pop_afetada

    Nov_casos = modelo_infeccao(Nd,Sg,Pp)[1]

    for dia in range(len(Nov_casos)):
        n_admissoes = Nov_casos[dia] * (taxa_internacao_hospitalar(Sg) + Aumento_taxa_Intern)
        Novas_admissoes.append(int(n_admissoes))


        # Demanda leitos Gerais - Geral: Leitos Gerais: Novas admissões
        lg_novas_admissoes = Novas_admissoes[dia] * Taxa_internacao_LG
        LG_novas_admissoes.append(lg_novas_admissoes)

        # Demanda leitos Gerais - Geral Sobrevive: Internações Leitos Gerais
        lgs_internacoes = LG_novas_admissoes[dia] * Taxa_sobrev_LG
        LGS_internacoes.append(lgs_internacoes)

        # Demanda leitos Gerais - Geral Sobrevive: Dia prov saída Leito Geral
        lgs_dia_saida = (dia+1) + Tempo_LG_sobrev - 1
        LGS_dia_saida.append(lgs_dia_saida)

        # Demanda leitos Gerais - Geral Óbito: Internações Leito Geral
        lgo_internacoes = LG_novas_admissoes[dia] * Taxa_obitos_LG
        LGO_internacoes.append(lgo_internacoes)

        # Demanda leitos Gerais - Geral Óbito: Dia prov saída Leito Geral
        lgo_dia_saida = (dia+1) + Tempo_LG_obito - 1
        LGO_dia_saida.append(lgo_dia_saida)

        # Demanda leitos UTI - UTI: UTI: Novas admissões
        uti_novas_admissoes = Novas_admissoes[dia] * Taxa_internacao_UTI
        UTI_novas_admissoes.append(uti_novas_admissoes)

        # Demanda leitos UTI - UTI Sobrevive: Internações para Leitos Gerais
        utis_internacoes_LGS = UTI_novas_admissoes[dia] * Taxa_sobrev_UTI
        UTIS_internacoes_LGS.append(utis_internacoes_LGS)

        # Demanda leitos UTI - UTI Sobrevive: Dia prov de transf para UTI
        utis_dia_transf_de_LG_para_UTI = (dia+1) + Tempo_UTI_sobrev_inicia_LG
        UTIS_dia_transf_de_LG_para_UTI.append(utis_dia_transf_de_LG_para_UTI)

        # Demanda leitos UTI - UTI Sobrevive: Dia prov de transf de volta para Geral
        utis_dia_transf_de_UTI_para_LG = UTIS_dia_transf_de_LG_para_UTI[dia] +\
                                         Tempo_UTI_sobrev_durante
        UTIS_dia_transf_de_UTI_para_LG.append(utis_dia_transf_de_UTI_para_LG)

        # Demanda leitos UTI - UTI Sobrevive: Dia prov saída Geral
        utis_dia_saida_LGS = UTIS_dia_transf_de_UTI_para_LG[dia] +\
                             Tempo_UTI_sobrev_fim_em_LG
        UTIS_dia_saida_LGS.append(utis_dia_saida_LGS)

        # Demanda leitos UTI - UTI Óbito: Admissões para Leitos Gerais
        utio_novas_admissoes = UTI_novas_admissoes[dia] * Taxa_obitos_UTI
        UTIO_novas_admissoes.append(utio_novas_admissoes)

        # Demanda leitos UTI - UTI Óbito: Dia prov transf para UTI
        utio_dia_transf_de_LG_para_UTI = (dia+1) + Tempo_UTI_obito_inicia_LG
        UTIO_dia_transf_de_LG_para_UTI.append(utio_dia_transf_de_LG_para_UTI)

        # Demanda leitos UTI - UTI Óbito: Dia prov saida da UTI
        utio_dia_saida = UTIO_dia_transf_de_LG_para_UTI[dia] +\
                         Tempo_UTI_obito_fim_em_LG
        UTIO_dia_saida.append(utio_dia_saida)

    # #############################################
    # USO DE leitos (Por dia)
    # #############################################
    for dia in range(len(Nov_casos)):
        uti = sum(UTIS_internacoes_LGS[i] \
                       for i in range(len(UTIS_internacoes_LGS)) \
                       if (UTIS_dia_transf_de_LG_para_UTI[i] <= (dia+1)) \
                       and (UTIS_dia_transf_de_UTI_para_LG[i] > (dia+1))) + \
              sum(UTIO_novas_admissoes[i] \
                       for i in range(len(UTIO_novas_admissoes)) \
                       if (UTIO_dia_transf_de_LG_para_UTI[i] <= (dia+1)) \
                       and (UTIO_dia_saida[i] >= (dia+1)))

        UTI.append(uti)

        geral = sum(LGS_internacoes[i] \
                       for i in range(len(LGS_internacoes)) \
                       if (LGS_dia_saida[i] >= (dia+1)) \
                       and (range(len(LGS_internacoes))[i]+1) <= (dia+1)) + \
                sum(LGO_internacoes[i] \
                       for i in range(len(LGO_internacoes)) \
                       if (LGO_dia_saida[i] >= (dia+1)) \
                       and (range(len(LGS_internacoes))[i]+1) <= (dia+1)) + \
                sum(UTIS_internacoes_LGS[i] \
                       for i in range(len(UTIS_internacoes_LGS)) \
                       if (UTIS_dia_transf_de_LG_para_UTI[i] > (dia+1)) \
                       and (range(len(LGS_internacoes))[i]+1) <= (dia+1)) + \
                sum(UTIS_internacoes_LGS[i] \
                       for i in range(len(UTIS_internacoes_LGS)) \
                       if (UTIS_dia_transf_de_UTI_para_LG[i] <= (dia+1)) \
                       and (UTIS_dia_saida_LGS[i] >= (dia+1))) + \
                sum(UTIO_novas_admissoes[i] \
                       for i in range(len(UTIO_novas_admissoes)) \
                       if (UTIO_dia_transf_de_LG_para_UTI[i] > (dia+1)) \
                       and (range(len(LGS_internacoes))[i]+1) <= (dia+1))

        GERAL.append(geral)

        # #############################################
        # Leitos disponíveis (por dia)
        # #############################################
        uti_disp = UTI_p_COVID19 - UTI[dia]
        lg__disp = LG_p_COVID19 - GERAL[dia]

        UTI_disp.append(uti_disp)
        LG__disp.append(lg__disp)

    with cd(IMAG_DIR):
        with open("chart_leitos_disp_por_dia.csv", 'w') as saida:
            print('dia,',
                  'leitos_disp_UTI_por_dia,',
                    'leitos_disp_LG_por_dia', file = saida)
            for i in range(len(UTI_disp)):
                    print('{0:d},{1:d},{2:d}'
                        .format(int(i),
                                int(UTI_disp[i]),
                                int(LG__disp[i])), file = saida)
    return UTI_disp, LG__disp
