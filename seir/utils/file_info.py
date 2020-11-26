#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

# print(f'importado! Módulo: {__name__}\tPacote: {__package__}')

import sys # Importa modulos do sistema
import os # Local (path) do sistema
import pathlib  # Local do diretório (pasta)
import time     # Para informações sobre o arquivo

def file_info():
    linha = 95*'='
    print("Informações do arquivo:")
    print(linha)
    nome = sys.argv[0]
    print("Nome: %s" % nome)
    print("Tamanho: %d" % os.path.getsize(nome))
    print("Criado: %s" % time.ctime(os.path.getctime(nome)))
    print("Modificado: %s" % time.ctime(os.path.getmtime(nome)))
    print("Acessado: %s" % time.ctime(os.path.getatime(nome)))
    print(linha)


if __name__ == '__main__':
    file_info()
