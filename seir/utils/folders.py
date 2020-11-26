# -*- coding: utf-8 -*-
#!/usr/bin/env python3

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
from contextlib import contextmanager

import shutil   # Endereços para apagar todos os arquivos da pasta


@contextmanager
def cd(newdir):
    # Mudar de diretório e retornar ao diretório original
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    # Diretório atual: Retorna ao diretório atual
    # print('>>> Local atual:', pathlib.Path().absolute())
    try:
        yield
    finally:
        os.chdir(prevdir)
        # Diretório atual: Retorna ao diretório atual
        # print('>>> Local atual:', pathlib.Path().absolute())


def del_previous_files(the_folder):
    folder = os.getcwd()+the_folder
    # print('folder:', folder)

    if not os.path.exists(folder):
        print(f'A pasta \"{folder}\" não existe!')
        # sys.exit(1)
        os.mkdir(folder)   # Cria diretorio
        print(f'A pasta \"{folder}\" foi criada!')
    else:
        with cd(folder):
            # Verifica se não há arquivo na pasta. Nem arquivo oculto
            if [f for f in os.listdir(folder) if not f.startswith('.')] == []:
                print('Ok. Sem arquivos na pasta {0:s}.'.format(the_folder))
            else:
                # dir = os.getcwd()
                for filename in os.listdir(folder):
                    file_path = os.path.join(folder, filename)
                    try:
                        if os.path.isfile(file_path) or os.path.islink(file_path):
                            os.unlink(file_path)
                        elif os.path.isdir(file_path):
                            shutil.rmtree(file_path)
                    except Exception as e:
                        print(f'Falha em deletar {filename}.')
                        print(f'Motivo: {e}.')
                # else:
                print('Todos os arquivos da pasta {0:s} foram eliminados.'
                      .format(str(the_folder)), end=' ')
                print("A pasta está limpa.")
