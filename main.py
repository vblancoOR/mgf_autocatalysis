# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 13:04:53 2024

@author: Trabajador
"""


import funciones
import scenario_generator as generator


def main():
    
    
    # Obsolete creationg of scenarios
    # # -------------------------------
    # numberSpecies = 10
    # numberReactions = 10
    # version = 69
    # createScenario(numberSpecies, numberReactions, version)
    # # -------------------------------
    
    
    # nameScenario = "s4_r4_v0"
    # nameScenario = "s4_r4_v1"
    # nameScenario = "s4_r4_v2"
    # nameScenario = "s4_r4_v3"
    # nameScenario = "s4_r4_v4"
    # nameScenario = "s4_r4_v69" # no autocatalytic
    # nameScenario = "s5_r5_v0"
    # nameScenario = "s6_r6_v0"
    # nameScenario = "s7_r7_v0"
    # nameScenario = "s8_r8_v0"
    # nameScenario = "s9_r9_v0"
    # nameScenario = "s10_r10_v0"
    nameScenario = "s10_r10_v1"
    # nameScenario = "s10_r10_v69" # no autocatalytic
    # nameScenario = "s20_r20_v0" # no autocatalytic
    # nameScenario = "formose" # no autocatalytic
    # nameScenario = "praful_1" # no autocatalytic
    # nameScenario = "praful_2" # no autocatalytic
    # nameScenario = "e_coli"

    number_periods = 3
    # funciones.tryGrowthRateGraphRecordData(nameScenario)
    # funciones.tryGrowthRateGraph(nameScenario)
    print('----------------------------------')
    # funciones.tryGrowthRateGraphAutocatalyticRecordData(nameScenario)
    # funciones.tryGrowthRateGraphAutocatalytic(nameScenario)
    print('----------------------------------')
    funciones.tryGrowthRateInSubGraphRecordData(nameScenario)
    funciones.tryGrowthRateInSubGraphRecordDataAlt1(nameScenario)
    funciones.tryGrowthRateInSubGraphRecordDataAlt2(nameScenario)

    # funciones.tryGrowthRateInSubGraph(nameScenario)
    print('----------------------------------')
    # funciones.tryGrowthRateInSubGraphWithTime(nameScenario, number_periods)
    print('----------------------------------')
    # funciones.tryGrowthRateInSubGraphFoodWaste(nameScenario, number_periods)
    print('----------------------------------')

    # Create scenarios
    # ---------------------------------------------------------------------
    # for numero_filas in [10, 25, 50]:
    #     for numero_columnas in [10, 25, 50]:
    #         for factor_de_densidad in [3, 4, 5]:
    #             for version in [0, 1, 2, 3, 4]:
    #                 for valor_maximo in [10, 100, 500]:
    #                     print(numero_filas, numero_columnas, factor_de_densidad, version, valor_maximo)
    #                     generator.generadorEscenario(numero_filas, numero_columnas, factor_de_densidad, version, valor_maximo)
    # ---------------------------------------------------------------------


    # Run experiments
    # ---------------------------------------------------------------------
    # for numero_filas in [10, 25, 50]:
    #     for numero_columnas in [10, 25, 50]:
    #         for factor_de_densidad in [3, 4, 5]:
    #             for valor_maximo in [10, 100, 500]:
    #                 for version in [0, 1, 2, 3, 4]:
    #                     nameScenario = "n" + str(numero_filas) + "m" + str(numero_columnas) + "d" + str(factor_de_densidad) + "max" + str(valor_maximo) + "v" + str(version)
    #                     funciones.tryGrowthRateGraphRecordData(nameScenario)
    #                     funciones.tryGrowthRateGraphAutocatalyticRecordData(nameScenario)
    #                     funciones.tryGrowthRateInSubGraphRecordData(nameScenario)
    # ---------------------------------------------------------------------



if __name__ == "__main__":
    main()







