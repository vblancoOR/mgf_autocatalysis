# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 01:03:40 2024

@author: Trabajador
"""

import algorithm_3 as algo
import auxiliary_functions as aux


def main():

    nameScenario = "formose"
    input_matrix, output_matrix = aux.readScenario(nameScenario)
    time_limit_iteration = 1000
    formose_species = ["C1a fomaldehyd", "C2a", "C2b", "C3a", "C3b", 
                        "C3c dihydroxy acetone", "C4a", "C4b", "C4c", "C5a", "C5b",
                        "C5c", "C5d", "C5e", "C6a", "C6b", "C6c", "C6d", "C6e",
                        "C7a", "C7b", "C7c", "C7d", "C7e", "C7f", "C8a", "C8b",
                        "C8c", "C8d"]
    algo.growthRateInSubgraphDefinitive(input_matrix, output_matrix, nameScenario, time_limit_iteration, formose_species)

if __name__ == "__main__":
    main()










