# -*- coding: utf-8 -*-

import algorithm_1 as algo
import auxiliary_functions as aux

def main():

    nameScenario = "n20m40d5max100v38"
    input_matrix, output_matrix = aux.readScenario(nameScenario)
    time_limit_iteration = 1000
    algo.tryGrowthRateGraph(input_matrix, output_matrix, nameScenario, time_limit_iteration, name_species="")


if __name__ == "__main__":
    main()
    

