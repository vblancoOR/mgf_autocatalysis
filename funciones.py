# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 14:52:03 2024

@author: Trabajador
"""


import scenario_generator as generator
import algoritmos as code
import matplotlib.pyplot as plt
import numpy as np
import drawing
import auxiliar as aux






# =============================================================================
def tryGrowthRateInSubGraphRecordDataAlt2(nameScenario):

    input_matrix, output_matrix = generator.readScenario(nameScenario)
    # Check autonomy
    autonomous = aux.checkAutonomy(input_matrix, output_matrix)
    #
    max_steps = 1000
    
    # Check if autonomous
    if autonomous[2] == False:
        # Remove rows and columns
        input_modified, output_modified, null_rows, null_columns = aux.removeNullRowsAndColumns(input_matrix, output_matrix, 0.1)
        # Mapping new rows and columns
        row_mapping, column_mapping = aux.mappingRowsAndColumns(input_matrix, null_rows, null_columns)
        # Run algorithm
        xx, alpha, step, alphaDict, a_sol, zz, time = code.growthRateinSubgraphAlt2(output_modified, input_modified, max_steps)
        cadena = []
        cadena = []
        cadena.append("Number of steps:\n")
        cadena.append(str(step) + '\n')
        cadena.append("Total time:\n")
        cadena.append(str(time) + "\n")
        cadena.append("average time per iteration:\n")
        if step != 0:   
            cadena.append(str(time/step) + "\n")
        else:
            cadena.append("---\n")
        cadena.append("S dimension:\n")
        cadena.append(str(input_matrix.shape) + "\n")
        cadena.append("Autonomy S:\n")
        cadena.append(str(aux.checkAutonomy(input_matrix, output_matrix)[2]) + "\n")
        cadena.append("S_modified dimension:\n")
        cadena.append(str(input_modified.shape) + "\n")
        cadena.append("Autonomy S_modified:\n")
        cadena.append(str(aux.checkAutonomy(input_modified, output_modified)[2]) + "\n")
        cadena.append("S_result dimension:\n")
        SSp = output_modified[:, zz]
        SSm = input_modified[:, zz]
        cadena.append(str((SSp - SSm).shape) + "\n")
        cadena.append("Autonomy S_result:\n")
        cadena.append(str(aux.checkAutonomy(SSm, SSp)[2]) + "\n")
        cadena.append("Number of autocatalytic species (a):\n")
        cadena.append(str(len(a_sol)) + "\n")
        cadena.append("Autocatalytic species (a):\n")
        cadena.append(" ".join(["s" + str(row_mapping[i] + 1) for i in a_sol]) + "\n")
        cadena.append("Growth factor:\n")
        cadena.append(str(alpha) + "\n")        
        cadena.append("S_result reactions:\n")
        
        new_column_mapping = dict(zip(range(len(zz)), zz))

        cadena_ayuda = aux.recordReactions(SSp, SSm, xx, row_mapping, new_column_mapping)
        # row_mapping = dict(zip(range(SSm.shape[0]), range(SSm.shape[0])))
        # column_mapping = dict(zip(range(SSm.shape[1]), range(SSm.shape[1])))
        # cadena_ayuda = aux.recordReactions(SSp, SSm, xx, row_mapping, column_mapping)

        
        
        for i in cadena_ayuda:
            cadena.append(i)
            
        with open(nameScenario + '_algorithm_3_alt2.txt', 'w') as f:
            for line in cadena:
                f.write(line)
    else:
        # Run algorithm
        xx, alpha, step, alphaDict, a_sol, zz, time = code.growthRateinSubgraphAlt2(output_matrix, input_matrix, max_steps)
        cadena = []
        cadena = []
        cadena.append("Number of steps:\n")
        cadena.append(str(step) + '\n')
        cadena.append("Total time:\n")
        cadena.append(str(time) + "\n")
        cadena.append("average time per iteration:\n")
        if step != 0:   
            cadena.append(str(time/step) + "\n")
        else:
            cadena.append("---\n")
        cadena.append("S dimension:\n")
        cadena.append(str(input_matrix.shape) + "\n")
        cadena.append("Autonomy S:\n")
        cadena.append(str(aux.checkAutonomy(input_matrix, output_matrix)[2]) + "\n")
        cadena.append("S_modified dimension:\n")
        cadena.append("---\n")
        cadena.append("Autonomy S_modified:\n")
        cadena.append("---\n")
        SSp = output_matrix[:, zz]
        SSm = input_matrix[:, zz]
        cadena.append(str((SSp - SSm).shape) + "\n")
        cadena.append("Autonomy S_result:\n")
        cadena.append(str(aux.checkAutonomy(SSm, SSp)[2]) + "\n")
        cadena.append("Number of autocatalytic species (a):\n")
        cadena.append(str(len(a_sol)) + "\n")
        cadena.append("Autocatalytic species (a):\n")
        cadena.append(" ".join(["s" + str(i + 1) for i in a_sol]) + "\n")
        cadena.append("Growth factor:\n")
        cadena.append(str(alpha) + "\n")        
        cadena.append("S_result reactions:\n")
        
        new_column_mapping = dict(zip(range(len(zz)), zz))
        row_mapping = dict(zip(range(SSm.shape[0]), range(SSm.shape[0])))
        cadena_ayuda = aux.recordReactions(SSp, SSm, xx, row_mapping, new_column_mapping)
        for i in cadena_ayuda:
            cadena.append(i)
    
        with open(nameScenario + '_algorithm_3_alt2.txt', 'w') as f:
            for line in cadena:
                f.write(line)
# =============================================================================





















# =============================================================================
def tryGrowthRateInSubGraphRecordDataAlt1(nameScenario):

    input_matrix, output_matrix = generator.readScenario(nameScenario)
    # Check autonomy
    autonomous = aux.checkAutonomy(input_matrix, output_matrix)
    #
    max_steps = 1000
    
    
    # Check if autonomous
    if autonomous[2] == False:
        # Remove rows and columns
        input_modified, output_modified, null_rows, null_columns = aux.removeNullRowsAndColumns(input_matrix, output_matrix, 0.1)
        # Mapping new rows and columns
        row_mapping, column_mapping = aux.mappingRowsAndColumns(input_matrix, null_rows, null_columns)
        # Run algorithm
        xx, alpha, step, alphaDict, a_sol, yy, zz, time = code.growthRateinSubgraphAlt1(output_modified, input_modified, max_steps)
        cadena = []
        cadena = []
        cadena.append("Number of steps:\n")
        cadena.append(str(step) + '\n')
        cadena.append("Total time:\n")
        cadena.append(str(time) + "\n")
        cadena.append("average time per iteration:\n")
        if step != 0:   
            cadena.append(str(time/step) + "\n")
        else:
            cadena.append("---\n")
        cadena.append("S dimension:\n")
        cadena.append(str(input_matrix.shape) + "\n")
        cadena.append("Autonomy S:\n")
        cadena.append(str(aux.checkAutonomy(input_matrix, output_matrix)[2]) + "\n")
        cadena.append("S_modified dimension:\n")
        cadena.append(str(input_modified.shape) + "\n")
        cadena.append("Autonomy S_modified:\n")
        cadena.append(str(aux.checkAutonomy(input_modified, output_modified)[2]) + "\n")
        cadena.append("S_result dimension:\n")
        SSp = output_modified[yy, :][:, zz]
        SSm = input_modified[yy, :][:, zz]
        new_row_mapping = dict(zip(range(len(yy)), yy))
        new_column_mapping = dict(zip(range(len(zz)), zz))
        final_row_mapping = {}
        for key, value in new_row_mapping.items():
            final_row_mapping[key] = row_mapping[value]

        cadena.append(str((SSp - SSm).shape) + "\n")
        cadena.append("Autonomy S_result:\n")
        cadena.append(str(aux.checkAutonomy(SSm, SSp)[2]) + "\n")
        cadena.append("Number of species:\n")
        cadena.append(str(len(yy)) + "\n")
        cadena.append("Species (y):\n")
        cadena.append(" ".join(["s" + str(row_mapping[i] + 1) for i in yy]) + "\n")
        cadena.append("Number of autocatalytic species (a):\n")
        cadena.append(str(len(a_sol)) + "\n")
        cadena.append("Autocatalytic species (a):\n")
        cadena.append(" ".join(["s" + str(row_mapping[i] + 1) for i in a_sol]) + "\n")
        cadena.append("Growth factor:\n")
        cadena.append(str(alpha) + "\n")        
        cadena.append("S_result reactions:\n")
        

        cadena_ayuda = aux.recordReactions(SSp, SSm, xx, final_row_mapping, new_column_mapping)
        for i in cadena_ayuda:
            cadena.append(i)
            
        with open(nameScenario + '_algorithm_3_alt1.txt', 'w') as f:
            for line in cadena:
                f.write(line)
    else:
        # Run algorithm
        xx, alpha, step, alphaDict, a_sol, yy, zz, time = code.growthRateinSubgraphAlt1(output_matrix, input_matrix, max_steps)
        cadena = []
        cadena = []
        cadena.append("Number of steps:\n")
        cadena.append(str(step) + '\n')
        cadena.append("Total time:\n")
        cadena.append(str(time) + "\n")
        cadena.append("average time per iteration:\n")
        if step != 0:   
            cadena.append(str(time/step) + "\n")
        else:
            cadena.append("---\n")
        cadena.append("S dimension:\n")
        cadena.append(str(input_matrix.shape) + "\n")
        cadena.append("Autonomy S:\n")
        cadena.append(str(aux.checkAutonomy(input_matrix, output_matrix)[2]) + "\n")
        cadena.append("S_modified dimension:\n")
        cadena.append("---\n")
        cadena.append("Autonomy S_modified:\n")
        cadena.append("---\n")
        SSp = output_matrix[yy, :][:, zz]
        SSm = input_matrix[yy, :][:, zz]
        cadena.append(str((SSp - SSm).shape) + "\n")
        cadena.append("Autonomy S_result:\n")
        cadena.append(str(aux.checkAutonomy(SSm, SSp)[2]) + "\n")
        cadena.append("Number of species:\n")
        cadena.append(str(len(yy)) + "\n")
        cadena.append("Species (y):\n")
        cadena.append(" ".join(["s" + str(i + 1) for i in yy]) + "\n")
        cadena.append("Number of autocatalytic species (a):\n")
        cadena.append(str(len(a_sol)) + "\n")
        cadena.append("Autocatalytic species (a):\n")
        cadena.append(" ".join(["s" + str(i + 1) for i in a_sol]) + "\n")
        cadena.append("Growth factor:\n")
        cadena.append(str(alpha) + "\n")        
        cadena.append("S_result reactions:\n")
        
        new_row_mapping = dict(zip(range(len(yy)), yy))
        new_column_mapping = dict(zip(range(len(zz)), zz))
        cadena_ayuda = aux.recordReactions(SSp, SSm, xx, new_row_mapping, new_column_mapping)
        for i in cadena_ayuda:
            cadena.append(i)
    
        with open(nameScenario + '_algorithm_3_alt1.txt', 'w') as f:
            for line in cadena:
                f.write(line)
# =============================================================================

































# =============================================================================
def tryGrowthRateInSubGraphRecordData(nameScenario):

    input_matrix, output_matrix = generator.readScenario(nameScenario)

    stoichiometric_matrix = output_matrix - input_matrix

    ### Construct Subnetwork with maximum Growth Factor:
    xx, alpha, step, alphaDict, a_sol, yy, zz, time = code.growthRateinSubgraph(output_matrix, input_matrix, 10000)

    cadena = []
    cadena.append("Number of steps:\n")
    cadena.append(str(step) + '\n')
    cadena.append("Total time:\n")
    cadena.append(str(time) + "\n")
    cadena.append("average time per iteration:\n")
    if step != 0:   
        cadena.append(str(time/step) + "\n")
    else:
        cadena.append("---\n")
    cadena.append("S dimension:\n")
    cadena.append(str(input_matrix.shape) + "\n")
    cadena.append("Autonomy S:\n")
    cadena.append(str(aux.checkAutonomy(input_matrix, output_matrix)[2]) + "\n")
    cadena.append("S_result dimension:\n")
    SSp = output_matrix[yy, :][:, zz]
    SSm = input_matrix[yy, :][:, zz]
    cadena.append(str((SSp - SSm).shape) + "\n")
    cadena.append("Autonomy S_result:\n")
    cadena.append(str(aux.checkAutonomy(SSm, SSp)[2]) + "\n")
    cadena.append("Number of species:\n")
    cadena.append(str(len(yy)) + "\n")
    cadena.append("Species (y):\n")
    cadena.append(" ".join(["s" + str(i + 1) for i in yy]) + "\n")
    cadena.append("Number of autocatalytic species (a):\n")
    cadena.append(str(len(a_sol)) + "\n")
    cadena.append("Autocatalytic species (a):\n")
    cadena.append(" ".join(["s" + str(i + 1) for i in a_sol]) + "\n")
    cadena.append("Growth factor:\n")
    cadena.append(str(alpha) + "\n")        
    cadena.append("S_result reactions:\n")
    number_species = stoichiometric_matrix.shape[0]
    number_reactions = stoichiometric_matrix.shape[1]
    especiesNoA = [s for s in range(number_species) if s not in yy]
    reaccionesNoA = [r for r in range(number_reactions) if r not in zz]
    row_mapping, column_mapping = aux.mappingRowsAndColumns(input_matrix, especiesNoA, reaccionesNoA)
    cadena_ayuda = aux.recordReactions(SSp, SSm, xx, row_mapping, column_mapping)
    for i in cadena_ayuda:
        cadena.append(i)
        
    with open(nameScenario + '_algorithm_3.txt', 'w') as f:
        for line in cadena:
            f.write(line)
# =============================================================================

# =============================================================================
def tryGrowthRateGraphAutocatalyticRecordData(nameScenario):

    # Read scenario
    input_matrix, output_matrix = generator.readScenario(nameScenario)
    # Check autonomy
    autonomous = aux.checkAutonomy(input_matrix, output_matrix)

    max_steps = 1000
    
    num_a = input_matrix.shape[0]
    
    # Check if autonomous
    if autonomous[2] == False:
        # Remove rows and columns
        input_modified, output_modified, null_rows, null_columns = aux.removeNullRowsAndColumns(input_matrix, output_matrix, 0.1)
        # Mapping new rows and columns
        row_mapping, column_mapping = aux.mappingRowsAndColumns(input_matrix, null_rows, null_columns)
        # Run algorithm
        x, alpha, step, alphaDict, time, a_sol = code.growthRateGraphAutocatalytic(output_modified, input_modified, max_steps, num_a)
        cadena = []
        cadena.append("Number of steps:\n")
        cadena.append(str(step) + '\n')
        cadena.append("Total time:\n")
        cadena.append(str(time) + "\n")
        cadena.append("average time per iteration:\n")
        if step != 0:   
            cadena.append(str(time/step) + "\n")
        else:
            cadena.append("---\n")
        cadena.append("S dimension:\n")
        cadena.append(str(input_matrix.shape) + "\n")
        cadena.append("Autonomy S:\n")
        cadena.append(str(aux.checkAutonomy(input_matrix, output_matrix)[2]) + "\n")
        cadena.append("S_modified dimension:\n")
        cadena.append(str(input_modified.shape) + "\n")
        cadena.append("Autonomy S_modified:\n")
        cadena.append(str(aux.checkAutonomy(input_modified, output_modified)[2]) + "\n")
        cadena.append("Number of autocatalytic species (a):\n")
        cadena.append(str(len(a_sol)) + "\n")
        cadena.append("Autocatalytic species (a):\n")
        cadena.append(" ".join(["s" + str(i + 1) for i in a_sol]) + "\n")
        cadena.append("Growth factor:\n")
        cadena.append(str(alpha) + "\n")        
        cadena.append("S_modified reactions:\n")
        column_mapping_ayuda = dict(zip(range(input_matrix.shape[1]), range(input_matrix.shape[1])))
        cadena_ayuda = aux.recordReactions(output_modified, input_modified, x, row_mapping, column_mapping_ayuda)
        for i in cadena_ayuda:
            cadena.append(i)
    else:
        # Run algorithm
        x, alpha, step, alphaDict, time, a_sol = code.growthRateGraphAutocatalytic(output_matrix, input_matrix, max_steps, num_a)
        cadena = []
        cadena.append("Number of steps:\n")
        cadena.append(str(step) + '\n')
        cadena.append("Total time:\n")
        cadena.append(str(time) + "\n")
        cadena.append("average time per iteration:\n")
        if step != 0:   
            cadena.append(str(time/step) + "\n")
        else:
            cadena.append("---\n")
        cadena.append("S dimension:\n")
        cadena.append(str(input_matrix.shape) + "\n")
        cadena.append("Autonomy S:\n")
        cadena.append(str(aux.checkAutonomy(input_matrix, output_matrix)[2]) + "\n")
        cadena.append("S_modified dimension:\n")
        cadena.append("---\n")
        cadena.append("Autonomy S_modified:\n")
        cadena.append("---\n")
        cadena.append("Number of autocatalytic species (a):\n")
        cadena.append(str(len(a_sol)) + "\n")
        cadena.append("Autocatalytic species (a):\n")
        cadena.append(" ".join(["s" + str(i + 1) for i in a_sol]) + "\n")
        cadena.append("Growth factor:\n")
        cadena.append(str(alpha) + "\n")
        row_mapping = dict(zip(range(input_matrix.shape[0]), range(input_matrix.shape[0])))
        column_mapping = dict(zip(range(input_matrix.shape[1]), range(input_matrix.shape[1])))
        cadena_ayuda = aux.recordReactions(output_matrix, input_matrix, x, row_mapping, column_mapping)
        cadena.append("S reactions:\n")
        for i in cadena_ayuda:
            cadena.append(i)

    with open(nameScenario + '_algorithm_2.txt', 'w') as f:
        for line in cadena:
            f.write(line)
# =============================================================================


# =============================================================================
def tryGrowthRateGraphRecordData(nameScenario):
     
    # Read scenario
    input_matrix, output_matrix = generator.readScenario(nameScenario)
    # Check autonomy
    autonomous = aux.checkAutonomy(input_matrix, output_matrix)
    maxIt = 1000
    
    # Check if autonomous
    if autonomous[2] == False:
        # Remove rows and columns
        input_modified, output_modified, null_rows, null_columns = aux.removeNullRowsAndColumns(input_matrix, output_matrix, 0.1)
        # Mapping new rows and columns
        row_mapping, column_mapping = aux.mappingRowsAndColumns(input_matrix, null_rows, null_columns)
        # Run algorithm
        x, alpha, step, alphaDict, time = code.growthRateGraph(output_modified, input_modified, maxIt)
        cadena = []
        cadena.append("Number of steps:\n")
        cadena.append(str(step) + '\n')
        cadena.append("Total time:\n")
        cadena.append(str(time) + "\n")
        cadena.append("average time per iteration:\n")
        if step != 0:   
            cadena.append(str(time/step) + "\n")
        else:
            cadena.append("---\n")
        cadena.append("S dimension:\n")
        cadena.append(str(input_matrix.shape) + "\n")
        cadena.append("Autonomy S:\n")
        cadena.append(str(aux.checkAutonomy(input_matrix, output_matrix)[2]) + "\n")
        cadena.append("S_modified dimension:\n")
        cadena.append(str(input_modified.shape) + "\n")
        cadena.append("Autonomy S_modified:\n")
        cadena.append(str(aux.checkAutonomy(input_modified, output_modified)[2]) + "\n")
        cadena.append("Growth factor:\n")
        cadena.append(str(alpha) + "\n")        
        cadena.append("S_modified reactions:\n")
        column_mapping_ayuda = dict(zip(range(input_matrix.shape[1]), range(input_matrix.shape[1])))
        cadena_ayuda = aux.recordReactions(output_modified, input_modified, x, row_mapping, column_mapping_ayuda)
        for i in cadena_ayuda:
            cadena.append(i)
    else:
        # Run algorithm
        x, alpha, step, alphaDict, time = code.growthRateGraph(output_matrix, input_matrix, maxIt)
        cadena = []
        cadena.append("Number of steps:\n")
        cadena.append(str(step) + '\n')
        cadena.append("Total time:\n")
        cadena.append(str(time) + "\n")
        cadena.append("average time per iteration:\n")
        if step != 0:   
            cadena.append(str(time/step) + "\n")
        else:
            cadena.append("---\n")
        cadena.append("S dimension:\n")
        cadena.append(str(input_matrix.shape) + "\n")
        cadena.append("Autonomy S:\n")
        cadena.append(str(aux.checkAutonomy(input_matrix, output_matrix)[2]) + "\n")
        cadena.append("S_modified dimension:\n")
        cadena.append("---\n")
        cadena.append("Autonomy S_modified:\n")
        cadena.append("---\n")
        cadena.append("Growth factor:\n")
        cadena.append(str(alpha) + "\n")
        row_mapping = dict(zip(range(input_matrix.shape[0]), range(input_matrix.shape[0])))
        column_mapping = dict(zip(range(input_matrix.shape[1]), range(input_matrix.shape[1])))
        cadena_ayuda = aux.recordReactions(output_matrix, input_matrix, x, row_mapping, column_mapping)
        cadena.append("S reactions:\n")
        for i in cadena_ayuda:
            cadena.append(i)

    with open(nameScenario + '_algorithm_1.txt', 'w') as f:
        for line in cadena:
            f.write(line)
# =============================================================================



# =============================================================================
def createScenario(numberSpecies, numberReactions, version):
    # Name scenario
    name = "s" + str(numberSpecies) + "_r" + str(numberReactions) + "_v" + str(version)
    # Generate matrizes
    inputMatrix, outputMatrix = generator.inputAndOutputMatrices(numberSpecies, numberReactions)
    # # Make sure they are autonomouses
    while generator.checkAutonomy(outputMatrix, inputMatrix) == False:
        inputMatrix, outputMatrix = generator.inputAndOutputMatrices(numberSpecies, numberReactions)
    generator.saveMatrices(inputMatrix, outputMatrix, name)
# =============================================================================


# =============================================================================
def tryGrowthRateGraph(nameScenario):
    
    # Read scenario
    input_matrix, output_matrix = generator.readScenario(nameScenario)
    # Check autonomy
    autonomous = aux.checkAutonomy(input_matrix, output_matrix)
    maxIt = 1000
    
    # Check if autonomous
    if autonomous[2] == False:
        # Remove rows and columns
        input_modified, output_modified, null_rows, null_columns = aux.removeNullRowsAndColumns(input_matrix, output_matrix, 0.1)
        # Mapping new rows and columns
        row_mapping, column_mapping = aux.mappingRowsAndColumns(input_matrix, null_rows, null_columns)
        # Run algorithm
        x, alpha, step, alphaDict, time = code.growthRateGraph(output_modified, input_modified, maxIt)
        print("Number of steps: ", step)
        print("Total time: ", time)
        print("S dimension: " + str(input_matrix.shape))
        print("Autonomy S: ", aux.checkAutonomy(input_matrix, output_matrix)[2])
        print("S_modified dimension: " + str(input_modified.shape))
        print("Autonomy S_modified: ", aux.checkAutonomy(input_modified, output_modified)[2])
        print("S_modified reactions: ")
        aux.printReactions(output_modified, input_modified, x, row_mapping, column_mapping)
        print("Growth factor:", alpha)
    else:
        # Run algorithm
        x, alpha, step, alphaDict, time = code.growthRateGraph(output_matrix, input_matrix, maxIt)
        print("Number of steps: ", step)
        print("Total time: ", time)
        print("S dimension: " + str(input_matrix.shape))
        print("Autonomy S: ", aux.checkAutonomy(input_matrix, output_matrix)[2])
        print("Reacciones S: ")
        row_mapping = dict(zip(range(input_matrix.shape[0]), range(input_matrix.shape[0])))
        column_mapping = dict(zip(range(input_matrix.shape[1]), range(input_matrix.shape[1])))
        aux.printReactions(output_matrix, input_matrix, x, row_mapping, column_mapping)
        print("Growth factor:", alpha)          
# =============================================================================


# =============================================================================
def tryGrowthRateGraphAutocatalytic(nameScenario):
    
    # Read scenario
    input_matrix, output_matrix = generator.readScenario(nameScenario)
    # Check autonomy
    autonomous = aux.checkAutonomy(input_matrix, output_matrix)

    max_steps = 1000
    
    num_a = input_matrix.shape[0]
    
    # Check if autonomous
    if autonomous[2] == False:
        # Remove rows and columns
        input_modified, output_modified, null_rows, null_columns = aux.removeNullRowsAndColumns(input_matrix, output_matrix, 0.1)
        # Mapping new rows and columns
        row_mapping, column_mapping = aux.mappingRowsAndColumns(input_matrix, null_rows, null_columns)
        # Run algorithm
        x, alpha, step, alphaDict, time, a_sol = code.growthRateGraphAutocatalytic(output_modified, input_modified, max_steps, num_a)
        print("Number of steps: ", step)
        print("Total time: ", time)
        print("S dimension: " + str(input_matrix.shape))
        print("Autonomy S: ", aux.checkAutonomy(input_matrix, output_matrix)[2])
        print("S_modified dimension: " + str(input_modified.shape))
        print("Autonomy S_modified: ", aux.checkAutonomy(input_modified, output_modified)[2])
        print("Autocatalytic species (a): ", end = "")
        [print("s" + str(i + 1), end = " ") for i in a_sol]
        print("")
        print("S_modified reactions: ")
        aux.printReactions(output_modified, input_modified, x, row_mapping, column_mapping)
        print("Growth factor:", alpha)
    else:
        # Run algorithm
        x, alpha, step, alphaDict, time, a_sol = code.growthRateGraphAutocatalytic(output_matrix, input_matrix, max_steps, num_a)
        print("Number of steps: ", step)
        print("Total time: ", time)
        print("S dimension: " + str(input_matrix.shape))
        print("Autonomy S: ", aux.checkAutonomy(input_matrix, output_matrix)[2])
        print("Autocatalytic species (a): ", end = "")
        [print("s" + str(i + 1), end = " ") for i in a_sol]
        print("")
        print("Reacciones S: ")
        row_mapping = dict(zip(range(input_matrix.shape[0]), range(input_matrix.shape[0])))
        column_mapping = dict(zip(range(input_matrix.shape[1]), range(input_matrix.shape[1])))
        aux.printReactions(output_matrix, input_matrix, x, row_mapping, column_mapping)
        print("Growth factor:", alpha)
# =============================================================================


# =============================================================================
def tryGrowthRateInSubGraph(nameScenario):

    input_matrix, output_matrix = generator.readScenario(nameScenario)

    stoichiometric_matrix = output_matrix - input_matrix

    ### Construct Subnetwork with maximum Growth Factor:
    xx, alpha, t, alphaDict, aa, yy, zz, time = code.growthRateinSubgraph(output_matrix, input_matrix, 10000)

    print("Number of steps: ", t)
    print("Total time: ", time)
    print("S dimension: " + str(stoichiometric_matrix.shape))
    print("Autonomy S: ", aux.checkAutonomy(input_matrix, output_matrix)[2])
    SSp = output_matrix[yy, :][:, zz]
    SSm = input_matrix[yy, :][:, zz]
    print("S_result dimension: " + str((SSp - SSm).shape))
    print("Autonomy S_result: ", aux.checkAutonomy(SSm, SSp)[2])
    print("Species (y): ", end = "")
    [print("s" + str(i + 1), end = " ") for i in yy]
    print("")
    print("Autocatalytic species (a): ", end = "")
    [print("s" + str(i + 1), end = " ") for i in aa]
    print("")
    print("S_result reactions: ")
    number_species = stoichiometric_matrix.shape[0]
    number_reactions = stoichiometric_matrix.shape[1]
    especiesNoA = [s for s in range(number_species) if s not in yy]
    reaccionesNoA = [r for r in range(number_reactions) if r not in zz]
    row_mapping, column_mapping = aux.mappingRowsAndColumns(input_matrix, especiesNoA, reaccionesNoA)
    aux.printReactions(SSp, SSm, xx, row_mapping, column_mapping)
    print("Growth factor:", alpha)
    generator.saveMatrices(SSm, SSp, "xxx")
# =============================================================================



# ============================================================================= 
def giveMeSpeciesAndReactios(inputMatrix, outputMatrix):
    nReactions = len(inputMatrix[0])
    inputMatrix = inputMatrix.transpose()
    outputMatrix = outputMatrix.transpose()
    nSpecies = len(inputMatrix[0])
    
    species = ["s" + str(i+1) for i in range(nSpecies)]
    
    reactions = []
    for row in range(nReactions):
        listLeft = [str(entry) + "*s" + str(idx + 1) 
                    for idx, entry in enumerate(inputMatrix[row]) 
                    if entry > 0]
        listLeft = [i[2:] if i[0:2] == "1*" else i for i in listLeft]
        stringLeft = ' + '.join(i for i in listLeft)
        listRight = [str(entry) + "*s" + str(idx + 1) 
                    for idx, entry in enumerate(outputMatrix[row]) 
                    if entry > 0]
        listRight = [i[2:] if i[0:2] == "1*" else i for i in listRight]
        stringRight = ' + '.join(i for i in listRight)
        reaction = stringLeft + " -> " + stringRight
        reactions.append(reaction)
    
    return species, reactions
# =============================================================================





# =============================================================================
def tryGrowthRateInSubGraphWithTime(nameScenario, number_periods):
    #
    input_matrix, output_matrix = generator.readScenario(nameScenario)
    #
    stoichiometric_matrix = output_matrix - input_matrix
    
    xx, alpha, step, alphaDict, alphabar, aa, yy, zz = code.growthRateWithTime(output_matrix, input_matrix, 10000, number_periods)
    
    # -- Dibujos --
    especies_sol = aux.cambiarFormatoEspecies(yy)
    species, reactions = giveMeSpeciesAndReactios(input_matrix, output_matrix)
    reacciones_sol = aux.cambiarFormatoReacciones(zz, reactions)
    plt.figure(1)
    drawing.dibujaPasos(especies_sol, reacciones_sol)
    plt.figure(2)
    drawing.dibuja(especies_sol, reacciones_sol)
    plt.figure(3)
    drawing.dibujaRedDePetriV0(especies_sol, reacciones_sol)
    plt.figure(4)
    drawing.dibujaRedDePetriV1(especies_sol, reacciones_sol)
    # -------------
    
    print("S-:")
    print(input_matrix)
    print("S+:")
    print(output_matrix)
    print("S:")
    print(stoichiometric_matrix)
    print("Alpha:", alpha)
    print("Species (y):")
    dictyy = aux.convierteDiccionarioDeTiempo(yy)
    print(dictyy)
    print("Autocatalytic species (a):")
    dictaa = aux.convierteDiccionarioDeTiempo(aa)
    print(dictaa)
    print("Reactions (z):")
    dictzz = aux.convierteDiccionarioDeTiempo(zz)
    print(dictzz)
    print("Flow (x) =", xx)

    print("t:", step)
    
    for t in dictyy.keys():
        ys = dictyy[t]
        zs = dictzz[t]
        SSp = output_matrix[ys, :][:, zs]
        SSm = input_matrix[ys, :][:, zs]
        print("Autonoma:", aux.checkAutonomy(SSp, SSm)[2])
        print("S:")
        print(SSp - SSm)
# =============================================================================


# =============================================================================
def tryGrowthRateInSubGraphFoodWaste(nameScenario, number_periods):
    
    input_matrix, output_matrix = generator.readScenario(nameScenario)

    stoichiometric_matrix = output_matrix - input_matrix
    print(stoichiometric_matrix.shape)

    xx, alpha, step, alphaDict, alphabar, aa, yy, zz, nn, ww = code.growthRateFoodWaste(output_matrix, input_matrix, 10000, number_periods)
    
    # -- Dibujos --
    especies_sol = aux.cambiarFormatoEspecies(yy)
    species, reactions = giveMeSpeciesAndReactios(input_matrix, output_matrix)
    reacciones_sol = aux.cambiarFormatoReacciones(zz, reactions)
    plt.figure(1)
    drawing.dibujaPasos(especies_sol, reacciones_sol)
    plt.figure(2)
    drawing.dibuja(especies_sol, reacciones_sol)
    plt.figure(3)
    drawing.dibujaRedDePetriV0(especies_sol, reacciones_sol)
    plt.figure(4)
    drawing.dibujaRedDePetriV1(especies_sol, reacciones_sol)
    # -------------
    print("S-:")
    print(input_matrix)
    print("S+:")
    print(output_matrix)
    print("S:")
    print(stoichiometric_matrix)
    print("Alpha:", alpha)
    print("Species (y):")
    dictyy = aux.convierteDiccionarioDeTiempoYNombre(yy, species)
    print(dictyy)
    print("Autocatalytic species (a):")
    dictaa = aux.convierteDiccionarioDeTiempoYNombre(aa, species)
    print(dictaa)
    print("Reactions (z):")
    dictzz = aux.convierteDiccionarioDeTiempoYNombre(zz, reactions)
    print(dictzz)
    print("Food (n) =")
    dictnn = aux.convierteDiccionarioDeTiempoYNombre(nn, species)
    print(dictnn)
    print("Waste (w) =")
    dictww = aux.convierteDiccionarioDeTiempoYNombre(ww, species)
    print(dictww)
    xxName = aux.modifyFlow(xx, reactions)
    print("Flow (x) =", xxName)
    print("t:", step)
    dictSimpleyy = aux.convierteDiccionarioDeTiempo(yy)
    dictSimplezz = aux.convierteDiccionarioDeTiempo(zz)
    for t in dictSimpleyy.keys():
        ys = dictSimpleyy[t]
        zs = dictSimplezz[t]
        SSp = output_matrix[ys, :][:, zs]
        SSm = input_matrix[ys, :][:, zs]
        print("\t Autonoma (t=%d):"%t, aux.checkAutonomy(SSp, SSm)[2])
        print("\t S:(t=%d):"%t, "\n", SSp - SSm)
# =============================================================================
