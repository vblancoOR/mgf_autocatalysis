# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 13:04:53 2024

@author: Trabajador
"""
import scenario_generator as generator
import codeModel as code
import matplotlib.pyplot as plt
import numpy as np
import drawing
import auxiliar as aux


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
    
    # =========================================================================
    def printSolGrowth(output_matrix, input_matrix, x, row_mapping):

        # Parameters input
        # ---------------------------
        # Number Species (int)
        number_species = output_matrix.shape[0]
        # Number Reactions (int)
        number_reactions = output_matrix.shape[1]
        # Species (list)
        species = range(number_species)
        # Reactions (list)
        reactions = range(number_reactions)
        # --------------------------------------

        for j in reactions:
            for i in species:
                if input_matrix[i,j] > 0.5:
                    coef = input_matrix[i,j] if input_matrix[i,j] > 1 else ''
                    sg = '+' if sum(ii for ii in species if ii > i and input_matrix[ii, j] > 0) else ''
                    print("%ss%d %s"%(coef, row_mapping[i] + 1, sg), end = " ")
            print ("->", end=" ")
            for i in species:
                if output_matrix[i, j] > 0.5:
                    coef = output_matrix[i, j] if output_matrix[i, j] > 1 else ''
                    sg = '+' if sum(ii for ii in species if ii > i and output_matrix[ii, j] > 0) else ''
                    print("%ss%d %s"%(coef, row_mapping[i] + 1, sg), end = " ")
            print("[%f] "%(x[j]))
    # =========================================================================
    
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
        printSolGrowth(output_modified, input_modified, x, row_mapping)
        print("Growth factor:", alpha)
        # alphaList = list(alphaDict.values())
        # plt.scatter(range(len(alphaList)), alphaList)
        generator.saveMatrices(input_modified, output_modified, "xxx")
    else:
        # Run algorithm
        x, alpha, step, alphaDict, time = code.growthRateGraph(output_matrix, input_matrix, maxIt)
        print("Number of steps: ", step)
        print("Total time: ", time)
        print("S dimension: " + str(input_matrix.shape))
        print("Autonomy S: ", aux.checkAutonomy(input_matrix, output_matrix)[2])
        print("Reacciones S: ")
        row_mapping = dict(zip(range(input_matrix.shape[0]), range(input_matrix.shape[0])))
        printSolGrowth(output_matrix, input_matrix, x, row_mapping)
        print("Growth factor:", alpha)            # alphaList = list(alphaDict.values())
        # plt.scatter(range(len(alphaList)), alphaList)
        # plt.show()
        # Define mapping
# =============================================================================


# =============================================================================
def tryGrowthRateGraphAutocatalytic(nameScenario):
    
    # =========================================================================
    def printSolGrowth(output_matrix, input_matrix, x, row_mapping):

        # Parameters input
        # ---------------------------
        # Number Species (int)
        number_species = output_matrix.shape[0]
        # Number Reactions (int)
        number_reactions = output_matrix.shape[1]
        # Species (list)
        species = range(number_species)
        # Reactions (list)
        reactions = range(number_reactions)
        # --------------------------------------

        for j in reactions:
            for i in species:
                if input_matrix[i,j] > 0.5:
                    coef = input_matrix[i,j] if input_matrix[i,j] > 1 else ''
                    sg = '+' if sum(ii for ii in species if ii > i and input_matrix[ii, j] > 0) else ''
                    print("%ss%d %s"%(coef, row_mapping[i] + 1, sg), end = " ")
            print ("->", end=" ")
            for i in species:
                if output_matrix[i, j] > 0.5:
                    coef = output_matrix[i, j] if output_matrix[i, j] > 1 else ''
                    sg = '+' if sum(ii for ii in species if ii > i and output_matrix[ii, j] > 0) else ''
                    print("%ss%d %s"%(coef, row_mapping[i] + 1, sg), end = " ")
            print("[%f] "%(x[j]))
    # =========================================================================
    
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
        printSolGrowth(output_modified, input_modified, x, row_mapping)
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
        printSolGrowth(output_matrix, input_matrix, x, row_mapping)
        print("Growth factor:", alpha)
# =============================================================================


# =============================================================================
def tryGrowthRateInSubGraph(nameScenario):
    
    # =========================================================================
    def printSolGrowth(output_matrix, input_matrix, x, row_mapping, column_mapping):

        # Parameters input
        # ---------------------------
        # Number Species (int)
        number_species = output_matrix.shape[0]
        # Number Reactions (int)
        number_reactions = output_matrix.shape[1]
        # Species (list)
        species = range(number_species)
        # Reactions (list)
        reactions = range(number_reactions)
        # --------------------------------------

        for j in reactions:
            for i in species:
                if input_matrix[i,j] > 0.5:
                    coef = input_matrix[i,j] if input_matrix[i,j] > 1 else ''
                    sg = '+' if sum(ii for ii in species if ii > i and input_matrix[ii, j] > 0) else ''
                    print("%ss%d %s"%(coef, row_mapping[i] + 1, sg), end = " ")
            print ("->", end=" ")
            for i in row_mapping.keys():
                if output_matrix[i, j] > 0.5:
                    coef = output_matrix[i, j] if output_matrix[i, j] > 1 else ''
                    sg = '+' if sum(ii for ii in species if ii > i and output_matrix[ii, j] > 0) else ''
                    print("%ss%d %s"%(coef, row_mapping[i] + 1, sg), end = " ")
            print("[%f] "%(x[column_mapping[j]]))
    # =========================================================================

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
    print("Autonomy S_result: ", aux.checkAutonomy(SSm, SSp))
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
    printSolGrowth(SSp, SSm, xx, row_mapping, column_mapping)
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
        print("Autonoma:", aux.checkAutonomy(SSp, SSm))
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
        print("\t Autonoma (t=%d):"%t, generator.checkAutonomy(SSp, SSm))
        print("\t S:(t=%d):"%t, "\n", SSp - SSm)
# =============================================================================




def main():
    
    
    
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
    # nameScenario = "s10_r10_v1"
    # nameScenario = "s10_r10_v69" # no autocatalytic
    # nameScenario = "s20_r20_v0" # no autocatalytic
    # nameScenario = "formose" # no autocatalytic
    # nameScenario = "praful_1" # no autocatalytic
    # nameScenario = "praful_2" # no autocatalytic
    nameScenario = "e_coli"

    number_periods = 4
    
    tryGrowthRateGraph(nameScenario)
    print('----------------------------------')
    # tryGrowthRateGraphAutocatalytic(nameScenario)
    print('----------------------------------')
    # tryGrowthRateInSubGraph(nameScenario)
    print('----------------------------------')
    # tryGrowthRateInSubGraphWithTime(nameScenario, number_periods)
    print('----------------------------------')
    # tryGrowthRateInSubGraphFoodWaste(nameScenario, number_periods)
    print('----------------------------------')






if __name__ == "__main__":
    main()







