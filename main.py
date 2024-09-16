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

def IsAutonomous(input_matrix, output_matrix, aut_species, aut_reactions):
    
    autonomous_species=True
    autonomous_reactions=True
    #n, m =input_matrix.shape
    for s in aut_species:
        if sum(input_matrix[s,r] for r in aut_reactions)<0.5 or sum(output_matrix[s,r] for r in aut_reactions)<0.5:
            #print(s, " is not s-autonomous")
            autonomous_species=False
    for r in aut_reactions:
        if sum(input_matrix[s,r] for s in aut_species)<0.5 or sum(output_matrix[s,r] for s in aut_species)<0.5:
            #print(r, " is not r-autonomous")
            autonomous_reactions=False

    return autonomous_species, autonomous_reactions 


def CleanData(SM):

    #SM=output_m-input_m
    n, m =SM.shape
    #print(n,m)
    del_c=100
    del_r=100
    col_delete=[]
    row_delete=[]
    while del_r>0 or del_c>0:
        del_c=0
        for  r in range(m):
            if r not in col_delete and (len([s for s in range(n) if SM[s,r]>=0.5])<0.001 or len([s for s in range(n) if SM[s,r]<=-0.5])<0.001):
                col_delete.append(r)
                SM[:, r]=0
                del_c+=1
        del_r=0
        for s in range(n):
            if s not in row_delete and (len([r for r in range(m) if SM[s,r]>=0.5])<0.001 or len([r for r in range(m) if  SM[s,r]<=-0.5])<0.001):
                row_delete.append(s)
                SM[s,:]=0
                del_r+=1

    n, m =SM.shape
    input_matrix=np.zeros((n,m))
    output_matrix=np.zeros((n,m))
    for s in range(n):
        for r in range(m):
            if SM[s,r]>0.5:
                output_matrix[s,r]=SM[s,r]
            if SM[s,r]<-0.5:
                input_matrix[s,r]=-SM[s,r]
    input_matrix=np.delete(input_matrix, row_delete, 0)
    output_matrix=np.delete(output_matrix, row_delete, 0)
    input_matrix=np.delete(input_matrix, col_delete, 1)
    output_matrix=np.delete(output_matrix, col_delete, 1)
    #n, m =input_matrix.shape
    #print(n,m)

    names_rows=[s for s in range(n) if s not in row_delete]
    names_cols=[r for r in range(m) if r not in col_delete]
    # print("len_del_c: ", len(col_delete), "len_del_r: ", len(row_delete))
    # print("len_c:", len(names_cols), "len_r: ", len(names_rows))

    return input_matrix, output_matrix, names_rows, names_cols

# =============================================================================
def tryGrowthRateGraph(nameScenario):
    
    # =========================================================================
    def printSolGrowth(output_matrix, input_matrix, x):

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
                    print("%s s%d %s "%(coef, i+1, sg), end = " ")
            print ("  -> ", end=" ")
            for i in species:
                if output_matrix[i, j]>0.5:
                        coef = output_matrix[i, j] if output_matrix[i, j] > 1 else ''
                        sg = '+' if sum(ii for ii in species if ii > i and output_matrix[ii, j] > 0) else ''
                        print("%s s%d %s "%(coef, i+1, sg), end = " ")
            print("[%f] "%(x[j]))
    # =========================================================================
    
    input_matrix, output_matrix = generator.readScenario(nameScenario)

    MaxIt = 1000
    x, alpha, step, alphaDict = code.growthRateGraph(output_matrix, input_matrix, MaxIt)
    if alpha < 1:
        print("Iterations: ", step, "alpha: ", alpha)
        alphaList = list(alphaDict.values())
        plt.scatter(range(len(alphaList)), alphaList)
        plt.show()
        printSolGrowth(output_matrix, input_matrix, x)
        print(np.dot(output_matrix, x)/np.dot(input_matrix, x))
        stop = 1
# =============================================================================


# =============================================================================
def tryGrowthRateInSubGraph(nameScenario):
    
    input_matrix, output_matrix = generator.readScenario(nameScenario)

    stoichiometric_matrix = output_matrix - input_matrix

    number_species = stoichiometric_matrix.shape[0]
    
    ### Construct Subnetwork with maximum Growth Factor:

    xx, alpha, t, alphaDict, aa, yy, zz = code.growthRateinSubgraph(output_matrix, input_matrix, 10000)

    print("S:")
    print(stoichiometric_matrix)
    print("Alfa:", alpha)
    print("Especies (y):", yy)
    print("E. autocataliticas (a):", aa)
    print("Reacciones (r):", zz)
    print("Flujo (x) =", xx)
    print("t:", t)

    SSp = output_matrix[yy, :][:, zz]
    SSm = input_matrix[yy, :][:, zz]
    print("Autonoma:", generator.checkAutonomy(SSp, SSm))
    print("S_final:")
    print(SSp - SSm)
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
        print("Autonoma:", generator.checkAutonomy(SSp, SSm))
        print("S:")
        print(SSp - SSm)

# =============================================================================


# =============================================================================
def tryGrowthRateInSubGraphFoodWaste(nameScenario, number_periods):
    
    input_matrix, output_matrix = generator.readScenario(nameScenario)

    stoichiometric_matrix = output_matrix - input_matrix

    xx, alpha, step, alphaDict, alphabar, aa, yy, zz, nn, ww = code.growthRateFoodWaste(output_matrix, input_matrix, 10000, number_periods)
    
    # -- Dibujos --
    especies_sol = aux.cambiarFormatoEspecies(yy)
    species, reactions = giveMeSpeciesAndReactios(input_matrix, output_matrix)
    reacciones_sol = aux.cambiarFormatoReacciones(zz, reactions)
    # plt.figure(1)
    # drawing.dibujaPasos(especies_sol, reacciones_sol)
    # plt.figure(2)
    # drawing.dibuja(especies_sol, reacciones_sol)
    # plt.figure(3)
    # drawing.dibujaRedDePetriV0(especies_sol, reacciones_sol)
    # plt.figure(4)
    # drawing.dibujaRedDePetriV1(especies_sol, reacciones_sol)
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
    xxName = aux.modifyFlow(xx, species)
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




# def main():
    
    
    
#     # # -------------------------------
#     # numberSpecies = 10
#     # numberReactions = 10
#     # version = 69
#     # createScenario(numberSpecies, numberReactions, version)
#     # # -------------------------------
    
    
    
#     # nameScenario = "formose"
#     #nameScenario = "s10_r10_v1"
#     nameScenario = "formose"

#     number_periods = 3
    
#     # tryGrowthRateGraph(nameScenario)
#     print('----------------------------------')
#     # tryGrowthRateInSubGraph(nameScenario)
#     print('----------------------------------')
#     #tryGrowthRateInSubGraphWithTime(nameScenario, number_periods)
#     print('----------------------------------')
#     tryGrowthRateInSubGraphFoodWaste(nameScenario, number_periods)
#     print('----------------------------------')

# if __name__ == "__main__":
#     main()










