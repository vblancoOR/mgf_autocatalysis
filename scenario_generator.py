# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 17:22:53 2024

@author: Trabajador
"""
import numpy as np
import random as rd


# https://cristianbastidas.com/my-blog/en/algorithms/python/challenge/guide/2021/06/06/parts-of-number.html
# function for getting all partitions of a number
# Example: 
# - input:3
# - output: [[3], [2, 1], [1, 1, 1]]
# =============================================================================
def listPartitions(x:int) -> list:
    # -------------------------------------------------------------------------
    def partitions(ones:list, x:int, origin:list=[]) -> list:
      total = []
      for i in range(ones.count(1), 1, -1):
        aux = ones[:ones.index(1)]
        aux.append(i)
        while sum(aux) < x:
          aux.append(1)
        if not sorted(aux) in origin:
          total.append(aux)
          origin.append(sorted(aux))
      for l in total:
        total = total + partitions(l,x,origin)
      return total
    # -------------------------------------------------------------------------
    ones = [1]*x
    parts = partitions(ones, x, [])
    parts.append(ones)
    return parts
# =============================================================================


# get list of stoichiometric vectors from list of partitions given max cluster
# size and highest reaction order
# =============================================================================
def turnIntoVectors(parts, n_max, rxn_order):
    # -------------------------------------------------------------------------
    def turnIntoVector(part, n_max):
        stoich_vec = np.zeros(n_max)
        for i in range(len(part)):
            stoich_vec[part[i]-1] += 1 
        return stoich_vec
    # -------------------------------------------------------------------------
    s_vec = []
    for part in parts:
        if len(part) <= rxn_order:
            if all(j <= n_max for j in part):
                s_vec.append(turnIntoVector(part, n_max))
    return np.array(s_vec)       
# =============================================================================


# get input-output matrix pair for the array of stoichiometric vectors
# =============================================================================
def inputAndOutputMatrices(numberSpecies, numberReactions):
    # -------------------------------------------------------------------------
    def giveMePossibleReactions(numberSpecies):
        parts = listPartitions(numberSpecies)
        rxn_order = len(parts)
        listVectors = turnIntoVectors(parts, numberSpecies, rxn_order)
        return listVectors
    # -------------------------------------------------------------------------

    listVectors = giveMePossibleReactions(numberSpecies)
    
    listPossible = []
    for i in range(len(listVectors)):
        for j in range(len(listVectors)):
            if i != j:
                listPossible.append((i, j))
                
    listSelected = rd.sample(listPossible, numberReactions)
    
    inputMatrix = []
    outputMatrix = []
    for i in listSelected:
        inputMatrix.append(listVectors[i[0]])
        outputMatrix.append(listVectors[i[1]])
        
    inputMatrix = np.array(inputMatrix)
    inputMatrix = np.transpose(inputMatrix)

    outputMatrix = np.array(outputMatrix)
    outputMatrix = np.transpose(outputMatrix)

    return inputMatrix, outputMatrix
# =============================================================================


# # =============================================================================        
# def saveMatrices(mMinus, mPlus, name):
#     f = open(name + ".txt", "w")
#     np.savetxt(f, mMinus, fmt='%d')
#     f.write("\n")
#     np.savetxt(f, mPlus, fmt='%d')
#     f.close()
# # =============================================================================
 

# =============================================================================        
def saveMatrices(mMinus, mPlus, name):
    np.savetxt(name + "_minus.txt", mMinus, fmt = '%i')
    np.savetxt(name + "_plus.txt", mPlus, fmt = '%i')
# =============================================================================


# # =============================================================================
# def checkAutonomy(output_matrix, input_matrix):

#     # Parameters input
#     # ---------------------------
#     #
#     stoichiometric_matrix = output_matrix - input_matrix
#     #
#     number_species = stoichiometric_matrix.shape[0]
#     #
#     number_reactions = stoichiometric_matrix.shape[1]
#     #
#     species = range(number_species)
#     #
#     reactions = range(number_reactions)
#     # --------------------------------------
    
#     autonomous = True

#     # Check autonomy in species
#     # --------------------------------------
#     for s in species:
#         pos = 0
#         neg = 0
#         for r in reactions:
#             if output_matrix[s, r] > 0.1:
#                 pos += 1
#             if input_matrix[s, r] > 0.1:
#                 neg += 1
#         if pos == 0 or neg == 0:
#             autonomous = False
#     # --------------------------------------
            
#     # Check autonomy in reactions
#     # --------------------------------------
#     for r in reactions:
#         pos = 0
#         neg = 0
#         for s in species:
#             if output_matrix[s, r] > 0.1:
#                 pos += 1
#             if input_matrix[s, r] > 0.1:
#                 neg += 1
#         if pos == 0 or neg == 0:
#             autonomous = False
#     # --------------------------------------

#     return autonomous 
# # =============================================================================


# =============================================================================
def readScenario(name):
    mMinus = np.loadtxt("./scenarios/" + name + '_minus.txt', dtype= int)
    mPlus = np.loadtxt("./scenarios/" + name + '_plus.txt', dtype= int)
    return mMinus, mPlus
# =============================================================================
           
            
# =============================================================================
def transformStoichiometricIntoInputAndOuput(path, name):
    
    mm = np.loadtxt(path, dtype= int)
    
    input_matrix = []
    output_matrix = []
    for i in mm:
        vector_input_matrix = []
        vector_output_matrix = []
        for j in i:
            if j > 0.5:
                vector_input_matrix.append(0)
                vector_output_matrix.append(1)
            elif j < -0.5:
                vector_input_matrix.append(1)
                vector_output_matrix.append(0)
            else:
                vector_input_matrix.append(0)
                vector_output_matrix.append(0)
        input_matrix.append(vector_input_matrix)
        output_matrix.append(vector_output_matrix)
        
    input_matrix = np.matrix(input_matrix)
    output_matrix = np.matrix(output_matrix)
    
    f = open(name + "_minus.txt", "w")
    np.savetxt(f, input_matrix, fmt='%d')
    
    q = open(name + "_plus.txt", "w")
    np.savetxt(q, output_matrix, fmt='%d')
# =============================================================================



# =============================================================================
def generadorEscenario(numero_filas, numero_columnas, factor_de_densidad, version, valor_maximo):
    # Función crea pares para matriz de inputs
    # -----------------------------------------------------------------------------
    def genPairs(numero_filas, numero_columnas, N):
        L = set()
        i, j = np.random.randint(0, numero_filas), np.random.randint(0, numero_columnas)
        L.add((i, j))
        while len(L) < N:
            while (i, j) in L:
                i, j = np.random.randint(0, numero_filas), np.random.randint(0, numero_columnas)
            L.add((i,j))
        return L
    # -----------------------------------------------------------------------------
    # Función crea pares para matriz de outputs
    # -----------------------------------------------------------------------------
    def genPairsNotInL(numero_filas, numero_columnas, N, LL):
        L = set()
        i, j = np.random.randint(0, numero_filas), np.random.randint(0, numero_columnas)
        if (i, j) not in LL:
            L.add((i, j))
        while len(L) < N:
            while (i, j) in L or (i,j) in LL:
                i, j = np.random.randint(0, numero_filas), np.random.randint(0, numero_columnas)
            L.add((i,j))
        return L
    # -----------------------------------------------------------------------------
    # Matriz aleatoria
    matrix = np.random.randint(0, valor_maximo, size=(numero_filas, numero_columnas))
    # Pares inputs
    pairs1 = genPairs(numero_filas, numero_columnas, numero_filas*numero_columnas/factor_de_densidad)
    # Pares outputs
    pairs2 = genPairsNotInL(numero_filas, numero_columnas, numero_filas*numero_columnas/factor_de_densidad, pairs1)
    # Inicialización de matrices
    output_matrix = np.zeros((numero_filas, numero_columnas))
    input_matrix = np.zeros((numero_filas, numero_columnas))
    # Rellenar las matrices
    for i in range(numero_filas):
        for j in range(numero_columnas):
            if (i, j) in pairs1:
                input_matrix[i, j] = matrix[i, j]
            elif (i, j) in pairs2:
                output_matrix[i, j] = matrix[i, j]
    # Guardar matrices
    np.savetxt("instances/left_n%dm%dd%d_%d.txt"%(numero_filas, numero_columnas, factor_de_densidad, version), input_matrix, fmt = "%d")
    np.savetxt("instances/right_n%dm%dd%d_%d.txt"%(numero_filas, numero_columnas, factor_de_densidad, version), output_matrix, fmt = "%d")
# =============================================================================



