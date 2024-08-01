# -*- coding: utf-8 -*-
"""
Created on Sat May 25 13:59:40 2024

@author: Trabajador
"""

import numpy as np
import gurobipy as gb


# =============================================================================
def growthRateGraph(output_matrix, input_matrix, max_steps):

    # Parameters input
    # ---------------------------
    # Stoichiometric Matrix
    stoichiometric_matrix = output_matrix - input_matrix
    # Number Species (int)
    number_species = stoichiometric_matrix.shape[0]
    # Number Reactions (int)
    number_reactions = stoichiometric_matrix.shape[1]
    # Species (list)
    species = range(number_species)
    # Reactions (list)
    reactions = range(number_reactions)
    # Alpha_0 (float)
    x_0 = np.ones(number_reactions)
    # alpha_0 = np.min([sum(output_matrix[s, r] * x_0[r] 
    #                      for r in reactions)
    #                  /
    #                  sum(input_matrix[s, r] * x_0[r] 
    #                      for r in reactions) 
    #                  for s in species])
    vector = []
    for s in species:
        numerador = sum(output_matrix[s, r] * x_0[r] 
                        for r in reactions)
        denominador = sum(input_matrix[s, r] * x_0[r] 
                          for r in reactions)
        if denominador == 0.0:
            denominador = 0.00000000000001
        ayuda = numerador/denominador
        vector.append(ayuda)
        
    alpha_0 = np.min(vector)
            

    # --------------------------------------

    # =========================================================================
    def modelGrowthRateFixed(previous_alpha):
    
        # Initialize model
        # --------------------------------------
        m = gb.Model("Growth_Rate_Model")
        # --------------------------------------
    
        # Variables
        # --------------------------------------
        # Flows
        x = m.addVars(number_reactions,
                      lb = 1,
                      ub = 100000,
                      name = "x")
        # Rate
        alpha = m.addVar(name = "alpha")
        # --------------------------------------
    
        # Objective function
        # --------------------------------------
        m.setObjective(alpha, gb.GRB.MAXIMIZE)
        # --------------------------------------
    
        # Constraints
        # --------------------------------------
        # 
        print(number_species)
        print([r for r in reactions])
        print([s for s in species])
        print(output_matrix.shape)
        
        m.addConstrs(
            (alpha <= gb.quicksum(output_matrix[s, r] * x[r] 
                                  for r in reactions) 
                     - previous_alpha * gb.quicksum(input_matrix[s, r] * x[r] 
                                            for r in reactions)
                     for s in species), "x")
        #
        m.addConstrs(
            (gb.quicksum(input_matrix[s, r] * x[r] 
                         for r in reactions) 
             >= 1 
             for s in species), "y")
        # --------------------------------------
    
        # Gurobi parameters
        # --------------------------------------
        # m.Params.LogToConsole = 0
        # m.Params.MIPGap = gap
        # m.Params.TimeLimit = timelimit
        m.Params.OutputFlag = 0
        # --------------------------------------
    
        # Run model
        # --------------------------------------
        m.optimize()
        # --------------------------------------
    
        # Result
        # --------------------------------------
        # print('xxxxxxxxxxxxxx')
        # print(m.status)
        # m.Params.DualReductions = 0
        # m.optimize()
        # print(m.status)

        # all_vars = m.getVars()
        # values = m.getAttr("X", all_vars)
        # names = m.getAttr("VarName", all_vars)
        
        # for name, val in zip(names, values):
        #     print(f"{name} = {val}")

        # print('xxxxxxxxxxxxxx')
        
        xsol = np.array([x[r].x for r in reactions])
        alphasol = alpha.x
        # --------------------------------------
    
        return xsol, alphasol
    # =========================================================================

    # Initialize algorithm
    # --------------------------------------
    stop = False
    step = 0
    alphaDict = {0: alpha_0}
    previous_alpha = alpha_0
    alpha_t = alpha_0
    # --------------------------------------


    while stop == False:
        # Solve model
        x_t, alphabar = modelGrowthRateFixed(previous_alpha)
        # In case alphabar too little or number of steps larger than maximum,
        # stop
        if (np.abs(alphabar) < 0.00000001 or
            step > max_steps):
            stop = True
            alpha_t = np.min([sum(output_matrix[s, r] * x_t[r]
                                 for r in reactions)
                             /
                             sum(input_matrix[s, r] * x_t[r]
                                 for r in reactions) 
                             for s in species])
            alphaDict[step] = alpha_t
            return x_t, alpha_t, step, alphaDict
        # Otherwise iterate
        else:
            alpha_t = np.min([sum(output_matrix[s, r] * x_t[r] 
                                 for r in reactions)
                             /
                             sum(input_matrix[s, r] * x_t[r] 
                                 for r in reactions) 
                             for s in species])
            alphaDict[step] = alpha_t
            # Initialize next step
            step += 1
            previous_alpha = alpha_t

# =============================================================================






# =============================================================================
def growthRateinSubgraph(output_matrix, input_matrix, t_max):
    
    # Parameters
    # ---------------------------
    # Stoichiometric Matrix
    stoichiometric_matrix = output_matrix - input_matrix
    # Number Species (int)
    number_species = stoichiometric_matrix.shape[0]
    # Number Reactions (int)
    number_reactions = stoichiometric_matrix.shape[1]
    # Species (list)
    species = range(number_species)
    # Reactions (list)
    reactions = range(number_reactions)
    # Alpha_0 (float)
    x_0 = np.ones(number_reactions)
    alpha_0 = np.min([sum(output_matrix[s, r] * x_0[r] 
                         for r in reactions)
                     /
                     sum(input_matrix[s, r] * x_0[r] 
                         for r in reactions) 
                     for s in species])
    # --------------------------------------
    
    # =========================================================================
    def modelGrowthRateFixed(alpha0):
        
        # Parameters
        # ---------------------------
        #
        bigM1 = 1000
        # 
        bigM2 = 10000
        # --------------------------------------
    
        # Initialize model
        # -------------------------------------
        m = gb.Model("GR_Model_SN")
        # -------------------------------------
    
        # Variables
        # --------------------------------------
        # Flows
        x = m.addVars(number_reactions,
                          ub = 100000, 
                          name = "x") 
        # Rate
        alpha = m.addVar(name = "alpha") 
        #
        y = m.addVars(number_species,
                          vtype = gb.GRB.BINARY,
                          name = "y")
        #
        a = m.addVars(number_species,
                          vtype = gb.GRB.BINARY,
                          name = "a")
        #
        z = m.addVars(number_reactions,
                          vtype = gb.GRB.BINARY,
                          name = "z")
        # --------------------------------------
    
        # Objective function
        # --------------------------------------      
        m.setObjective(alpha, gb.GRB.MAXIMIZE)
        # --------------------------------------      
    
        # Constraints
        # --------------------------------------
        #
        m.addConstrs(
            (
            alpha <= gb.quicksum(output_matrix[s, r] * x[r] 
                                 for r in reactions)
                    - alpha0 * gb.quicksum(input_matrix[s, r] * x[r] 
                                           for r in reactions) 
                    for s in species),
            name = "name1")
        
        #
        m.addConstrs(
            (gb.quicksum(input_matrix[s, r] * x[r] 
                         for r in reactions) 
            >= 1
            for s in species),
            name = "name2")
        #
        m.addConstrs(
            (y[s] <= gb.quicksum(z[r] 
                                 for r in reactions 
                                 if (output_matrix[s, r] > 0 
                                     or input_matrix[s, r] > 0)) 
             for s in species), 
            name = "name3")
        #
        m.addConstrs(
            (a[s] <= gb.quicksum(z[r] 
                                 for r in reactions 
                                 if output_matrix[s,r] > 0) 
             for s in species), 
            name = "name4")
        #
        m.addConstrs(
            (a[s] <= gb.quicksum(z[r] 
                                 for r in reactions
                                 if input_matrix[s,r] > 0) 
             for s in species), 
            name = "name5")
        #
        m.addConstrs(
            (a[s] <= y[s] 
             for s in species), 
            name = "name6")
        #
        m.addConstrs(
            (z[r] <= gb.quicksum(a[s] for s in species 
                                 if output_matrix[s, r] > 0) 
             for r in reactions), 
            name = "name7")
        #
        m.addConstrs(
            (z[r] <= gb.quicksum(a[s] for s in species 
                                 if input_matrix[s, r] > 0) 
             for r in reactions), 
            name = "name8")
        #
        m.addConstrs(
            (x[r] <= bigM2 * z[r]
             for r in reactions), 
            name = "name9")
        #
        m.addConstrs(
            (z[r] <= x[r] 
             for r in reactions), 
            name = "name10")
        #   
        m.addConstr(gb.quicksum(a[s] 
                                for s in species) >= 2, name = "name11")
        # --------------------------------------      
    
        # Gurobi parameters
        # --------------------------------------
        m.Params.OutputFlag = 0
        # --------------------------------------
        
        # Run model
        # --------------------------------------
        m.optimize()
        # --------------------------------------
    
        # Result
        # --------------------------------------
        if m.status != gb.GRB.OPTIMAL:
            # st=0
            model.computeIIS()
            
            IISfile="inf2.ILP"
            model.write(IISfile)
                
            print("INFEASIBLE!!!!!")
            with open(IISfile) as f: 
                for line in f: 
                    print(line.strip())
    
            return [], 0, [], [], []
        else:
            xsol = np.array([x[r].x for r in reactions])
            alphasol = alpha.x
            asol = [s for s in species
                    if a[s].x > 0]
            ysol = [s for s in species 
                    if y[s].x > 0]
            zsol = [r for r in reactions 
                    if z[r].x > 0]
    
            return xsol, alphasol, asol, ysol, zsol
        # --------------------------------------
    # =========================================================================

    # x0, a0, y0, z0 = modelGrowthRateFixedCase0()
    
    # if len(x0) >= 1:
    #     alpha0 = np.min([sum(output_matrix[s, r] * x0[r] 
    #                          for r in reactions)/
    #                      sum(input_matrix[s, r] * x0[r] 
    #                          for r in reactions) 
    #                      for s in a0])
        
    stop = False
    step = 0
    alphaDict = {}
    # alphabar = 10000
    alpha = alpha_0
    while stop == False:
            

        xx, alphabar, aa, yy, zz = modelGrowthRateFixed(alpha)
        alpha = np.min([sum(output_matrix[s, r] * xx[r] 
                        for r in reactions)
                        /
                        sum(input_matrix[s, r] * xx[r] 
                        for r in reactions) 
                    for s in aa])
            
        if (len(aa) < 1 or 
            (np.abs(alphabar) < 0.00000001 or 
             step > t_max)):
            stop = True
            alphaDict[step] = alpha
            return xx, alpha, step, alphaDict, aa, yy, zz
        else:
            alphaDict[step] = alpha
            step += 1
# =============================================================================





# =============================================================================
def growthRateWithTime(output_matrix, input_matrix, t_max, number_periods):
    
    # Parameters
    # ---------------------------
    # Stoichiometric Matrix
    stoichiometric_matrix = output_matrix - input_matrix
    # Number Species (int)
    number_species = stoichiometric_matrix.shape[0]
    # Number Reactions (int)
    number_reactions = stoichiometric_matrix.shape[1]
    # Species (list)
    species = range(number_species)
    # Reactions (list)
    reactions = range(number_reactions)
    # Periods (list)
    periods = range(number_periods)
    # Alpha_0 (float)
    x_0 = np.ones(number_reactions)
    alpha_0 = np.min([sum(output_matrix[s, r] * x_0[r] 
                         for r in reactions)
                     /
                     sum(input_matrix[s, r] * x_0[r] 
                         for r in reactions) 
                     for s in species])
    # alpha_0 = [alpha_0 for i in periods]
    # --------------------------------------
    
    
    # =========================================================================
    def growthRateWithTimeFixed(alpha_0):
        
        # Parameters
        # ---------------------------
        #
        bigM1 = 100
        # 
        bigM2 = 200
        # 
        tt = periods[-1]
        # --------------------------------------
    
        # Initialize model
        # -------------------------------------
        m = gb.Model("GR_Model_Time")
        # -------------------------------------
    
        # Variables
        # --------------------------------------
        # Rate
        alpha = m.addVar(name = "alpha") 
        # Flows
        x = m.addVars(reactions, 
                      periods,
                      ub = 100000, 
                      name = "x")
        #
        z = m.addVars(reactions,
                      periods, 
                      vtype = gb.GRB.BINARY,
                      name = "z")
        #
        y = m.addVars(species,
                      periods,
                      vtype = gb.GRB.BINARY,
                      name = "y")
        # autocatalitic
        a = m.addVars(species,
                      periods,
                      vtype = gb.GRB.BINARY,
                      name = "a")
        # --------------------------------------
    
        # Objective function
        # --------------------------------------      
        m.setObjective(alpha, gb.GRB.MAXIMIZE)
        # --------------------------------------
        
        # Constraints
        # --------------------------------------
        # --------- 2
        m.addConstrs(
            (
            alpha <= gb.quicksum(output_matrix[s, r] * x[r, tt] 
                                  for r in reactions)
                    - alpha_0 * gb.quicksum(input_matrix[s, r] * x[r, tt] 
                                            for r in reactions) 
                    for s in species),
            name = "name2")
        # --------- 3
        print("*****", reactions)
        m.addConstrs(
            (gb.quicksum(input_matrix[s, r] * x[r, tt] 
                         for r in reactions) 
            >= 1
            for s in species),
            name = "name3")
        # --------- 4
        m.addConstrs(
            (y[s, t] <= gb.quicksum(z[r, t] 
                                 for r in reactions 
                                 if (output_matrix[s, r] > 0 
                                     or input_matrix[s, r] > 0)) 
             for s in species 
             for t in periods), 
            name = "name4")
        # --------- 5
        m.addConstrs(
            (a[s, t] <= gb.quicksum(z[r, t] 
                                 for r in reactions 
                                 if output_matrix[s, r] > 0) 
             for s in species
             for t in periods), 
            name = "name5")
        # --------- 6
        m.addConstrs(
            (a[s, t] <= gb.quicksum(z[r, t] 
                                 for r in reactions
                                 if input_matrix[s, r] > 0) 
             for s in species
             for t in periods), 
            name = "name6")
        # --------- 7
        m.addConstrs(
            (a[s, t] <= y[s, t] 
             for s in species
             for t in periods), 
            name = "name7")
        # --------- 8
        m.addConstrs(
            (z[r, t] <= gb.quicksum(a[s, t] for s in species 
                                 if output_matrix[s, r] > 0) 
             for r in reactions
             for t in periods), 
            name = "name8")
        # --------- 9
        m.addConstrs(
            (z[r, t] <= gb.quicksum(a[s, t] for s in species 
                                 if input_matrix[s, r] > 0) 
             for r in reactions
             for t in periods), 
            name = "name9")
        # --------- 10
        m.addConstrs(
            (x[r, t] <= bigM2 * z[r, t]
             for r in reactions
             for t in periods), 
            name = "name10")
        # --------- 11
        m.addConstrs(
            (z[r, t] <= x[r, t] 
             for r in reactions
             for t in periods), 
            name = "name11")
        # --------- 12
        m.addConstrs(
            (gb.quicksum(a[s, t] 
                          for s in species) 
              >= 2 
              for t in periods), 
            name = "name12")
        # --------- 13
        m.addConstrs((a[s, t + 1] >= a[s, t]
                      for s in species
                      for t in periods[:-1]), 
                      name = "Persistence_of_autocatalytic_species") 
        # # --------- 16
        m.addConstrs(
            (z[r, t] >= z[r, t - 1]
              for r in reactions     
              for t in periods[1:]), 
              name = "Persistence_of_reactions_1")
        # # --------- 17
        m.addConstrs(
            (gb.quicksum(z[r, t] 
                          for r in reactions) 
              - gb.quicksum(z[r, t - 1] 
                            for r in reactions) 
              >= 1
              for t in periods[1:]), 
            name = "Persistence_of_reactions_2")
        # --------------------------------------
    
        # Gurobi parameters
        # --------------------------------------
        # m.Params.MIPGap = gap
        # m.Params.TimeLimit = timelimit
        m.Params.OutputFlag = 0
        # --------------------------------------
        
        # Run model
        # --------------------------------------
        m.optimize()
        # --------------------------------------
    
        # Result
        # --------------------------------------
        if m.status != gb.GRB.OPTIMAL:
            # st=0
            m.computeIIS()
            
            IISfile="inf2.ILP"
            m.write(IISfile)
                
            print("INFEASIBLE!!!!!")
            with open(IISfile) as f: 
                for line in f: 
                    print(line.strip())
    
            return [], [], [], [], []
        else:
            xsol = {}
            for r in reactions:
                for t in periods:
                    xsol[(r, t)] = x[r, t].x 

            alphasol = alpha.x

            asol = [(s, t) 
                    for s in species
                    for t in periods
                    if a[s, t].x > 0]
            
            # imprime = [a[s, t].x 
            #         for s in species
            #         for t in periods]
            
            # print(imprime)
            
            ysol = [(s, t) 
                    for s in species
                    for t in periods
                    if y[s, t].x > 0]
            zsol = [(r, t)
                    for r in reactions
                    for t in periods
                    if z[r, t].x > 0]

            return xsol, alphasol, asol, ysol, zsol
        # --------------------------------------
    # =========================================================================

    stop = False
    step = 0
    alphaDict = {}
    # alphabar = 10000
    alpha = alpha_0
    while stop == False:
        xx, alphabar, aa, yy, zz = growthRateWithTimeFixed(alpha)
        alpha = np.min([sum(output_matrix[s, r] * xx[r, periods[-1]] 
                        for r in reactions)
                        /
                        sum(input_matrix[s, r] * xx[r, periods[-1]] 
                        for r in reactions) 
                    for s in aa])
        # print(aa)
        # print('-------------')
        # print(alphabar)
        if (len(aa) < 1 or 
            (np.abs(alphabar) < 0.00000001 or 
             step > t_max)):
            stop = True
            alphaDict[step] = alpha
            return xx, alpha, step, alphaDict, alphabar, aa, yy, zz
        else:
            alphaDict[step] = alpha
            step += 1
# =============================================================================
    

# =============================================================================
def growthRateFoodWaste(output_matrix, input_matrix, t_max, number_periods):
    
    # Parameters
    # ---------------------------
    # Stoichiometric Matrix
    stoichiometric_matrix = output_matrix - input_matrix
    # Number Species (int)
    number_species = stoichiometric_matrix.shape[0]
    # Number Reactions (int)
    number_reactions = stoichiometric_matrix.shape[1]
    # Species (list)
    species = range(number_species)
    # Reactions (list)
    reactions = range(number_reactions)
    # Periods (list)
    periods = range(number_periods)
    # Alpha_0 (float)
    x_0 = np.ones(number_reactions)
    alpha_0 = np.min([sum(output_matrix[s, r] * x_0[r] 
                         for r in reactions)
                     /
                     sum(input_matrix[s, r] * x_0[r] 
                         for r in reactions) 
                     for s in species])
    # --------------------------------------
    
    
    # =========================================================================
    def growthRateFoodWasteFixed(alpha_0):
        
        # Parameters
        # ---------------------------
        #
        bigM2 = 200
        # 
        tt = periods[-1]
        # --------------------------------------
    
        # Initialize model
        # -------------------------------------
        m = gb.Model("GR_Model_Food_Waste")
        # -------------------------------------
    
        # Variables
        # --------------------------------------
        # Rate
        alpha = m.addVar(name = "alpha") 
        # Flows
        x = m.addVars(reactions, 
                      periods,
                      ub = 100000, 
                      name = "x")
        #
        z = m.addVars(reactions,
                      periods, 
                      vtype = gb.GRB.BINARY,
                      name = "z")
        #
        y = m.addVars(species,
                      periods,
                      vtype = gb.GRB.BINARY,
                      name = "y")
        # autocatalitic
        a = m.addVars(species,
                      periods,
                      vtype = gb.GRB.BINARY,
                      name = "a")
        # food
        n = m.addVars(species,
                      periods,
                      vtype = gb.GRB.BINARY,
                      name = "n")
        # waste
        w = m.addVars(species,
                      periods,
                      vtype = gb.GRB.BINARY,
                      name = "w")
        # If food waste is transformed into autoc.
        h = m.addVars(species, periods[:-1], 
                      vtype = gb.GRB.BINARY, 
                      name = "h") 
        # --------------------------------------
    
        # Objective function
        # --------------------------------------      
        m.setObjective(alpha, gb.GRB.MAXIMIZE)
        # --------------------------------------
        
        # Constraints
        # --------------------------------------
        # --------- 2
        m.addConstrs(
            (
            alpha <= gb.quicksum(output_matrix[s, r] * x[r, tt] 
                                  for r in reactions)
                    - alpha_0 * gb.quicksum(input_matrix[s, r] * x[r, tt] 
                                            for r in reactions) 
                    for s in species),
            name = "name2")
        # --------- 3
        m.addConstrs(
            (gb.quicksum(input_matrix[s, r] * x[r, tt] 
                         for r in reactions) 
            >= 1
            for s in species),
            name = "name3")
        # --------- 4
        m.addConstrs(
            (y[s, t] <= gb.quicksum(z[r, t] 
                                 for r in reactions 
                                 if (output_matrix[s, r] > 0 
                                     or input_matrix[s, r] > 0)) 
             for s in species 
             for t in periods), 
            name = "name4")
        # --------- 5
        m.addConstrs(
            (a[s, t] <= gb.quicksum(z[r, t] 
                                 for r in reactions 
                                 if output_matrix[s, r] > 0) 
             for s in species
             for t in periods), 
            name = "name5")
        # --------- 6
        m.addConstrs(
            (a[s, t] <= gb.quicksum(z[r, t] 
                                 for r in reactions
                                 if input_matrix[s, r] > 0) 
             for s in species
             for t in periods), 
            name = "name6")
        # --------- 7
        m.addConstrs(
            (a[s, t] <= y[s, t] 
             for s in species
             for t in periods), 
            name = "name7")
        # --------- 8
        m.addConstrs(
            (z[r, t] <= gb.quicksum(a[s, t] for s in species 
                                 if output_matrix[s, r] > 0) 
             for r in reactions
             for t in periods), 
            name = "name8")
        # --------- 9
        m.addConstrs(
            (z[r, t] <= gb.quicksum(a[s, t] for s in species 
                                 if input_matrix[s, r] > 0) 
             for r in reactions
             for t in periods), 
            name = "name9")
        # --------- 10
        m.addConstrs(
            (x[r, t] <= bigM2 * z[r, t]
             for r in reactions
             for t in periods), 
            name = "name10")
        # --------- 11
        m.addConstrs(
            (z[r, t] <= x[r, t] 
             for r in reactions
             for t in periods), 
            name = "name11")
        # --------- 12
        m.addConstrs(
            (gb.quicksum(a[s, t] 
                          for s in species) 
              >= 2 
              for t in periods), 
            name = "name12")
        # --------- 13
        m.addConstrs(
            (a[s, t] + n[s, t] + w[s, t] <= 1 
              for s in species 
              for t in periods), 
              name = "Types_of_species")
        # --------- 14
        m.addConstrs((z[r, t] <= a[s, t] + n[s, t] + w[s, t]
                      for s in species
                      for r in reactions
                      for t in periods
                      if output_matrix[s][r] > 0 
                      and input_matrix[s][r] > 0), 
                        name = "Activation_of_reactions") 
        # ---------
        m.addConstrs((a[s, t + 1] >= a[s, t]
                      for s in species
                      for t in periods[:-1]), 
                      name = "Persistence_of_autocatalytic_species") 
        # # ---------
        m.addConstrs(
            (z[r, t] >= z[r, t - 1]
              for r in reactions     
              for t in periods[1:]), 
              name = "Persistence_of_reactions_1")
        # # ---------
        m.addConstrs(
            (gb.quicksum(z[r, t] 
                          for r in reactions) 
              - gb.quicksum(z[r, t - 1] 
                            for r in reactions) 
              >= 1
              for t in periods[1:]), 
            name = "Persistence_of_reactions_2")
        # --------- 18
        m.addConstrs((
            h[s, t] <= a[s, t + 1] + 1 - (n[s, tt] + w[s, tt])
            for s in species
            for r in reactions
            for t in periods[:-1]
            for tt in periods
            if tt <= t
            ),
            name = "Transformation_of_Species_1")
        # --------- 19
        m.addConstrs((h[s, t] <= gb.quicksum(n[s, tt] + w[s, tt] 
                                            for tt in periods
                                            if tt <= t)
                      for s in species
                      for r in reactions
                      for t in periods[:-1]),
                        name = "Persistence_of_reactions_2")    
        # --------- 20
        m.addConstrs(
            (h[s, t] <= 1 - a[s, t]
                      for s in species
                      for t in periods[:-1]), 
            name = "Transformation_of_Species_3")
        # --------------------------------------
        
        # Gurobi parameters
        # --------------------------------------
        # m.Params.MIPGap = gap
        # m.Params.TimeLimit = timelimit
        m.Params.OutputFlag = 0
        # --------------------------------------
        
        # Run model
        # --------------------------------------
        m.optimize()
        # --------------------------------------
    
        # Result
        # --------------------------------------
        if m.status != gb.GRB.OPTIMAL:
            # st=0
            m.computeIIS()
            
            IISfile="inf2.ILP"
            m.write(IISfile)
                
            print("INFEASIBLE!!!!!")
            with open(IISfile) as f: 
                for line in f: 
                    print(line.strip())
    
            return [], [], [], [], [], [], [], []
        else:
            xsol = {}
            for r in reactions:
                for t in periods:
                    xsol[(r, t)] = x[r, t].x 

            alphasol = alpha.x

            asol = [(s, t) 
                    for s in species
                    for t in periods
                    if a[s, t].x > 0]
            
            nsol = [(s, t) 
                    for s in species
                    for t in periods
                    if n[s, t].x > 0]
            wsol = [(s, t) 
                    for s in species
                    for t in periods
                    if w[s, t].x > 0]
            hsol = [(s, t) 
                    for s in species
                    for t in periods[:-1]
                    if w[s, t].x > 0]
            
            ysol = [(s, t) 
                    for s in species
                    for t in periods
                    if y[s, t].x > 0]
            zsol = [(r, t)
                    for r in reactions
                    for t in periods
                    if z[r, t].x > 0]

            return xsol, alphasol, asol, ysol, nsol, wsol, hsol, zsol
        # --------------------------------------
    # =========================================================================

    stop = False
    step = 0
    alphaDict = {}
    # alphabar = 10000
    alpha = alpha_0
    while stop == False:
        xx, alphabar, aa, yy, nn, ww, hh, zz = growthRateFoodWasteFixed(alpha)
        alpha = np.min([sum(output_matrix[s, r] * xx[r, periods[-1]] 
                        for r in reactions)
                        /
                        sum(input_matrix[s, r] * xx[r, periods[-1]] 
                        for r in reactions) 
                    for s in aa])
        # print(aa)
        # print('-------------')
        # print(alphabar)
        if (len(aa) < 1 or 
            (np.abs(alphabar) < 0.00000001 or 
             step > t_max)):
            stop = True
            alphaDict[step] = alpha
            return xx, alpha, step, alphaDict, alphabar, aa, yy, zz, nn, ww
        else:
            alphaDict[step] = alpha
            step += 1
# =============================================================================
    

    

