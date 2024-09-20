# -*- coding: utf-8 -*-
"""
Created on Sat May 25 13:59:40 2024

@author: Trabajador
"""

import numpy as np
import gurobipy as gb
import time

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
                      lb = 1, ub=1000,
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
        # print(number_species)
        # print([r for r in reactions])
        # print([s for s in species])
        # print(output_matrix.shape)
        
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
        if m.status != gb.GRB.OPTIMAL:
            # st=0
            # m.computeIIS()
            
            # IISfile="inf2.ILP"
            # m.write(IISfile)
                
            # print("INFEASIBLE!!!!!")
            # with open(IISfile) as f: 
            #     for line in f: 
            #         print(line.strip())
    
            return [], 0
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
    alpha_old=0
    # --------------------------------------

    start = time.time()
    while stop == False:
        # Solve model
        x_t, alphabar = modelGrowthRateFixed(previous_alpha)
        #print("[%d] alphabar: %.4f, alpha: %.8f"%(step,alphabar,alpha_t))
        if len(x_t)==0:
            return [], -1, step, alphaDict, time.time()-start
        else:
            # In case alphabar too little or number of steps larger than maximum,
            # stop
            alpha_t = np.min([sum(output_matrix[s, r] * x_t[r]
                                        for r in reactions)
                                    /
                                    sum(input_matrix[s, r] * x_t[r]
                                        for r in reactions) 
                                    for s in species])
            if (np.abs(alphabar) < 0.00000001 or
                step > max_steps or np.abs(alpha_old-alpha_t)<0.0000001):
                stop = True
                #print(np.abs(alphabar) < 0.00000001, step > max_steps, np.abs(alpha_old-alpha_t)<0.0000001)
                alphaDict[step] = alpha_t
                return x_t, alpha_t, step, alphaDict, time.time()-start
            # Otherwise iterate
            else:
                alphaDict[step] = alpha_t
                alpha_old=alpha_t
                # Initialize next step
                step += 1
                previous_alpha = alpha_t

# =============================================================================


# =============================================================================
def growthRateGraph_a(output_matrix, input_matrix, max_steps, num_a):

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
    x_0 = 10*np.ones(number_reactions)
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
    #print(vector)
    alpha_0 = np.min(vector)
    
            

    # --------------------------------------

    # =========================================================================
    def modelGrowthRateFixeda(previous_alpha):
    
        # Initialize model
        # --------------------------------------
        m = gb.Model("Growth_Rate_Model")
        # --------------------------------------
    
        # Variables
        # --------------------------------------
        # Flows
        UB=1000
        x = m.addVars(number_reactions,
                      lb = 1,
                      ub = UB,
                      name = "x")
        # Rate
        alpha = m.addVar(name = "alpha")
        a = m.addVars(number_species, vtype=gb.GRB.BINARY, name="a")
        # for s in species:
        #     if s in [0, 2, 7, 8, 9, 13, 17]:
        #         a[s].lb=1
        #     else:
        #         a[s].ub=0
        # --------------------------------------
    
        # Objective function
        # --------------------------------------
        m.setObjective(alpha, gb.GRB.MAXIMIZE)

        #m.addConstr(alpha>=0.1)
        # --------------------------------------
    
        # Constraints
        # --------------------------------------
        # 
        # print(number_species)
        # print([r for r in reactions])
        # print([s for s in species])
        # print(output_matrix.shape)
        
        m.addConstrs(
            (alpha <= gb.quicksum(output_matrix[s, r] * x[r] 
                                  for r in reactions) 
                     - a[s]*previous_alpha * gb.quicksum(input_matrix[s, r] * x[r] for r in reactions) + UB*(1-a[s])
                     for s in species), "x")
        
        #
        m.addConstrs(
            (gb.quicksum(input_matrix[s, r] * x[r] 
                         for r in reactions) 
             >= a[s] 
             for s in species), "y")
        m.addConstr(gb.quicksum(a[s] for s in species)>=2, name="suma")

        for s in species:
            if sum(output_matrix[s,r] for r in reactions)<0.5 or sum(input_matrix[s,r] for r in reactions)<0.5:
                a[s].ub=0
        for r in reactions:
            m.addConstr(gb.quicksum(a[s] for s in species if input_matrix[s,r]>0.5)>=1, name="aut_1[%d]"%r)
            m.addConstr(gb.quicksum(a[s] for s in species if output_matrix[s,r]>0.5)>=1, name="aut_2[%d]"%r)

        #m.addConstr(gb.quicksum(a[s] for s in species) <= num_a)
        # --------------------------------------
    
        # Gurobi parameters
        # --------------------------------------
        # m.Params.LogToConsole = 0
        # m.Params.MIPGap = gap
        # m.Params.TimeLimit = timelimit
        m.Params.OutputFlag = 0
        m.Params.DualReductions=0
        #m.write("graph_a.lp")
        # --------------------------------------
    
        # Run model
        # --------------------------------------
        m.optimize()
        if m.status != gb.GRB.OPTIMAL:
            print("Status: ", m.status)
            m.computeIIS()
            
            IISfile="infa.ILP"
            m.write(IISfile)
                
            print("INFEASIBLE!!!!!")
            with open(IISfile) as f: 
                for line in f: 
                    print(line.strip())
    
            return [], 0, []
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
        asol = np.array([s for s in species if a[s].x>0.5])
        # --------------------------------------
    
        return xsol, alphasol, asol
    # =========================================================================

    # Initialize algorithm
    # --------------------------------------
    stop = False
    step = 0
    alphaDict = {0: alpha_0}
    previous_alpha = alpha_0
    alpha_t = alpha_0
    #print("alpha_0: ", alpha_t)
    alpha_old=0
    # --------------------------------------

    start = time.time()
    while stop == False:
        # Solve model
        x_t, alphabar, a_t = modelGrowthRateFixeda(previous_alpha)
        #print("alphabarA: ", alphabar, "alpha: ", alpha_t, "|a|=", len(a_t))
        if len(x_t)==0:
            return [], -1, step, alphaDict, time.time()-start, []
        else:
            # In case alphabar too little or number of steps larger than maximum,
            # stop
            alpha_t = np.min([sum(output_matrix[s, r] * x_t[r] 
                                    for r in reactions)
                                /
                                sum(input_matrix[s, r] * x_t[r] 
                                    for r in reactions) 
                                for s in a_t])
            if (np.abs(alphabar) < 0.00001 or
                step > max_steps or np.abs(alpha_old-alpha_t)<0.00000001):
                stop = True
                #print("alpha_t[%d])%.4f"%(step, alpha_t))
                alphaDict[step] = alpha_t
                return x_t, alpha_t, step, alphaDict, time.time()-start, a_t
            # Otherwise iterate
            else:
                alphaDict[step] = alpha_t
                # Initialize next step
                alpha_old=alpha_t
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
    # Potential autocatalytic:
    potential_aut=[s for s in species if sum(output_matrix[s,r] for r in reactions)>0.5 and sum(input_matrix[s,r] for r in reactions)>0.5]
    # Alpha_0 (float)
    x_0 = np.ones(number_reactions)
    alpha_0 = np.min([sum(output_matrix[s, r] * x_0[r] 
                         for r in reactions)
                     /
                     sum(input_matrix[s, r] * x_0[r] 
                         for r in reactions) 
                     for s in potential_aut])
    if alpha_0<0.00001:
        x_0 = np.random.randint(1,100, size=number_reactions)
        #print(output_matrix)
        alpha_0 = np.min([sum(output_matrix[s, r] * x_0[r] 
                         for r in reactions)
                     /
                     sum(input_matrix[s, r] * x_0[r] 
                         for r in reactions) 
                     for s in potential_aut])

    # --------------------------------------
    #print("alpha0: ", alpha_0)
    # =========================================================================
    def modelGrowthRateFixed(alpha0):
        
        # Parameters
        # ---------------------------
        #
        bigM1 = 1000
        # 
        bigM2 = 1000
        # --------------------------------------
    
        # Initialize model
        # -------------------------------------
        m = gb.Model("GR_Model_SN")
        # -------------------------------------
        print("nR: ", number_reactions)
        # Variables
        # --------------------------------------
        # Flows
        x = m.addVars(number_reactions,
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
        for s in species:
            if s not in potential_aut:
                a[s].ub=0

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
                                           for r in reactions) + bigM1*(1-a[s])
                    for s in species),
            name = "name1")
        
        #
        m.addConstrs(
            (gb.quicksum(input_matrix[s, r] * x[r] 
                         for r in reactions) 
            >= a[s]
            for s in species if sum(input_matrix[s,r] for r in reactions)>=1),
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
                                 if output_matrix[s,r] > 0.5) 
             for s in potential_aut), 
            name = "name4")
        #
        m.addConstrs(
            (a[s] <= gb.quicksum(z[r] 
                                 for r in reactions
                                 if input_matrix[s,r] > 0.5) 
             for s in potential_aut), 
            name = "name5")
        #
        m.addConstrs(
            (a[s] <= y[s] 
             for s in species), 
            name = "name6")
        #
        m.addConstrs(
            (z[r] <= gb.quicksum(a[s] for s in potential_aut
                                 if output_matrix[s, r] > 0.5) 
             for r in reactions), 
            name = "name7")
        #
        m.addConstrs(
            (z[r] <= gb.quicksum(a[s] for s in potential_aut 
                                 if input_matrix[s, r] > 0.5) 
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
        m.addConstr(gb.quicksum(a[s] for s in potential_aut) >= 1, name = "name11")
        m.addConstr(gb.quicksum(z[r] for r in reactions) >= 1, name = "name12")
        # --------------------------------------      
    
        # Gurobi parameters
        # --------------------------------------
        m.Params.OutputFlag = 0
        # --------------------------------------
        m.write("model_SN.lp")
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
    alpha_old=0#gb.GRB.INFINITY
    start=time.time()
    while stop == False:
        xx, alphabar, aa, yy, zz = modelGrowthRateFixed(alpha)
        alpha = np.min([sum(output_matrix[s, r] * xx[r]
                        for r in reactions)
                        /
                        sum(input_matrix[s, r] * xx[r] 
                        for r in reactions) 
                    for s in aa])
        #print("**[%d] alphabar: "%step, alphabar, "alpha: ", alpha, "|a|=", len(aa), "|z|=", len(zz), end="") 
        # if alpha>1:
        #     print("  AUTOCATALYTIC!")
        # else:
        #     print()
        if (len(aa) < 1 or 
            (np.abs(alphabar) < 0.0001 or 
             step > t_max or np.abs(alpha-alpha_old)<0.00001)):
            stop = True
            alphaDict[step] = alpha
            return xx, alpha, step, alphaDict, aa, yy, zz, time.time()-start
        else:
            alphaDict[step] = alpha
            alpha_old=alpha
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
                     for s in species if sum(input_matrix[s, r] * x_0[r] 
                         for r in reactions)>0.0001])
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
        for s in species:
            print("sp ", s, "--",input_matrix[s,:])
        m.addConstrs(
            (gb.quicksum(input_matrix[s, r] * x[r, tt] 
                         for r in reactions) 
            >= 1
            for s in species if sum(abs(input_matrix[s,r]) for r in reactions)>0.01),
            name = "name3a")
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
            # m.computeIIS()
            
            # IISfile="inf2.ILP"
            # m.write(IISfile)
                
            # print("INFEASIBLE!!!!!")
            # with open(IISfile) as f: 
            #     for line in f: 
            #         print(line.strip())
    
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
                     for s in species if sum(input_matrix[s, r] * x_0[r] 
                         for r in reactions)>0.0001])
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
        for s in species:
            print("sp ", s, "--",input_matrix[s,:])
        m.addConstrs(
            (gb.quicksum(input_matrix[s, r] * x[r, tt] 
                         for r in reactions) 
            >= 1
            for s in species if sum(abs(input_matrix[s,r]) for r in reactions)>0.01),
            name = "name3b")
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
    

