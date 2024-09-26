from gurobipy import Model, GRB, LinExpr
import gurobipy as gp
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix
import time
import os

def decomp_model(A, sizeA, con, nonzeros, vars, instance, nBlocks):
    print("Create a Row-Net Hypergraph.")
    ## Create Row-Net Hypergraph
    rNetHg = []
    for i in range(sizeA[1]):                       # sizeA[0]: number of constraints
        rNetHg.append(A.getcol(i).nonzero()[0])     # for each constraint, add the index of its non-zero element into rNetHg
    nNodes = sizeA[0]+round(len(nonzeros)*0.2)
    nHedges = sizeA[1]
    wrtStr = str(nHedges)+'\t'+str(nNodes)+'\n'     # Convert the number of hyperedges nHedges and the number of nodes nNodes into strings, separated by tabs and newlines
    for i in range(nHedges):
        for k in range(len(rNetHg[i])):
            wrtStr += str(rNetHg[i][k]+1)
            if k < len(rNetHg[i])-1:                # add a tab after the node index if it is not the last node of the hyperedge
                wrtStr += '\t'
        wrtStr += '\n'
    # 'wrtStr' stores the information of the hypergraph
    ## each row represents a hyperedge
    ## each hyperedge is followed by the index of the node it contains; 
    ## Hyperedges are separated by newlines
        
    print(os.getcwd())
    os.chdir('/Users/aaron/Decomposition/Tests/readMPS')
    f = open('HGraphFiles/'+instance+'rHG_Benders','w')
    f.write(wrtStr)
    f.close


    # Use hmetis to decompose the constraint matrix
    os.chdir('/Users/aaron/Decomposition/hmetis-1.5-osx-i686')

    open('../Tests/readMPS/HGraphFiles/'+instance+'rHG_Benders')
    print('shmetis ../Tests/readMPS/HGraphFiles/'+instance+'rHG_Benders '+str(nBlocks)+' 1')
    print("Current directory:", os.getcwd())
    #os.system('/bin/bash shmetis ../Tests/readMPS/HGraphFiles/'+instance+'rHG_Benders '+str(nBlocks)+' 1')
    # Use the HMetis tool to decompose the hypergraph
    if not os.path.exists('HGraphFiles/'+instance+'rHG_Benders.part.'+str(nBlocks)): 
        # If the decomposition file does not exist, call the shmetis command to decompose the hypergraph, and store the result in a file under the HGraphFiles directory
        os.system('./shmetis ../Tests/readMPS/HGraphFiles/'+instance+'rHG_Benders '+str(nBlocks)+' 1')
    else:
        print("Hypergraph exist!!")

    # Read and plot the reordered constraint matrix
    # os.chdir('../')
    os.chdir('../Tests/readMPS/')
    print(os.getcwd())
    f = open('HGraphFiles/'+instance+'rHG_Benders.part.'+str(nBlocks))
    xx = f.readlines() # Read the above decomposed file content into the list 'xx'
    f.close

    # rowgroup: re-group the constraints corresponding to hypergraph partition
    rowgroup = {}
    colToGroup = {}
    for i in range(nBlocks+1):
        rowgroup[i] = [] # each group has an empty set
        # rowgroup[i] contains all constraints in block i

    for i in range(sizeA[0]):
        # xx[i] is the group number of constraint i, i.e., xx[i] = 1,2,...,nBlocks
        rowgroup[int(xx[i].strip('\n'))].append(i) # 'xx' contains the group info from hypergraph partition
        colToGroup[i] = int(xx[i].strip('\n'))     # constriant i belongs to which group

    constmap = {}
    ind = 0
    for i in range(nBlocks):
        for k in rowgroup[i]:
            constmap[k] = ind
            ind += 1

    colgroup = {}
    for i in range(nBlocks+1):
        colgroup[i] = []

    for i in range(sizeA[1]):
        groupind = colToGroup[A.getcol(i).nonzero()[0][0]]
        if all([colToGroup[A.getcol(i).nonzero()[0][k]] == groupind for k in range(len(A.getcol(i).nonzero()[0]))]):
            # if all variables in this row belong to the same group (groupind), then this constraint belongs to this group
            colgroup[groupind].append(i)
        else:
            # otherwise, it is a linking constraint 
            colgroup[0].append(i)

    varmap = {}
    ind = 0
    for i in range(nBlocks):
        for k in colgroup[i]:
            varmap[k] = ind
            ind += 1
            
    A_coo = A.tocoo(copy=True) # COO(Coordinate Format) is a format for sparse matrix --> (row, col, data); copy = True/False
    A_row = A_coo.row
    A_col = A_coo.col
    A_data = A_coo.data
    A_reord_col = [varmap[i] for i in A_col]
    A_reord_row = [constmap[i] for i in A_row]
    # A_reord is only for plotting, its corresponding data is incorrect
    A_reord = csr_matrix((A_data,(A_reord_row,A_reord_col))) # create a Compressed Sparse Row format (CSR) sparse matrix
    # plot
    plt.spy(A_reord,markersize=0.5)
    plt.show()
    os.chdir('../../')

    # Print the names of the variables
    #for v in vars:
        #print(v.varName)

    constraint_group = [x.strip() for x in xx]

    counter = {}

    for element in constraint_group:
        counter[element] = counter.get(element, 0) + 1

    element_counts = list(counter.items())
    element_counts.sort(key=lambda x: x[1])

    lowest_count = element_counts[0][1]
    highest_count = element_counts[-1][1]

    #get the row need to add slack variable
    add_slack = {}
    num_linking_cons = 0
    for i in range(nBlocks):
        for row in rowgroup[i]:
            add_slack[row] = True

    for row in rowgroup[nBlocks]:
        add_slack[row] = False
        num_linking_cons += 1
    
    return add_slack, num_linking_cons, lowest_count, highest_count
    



































