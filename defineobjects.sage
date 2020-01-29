# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 08:14:13 2019

@author: samis
"""

    # ========================================================================
    # --- Defining some Matrcies 
    # ========================================================================
myring = GF(2)

for i in range(1):
    matrix_A = Matrix(myring,
            [[1,1,1],
             [0,1,0],
             [0,0,0],
             [0,0,0]
            ])
    matrix_B = Matrix(myring,
            [[0,0,1],
             [1,1,1],
             [1,0,1],
             [0,0,0]
            ])
        
    matrix_C = Matrix(myring,
            [[1,1,1,0],
             [0,1,1,1],
             [0,0,0,0],
             [0,0,0,0],
             [0,0,0,0]
            ])
    matrix_D = Matrix(myring,
            [[0,0],
             [1,1],
             [1,0],
             [0,0]
            ])
    matrix_E = Matrix(myring,
            [[0,0],
             [1,1],
             [0,0],
             [0,0]
            ])
    matrix_F = Matrix(myring,
            [[0,1,0,0],
             [0,1,1,1],
             [0,0,0,0],
             [1,0,0,1],
             [0,0,0,0]
            ])
    matrix_G = Matrix(myring,
            [[4,4,0],
             [1,1,0],
             [0,0,1],
             [0,0,1]
            ]) 
    #---------------------------------------------------------------
    matrix_B1 = Matrix(myring,
            [[1,0,1,1],
             [0,1,1,1],
             [0,0,0,1],
             [0,0,0,0]
            ])
    matrix_B2 = Matrix(myring,
            [[1,0,1,0],
             [0,0,0,0]
            ])
    matrix_B3 = Matrix(myring,
            [[1,1,1,1],
             [0,1,0,1],
             [0,1,0,1]
            ])
    matrix_B4 = Matrix(myring,
            [[1,0,1],
             [0,1,1]
            ])           

    #---------------------------------------------------------------
    matrix_C1 = Matrix(myring,
            [[0,0,0],
             [1,1,0],
             [1,1,1]
            ])
    matrix_C2 = Matrix(myring,
            [[1,1,1],
             [0,1,0],
             [0,1,0]
            ])
    matrix_C3 = Matrix(myring,
            [[1,0,0],
             [1,1,1],
             [1,0,1]
            ])
    matrix_C4 = Matrix(myring,
            [[0,1,0],
             [0,0,0],
             [1,1,1]
            ])
            
    #---------------------------------------------------------------
    matrix_D1 = Matrix(myring,
            [[1,0,0],
             [1,1,0],
             [1,1,0],
             [0,1,0],
             [0,0,1]
            ])
    matrix_D2 = Matrix(myring,
            [[1,1,0,0,0],
             [1,1,1,0,0],
             [0,0,1,1,0],
             [0,0,1,1,1],
             [0,0,0,0,0]
            ])
    matrix_D3 = Matrix(myring,
            [[1,0,0,0,0],
             [1,1,0,0,0]
            ])
    matrix_D4 = Matrix(myring,
            [[0,1,0],
             [0,0,0]
            ])
#---------------------------------------------------------------
# This is an example that has object-excess. 
# The kernals of the maps 0->1, 0->2, and 0->3 are different 1-dim subspaces of R^2,
# so the multiflag on object 0 has the shape of the poset { {}, {a}, {b}, {c}, {a,b,c}  }, which is NOT in general position
# This can be repeated arbitarily to get as large an object-excess as you want
    matrix_E1 = Matrix(myring,
            [[1,0],
             [0,0],
             [0,0]
            ])
    matrix_E2 = Matrix(myring,
            [[1,1],
             [1,1],
             [1,1]
            ])
    matrix_E3 = Matrix(myring,
            [[0,0],
             [0,0],
             [0,1]
            ])
    matrix_E4 = Matrix(myring,
            [[0,1,1],
             [0,1,1]
            ])
    matrix_E5 = Matrix(myring,
            [[1,-1,0],
             [-1,1,0]
            ])
    matrix_E6 = Matrix(myring,
            [[1,1,0],
             [1,1,0]
            ])
# -------------------------------------------------------------
# Now an attempt to get morphism excess without object excess.
    matrix_F1 = Matrix(myring,
            [[0,1,0],
             [0,0,1],
             [0,0,0]
            ])
    matrix_F2 = Matrix(myring,
            [[1,1,0],
             [0,0,1],
             [0,0,0]
            ])
    matrix_F3 = Matrix(myring,
            [[1,0,0],
             [0,0,1],
             [0,0,0]
            ])
    matrix_F4 = Matrix(myring,
            [[0,0,1],
             [0,0,1],
             [0,0,1]
            ])
    matrix_F5 = Matrix(myring,
            [[0,0,1],
             [0,0,1],
             [0,0,1]
            ])
    matrix_F6 = Matrix(myring,
            [[0,0,1],
             [0,0,1],
             [0,0,1]
            ])
# -------------------------------------------------------------
# A 2x2 digram with spaces of dimension 6
    matrix_G1 = Matrix(myring,
            [[0 ,1 ,1 ,0 ,1 ,0],
            [1 ,1 ,1 ,0 ,0 ,0],
            [0 ,1 ,0 ,0 ,0 ,1],
            [1 ,1 ,1 ,1 ,0 ,1],
            [1 ,0 ,1 ,0 ,0 ,1],
            [0 ,0 ,1 ,0 ,1 ,1]
            ])
    matrix_G3 = Matrix(myring,
            [[0 ,0 ,1 ,1 ,1 ,1],
            [1 ,0 ,0 ,0 ,1 ,0],
            [0 ,0 ,0 ,1 ,1 ,1],
            [1 ,1 ,1 ,1 ,0 ,1],
            [1 ,0 ,0 ,0 ,1 ,1],
            [1 ,1 ,0 ,0 ,0 ,1]
            ])
    matrix_G2 = Matrix(myring,
            [[0, 0, 1, 1, 1, 0],
             [1, 1, 0, 0, 1, 1],
             [0, 1, 1, 1, 1, 1],
             [0, 0, 0, 1, 0, 1],
             [0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0]
             ])
    matrix_G4 = Matrix(myring,
            [[1, 0, 0, 0, 0, 0],
             [0, 1, 0, 0, 0, 0],
             [0, 0, 1, 0, 0, 0],
             [0, 0, 0, 1, 0, 0],
             [1, 1, 0, 1, 0, 0],
             [0, 1, 1, 1, 0, 0]
             ])

# -------------------------------------------------------------
# Trying to get a 2x2 diagram which takes more than 2 LS2 steps to stabilize.
    array_H1 = list_transpose([[0 ,1 ,1 ,0 ,0 ,1 ,1 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,1 ,0],[0 ,1 ,1 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1],[0 ,1 ,1 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1],[0 ,1 ,1 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1],[0 ,1 ,1 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1],[1 ,0 ,0 ,1 ,1 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,1 ,0],[1 ,1 ,1 ,1 ,1 ,0 ,1 ,1 ,0 ,1 ,1 ,0 ,1 ,0 ,1 ,0],[0 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,1 ,1 ,0 ,0 ,1 ,1 ,0],[0 ,0 ,0 ,0 ,1 ,0 ,0 ,1 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,1],[0 ,1 ,1 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1],[0 ,1 ,1 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1],[0 ,1 ,1 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1],[0 ,1 ,1 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1],[0 ,1 ,1 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1],[0 ,1 ,0 ,1 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,1],[1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,1 ,0
            ]])
    shuffle(array_H1)
    array_H1 = list_transpose(array_H1)
    matrix_H1 = Matrix(myring,
            array_H1
            )
    matrix_H2 = Matrix(myring,
            [[
            1 ,0 ,0 ,1 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,1 ,1],[0 ,0 ,0 ,1 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1 ,1],[1 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,1 ,1 ,0 ,1],[1 ,0 ,1 ,0 ,1 ,1 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,1 ,1 ,0],[1 ,1 ,1 ,0 ,1 ,1 ,0 ,1 ,0 ,1 ,1 ,0 ,0 ,0 ,1 ,1],[1 ,0 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0],[0 ,1 ,0 ,0 ,1 ,1 ,0 ,1 ,1 ,1 ,1 ,0 ,0 ,1 ,0 ,1],[1 ,0 ,0 ,1 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,1 ,1],[1 ,0 ,0 ,1 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,1 ,1],[1 ,0 ,0 ,1 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,1 ,1],[1 ,0 ,0 ,1 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,1 ,1],[0 ,1 ,1 ,1 ,1 ,0 ,1 ,0 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,1 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1 ,1],[0 ,0 ,0 ,1 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1 ,1],[0 ,0 ,0 ,1 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1 ,1],[0 ,0 ,0 ,1 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,1 ,1
            ]])
    matrix_H3 = Matrix(myring,
            [[
            0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0
            ]])
    matrix_H4 = Matrix(myring,
            [[
            0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0],[0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0
            ]])

    # -------------------------------------------------------------
#    basis_A = matrix_to_columns(
#            Matrix(myring,
#             [[Rational(1,2),Rational(1,3),Rational(1,5)],
#              [0,            Rational(1,15),Rational(1,6)],
#              [0,            0,          Rational(1,10)]
#             ]))
#    basis_B = matrix_to_columns(
#            Matrix(myring,
#             [[Rational(1,2),Rational(1,3),0],
#              [0,            Rational(1,15),0],
#              [0,            0,          0]
#             ]))
#    basis_C = matrix_to_columns(
#            Matrix(myring,
#             [[0,0,0],
#              [0,1,0],
#              [0,0,1]
#             ]))
    vector_A = Matrix(myring,[1,2,3])
    vector_B = Matrix(myring,
             [[0],
              [0],
              [0],
              [1]
             ])
    #-------------------------------------------------------------------------
    moduleA_size = 6
    mapsA = [[0 for x in range(moduleA_size)] for y in range(moduleA_size)]
    mapsA[1][2] = linear_transformation(matrix_A, side = 'right')
    mapsA[1][4] = linear_transformation(matrix_B, side = 'right')
    mapsA[2][5] = linear_transformation(matrix_C, side = 'right')
    mapsA[3][2] = linear_transformation(matrix_D, side = 'right')
    mapsA[3][4] = linear_transformation(matrix_E, side = 'right')
    mapsA[4][5] = linear_transformation(matrix_F, side = 'right')
    
    moduleB_size = 4
    mapsB = [[0 for x in range(moduleB_size)] for y in range(moduleB_size)]
    mapsB[0][1] = linear_transformation(matrix_B1, side = 'right')
    mapsB[0][2] = linear_transformation(matrix_B3, side = 'right')
    mapsB[1][3] = linear_transformation(matrix_B2, side = 'right')
    mapsB[2][3] = linear_transformation(matrix_B4, side = 'right')

    moduleC_size = 5
    mapsC = [[0 for x in range(moduleC_size)] for y in range(moduleC_size)]
    mapsC[0][1] = linear_transformation(matrix_C1, side = 'right')
    mapsC[1][2] = linear_transformation(matrix_C2, side = 'right')
    mapsC[2][3] = linear_transformation(matrix_C3, side = 'right')
    mapsC[3][4] = linear_transformation(matrix_C4, side = 'right')
    
    moduleD_size = 5
    mapsD = [[0 for x in range(moduleD_size)] for y in range(moduleD_size)]
    mapsD[0][1] = linear_transformation(matrix_D1, side = 'right')
    mapsD[2][1] = linear_transformation(matrix_D2, side = 'right')
    mapsD[2][3] = linear_transformation(matrix_D3, side = 'right')
    mapsD[4][3] = linear_transformation(matrix_D4, side = 'right')
    
    moduleE_size = 5
    mapsE = [[0 for x in range(moduleE_size)] for y in range(moduleE_size)]
    mapsE[0][1] = linear_transformation(matrix_E1, side = 'right')
    mapsE[0][2] = linear_transformation(matrix_E2, side = 'right')
    mapsE[0][3] = linear_transformation(matrix_E3, side = 'right')
    mapsE[1][4] = linear_transformation(matrix_E4, side = 'right')
    mapsE[2][4] = linear_transformation(matrix_E5, side = 'right')
    mapsE[3][4] = linear_transformation(matrix_E6, side = 'right')
    
    moduleF_size = 5
    mapsF = [[0 for x in range(moduleF_size)] for y in range(moduleF_size)]
    mapsF[0][1] = linear_transformation(matrix_F1, side = 'right')
    mapsF[0][2] = linear_transformation(matrix_F2, side = 'right')
    mapsF[0][3] = linear_transformation(matrix_F3, side = 'right')
    mapsF[1][4] = linear_transformation(matrix_F4, side = 'right')
    mapsF[2][4] = linear_transformation(matrix_F5, side = 'right')
    mapsF[3][4] = linear_transformation(matrix_F6, side = 'right')
    
    moduleG_size = 4
    mapsG = [[0 for x in range(moduleG_size)] for y in range(moduleG_size)]
    mapsG[0][1] = linear_transformation(matrix_G1, side = 'right')
    mapsG[0][2] = linear_transformation(matrix_G2, side = 'right')
    mapsG[1][3] = linear_transformation(matrix_G3, side = 'right')
    mapsG[2][3] = linear_transformation(matrix_G4, side = 'right')
    
    moduleH_size = 4
    mapsH = [[0 for x in range(moduleH_size)] for y in range(moduleH_size)]
    mapsH[0][1] = linear_transformation(matrix_H1, side = 'right')
    mapsH[0][2] = linear_transformation(matrix_H2, side = 'right')
    mapsH[1][3] = linear_transformation(matrix_H3, side = 'right')
    mapsH[2][3] = linear_transformation(matrix_H4, side = 'right')

#except:
    #print "There was a problem defining some objects"
    

    
    