# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 08:14:13 2019

@author: samis
"""

from PIL import Image
import networkx as nx
import matplotlib.pyplot as plt
from sage.plot.arrow import Arrow
from random import shuffle

    # ========================================================================
    # --- Load some Matrcies 
    # ========================================================================

myring = GF(2)
load("defineobjects.sage")

    # ===========================================================
    #------------  More variables you have to define ----------
    # ===========================================================

#space_dimensions = [0,3,4,2,4,5]
#module_size = len(space_dimensions)

#maps = [[0 for x in range(module_size)] for y in range(module_size)]
#maps[1][2] = linear_transformation(matrix_A, side = 'right')
#maps[1][4] = linear_transformation(matrix_B, side = 'right')
#maps[2][5] = linear_transformation(matrix_C, side = 'right')
#maps[3][2] = linear_transformation(matrix_D, side = 'right')
#maps[3][4] = linear_transformation(matrix_E, side = 'right')
#maps[4][5] = linear_transformation(matrix_F, side = 'right')



    # ========================================================================
    # ---
    # ========================================================================


def space_dimensions_from_map_array(map_array):
    module_size = len(map_array)
    space_dimensions           = [-1]*module_size
    
    for i in range(module_size):
        for j in range(module_size):
            if map_array[i][j] != 0:
                morphism = map_array[i][j]
                domain_dimension = morphism.domain().rank()
                range_dimension = morphism.codomain().rank()
                
                if space_dimensions[i] == -1:
                    space_dimensions[i] = domain_dimension
                elif space_dimensions[i] != domain_dimension:
                    raise ValueError("matricies in maps are of inconsistent size. See the function space_dimensions_from_map_array")
                    
                if space_dimensions[j] == -1:
                    space_dimensions[j] = range_dimension
                elif space_dimensions[j] != range_dimension:
                    raise ValueError("matricies in maps are of inconsistent size. See the function space_dimensions_from_map_array")
    if -1 in space_dimensions:
        raise ValueError("The module size is wrong ")
        
    return space_dimensions
    
def trivial_module(space_dimensions_list, ring=ZZ):
    module_size = len(space_dimensions_list)
    module = []
    
    for i in range(module_size):
        module.append([])
        space = ring^space_dimensions_list[i]
        module[i].append(space)
        module[i].append(space.zero_subspace())
        #print(module[i])
    
    return module

    # ========================================================================
    # --- Custom printing functions
    # ========================================================================

def myprint(array):
    for i in range( len(array) ):
        show(i, '.    .', array[i])
        
def list_transpose(array):
    return [list(x) for x in zip(*array)]

def subspace_to_basis_vector_array(my_subspace):
    #print '---------------------------'
    #show('==== started subspace to basis vector using the subspace ====', my_subspace)
    #print 'printed, the subspace is', my_subspace
    mybasis=my_subspace.basis()
    #show('mybasis is', mybasis)
    vector_array=[]
    for my_vector in mybasis:
        #show( 'my_vector is', my_vector, 'e has type', type(my_vector))
        #print 'printed my vector is', my_vector, 'parent', my_vector.parent()
        #show('my_vector has terms:  \t ', my_vector.terms())
        try: 
            vector_array.append(matrix(myring, my_vector).T)
        except:
            vector_array.append(matrix(myring, my_vector.list()).T)
        #show('added', my_vector, 'to vector array, which is now', vector_array)
    
    return vector_array

        
    # ========================================================================
    # --- Step 1
    # ========================================================================

def step1(module, maps):
    module_copy = deepcopy(module)
    for i in range(len(module)):
        for j in range(len(module)):
            morphism = maps[i][j]
            if morphism!=0:
                nullspace = morphism.kernel()
                domain_flag = module_copy[i]
                if nullspace not in domain_flag:
                    domain_flag.append(nullspace)
                    #show('got', nullspace, 'from the kernal of', morphism)
                    for subspace in domain_flag:
                        if nullspace.intersection(subspace) not in domain_flag:
                            domain_flag.append( nullspace.intersection(subspace) )
                            #show("We got an an intersection we didn't have before: \t (i=\t",i, " j=\t", j, ") \t \t", nullspace, subspace )
                        
                    
                image = morphism.image()
                codomain_flag = module_copy[j]
                if image not in codomain_flag:
                    codomain_flag.append(image)
                    #show('got', image, 'from the image of', morphism)
                    for cosubspace in codomain_flag:
                        if image.intersection(cosubspace) not in codomain_flag:
                            codomain_flag.append( image.intersection(cosubspace) )
                            #show("We got an an intersection we didn't have before: \t (i=\t",i, " j=\t", j, ") \t \t" , image, cosubspace )
    return module_copy
    
    # ========================================================================
    # --- Step 2
    # ========================================================================

def step2(module,maparray):
    module_copy = deepcopy(module)
    for i in range(len(module)):
        for j in range(len(module)):
            morphism = maparray[i][j]
            if morphism == 0:
                continue
            domain_flag = module_copy[i]
            codomain_flag = module_copy[j]
            V = morphism.domain()
            W = morphism.codomain()
            
            for cosubspace in module[j]:
                #show( "---------------", 'i= \t ' , i, '\t j=\t ', j,'\t --- \t V is \t' , V, '\t W is:\t ', W, "---------------")
                #show("morphism ", matrix(morphism).T)  
                #show("cosubspace",  subspace_to_basis_vector_array(cosubspace), '\t index: \t ', codomain_flag.index(cosubspace))
                inverse_image = morphism.inverse_image(cosubspace)
                #show("inverse_image \t ",  subspace_to_basis_vector_array(inverse_image))
                
                if inverse_image not in domain_flag:
                    #show('added the inverse image\t ', inverse_image, '\t of the cosubspace\t', cosubspace, '\t under the morphism\t', Matrix(morphism)) #, '\t to\t ', domain_flag)
                    domain_flag.append(inverse_image)
                    for subspace in domain_flag:
                        if inverse_image.intersection(subspace) not in domain_flag:
                            domain_flag.append( inverse_image.intersection(subspace) )
                            #show("We got a new intersection from inverse images: \t", inverse_image.intersection(subspace), "\t ,  (i=\t",i, " j=\t", j, ") \t \t", inverse_image, subspace )
                #else:
                    #show('the subspace: \t ', inverse_image, ' \t is already in \t ', domain_flag)
            
            for my_subspace in module[i]:
                #show( "---------------", 'i= \t ' , i, '\t j=\t ', j,'\t --- \t V is \t' , V, '\t W is:\t ', W, "---------------")
                #show("morphism ", matrix(morphism).T, "my_subspace", my_subspace)  
                subspace_basis = my_subspace.basis()
                #show("subspace",  subspace_to_basis_vector_array(my_subspace), '\t basis \t', subspace_basis, '\t index: \t ', domain_flag.index(my_subspace))
                    
                image_basis = morphism( my_subspace.ambient_vector_space().subspace( my_subspace.basis_matrix()   ) ).basis()
                #show('image_basis: \t' , image_basis)
                image_subspace = W.subspace(image_basis)
                #show('image_subspace: \t' , image_subspace)
                
                if image_subspace not in codomain_flag:
                    #show('added the image \t ', image_subspace) #, '\t to\t ', codomain_flag)
                    codomain_flag.append(image_subspace)
                    for cosubspace in codomain_flag:
                        if cosubspace.intersection(image_subspace) not in codomain_flag:
                            codomain_flag.append( cosubspace.intersection(image_subspace) )
                            #show("We got a new intersection from images: \t", cosubspace.intersection(image_subspace),  "(i=\t",i, " j=\t", j, ") \t \t", image_subspace, cosubspace )
                #else:
                    #show('the subspace: \t ', image_subspace, ' \t is already in \t ', codomain_flag)
                
                
    return module_copy
 
    # ========================================================================
    # --- Some Functions you need to draw multiflags as digraphs:
    # ========================================================================
   
def subspace_to_basis_as_tuples(mysubspace):
    mybasis = mysubspace.basis()
    
    if len(mybasis)==0:
        return tuple([tuple([0]*mysubspace.degree())])
    
    new_basis=[]
    for myvector in mybasis:
        myvector_as_tuple = tuple(myvector.list())
        new_basis.append(myvector_as_tuple )
    return tuple(new_basis)
    
def multiflag_using_subspaces_to_flag_using_tuple_of_basis_tuples(mymultiflag):
    new_multiflag = []
    for mysubspace in mymultiflag:
        new_multiflag.append(subspace_to_basis_as_tuples(mysubspace))
    return tuple(new_multiflag)
    
def basis_as_tuples_to_subspace(my_tuple_tuple):
    sample_vector = my_tuple_tuple[0]
    dimension = len(sample_vector)
    ambient_space = VectorSpace(myring, dimension)
    
    my_basis =[]
    for tuple_vector in my_tuple_tuple:
        true_vector = ambient_space(tuple_vector)
        my_basis.append(true_vector)
    
    true_subspace = ambient_space.subspace(my_basis)
    return true_subspace   

def subspace_to_column_matrix(my_subspace):
    return matrix(myring, subspace_to_basis_as_tuples( my_subspace ) ).T

def multiflag_using_subspaces_to_multiflag_using_column_matrices(my_multiflag):
    new_multiflag = []
    for my_subspace in my_multiflag:
        subspace_matrix =subspace_to_column_matrix(my_subspace)
        subspace_matrix.set_immutable()
        new_multiflag.append(subspace_matrix)
    return new_multiflag
    
def column_span_of_matrix_is_subspace(my_matrix_1, my_matrix_2):
    return span(my_matrix_1.T).is_subspace(span(my_matrix_2.T))
    
def multiflag_to_poset(multiflag):
    multiflag_using_column_matricies = multiflag_using_subspaces_to_multiflag_using_column_matrices(multiflag)
    my_poset = Poset([multiflag_using_column_matricies, column_span_of_matrix_is_subspace])
    
    return my_poset      

def draw_multiflag(my_multiflag):
    multiflag_using_column_matricies = multiflag_using_subspaces_to_multiflag_using_column_matrices(my_multiflag)
    #show(multiflag_using_column_matricies)
   
    space_dimension =  multiflag_using_column_matricies[0].nrows()
    space_codimension =  multiflag_using_column_matricies[0].ncols()
    my_poset = Poset([multiflag_using_column_matricies, column_span_of_matrix_is_subspace])
    
    cell_height = 0.2 * space_dimension + 0.1
    cell_width =  0.2 * space_codimension + 0.1
    figure_height = max(1 + my_poset.height() * cell_height,2)
    figure_width  = max(1 + my_poset.width()  * cell_width,2)
    
    #show(my_poset, vertex_size = 500+500*space_dimension, edge_thickness=30,vertex_color='white', edge_color='gray', vertex_shape=',', figsize=(figure_width, figure_height)) 
    return my_poset.plot(vertex_size = 500+500*space_dimension, edge_thickness=30,vertex_color='white', edge_color='gray', vertex_shape=',', figsize=(figure_width, figure_height)) 

    # ========================================================================
    # --- Some different functions for turning multiflags into posets
    # ========================================================================    

def list_of_relations_from_multiflag(multiflag):
    multiflag_size = len(multiflag)
    base_list = range(multiflag_size)
    product_list = [[x,y] for x in base_list for y in base_list  ]
    
    relation_list = []
    for i in product_list:
        if multiflag[i[0]].is_subspace( multiflag[ i[1] ] ):
            relation_list.append(i)
    
    return relation_list
    
def poset_of_multiflag_with_integer_labels( multiflag ):
    myposet = Poset([ [],   list_of_relations_from_multiflag(multiflag)])
    
    return myposet
    
def draw_multiflag_with_integer_labels(my_multiflag):
   
    space_dimension   = my_multiflag[0].dimension()
    space_codimension = space_dimension
    my_poset = poset_of_multiflag_with_integer_labels( my_multiflag )
    
    cell_height = 0.2 * space_dimension + 0.1
    if space_dimension >8: cell_height +=5 
    
    cell_width =  0.2 * space_codimension + 0.1
    if space_codimension >8: cell_width +=5 
    
    figure_height = max(1 + my_poset.height() * cell_height, 2)
    figure_width  = max(1 + my_poset.width()  * cell_width, 2)
    
    #show(my_poset, vertex_size = 500+500*space_dimension, edge_thickness=30,vertex_color='white', edge_color='gray', vertex_shape=',', figsize=(figure_width, figure_height)) 
    return my_poset.plot(vertex_size = 500+500*space_dimension, edge_thickness=30,vertex_color='white', edge_color='gray', vertex_shape=',', figsize=(figure_width, figure_height)) 

def get_xy_coords_of_subspaces_in_multiflag( multiflag ):
    
    mygraphic = draw_multiflag_with_integer_labels(multiflag)
    labeled_coords = get_vertex_coords_from_graphics_object( mygraphic )
    
    coord_list = [0]*len(labeled_coords)
    
    for vertex in labeled_coords:
        vertex_index = int( vertex[0] )
        vertex_coord = vertex[1]
        
        coord_list[vertex_index] = vertex_coord
    
    return coord_list
    
def local_xy_coords_for_all_subspaces_in_module( mymodule ):
    module_size = len(mymodule)
    coord_list = [0]*module_size
    
    for i in range(module_size):
        coord_list[i] = get_xy_coords_of_subspaces_in_multiflag( mymodule[i] )
        
    return coord_list

    # ========================================================================
    # --- Drawing Modules ---------
    # ======================================================================== 

def map_array_to_dictionary(map_array):
    module_size = len(map_array)
    my_dict = {}
    for i in range(len(map_array)):
        my_dict[i]={}
        for j in range(len(map_array)):
            if map_array[i][j] != 0:
                my_dict[i][j]=Matrix(map_array[i][j]).T
    
    return my_dict

def digraph_of_map_array(map_array):
    return DiGraph( map_array_to_dictionary(map_array) )

#def draw_module_outline_from_map_array(map_array): # For some reason, this doesn't plot anything in the Jupyter notebook, but run manually, the code works.
#    print 'I ran'
#    my_digraph = digraph_of_map_array(map_array)
#    print 'got the digraph it has sink:', my_digraph.sinks()
#    my_digraph.plot(edge_labels=True)

    # ========================================================================
    # --- Checking that a map array commutes ---------
    # ======================================================================== 

def morphism_from_path_and_map_array(path, map_array):
    start_vertex = path[0]
    end_vertex = path[1]
    new_morphism = map_array[start_vertex][end_vertex]
    if new_morphism == 0:
        raise ValueError("Could not compute the composition for the path", path, "There is not a morphism from ", start_vertex, "to", end_vertex )
    
    if len(path) ==2:
        return new_morphism
    
    for i in range(1, len(path) -1):
        old_morphism = new_morphism
        start_vertex = path[i]
        end_vertex = path[i+1]
        
        if map_array[start_vertex][end_vertex] == 0:
            raise ValueError("Could not compute the composition for the path", path, "There is not a morphism from ", start_vertex, "to", end_vertex )
        new_morphism = map_array[start_vertex][end_vertex] * old_morphism
    
    return new_morphism
    
def does_map_array_commute(map_array):
    module_size = len(map_array)
    module_digraph = digraph_of_map_array(map_array)
    
    for start_vertex in range(module_size):
        paths = module_digraph.all_simple_paths( [start_vertex] )
        #print paths
        
        for end_vertex in range(module_size):
            good_paths = []
            for path in paths:
                if path[-1] == end_vertex:
                    #print "added ", path, " which ended in ", path[-1], "for the vertex", end_vertex
                    good_paths.append(path)
            #print "The paths from ", start_vertex, " to ", end_vertex, "are", good_paths
            if len(good_paths)>1:
                #print "There is more than one path in ", good_paths
                for i in range(len(good_paths) -1):
                    if morphism_from_path_and_map_array(good_paths[i], map_array)  != morphism_from_path_and_map_array(good_paths[i+1], map_array):
                        return False
    
    return True
    
    # ========================================================================
    # --- Stabilizing local structure
    # ======================================================================== 

def get_multiflag_sizes_of_module(my_module):
    module_size = len(my_module)
    size_list = [0]*module_size
    for i in range(module_size):
        size_list[i] = len(my_module[i])
    return size_list
    
def get_stable_local_structure(my_map_array, max_iterations=20, starting_module=0):
    if not does_map_array_commute(my_map_array):
        print '==============================================================\n' +        '-------------WARNING: MAP ARRAY DOES NOT COMMUTE-------------\n' +        '=============================================================='
    
    if starting_module == 0:
        my_space_dimensions = space_dimensions_from_map_array(my_map_array)
        my_trivial_module = trivial_module(my_space_dimensions, myring)

    else:
        my_trivial_module = starting_module
    
    
    old_size_list = get_multiflag_sizes_of_module(my_trivial_module)
    old_module=my_trivial_module
    iterations=0
    
    new_module = step1(my_trivial_module, my_map_array)
    new_size_list = get_multiflag_sizes_of_module(new_module)
    if old_size_list == new_size_list:
        return( [1, first_module] )
    iterations =1
    
    for i in range(max_iterations):
        print i, "'th iteration"
        old_module = new_module
        old_size_list = new_size_list
        
        new_module = step2(old_module, my_map_array)
        new_size_list = get_multiflag_sizes_of_module(new_module)
        iterations +=1
        if old_size_list == new_size_list:
            return( [iterations, new_module] )
            
    return [False, new_module]

def local_structure_from_map_array(map_array, max_iterations=20):
    my_space_dimensions_list = space_dimensions_from_map_array(map_array)
    my_trivial_module = trivial_module(my_space_dimensions_list, myring)
    return get_stable_local_structure(my_trivial_module, map_array, max_iterations)


    # ========================================================================
    # --- Drawing the local structure
    # ======================================================================== 

    
def image_of_local_structure_from_map_array(map_array, multiflag_image_name = "multiflag_"):
    stable_module = get_stable_local_structure(map_array)[1]
    
    n = len(stable_module)
    size_list = []
    
    for i in range(n):
        multiflag = stable_module[i]
        multiflag_graphics_object = draw_multiflag(multiflag)
        file_name = multiflag_image_name + str(i) + ".png"
        multiflag_graphics_object.save( file_name )
    
    list_of_points = object_positions_from_map_array(map_array) 
    
    #INCOMPLETE. After I write object_positions_from_map_array, I need to draw the multiflags at those positions
    
def object_positions_from_map_array(map_array):
    my_poset = Poset( digraph_of_map_array(map_array) )
    #my_poset.show()
    my_graphic = my_poset.plot()
    description_string = my_graphic.description()
    
    first_split = description_string.split('[')[1]
    second_split = first_split.split(']')[0]
    third_split = second_split.replace('),', ')//')
    list_of_point_strings = third_split.split('//')
    
    list_of_points = []
    for point_string in list_of_point_strings:
        list_of_points.append(  eval(point_string) )
            
    return list_of_points # The co-ordinates are given as mathematicians would use, not the way they are used for computer graphics. i.e. (0,1) is above (0,0) not below.

def transpose_xy_positions(point_list):
    for i in range(len(point_list)):
        point = list(point_list[i])
        point.reverse()
        print point
        point_list[i] = point
    
    return point_list


def draw_skeleton_of_module( map_array ):    
    my_poset = Poset( digraph_of_map_array(map_array) )
    
    point_list = transpose_xy_positions( object_positions_from_map_array(map_array) )
    n_points = len( point_list)
    
    mygraph = Graphics()
    
    for i in range(n_points):
        
        for j in range(n_points):
            if my_poset.covers(i,j):
                start_coord = point_list[i]
                end_coord = point_list[j]
                mycolor = 'gray'
                new_arrow = arrow(  start_coord, end_coord, color=mycolor, arrowsize =3, width = 1  )[0]
                mygraph.add_primitive( new_arrow )
                
                matrix_string = str( map_array[i][j].matrix().transpose() )
                
                mid_x = (end_coord[0] + start_coord[0])/2
                mid_y = (end_coord[1] + start_coord[1])/2
                
                new_text = text( matrix_string , (mid_x, mid_y)  )
                mygraph = mygraph + new_text
        #------
        new_elt = text( str(i) , point_list[i]  )
        mygraph = mygraph + new_elt
        
    
    return mygraph

   
    # ========================================================================
    # --- Computing Associated Gradeds
    # ======================================================================== 


def repack(mysubspace):
    ambient_ring = mysubspace.base_field()
    my_ambient_space = VectorSpace( ambient_ring, mysubspace.ambient_vector_space().dimension() )
    repacked_subspace = my_ambient_space.span( [  my_ambient_space( x.list() ) for x in mysubspace.basis()  ]  )
    
    return repacked_subspace

def list_of_sub_sub_spaces(mysubspace, multiflag):
    return_list = [] # Adding the zero subspace. I may need to fix this later
    for subspace in multiflag:
        if (subspace.is_subspace(mysubspace) and subspace != mysubspace) or subspace == multiflag[1]:
            return_list.append(subspace)       
    
    #print('the sub sub space list is')
    #myprint(return_list)
    return return_list

def sum_of_subspace_list(subspace_list):
    
    ambient_space = subspace_list[0].ambient_vector_space()
    
    spanning_set = []
    for subspace in subspace_list:
        subspace_basis = subspace.basis()
        for basis_vector in subspace_basis:
            spanning_set += [basis_vector.list()]
    
    #show('the sum is', ambient_space.subspace(spanning_set))
    return ambient_space.subspace(spanning_set)
    
    #subspace_sum = span(spanning_set)
    #show('the sum is', subspace_sum)
    #return subspace_sum
    
def get_sub_subspace_sum_from_multiflag(mysubspace, multiflag):
    subspace_list = list_of_sub_sub_spaces(mysubspace, multiflag)
    sum_of_subspaces = sum_of_subspace_list(subspace_list)
    
    return sum_of_subspaces

def list_of_subspace_sums(multiflag):
    multiflag_size = len( multiflag )
    sum_list = [0] * multiflag_size
    
    for i in range(multiflag_size):
        sum_list[i] = get_sub_subspace_sum_from_multiflag(multiflag[i], multiflag)
    
    return sum_list
   
def sub_quotient(mysubspace, multiflag):
    sum_of_subspaces = get_sub_subspace_sum_from_multiflag(mysubspace, multiflag)
    
    if sum_of_subspaces.dimension() > 0:
        mysubquotient = repack( mysubspace ).quotient( repack( sum_of_subspaces ) )
        return mysubquotient
    else:
        return mysubspace
    
def associated_graded(multiflag):
    zero_space = multiflag[1]
    
    graded_list = [0]*len(multiflag)
    for i in range( len(multiflag) ):
        my_sub_quotient = sub_quotient(multiflag[i], multiflag)
        # if multiflag[i].dimension() == 0:
            # graded_list[i] = zero_space
        #else:
            # graded_list[i] = sub_quotient(multiflag[i], multiflag)
        
        if my_sub_quotient not in graded_list[0:i]:
            graded_list[i] = my_sub_quotient
        else:
            graded_list[i] = zero_space
        
    return graded_list
        
def object_excess(multiflag):
    my_associated_graded = associated_graded(multiflag)
    total_graded_dimension = 0
    for i in range( len(multiflag) ):
        total_graded_dimension += my_associated_graded[i].dimension()
    
    ambient_dimension = multiflag[0].dimension()    
    
    return total_graded_dimension - ambient_dimension

    
    # ========================================================================
    # --- More Graphics for drawing modules with their local structure
    # ======================================================================== 

def get_vertex_coords_from_graphics_object(mygraphics): # This get both the **label** (as a string)  and the coordinate
    primative_list = list( mygraphics )
    vertex_list = []
    for prim in primative_list:
        if type(prim) == sage.plot.text.Text:
            vertex_list.append( [ prim.string, (prim.x, prim.y)] ) # The string prim.string is a unicode string. Hopefully that doesn't cause problems.
    
    return vertex_list

def get_labels_and_coords_of_module_digraph( map_array ):
    mydigraph  = digraph_of_map_array(map_array)
    mygraphics = Poset( mydigraph ).plot()
    
    vertical_coord_list =  get_vertex_coords_from_graphics_object(mygraphics) # Returns [label, (x,y)]
    
    horizontal_coord_list = [0]*len( vertical_coord_list )
    for i in range( len(vertical_coord_list)  ):
        horizontal_coord_list[i] = [vertical_coord_list[i][0], (vertical_coord_list[i][1][1], vertical_coord_list[i][1][0]  )  ]
    
    return horizontal_coord_list # looks like [[u'0', (0, 1)], [u'1', (1, 0.6)], [u'2', (1, 1.3)], [u'3', (2, 1)]]

def digraph_of_module_from_local_structure(mymodule):
    # -------- Define a set of verticies ------------
    module_size = len(mymodule)
    vertex_set = []
    
    for i in range(module_size):
        for j in range( len( mymodule[i])):
            vertex_set.append( (i,j) )
            
    # -------- Add the Verticies to the graph ------------
    mygraph = DiGraph()
    mygraph.add_vertices(vertex_set)
    mygraph.plot()
    
    # -------- Define a list of relations ------------
    relation_list = []
    for v in vertex_set:
        for w in vertex_set:
            if v!=w and v[0] == w[0] and mymodule[v[0]][v[1]].is_subspace( mymodule[w[0]][w[1]]):
                relation_list.append((v,w))
    
    # -------- Define a poset ------------
    #print(vertex_set)
    #print(relation_list)
    myposet = Poset( (vertex_set, relation_list) )
        
    return myposet

def get_max_xy_in_graphic(mygraphics):
    labeled_coord_list = get_vertex_coords_from_graphics_object(mygraphics)
    
    x_coord_list = [0]*len( labeled_coord_list )
    y_coord_list = [0]*len( labeled_coord_list )
    for i in range( len(labeled_coord_list)):
        x_coord_list[i] = labeled_coord_list[i][1][0]
        y_coord_list[i] = labeled_coord_list[i][1][1]
    
    x_max = max(x_coord_list)
    y_max = max(y_coord_list)
        
    return [x_max,y_max]
    
def get_min_xy_in_graphic(mygraphics):
    labeled_coord_list = get_vertex_coords_from_graphics_object(mygraphics)
    
    x_coord_list = [0]*len( labeled_coord_list )
    y_coord_list = [0]*len( labeled_coord_list )
    for i in range( len(labeled_coord_list)):
        x_coord_list[i] = labeled_coord_list[i][1][0]
        y_coord_list[i] = labeled_coord_list[i][1][1]
    
    x_max = min(x_coord_list)
    y_max = min(y_coord_list)
        
    return [x_max,y_max]

def get_max_xy_in_coord_list(coord_list):
    transposed = map(list, zip(*coord_list ))
    return [max( transposed[0] ), max( transposed[1] )  ]

def get_min_xy_in_coord_list(coord_list):
    transposed = map(list, zip(*coord_list ))
    return [min( transposed[0] ), min( transposed[1] )  ]
    
def get_plot_sizes_of_multiflags(mymodule):
    
    size_list = []
    
    for multiflag in mymodule:
        multiflag_plot = draw_multiflag( multiflag )
        multiflag_xy = get_max_xy_in_graphic( multiflag_plot  )
        size_list.append(multiflag_xy)
    
    return size_list

def get_max_plot_size_of_multiflag( mymodule ):
    size_list = get_plot_sizes_of_multiflags(mymodule)
    n_coords = len(size_list)
    
    x_coord_list = [0]*n_coords
    y_coord_list = [0]*n_coords
    for i in range( n_coords):
        x_coord_list[i] = size_list[i][0]
        y_coord_list[i] = size_list[i][1]
    
    x_max = max(x_coord_list)
    y_max = max(y_coord_list)
        
    return [x_max,y_max]

def generate_display_coords_for_modules_with_local_structure(mymodule, map_array, x_space = 'default', y_space = 'default'): # should return a list like [ [[0,0], (x,y)], [[0,1], (x,y)], ...   ]. i.e. specifying to coords of each vector subspace
    module_size = len(map_array)
    
    box_for_each_multiflag = get_max_plot_size_of_multiflag(mymodule)
    x_for_each_multi_flag  = box_for_each_multiflag[0]
    y_for_each_multi_flag  = box_for_each_multiflag[1]
    
    #It's better to make the figsize big than to add a lot of margin
    if x_space == 'default':
        x_space = -1 #x_for_each_multi_flag/4    # The margin between each multiflag  
    if y_space == 'default':
        y_space = y_for_each_multi_flag +10
    
    total_x = 0     # This will be increased below
    total_y = 0
    
    
    #------ Add a margin between each multiflag ----------
    module_poset = Poset( digraph_of_map_array(map_array) )
    nrows = module_poset.width()
    ncols = module_poset.height()
    
    total_x +=(ncols-1)*x_space
    total_y +=(nrows-1)*y_space
    
    #------ Add space for each multiflag ----------
    
    total_x += ncols * x_for_each_multi_flag
    total_y += nrows * y_for_each_multi_flag
    
    #------ Get the xy size of the digraph for the map array ----------
    
    base_box = get_max_xy_in_graphic( Poset( digraph_of_map_array(map_array) ).plot() )
    
    x_scale = total_x/base_box[0]
    y_scale = total_y/base_box[1]
    
    # ------ Start assigning coords to the vector subspaces ---------
    
    skeleton = get_labels_and_coords_of_module_digraph(map_array)
    unlabeled_positions = [0]*module_size
    relative_coords = local_xy_coords_for_all_subspaces_in_module( mymodule )
    
    for i in range( module_size ):
        multiflag_size = len( mymodule[i] )
        unlabeled_positions[i] = [0]*multiflag_size
    
        
        for j in range( multiflag_size ):
            unlabeled_positions[i][j] = list( skeleton[i][1] )                  # Sets coords for the multiflag the subspace is in
            unlabeled_positions[i][j][0] = unlabeled_positions[i][j][0]*x_scale # Spaces out the the multiflags
            unlabeled_positions[i][j][1] = unlabeled_positions[i][j][1]*y_scale # Spaces out the the multiflags
            unlabeled_positions[i][j][0] += relative_coords[i][j][0] +random()/20    # Adds the x-offset within the multiflag and a random factor
            unlabeled_positions[i][j][1] += relative_coords[i][j][1] +random()/20   # Adds the y-offset within the multiflag and a random factor
    
            # If there is a problem with the labels over lapping, then I think the most likely way to fix it is to increase the spacing in relative_coords()
            
            
            #print( i,j, 
            #round(unlabeled_positions[i][j][0],3), 
            #round(unlabeled_positions[i][j][1],3), 
            #'skel', 
            #round(skeleton[i][1][0],3), round(skeleton[i][1][1],3), 
            #'scaled', 
            #round(skeleton[i][1][0]*x_scale,3), 
            #round(skeleton[i][1][1]*y_scale,3), 
            #'relative', 
            #relative_coords[i][j]
            #)
        
    
    # --------- The graphic has shifted away from (0,0), so I should shift it back to (1,1)
    
    raw_positions = []
    for i in range(module_size):
        raw_positions += unlabeled_positions[i]
    
    shift = get_min_xy_in_coord_list( raw_positions  )
    shift[0]-=1
    shift[1]-=1
    
    for i in range( module_size ):
        multiflag_size = len( mymodule[i] )        
        for j in range( multiflag_size ):
            unlabeled_positions[i][j][0] -= shift[0]
            unlabeled_positions[i][j][1] -= shift[1]
    
    # --------- Purely for debugging
    test_graphic = Graphics()
    
    for i in range( module_size ):
        multiflag_size = len( mymodule[i] )
        for j in range( multiflag_size ):
            new_elt = text( str(i) + str(j), unlabeled_positions[i][j]  )
            test_graphic = test_graphic + new_elt
    
    #test_graphic.show(figsize=[9,9])
    
    #return test_graphic
    return unlabeled_positions
    
def get_string_rep_for_subspace_column_basis( mysubspace ):    
    # This whole part is just so that I can use the same function for subspaces and for lists of basis vectors.
    subspace_copy = mysubspace
    if type(mysubspace) == list:
        if mysubspace !=[]:
            ambient_space = VectorSpace(myring, len(mysubspace[0]))
            subspace_copy = ambient_space.span( mysubspace )
        else:
            return "[]"
    
    # For actual subspaces, this is the only code that is needed
    return str( subspace_copy.basis_matrix().transpose() )
    
def get_string_reps_for_all_subspace_column_basis( mymodule ):
    module_size = len(mymodule)
    string_reps = [0]*module_size
    
    for i in range(module_size):
        multiflag_size = len(mymodule[i])
        string_reps[i] = [0]*multiflag_size
        
        for j in range(multiflag_size):
            string_reps[i][j] = get_string_rep_for_subspace_column_basis( mymodule[i][j])
            
    return string_reps

def get_string_rep_for_indicies_of_subspaces(mymodule):
    module_size = len(mymodule)
    string_reps = [0]*module_size
    
    for i in range(module_size):
        multiflag_size = len(mymodule[i])
        string_reps[i] = [0]*multiflag_size
        
        for j in range(multiflag_size):
            string_reps[i][j] = str( (i,j) )
            
    return string_reps    
  
def default_color_list( pos_list, color = 'black' ):
    module_size = len( pos_list )
    color_list = [0]* module_size
    
    for i in range(module_size):
        multiflag_size = len(pos_list[i])
        color_list[i] = [0]*multiflag_size
        
        for j in range(multiflag_size):
            color_list[i][j] = color
    return color_list

def arrow_list_for_subspace_inclusions_within_multiflags( mymodule, color = (240, 240, 245) ):
    module_as_digraph = digraph_of_module_from_local_structure( mymodule )
    arrow_list = module_as_digraph.cover_relations()
    
    for i in range(len(arrow_list)):
        arrow_list[i].append(color)
            
    return arrow_list  # Looks like a list of pairs of tuples [ [(3, 1), (3, 2)],   [(3, 2), (3, 0)]  ]

def arrow_list_for_subspace_images( mymodule, map_array, color = 'gray' ):
    arrow_list=[]
    
    module_size = len( mymodule)
    for i in range(module_size):
        multiflag_size = len(mymodule[i])
        
        for j in range(multiflag_size):
            start_index = [i,j]
            
            for i2 in range(module_size):
                multiflag_size2 = len(mymodule[i2])
                
                for j2 in range(multiflag_size2):
                    end_index = [i2,j2]
                    start_subspace = mymodule[i][j]
                    end_subspace = mymodule[i2][j2]
                    morphism = map_array[i][i2]
                    
                    if morphism!=0 and morphism(repack(start_subspace)) == end_subspace:
                        arrow_list.append([[i,j], [i2,j2], color  ])
    
    return arrow_list

def arrow_list_for_subspace_inverse_images( mymodule, map_array, color = 'gray' ):
    arrow_list=[]
    
    module_size = len( mymodule)
    for i in range(module_size):
        multiflag_size = len(mymodule[i])
        
        for j in range(multiflag_size):
            start_index = [i,j]
            
            for i2 in range(module_size):
                multiflag_size2 = len(mymodule[i2])
                
                for j2 in range(multiflag_size2):
                    end_index = [i2,j2]
                    start_subspace = mymodule[i][j]
                    end_subspace = mymodule[i2][j2]
                    morphism = map_array[i][i2]
                    
                    if morphism!=0 and morphism.inverse_image(end_subspace) == start_subspace:
                        arrow_list.append([ [i2,j2], [i,j], color  ])
    
    return arrow_list
    
def arrow_list_for_blocks( mymodule, map_array, color = 'gray' ):
    old_arrow_list = arrow_list_for_subspace_images( mymodule, map_array)
    
    mymodule = basis_extensions_for_all_multflags_in_module( mymodule )
    new_arrow_list = []
    for arrow in old_arrow_list:
        arrow_start = arrow[0]
        arrow_end = arrow[1]
        
        if (mymodule[arrow_start[0]][arrow_start[1]] !=[]) and (mymodule[arrow_end[0]][arrow_end[1]] !=[]) :
            new_arrow_list.append(arrow)
    
    return new_arrow_list
 
def graphic_from_positions_and_labels(pos_list, label_list, color_list =[], arrow_list = []): # The arrow list should contain the colors
    if color_list == []:
        color_list = default_color_list(pos_list)
    
    module_size = len( pos_list )
    mygraphic = Graphics()
    
    for mylist in arrow_list:
        
        start_coord = list (  pos_list[mylist[0][0]][mylist[0][1]] )
        end_coord = list( pos_list[mylist[1][0]][mylist[1][1]] )
        
        vector_x = end_coord[0] - start_coord[0]
        vector_y = end_coord[1] - start_coord[1]
        unit_x = vector_x /sqrt( vector_x**2 +  vector_y**2)
        unit_y = vector_y /sqrt( vector_x**2 +  vector_y**2)
        #print 'x', round(unit_x, 3), 'y', round(unit_y, 3)
        
        start_coord[0] += unit_x * .3
        start_coord[1] += unit_y * .3
        end_coord[0] -= unit_x *.4
        end_coord[1] -= unit_y * .4
        
        new_arrow = arrow(  start_coord, end_coord, color = mylist[2], arrowsize =3, width = 1  )[0]
        mygraphic.add_primitive( new_arrow )
    
    for i in range(module_size):
        multiflag_size = len(pos_list[i])
        
        for j in range(multiflag_size):
            # Add labels
            new_label = text(label_list[i][j], pos_list[i][j], color = color_list[i][j])
            mygraphic += new_label
            
            # Add arrows
    
    return mygraphic

# ========================================================================
# ========================================================================
# --- The graphics Suite
# ========================================================================     
# ========================================================================

# For some reason this is giving me the error """ got an unexpected keyword argument 'arrow_type' """
def draw_local_subspace_inclusions_of_module_from_map_array(map_array, arrow_type = '', label_type=''):
    stablemodule = get_stable_local_structure(map_array, max_iterations=20)[1] 
    display_coords = generate_display_coords_for_modules_with_local_structure(stablemodule, map_array)
    
    # --- Labels ----
    if label_type == '':
        label_list = get_string_reps_for_all_subspace_column_basis( stablemodule)
    if label_type == 'subquotient':
        new_module = vectors_to_extend_basis(stablemoduleB[1][0], stablemoduleB[1][0], stablemoduleB[1][0])
        label_list = get_string_reps_for_all_subspace_column_basis( new_module)
    
    # --- Arrows ---
    if arrow_type == '':
        arrow_list = arrow_list_for_subspace_inclusions_within_multiflags( stablemodule )
    
    elif arrow_type =='images':
        arrow_list = arrow_list_for_subspace_images( stablemodule, map_array)
    
    elif arrow_type =='inverses':
        arrow_list = arrow_list_for_subspace_inverse_images( stablemodule, map_array)
    
    # --- Graph ---
    graph  = graphic_from_positions_and_labels(display_coords, label_list, arrow_list = arrow_list)
    
    return graph


# ----- The Functions below here should be obsolete, since the one above handles all their cases ---------

def draw_local_subspace_inclusions_of_module_from_map_array(map_array):
    stablemodule = get_stable_local_structure(map_array, max_iterations=20)[1] 
    display_coords = generate_display_coords_for_modules_with_local_structure(stablemodule, map_array)
    label_list = get_string_reps_for_all_subspace_column_basis( stablemodule)
    
    arrow_list = arrow_list_for_subspace_inclusions_within_multiflags( stablemodule )
    graph  = graphic_from_positions_and_labels(display_coords, label_list, arrow_list = arrow_list)
    
    return graph
    
def draw_subspace_images_of_module_from_map_array(map_array):
    stablemodule = get_stable_local_structure(map_array, max_iterations=20)[1] 
    display_coords = generate_display_coords_for_modules_with_local_structure(stablemodule, map_array)
    label_list = get_string_reps_for_all_subspace_column_basis( stablemodule)
    
    arrow_list = arrow_list_for_subspace_images( stablemodule, map_array)
    graph  = graphic_from_positions_and_labels(display_coords, label_list, arrow_list = arrow_list)
    
    return graph
    
def draw_subspace_inverse_images_of_module_from_map_array(map_array):
    stablemodule = get_stable_local_structure(map_array, max_iterations=20)[1] 
    display_coords = generate_display_coords_for_modules_with_local_structure(stablemodule, map_array)
    label_list = get_string_reps_for_all_subspace_column_basis( stablemodule)
    
    arrow_list = arrow_list_for_subspace_inverse_images( stablemodule, map_array)
    graph  = graphic_from_positions_and_labels(display_coords, label_list, arrow_list = arrow_list)
    
    return graph

def draw_subspace_indicies_of_module_from_map_array(map_array):
    stablemodule = get_stable_local_structure(map_array, max_iterations=20)[1] 
    display_coords = generate_display_coords_for_modules_with_local_structure(stablemodule, map_array)
    
    label_list = get_string_rep_for_indicies_of_subspaces( stablemodule )
    arrow_list = arrow_list_for_subspace_inclusions_within_multiflags( stablemodule )
    graph  = graphic_from_positions_and_labels(display_coords, label_list, arrow_list = arrow_list)
    
    return graph


def draw_blocks(map_array):
    stablemodule = get_stable_local_structure(map_array, max_iterations=20)[1]
    block_subspaces = basis_extensions_for_all_multflags_in_module( stablemodule ) 
    display_coords = generate_display_coords_for_modules_with_local_structure(stablemodule, map_array)
    label_list = get_string_reps_for_all_subspace_column_basis( block_subspaces )
    
    arrow_list = arrow_list_for_blocks( stablemodule, map_array)
    graph  = graphic_from_positions_and_labels(display_coords, label_list, arrow_list = arrow_list)
    
    return graph

# ========================================================================
# --- Computing Blocks -----
# ========================================================================

def vectors_to_extend_basis(small_space, large_space, ambient_space):
    if not small_space.is_subspace( large_space ):
        raise ValueError( small_space, "\n\n is not a subspace of \n\n", large_space  )
    
    small_basis = list( repack( small_space ).basis() )
    large_basis = repack( large_space ).basis()
    
    #dimension_difference = large_space.dimension() - small_space.dimension()
    
    new_vector_list = []
    new_basis = small_basis     # initialing for the for loop
    
    for i in large_basis: # Check if the vector i is already in the space built up so far. i.e. if i should be included in the extension for the basis    
        if ambient_space.span( new_basis + [i] ) != ambient_space.span( new_basis ):
            #print('got here', i)
            new_vector_list.append(i)
            #print('got past the first append', i)
            new_basis.append(i)
        #print 'too the end of for', i
            
    return( new_vector_list )

def basis_extensions_for_multiflag( multiflag ):
    multiflag_size = len( multiflag )
    extension_list =[[]]*multiflag_size
    
    for i in range(multiflag_size):
        subspace = multiflag[i]
        subspace_sum = get_sub_subspace_sum_from_multiflag(subspace, multiflag)
        
        extension_list[i] = vectors_to_extend_basis(subspace_sum, subspace, multiflag[0])
    
    return extension_list

def basis_extensions_for_all_multflags_in_module( mymodule ):
    extension_list = []
    
    for i in range( len( mymodule )):
        extension_list.append(basis_extensions_for_multiflag( mymodule[i] ))
    
    return extension_list
    