from snappy import *
from slopes import *
from turaev import *
from copy import *
import ast
import pandas
#import turaev

qht_final = pandas.read_csv('./src/QHSolidTori.csv.bz2')
#This was useful but does not work with old pandas
#pandas.options.display.max_colwidth= None

class man_inv:
    def __init__(self,Manifold):
        self.volume=Manifold.volume()
        self.homology=Manifold.homology()
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return ( abs(self.volume - other.volume) < 0.000001 and self.homology == other.homology)
        else:
            return False
    def __repr__(self):
        return str(self.volume)+", "+str(self.homology)
    def __str__(self):
        return str(self.volume)+", "+str(self.homology)
    

def inside_man_inv(element, lista):
    for x in lista:
        if ( (abs(x.volume - element.volume) < 0.000001) and x.homology == element.homology):
            return True
    return False



def is_hyperbolic(M):
    #Returns true if M (maybe canonized) admits a good solution_type and the triangulation that supports such solution

    #We do this so M does not change outside the function
    M=M.copy()

    if M.solution_type(enum=True) not in allowed_solution_type or M.volume()<0.94:
        #Sometimes we need several tries to find a good solution_type
        count=0
        while M.solution_type(enum=True) not in allowed_solution_type or M.volume()<0.94:
            try:
                M.canonize()
            except:
                pass
            count=count+1
            if count >50:
                return [False, M]
    return [True, M]

def my_boolean(value):
    if isinstance(value, bool):
        return bool
    elif isinstance(value, int):
        if value > 0:
            return True
        else:
            return False
    else:
        raise Exception("I should convert " + str(value) +" into a boolean, i do not know how.")
        
def Slope_valuation(stri):
    if stri=="None":
        return None
    elif stri.startswith("SlopeCone"):
        x=stri.split("SlopeCone")
        extr=ast.literal_eval(x[1])
        if isinstance(extr[0],tuple):
            S=SlopeCone(extr[0],extr[1])
            return S
        else:
            S=SlopeCone(extr)
            return S
    elif stri.startswith("SingleSlope"):
        x=stri.split("SingleSlope")
        extr=ast.literal_eval(x[1])
        S=SingleSlope(extr)
        return S
    else:
        raise Exception("I could not evaluate this string as a slope")

def search_in_census_if_L_space(string):
    identity=string
    if not identity == []:
        x=str(identity[0])
        y=str(x).split('(')
        name_qht=y[0]
        filling_qht_str='('+y[1]
        filling_qht=ast.literal_eval(filling_qht_str)
        #Here we check if the filling is already known
        y=qht_final[qht_final['name']==name_qht]['L_space_fillings'].values[0]
        dehn_fillings_L_spaces=ast.literal_eval(y)
        if filling_qht in dehn_fillings_L_spaces or (-filling_qht[0], -filling_qht[1]) in dehn_fillings_L_spaces:
            return 1
        y=qht_final[qht_final['name']==name_qht]['non_L_space_fillings'].values[0]
        dehn_fillings_non_L_spaces=ast.literal_eval(y)
        if filling_qht in dehn_fillings_non_L_spaces or (-filling_qht[0], -filling_qht[1]) in dehn_fillings_non_L_spaces:
            return -1
        y=qht_final[qht_final['name']==name_qht]['non_L_cone'].values[0]
        try:
            S=Slope_valuation(y)
            if Slope(filling_qht) in S:
                return -1
        except:
            pass
        y=qht_final[qht_final['name']==name_qht]['floer_simple'].values[0]
        if y==1:
            y=qht_final[qht_final['name']==name_qht]['non_L_cone'].values[0]
            try:
                S=Slope_valuation(y)
                if Slope(filling_qht) not in S:
                    return 1
            except:
                pass
    print("We found "+ str(identity[0]) +", a filling of one of the 0.2% of QHT in Dunfield census whose L-space fillings are not known. We move on; in principle, it should be possible to check its L-space status quite easily.")
    raise Exception('I was not able to find the value in the census, this is very strange.')

    
#We do not use this function; it will maybe useful to optimize the research
def second_minimal_vol(M, max_coeff):
    N=deepcopy(M)
    volumes=[]
    M.simplify()
    if M.num_tetrahedra() < 8:
        try:
            tau=TuraevTorsion(M)
            if tau.could_be_floer_simple():
                D=IotaInverseDtau(tau)
                X=D.possible_non_L_space_cones((1,0))
                if len(X)==2:
                    return(1000)
            else:
                return(2000)
        except:
            pass
    for h in range(-max_coeff, max_coeff+1):
        for k in range (0, max_coeff+1):
            if gcd(h,k)==1 and h!=-1:
                if M.num_tetrahedra() > 7:
                    M=deepcopy(N)
                    M.dehn_fill((h,k))
                    if M.solution_type(enum=True) in allowed_solution_type and M.volume()>0.94:
                        volumes.append(M.volume())
                else:
                    if Slope(h,k) not in X[0]:
                        M=deepcopy(N)
                        M.dehn_fill((h,k))
                        if M.solution_type(enum=True) in allowed_solution_type and M.volume()>0.94:
                            volumes.append(M.volume())
                    else:
                        volumes.append(1000)
    volumes.sort()
    if len(volumes)<2:
        print(M)
    return volumes[1]
    


def order_curves_by_volume(curves, M, init_string=''):
    #Drills the given curves, and sort them by increasing hyperbolic volume of the drilled manifold.
    #It also tries to identify the drilled manifolds (and using that, the L-space value of M).

    #We do this so M does not change outside the function
    M=M.copy()
    [a,M]=is_hyperbolic(M)
    
    volumes_drilled=[]
    to_pop=[]
    for k in range(0,len(curves)):
        curve=curves[k]
        N=M.drill(curve)
        N=N.filled_triangulation()
        #We try to identify the drilled manifold
        try:
            ids=N.identify()
            if not ids == []:
                if not str(ids[0]).startswith("ocube") and  not str(ids[0]).startswith("odod") and not str(ids[0]).startswith("oicocl"):
                    X=Manifold(ids[0])
                    B=N.is_isometric_to(X, return_isometries=True)
                    new_filling=tuple(B[0].cusp_maps()[0]*vector((1,0)))
                    Manifold_found=ids[0]
                    Manifold_found.dehn_fill(new_filling)
                    is_L_space=search_in_census_if_L_space(ids)
                    print (init_string+': '+M.name()+' is '+str(ids[0])+', whose L-space value is known to be '+str(is_L_space))
                    return int(is_L_space)
        except:
            pass
        [a,N]=is_hyperbolic(N)
        if not a:
            to_pop.append(k)
        else:
            volumes_drilled.append(N.volume())
    #We do not consider the drillings that are not hyperbolic
    if to_pop!=[]:
        to_pop.reverse()
        for h in to_pop:
            curves.pop(h)
    h=list(range(0,len(volumes_drilled)))
    orders=[x for y, x in sorted(zip(volumes_drilled, h))]
    return [x for y, x in sorted(zip(h, curves))]


def search_for_minimal_volume_fillings(M, non_L_space_interval, max_coefficient=15, init_string="", already_found_inv_fill=[]):
#This function takes a manifold M with one cusp and a candidate interval to be the
#L-space interval (the complementary of non_L_space_interval) and searches for
#the two dehn fillings in that interval that minimizes the volume.

    #We do this so M does not change outside the function
    M=M.copy()
    
    [a,M]=is_hyperbolic(M)

    assert a

    N=M.copy()
    vol_min=100000
    vol_min_2=100000
    N_1=N.copy()
    N_2=N.copy()
    minimizing_fillings=[(0,0),(0,0)]
    found_1_L_space=0
    found_2_L_space=0

    #We now look for minimal volume fillings
    for h in range(-max_coefficient, max_coefficient+1):
        for k in range (0, max_coefficient+1):
            if gcd(h,k)==1 and (h,k)!=(-1,0):
                slope=Slope(h,k)
                if slope not in non_L_space_interval:
                    M=N.copy()
                    M.dehn_fill((h,k),0)
                    [a,M]=is_hyperbolic(M)
                    if a:
                        if not inside_man_inv( man_inv(M), already_found_inv_fill ):
                            #Here we try to look in the census if the filling was already known
                            try:
                                ids=M.identify()
                                if not ids == []:
                                    if not str(ids[0]).startswith("ocube") and  not str(ids[0]).startswith("odod") and not str(ids[0]).startswith("oicocl"):
                                        if found_1_L_space==0:
                                            x=search_in_census_if_L_space(ids)
                                            print(init_string+"1: The manifold " + M.name() +" filled with " + str((h,k)) + " is " +str(ids) + ", its L-space value is " + str(x))
                                            found_1_L_space=int(x)
                                            M_vol=0.1
                                        elif found_2_L_space==0:
                                            x=search_in_census_if_L_space(ids)
                                            print(init_string+"2: Found in fillings another one: the manifold " + M.name() +" filled with " + str((h,k)) + " is "+str(ids)+ ": its L-space value is " + str(x))
                                            found_2_L_space=int(x)
                                            M_vol=0.2
                                else:
                                    M_vol=M.volume()
                            except:
                                M_vol=M.volume()

                            if M_vol < vol_min:
                                vol_min_2=vol_min
                                vol_min=M_vol
                                minimizing_fillings[1]=minimizing_fillings[0]
                                minimizing_fillings[0]=(h,k)
                                N_2=N_1.copy()
                                N_1=M.copy()
                            elif M_vol < vol_min_2:
                                vol_min_2=M_vol
                                minimizing_fillings[1]=(h,k)
                                N_2=M.copy()
    if minimizing_fillings[1]==(0,0):
        raise Exception("I could not find two fillings in the interval. Try raising max_coefficient.")
    N_1.set_name(N.name() + str(minimizing_fillings[0]) )
    N_2.set_name(N.name() + str(minimizing_fillings[1]) )
    return([N_1, N_2, found_1_L_space, found_2_L_space])



def search_for_minimal_volume_drillings_floer_simple(M, curves=None, init_string='', max_drills=10, max_segms=6, already_found_inv_dr=[]):
#This function searches for drillings that minimize volume

    #We do this so M does not change outside the function
    M=M.copy()

    #If dual curves were not given, we compute them
    if curves==None:
        curves=M.dual_curves(max_segments=max_segms)

    #We check if M is a hyperbolic manifold
    assert is_hyperbolic(M)[0]

    #We now search for a drilling that is turaev simple and was not already used. Curves should be given in volume order (so the first turaev simple that we find is the one we are looking for)
    minimizing_drilling=-1
    dictionary_volumes={}
    found_one=0
    for ind in range(0, min(max_drills, len(curves))):
        if found_one==0:
            N=M.drill(curves[ind])
            N=N.filled_triangulation()
            [a,N]=is_hyperbolic(N)
            if a:
                if not( inside_man_inv( man_inv(N), already_found_inv_dr ) ) :
                    try:
                        print(init_string+": Computing Turaev torsion drilling...")
                        tau=TuraevTorsion(N)
                        if tau.could_be_floer_simple():
                            found_one=1
                    except:
                        pass
    if found_one==0:
        M.dehn_fill((0,0))
        #print(M.isometry_signature())
        #print([M, curves, init_string, max_drills, max_segms, already_found_inv_dr])
        raise Exception("We could not find a drilling floer simple, hyperbolic and not already used. Try raising max_drills.")
    assert tau.could_be_floer_simple()

    already_found_inv_dr.append(man_inv(N))

    return([N, tau, already_found_inv_dr])



def is_certified_L_space(Man, num_iter=0, max_iter=25, max_coefficient=17, max_segms=6, max_drills=10, init_string="", already_found_inv=[], T=None, tau_T=None, which_interval=2, only_true_hyperbolic_structures=False, save_QHT=False, path_save_QHT="./proofs/", curves_to_avoid=-1):
    #The algorithm takes a QHS Man and tries to prove that it is an L-space. The answer True is rigorous, while the answer False is not.
        
    #We need this because "list" is a mutable type, hence its value is mantained through multiple calls of the function
    if num_iter==0:
        print("Inizializing...")
        if Man.homology().betti_number()!=0:
            raise Exception('The given manifold is not a rational homology sphere.')
        if not Man.is_orientable():
            raise Exception('M is not orientable.')
        try:
            ids=Man.identify()
            if not ids == []:
                if not str(ids[0]).startswith("ocube") and  not str(ids[0]).startswith("odod") and not str(ids[0]).startswith("oicocl"):
                    x=search_in_census_if_L_space(ids)
                    print(init_string+"The manifold is in the census, its L-space value is " + str(x))
                    return(my_boolean(int(x)))
        except:
            pass
        already_found_inv=[]
        global allowed_solution_type
    
        if only_true_hyperbolic_structures==True:
            allowed_solution_type=[1]
        else:
            allowed_solution_type=[1,2]
        
        
        
        Man.set_name('M')
        
    #We add the invariants of Man to the list, in order to avoid Man in the subsequent calls of the function
    already_found_inv.append(man_inv(Man))
        
    #We need this in order to let Man invariant outside the function
    Man=Man.copy()
    
    #We recompute the hyperbolic structure on Man, this solves some issues
    if Man.solution_type(enum=True) not  in allowed_solution_type:
        v1=Man.volume()
        [a, Man]=is_hyperbolic(Man)
        if not a:
            raise Exception(Man.name() + "is not hyperbolic.")
        if abs(v1-Man.volume()) > 0.0001:
            print("We recomputed the hyperbolic structure.")
    
    #If which_interval is not 2, we already have T and tau_T
    if which_interval==2:
        #We check that the given manifold is a QHS
        if Man.homology().betti_number()!=0:
            raise Exception('The given manifold is not a rational homology sphere.')

        print(init_string+": "+Man.name() +" has volume " + Man.volume().str(digits=6) + "... and homology " + str(Man.homology()))
        
        #We interrupt the algorithm if it is taking too long
        if num_iter > max_iter:
            print("Maximal depth reached; try raising max_iter")
            return False


        #We look for a curve to drill that gives a low-volume manifold and is turaev simple
        curves_Man = Man.dual_curves(max_segments=max_segms)
        
        #If curves to avoid is different from -1, we want to avoid some curves. This is usually done to "change the path" we are taking, to prove hard that a manifold is an L-space
        if curves_to_avoid!=-1:
            curves_Man=curves_Man[curves_to_avoid:]
        
        curves_Man = order_curves_by_volume(curves_Man, Man, init_string)
        #Here we check if some drilled manifold was identified
        if isinstance(curves_Man, int):
            if curves_Man==1:
                return True
            elif curves_Man==-1:
                return False
            else:
                raise Exception("The L-space value in the census is "+str(curves_Man)+"; this is not expected.")

        [T, tau_T, already_found_inv]= search_for_minimal_volume_drillings_floer_simple(Man, curves=curves_Man, init_string=init_string, max_drills=max_drills, already_found_inv_dr=already_found_inv)
        #TODO: if M is "simple" (in particular, if it is very fast to compute the turaev_torsion of the drilled manifolds), it could be a good idea
        #to sort the curves by the second minimal volume of the filling that we are going to use.

        T.set_name("T"+init_string)
        if save_QHT==True:
            T.save(path_save_QHT+T.name()+".tri")
            
        assert len(T.cusp_info())==1

        #We now look at minimal volume fillings:
        if not tau_T.could_be_floer_simple():
            raise Exception("Something that should have been floer simple was not floer simple.")

    #We now search for the minimal volume fillings of T in the same possible L-space interval of M

    #We define the possible L-space interval to look at. If there is only one, it is easy to do. Otherwise, if which_interval is not 2, we know which one we want to choose.
    D=IotaInverseDtau(tau_T)
    A=D.possible_non_L_space_cones((1,0))
    #If we have two possible intervals, we try both the possibilities
    if len(A)==2:
        if which_interval==0:
            non_L_sp_interval=A[0]
        elif which_interval==1:
            non_L_sp_interval=A[1]
        else:
            print(init_string+": Double interval..")
            return is_certified_L_space(Man, num_iter=num_iter+1, max_iter=max_iter, max_coefficient=max_coefficient, max_segms=max_segms, max_drills=max_drills, init_string=init_string+"A", already_found_inv=already_found_inv, T=T, tau_T=tau_T, which_interval=0, save_QHT=save_QHT, path_save_QHT=path_save_QHT) or is_certified_L_space(Man, num_iter=num_iter+1, max_iter=max_iter, max_coefficient=max_coefficient, max_segms=max_segms, max_drills=max_drills, init_string=init_string+"B", already_found_inv=already_found_inv, T=T, tau_T=tau_T, which_interval=1, save_QHT=save_QHT, path_save_QHT=path_save_QHT)
    else:
        non_L_sp_interval=A[0]

    #We get M_1 and M_2, the two fillings on T with lower volume that, if are L-spaces, prove that M is an L-space.
    [M_1, M_2, found_1_L_space, found_2_L_space]=search_for_minimal_volume_fillings(T, non_L_space_interval=non_L_sp_interval, max_coefficient=max_coefficient, init_string=init_string, already_found_inv_fill=already_found_inv)

    #If we identified some dehn filling, our work is simpler:
    if found_2_L_space!=0:
        return my_boolean(found_2_L_space) and my_boolean(found_1_L_space)
    elif found_1_L_space!=0:
        return my_boolean(found_1_L_space) and is_certified_L_space(M_2, num_iter=num_iter+1, max_iter=max_iter, max_coefficient=max_coefficient, max_segms=max_segms, max_drills=max_drills, init_string=init_string+"2", already_found_inv=already_found_inv, save_QHT=save_QHT, path_save_QHT=path_save_QHT)
    #Otherwise we need to check both M_1 and M_2
    else:
        return is_certified_L_space(M_2, num_iter=num_iter+1, max_iter=max_iter, max_coefficient=max_coefficient, max_segms=max_segms, max_drills=max_drills, init_string=init_string+"1", already_found_inv=already_found_inv, save_QHT=save_QHT, path_save_QHT=path_save_QHT) and is_certified_L_space(M_1, num_iter=num_iter+1, max_iter=max_iter, max_coefficient=max_coefficient, max_segms=max_segms, max_drills=max_drills, init_string=init_string+"2", already_found_inv=already_found_inv, save_QHT=save_QHT, path_save_QHT=path_save_QHT)

#The following function was used to check the correctness of the algorithm, it can be ignored
def is_coherent_possible_L_space_cone(T, fill1, fill2):
    tau=TuraevTorsion(T)
    D=IotaInverseDtau(tau)
    A=D.non_L_space_cone([fill1,fill2])
    if Slope((1,0)) in A:
        raise Exception("Check what you are doing.")
    else:
        return True
    