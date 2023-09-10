module JKPM
using LinearAlgebra
using Metis
using PyCall
export Vr_exponential_model
export Atoms
export NeighborList
export bulk
export ConstructHamiltonian
export Hamiltonian,dot_product,dot
export ConstructHamiltonian
using LinearAlgebra
using SparseMatricesCSR
using Metis
using StatsBase
pyase = pyimport("ase")
py"""
from ase import neighborlist
from ase.build import bulk
"""
const Pyneighborlist=py"neighborlist"
const PyAtom=pyase.Atoms
const Pybulk=py"bulk"
function Atoms(;kwargs...)
    return PyAtom(;kwargs...)
end
function NeighborList(cutoffs;kwargs...)
    return Pyneighborlist.NeighborList(cutoffs;kwargs...)
end 
function bulk(name;kwargs...)
    return Pybulk(name;kwargs...)
end 
function Vr_exponential_model(
            ss_sig1::Float64, pp_sig1::Float64, pp_pi1::Float64, dd_sig1::Float64, dd_pi1::Float64, dd_del1::Float64, sp_sig1::Float64, sd_sig1::Float64, pd_sig1::Float64, pd_pi1::Float64,
            a_ss_sig1::Float64, a_pp_sig1::Float64, a_pp_pi1::Float64, a_dd_sig1::Float64, a_dd_pi1::Float64, a_dd_del1::Float64, a_sp_sig1::Float64, a_sd_sig1::Float64, a_pd_sig1::Float64, a_pd_pi1::Float64,
            ss_sig2::Float64, pp_sig2::Float64, pp_pi2::Float64, dd_sig2::Float64, dd_pi2::Float64, dd_del2::Float64, sp_sig2::Float64, sd_sig2::Float64, pd_sig2::Float64, pd_pi2::Float64,
            a_ss_sig2::Float64, a_pp_sig2::Float64, a_pp_pi2::Float64, a_dd_sig2::Float64, a_dd_pi2::Float64, a_dd_del2::Float64, a_sp_sig2::Float64, a_sd_sig2::Float64, a_pd_sig2::Float64, a_pd_pi2::Float64,
            R0::Float64, R::Array{Float64,1})
    function getHoppingParameters(R::Array{Float64,1})
        r = norm(R)
        if r < (1/2)^0.5 * R0
            ss_sig = ss_sig1 * exp(-a_ss_sig1 * (r - R0))
            pp_sig = pp_sig1 * exp(-a_pp_sig1 * (r - R0))
            pp_pi = pp_pi1 * exp(-a_pp_pi1 * (r - R0))
            dd_sig = dd_sig1 * exp(-a_dd_sig1 * (r - R0))
            dd_pi = dd_pi1 * exp(-a_dd_pi1 * (r - R0))
            dd_del = dd_del1 * exp(-a_dd_del1 * (r - R0))
            sp_sig = sp_sig1 * exp(-a_sp_sig1 * (r - R0))
            sd_sig = sd_sig1 * exp(-a_sd_sig1 * (r - R0))
            pd_sig = pd_sig1 * exp(-a_pd_sig1 * (r - R0))
            pd_pi = pd_pi1 * exp(-a_pd_pi1 * (r - R0))
        else
            ss_sig = ss_sig2 * exp(-a_ss_sig2 * (r - R0))
            pp_sig = pp_sig2 * exp(-a_pp_sig2 * (r - R0))
            pp_pi = pp_pi2 * exp(-a_pp_pi2 * (r - R0))
            dd_sig = dd_sig2 * exp(-a_dd_sig2 * (r - R0))
            dd_pi = dd_pi2 * exp(-a_dd_pi2 * (r - R0))
            dd_del = dd_del2 * exp(-a_dd_del2 * (r - R0))
            sp_sig = sp_sig2 * exp(-a_sp_sig2 * (r - R0))
            sd_sig = sd_sig2 * exp(-a_sd_sig2 * (r - R0))
            pd_sig = pd_sig2 * exp(-a_pd_sig2 * (r - R0))
            pd_pi = pd_pi2 * exp(-a_pd_pi2 * (r - R0))
        end
        return sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del
    end
    Vr=zeros(9,9)
    r=norm(R)
    l,m,n=R/r
    sp_sig,ss_sig,pp_sig,pp_pi,sd_sig,pd_sig,pd_pi,dd_sig,dd_pi,dd_del=getHoppingParameters(norm(R))
    #s and d interaction
    #Same orbital interaction
    Vr[0,0]=ss_sig
    Vr[1,1]=( 3*l^2*m^2*dd_sig + ( l^2 + m^2 - 4*l^2*m^2)*dd_pi +(n^2+l^2*m^2)*dd_del)
    Vr[2,2]=( 3*l^2*n^2*dd_sig + ( l^2 + n^2 - 4*l^2*n^2)*dd_pi +(m^2+l^2*n^2)*dd_del)
    Vr[3,3]=( 3*m^2*n^2*dd_sig + ( m^2 + n^2 - 4*m^2*n^2)*dd_pi +(l^2+m^2*n^2)*dd_del)
    Vr[4,4]=(3/4*(l^2-m^2)^2*dd_sig+(l^2+m^2-(l^2-m^2)^2)*dd_pi+(n^2+(l^2-m^2)^2/4)*dd_del)
    Vr[5,5]=((n^2-.5*(l^2+m^2))^2*dd_sig+3*n^2*(l^2+m^2)*dd_pi +3/4*(l^2+m^2)^2*dd_del)
    Vr[6,6]=(l^2*pp_sig+(1-l^2)*pp_pi)
    Vr[7,7]=(m^2*pp_sig+(1-m^2)*pp_pi)
    Vr[8,8]=(n^2*pp_sig+(1-n^2)*pp_pi)
    #p-p interaction
    Vr[6,7]=(l*m*pp_sig-l*m*pp_pi)
    Vr[6,8]=(l*n*pp_sig-l*n*pp_pi)
    Vr[7,8]=(m*n*pp_sig-m*n*pp_pi)
    #dd interaction
    Vr[1,2]=(3*l^2*m*n*dd_sig+m*n*(1-4*l^2)*dd_pi+m*n*(l^2-1)*dd_del) #xy-xz
    Vr[1,3]=(3*l*m^2*n*dd_sig+l*n*(1-4*m^2)*dd_pi+l*n*(m^2-1)*dd_del) #xy-yz
    Vr[1,4]=(1.5*l*m*(l^2-m^2)*dd_sig+2*l*m*(m^2-l^2)*dd_pi+.5*l*m*(l^2-m^2)*dd_del) #xy - x2-y2
    Vr[1,5]=(3^.5*l*m*(n^2-.5*(l^2+m^2))*dd_sig-2*3^.5*l*m*n^2*dd_pi+.5*3^.5*l*m*(1+n^2)*dd_del) #xy - 3z^2-r^2
    Vr[2,3]=(3*n^2*m*l*dd_sig+m*l*(1-4*n^2)*dd_pi+l*m*(n^2-1)*dd_del) #xz-yz
    Vr[2,4]=(1.5*n*l*(l^2-m^2)*dd_sig+n*l*(1-2*(l^2-m^2))*dd_pi-n*l*(1-.5*(l^2-m^2))*dd_del) #xz ->x^2-y^2
    Vr[2,5]=(3^.5*l*n*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*l*n*(l^2+m^2-n^2)*dd_pi-.5*3^.5*l*n*(l^2+m^2)*dd_del) #xz ->z^2
    Vr[3,4]=(1.5*m*n*(l^2-m^2)*dd_sig-m*n*(1+2*(l^2-m^2))*dd_pi+m*n*(1+(l^2-m^2)/2)*dd_del) #yz ->x^2-y^2
    Vr[3,5]=(3^.5*m*n*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*m*n*(l^2+m^2-n^2)*dd_pi-.5*3^.5*m*n*(l^2+m^2)*dd_del)#yz->z^2
    Vr[4,5]=(.5*3^.5*(l^2-m^2)*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*n^2*(m^2-l^2)*dd_pi+3^.5*(1+n^2)*(l^2-m^2)/4*dd_del)
    #s-d interaction
    Vr[0,1]=3^.5 * l*m*sd_sig
    Vr[0,2]=3^.5*l*n*sd_sig
    Vr[0,3]=3^.5*n*m*sd_sig
    Vr[0,4]=3^.5/2*(l^2-m^2)*sd_sig
    Vr[0,5]=(n^2-.5*(l^2+m^2))*sd_sig
    #sp interaction
    Vr[0,6]=(l*sp_sig)
    Vr[0,7]=(m*sp_sig)
    Vr[0,8]=(n*sp_sig)
    #pd interaction
    Vr[6,1]=(3^.5*l^2*m*pd_sig+m*(1-2*l^2)*pd_pi)#x->xy 
    Vr[7,1]=(3^0.5*m^2*l*pd_sig+l*(1-2*m^2)*pd_pi)#y->xy
    Vr[8,1]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi)#z->xy
    Vr[6,2]=(3^.5*l^2*n*pd_sig+n*(1-2*l^2)*pd_pi) #x->xz
    Vr[7,2]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #xz->y
    Vr[8,2]=(3^0.5*n^2*l*pd_sig+l*(1-2*n^2)*pd_pi) #xz->z
    Vr[6,3]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #yz->x
    Vr[7,3]=(3^0.5*m^2*n*pd_sig+n*(1-2*m^2)*pd_pi) #yz->y
    Vr[8,3]=(3^0.5*n^2*m*pd_sig+m*(1-2*n^2)*pd_pi) #yz->z
    Vr[6,4]=(3^0.5/2*l*(l^2-m^2)*pd_sig+l*(1-l^2+m^2)*pd_pi) #x^2-y^2->x
    Vr[7,4]=(3^0.5/2*m*(l^2-m^2)*pd_sig-m*(1+l^2-m^2)*pd_pi) #x^2-y^2->y
    Vr[8,4]=(3^0.5/2*n*(l^2-m^2)*pd_sig-n*(l^2-m^2)*pd_pi) #x^2-y^2->z
    Vr[6,5]=(l*(n^2-(l^2+m^2)/2)*pd_sig-3^0.5*l*n^2*pd_pi) #3z^2-r^2->x
    Vr[7,5]=(m*(n^2-(l^2+m^2)/2)*pd_sig-3^0.5*m*n^2*pd_pi) #3z^2-r^2->y
    Vr[8,5]=(n*(n^2-(l^2+m^2)/2)*pd_sig+3^0.5*n*(l^2+m^2)*pd_pi)#3z^2-r^2->z
    #For those hopping not listed in SK table
    #Compute the opposite hopping and use Hermiticity of Hamiltonian to 
    #derive these elements. 
    #Now flip the direction:
    l,m,n=-R/norm(R)
    #dd interaction
    Vr[2,1]=(3*l^2*m*n*dd_sig+m*n*(1-4*l^2)*dd_pi+m*n*(l^2-1)*dd_del) #xy-xz
    Vr[3,1]=(3*l*m^2*n*dd_sig+l*n*(1-4*m^2)*dd_pi+l*n*(m^2-1)*dd_del) #xy-yz
    Vr[4,1]=(1.5*l*m*(l^2-m^2)*dd_sig+2*l*m*(m^2-l^2)*dd_pi+.5*l*m*(l^2-m^2)*dd_del) #xy - x2-y2
    Vr[5,1]=(3^.5*l*m*(n^2-.5*(l^2+m^2))*dd_sig-2*3^.5*l*m*n^2*dd_pi+.5*3^.5*l*m*(1+n^2)*dd_del) #xy - 3z^2-r^2
    Vr[3,2]=(3*n^2*m*l*dd_sig+m*l*(1-4*n^2)*dd_pi+l*m*(n^2-1)*dd_del) #xz-yz
    Vr[4,2]=(1.5*n*l*(l^2-m^2)*dd_sig+n*l*(1-2*(l^2-m^2))*dd_pi-n*l*(1-.5*(l^2-m^2))*dd_del) #xz ->x^2-y^2
    Vr[5,2]=(3^.5*l*n*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*l*n*(l^2+m^2-n^2)*dd_pi-.5*3^.5*l*n*(l^2+m^2)*dd_del) #xz ->z^2
    Vr[4,3]=(1.5*m*n*(l^2-m^2)*dd_sig-m*n*(1+2*(l^2-m^2))*dd_pi+m*n*(1+(l^2-m^2)/2)*dd_del) #yz ->x^2-y^2
    Vr[5,3]=(3^.5*m*n*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*m*n*(l^2+m^2-n^2)*dd_pi-.5*3^.5*m*n*(l^2+m^2)*dd_del)#yz->z^2
    Vr[5,4]=(.5*3^.5*(l^2-m^2)*(n^2-.5*(l^2+m^2))*dd_sig+3^.5*n^2*(m^2-l^2)*dd_pi+3^.5*(1+n^2)*(l^2-m^2)/4*dd_del)
    #p-p interaction
    Vr[7,6]=(l*m*pp_sig-l*m*pp_pi)
    Vr[8,6]=(l*n*pp_sig-l*n*pp_pi)
    Vr[8,7]=(m*n*pp_sig-m*n*pp_pi)
    #ds interaction
    Vr[1,0]=3^.5 * l*m*sd_sig
    Vr[2,0]=3^.5*l*n*sd_sig
    Vr[3,0]=3^.5*n*m*sd_sig  
    Vr[4,0]=3^.5/2*(l^2-m^2)*sd_sig
    Vr[5,0]=(n^2-.5*(l^2+m^2))*sd_sig
    #ps interaction
    Vr[6,0]=(l*sp_sig)
    Vr[7,0]=(m*sp_sig)
    Vr[8,0]=(n*sp_sig)
    #dp interaction
    Vr[1,6]=(3^.5*l^2*m*pd_sig+m*(1-2*l^2)*pd_pi)#x->xy 
    Vr[1,7]=(3^0.5*m^2*l*pd_sig+l*(1-2*m^2)*pd_pi)#y->xy
    Vr[1,8]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi)#z->xy
    Vr[2,6]=(3^.5*l^2*n*pd_sig+n*(1-2*l^2)*pd_pi) #x->xz
    Vr[2,7]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #xz->y
    Vr[2,8]=(3^0.5*n^2*l*pd_sig+l*(1-2*n^2)*pd_pi) #xz->z
    Vr[3,6]=(3^0.5*l*m*n*pd_sig-2*l*m*n*pd_pi) #yz->x
    Vr[3,7]=(3^0.5*m^2*n*pd_sig+n*(1-2*m^2)*pd_pi) #yz->y
    Vr[3,8]=(3^0.5*n^2*m*pd_sig+m*(1-2*n^2)*pd_pi) #yz->z
    Vr[4,6]=(3^0.5/2*l*(l^2-m^2)*pd_sig+l*(1-l^2+m^2)*pd_pi) #x^2-y^2->x
    Vr[4,7]=(3^0.5/2*m*(l^2-m^2)*pd_sig-m*(1+l^2-m^2)*pd_pi) #x^2-y^2->y
    Vr[4,8]=(3^0.5/2*n*(l^2-m^2)*pd_sig-n*(l^2-m^2)*pd_pi) #x^2-y^2->z
    Vr[5,6]=(l*(n^2-(l^2+m^2)/2)*pd_sig-3^0.5*l*n^2*pd_pi) #3z^2-r^2->x
    Vr[5,7]=(m*(n^2-(l^2+m^2)/2)*pd_sig-3^0.5*m*n^2*pd_pi) #3z^2-r^2->y
    Vr[5,8]=(n*(n^2-(l^2+m^2)/2)*pd_sig+3^0.5*n*(l^2+m^2)*pd_pi)#3z^2-r^2->z
    return Vr
end
struct Hamiltonian
    #=
    Hlist: a list of each Hamiltonian, if we partition the original graph to N subgraph,
    then Hlist has the size of N+1, where Hlist[:N] correspond to the Hamiltonian of all subgraphes, 
    Hlist[N] correspond to the hopping between subgraphes. 

    Hlist[i], the Hamiltonian of a particular part of the Nanoparticle, 

    It is multiplied with v in the following way, Hlist[i]*v[istart:iend]
    =#
    Hlist::Vector{SparseMatrixCSR{Float64, Int}}
    start_list::Vector{Int}
    end_list::Vector{Int}
end
function dot_product!(H::SparseMatrixCSR{Float64, Int}, v::Vector{Float64})
    return H * v
end
function dot(H::Hamiltonian, v::Vector{Float64})::Vector{Float64}
    # Extract sub-vectors based on the shape of sub-Hamiltonians
    sub_vectors = Vector{Vector{Float64}}()
    #=
    There is no need for this
    start = 1
    for i in 1:(length(H.Hlist) - 1)
        size = size(H.Hlist[i], 2)  # Number of columns corresponds to the size of sub-vector
        push!(sub_vectors, v[start:(start+size-1)])
        start += size
    end=#

    for i in 1:(length(H.Hlist) - 1)
        push!(sub_vectors, v[start_list[i]:end_list[i]])
    end

    # Use parallel computation for the dot products
    hopping_result = H.Hlist[end] * v
    Threads.@threads for i = 1:(length(H.Hlist) - 1)        
        ans_part=dot_product!(Hlist[i],sub_vectors[i])
        hopping_result[start_list[i]:end_list[i]]=hopping_result[start_list[i]:end_list[i]]+ans_part
    end
    #tasks = [Threads.@spawn dot_product!(H.Hlist[i], sub_vectors[i]) for i in 1:(length(H.Hlist) - 1)]
    #map(fetch, tasks)  # Wait for all tasks to complete
    return hopping_result
end
function compute_subgraph_hamiltonian(atoms, nl, start_idx::Int, end_idx::Int)
    rows, cols, values = Int[], Int[], Float64[]

    for i in start_idx:end_idx
        neighbors = nl.get_neighbors(i-1)[1]  # The indices of the neighboring atoms
        neighbors=neighbors.+1
        for neighbor in neighbors
            # Exclude interactions between different subgraphs
            if start_idx <= neighbor <= end_idx
                R = atoms.positions[neighbor, :] - atoms.positions[i, :]
                V_matrix = getV(R)
                for p in 1:9
                    for q in 1:9
                        push!(rows, (i-start_idx)*9 + p)
                        push!(cols, (neighbor-start_idx)*9 + q)
                        push!(values, V_matrix[p, q])
                    end
                end
            end
        end
    end

    return sparse(rows, cols, values)
end
function compute_inter_subgraph_hamiltonian(atoms, nl, Hsizes)
    rows, cols, values = Int[], Int[], Float64[]

    start_idx = 1
    for size in Hsizes
        end_idx = start_idx + size - 1

        for i in start_idx:end_idx
            neighbors = nl.get_neighbors(i-1)[1]
            neighbors=neighbors.+1
            for neighbor in neighbors
                # Only consider hopping between different subgraphs
                if !(start_idx <= neighbor <= end_idx)
                    R = atoms.positions[neighbor, :] - atoms.positions[i, :]
                    V_matrix = getV(R)
                    for p in 1:9
                        for q in 1:9
                            push!(rows, (i-1)*9 + p)
                            push!(cols, (neighbor-1)*9 + q)
                            push!(values, V_matrix[p, q])
                        end
                    end
                end
            end
        end

        start_idx += size
    end

    return sparse(rows, cols, values)
end
function ConstructHamiltonian(atoms::PyObject,cutoffs,TB_params::Dict;ncores=4)
    #setup the neighborlist 
    nl=Neighborlist(cutoffs,skin=0,self_interaction=false,bothways=true)#Need to set the skin etc. 
    nl.update(atoms)
    #Construct the sparse neighboring matrix. 
    rows=[]
    cols=[]
    vals=[]
    for iatom in 0:(length(atoms)-1)#Remember that atoms is a python object, so it is 0-indexed
        neighbors, _ = nl.get_neighbors(iatom)
        for jatom in neighbors
            push!(rows,iatom+1)#Add 1 to recover 1-index
            push!(cols,jatom+1)#Add 1 to recover 1-index 
            push!(vals,1)
        end
    end
    S=sparse(rows,cols,vals)
    G=Metis.graph(S)
    partition=Metis.partition(G,ncores)
    sorted_indices = sortperm(partition)
    new_atoms=atoms[sorted_indices.-1]
    group_counts = StatsBase.countmap(partition)
    nl.update(new_atoms)
    #Compute start_list and end_list
    start_list=[1]
    end_list=[]
    for i=1:ncores 
        Natom=group_counts[i]
        end_list[i]=start_list+Natom-1
        if i<ncores 
            push!(start_list,start_list[i-1]+Natom)
        end
    end
    #now compute the list of Hamiltonian 
    Hlist=[]
    for k=1:lastindex(start_list)
        push!(Hlist,compute_subgraph_hamiltonian!(new_atoms,nl,start_list[k],end_list[k]))
    end
    push!(Hlist,compute_inter_subgraph_hamiltonian(new_atoms,nl,end_list[end]))
    return Hamiltonian(Hlist,start_list,end_list)
end
end


