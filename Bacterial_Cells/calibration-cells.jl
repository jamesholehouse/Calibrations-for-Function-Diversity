# import the relevant packages
using DataFrames, CSV, BlackBoxOptim, Distances, StatsBase, Printf, DelimitedFiles, LinearAlgebra, Interpolations
include("<path to diversity simulator>")
using .dsim

data_ind = ARGS[1];
# data_ind = 1 <-- using cell_rf_data_1.csv

"""
Set-up: import and select data.
"""

# import the cell data
read_data = readdlm("./data/cell_rf_data_$(data_ind).csv");
import_data = [deleteat!(read_data[i,:],1) for i in 1:length(read_data[:,1])];
[filter!(e -> e ≠ "",import_data[i]) for i in 1:length(read_data[:,1])];

Ns = floor.(sum.(import_data));
Ds = length.(import_data);
    
# pick out the smallest cell from the entire list
all_data_raw = readdlm("./data/cell_rf_data.csv");
all_data = [deleteat!(all_data_raw[i,:],1) for i in 1:length(all_data_raw[:,1])]
[filter!(e -> e ≠ "",all_data[i]) for i in 1:length(all_data[:,1])]
small_cell = filter(e->e!=0.0,map(e->Int64(floor(e)), all_data[21])); 
# min size at 21 which is argmin(Ns)

println("Data imported.")
    
"""
Set up parallel prog
"""
# take the number of threads
nth = Threads.nthreads()
# find the number of sims to run on each thread, try and equate
runs_per_th = [floor(Int,length(import_data)/nth) for i in 1:nth]
if sum(runs_per_th) != length(import_data)+1
    remn = length(import_data) % sum(runs_per_th)
    for i in 1:remn
        ind = i%length(runs_per_th)
        runs_per_th[ind] += 1
    end
end
            
println("Parallel prog set-up.")
            

"""
The actual calibration: set up BBO and run.
"""

# method to calculate the transformed diversity
Dn(θ,cD,cN,starD,starN,N) = cD + (starD - cD)*(N^((1-θ)/(2-θ))-cN^((1-θ)/(2-θ)))/(starN^((1-θ)/(2-θ))-cN^((1-θ)/(2-θ)));

# padding function
function pad_rdist(d1, d2)
    min_not_0 = minimum(filter!(x->x!=0,vcat(d1, d2))); # so that the padded values are not exactly zero
    n = max(length(d1), length(d2));
    d1_padded = vcat(d1, ones(n - length(d1))*min_not_0/10); # divide by 10, so as not to punish the distance too badly (could make things too stochastic).
    d2_padded = vcat(d2, ones(n - length(d2))*min_not_0/10);
    return d1_padded, d2_padded, abs(length(d1) - length(d2));
end;

# the vectors that will store the required values of θ, γ
θ_vec_all = Vector{Union{Missing, Float64}}(undef, length(import_data)); θ_vec_all .= 0.0;
γ_vec_all = Vector{Union{Missing, Float64}}(undef, length(import_data)); γ_vec_all .= 0.0;
fitness_all = Vector{Union{Missing, Float64}}(undef, length(import_data)); fitness_all .= 0.0;
            
println("Storage vector initialized.")
            

# Loop over each agency in the selected part of the data, and store the parameters found via the optimisation into the vectors.
# Each thread runs a set number of calibrations specified in runs_per_th
Threads.@threads for iter in 1:nth
    println("Running thread $iter.")
    for j in 1:runs_per_th[iter]
        # find the index of the initial rf data array
        if iter == 1
            i = j
        else
            i = sum(runs_per_th[1:iter-1])+j
        end
        
        # Set up some things for the eventual calibration
        data_rf = import_data[i]
        initial_val_d = length(small_cell)
        initial_val_n = sum(small_cell)
                        
        # if the cell is the smallest one used as the initial condition then output is missing.
        if sum(data_rf) == 26925
            θ_vec_all[i] = missing; γ_vec_all[i] = missing; fitness_all[i] = missing;
            break
        end

        function ibb_fn(S::Array{Float64,1})
            # assign values of model pars
            θ, γ = S[1:2];
            # define the maximum value of simulated org size
            reduced_n = 5e5;
            # find out if the size is greater than the "reduced size" (if > red. size then Qred=1)
            Qred = 0; if Ns[i]>5e5 Qred = 1 end

            if Qred == 1
                # first need to back track the diversity
                chosen_div = length(data_rf);
                chosen_size = Ns[i];
                vals = [θ,initial_val_d,initial_val_n,chosen_div,chosen_size];
                reduced_d = Dn(vals...,reduced_n);

                # find the relative ranks of the full to use in the reduced
                chosen_rel_ranks = collect(0:1/chosen_div:1-1/chosen_div);
                mod_vec = zeros(floor(Int,reduced_d));
                mod_rel_ranks = collect(0:1/length(mod_vec):1-1/length(mod_vec))
                mod_closest_ind = [findmin(abs.(chosen_rel_ranks .- b))[2] for b in mod_rel_ranks]
                mod_vec = filter!(x->x!=0.0,round.(reduced_n .* data_rf[mod_closest_ind] ./sum(data_rf[mod_closest_ind])));

                # now simulate the organization with the proposed params
                osim = dsim.create_cell(θ,γ,35.9,convert(Int64,reduced_n), small_cell, :pdf)
                mod_vec_norm = mod_vec/sum(mod_vec)
                osim_norm = osim/sum(osim)

                # need to pad the distributions 
                mv_norm_pad, osim_norm_pad = pad_rdist(mod_vec_norm,osim_norm)
                # find the logarithm of elements for the distance calculation
                mv_norm_pad_log, osim_norm_pad_log = log.(mv_norm_pad), log.(osim_norm_pad)

                # calculate the distance
                euc_log = euclidean(mv_norm_pad_log, osim_norm_pad_log)

                return euc_log
            else #Qred = 0
                # simulate the org with proposed params
                osim = dsim.create_cell(θ,γ,35.9,convert(Int64,Ns[i]), small_cell, :pdf)
                data_rf_norm = data_rf/sum(data_rf)
                osim_norm = osim/sum(osim)

                # need to pad the distributions 
                data_rf_norm_pad, osim_norm_pad = pad_rdist(data_rf_norm,osim_norm)
                # find the logarithm of elements for the distance calculation
                data_rf_norm_pad_log, osim_norm_pad_log = log.(data_rf_norm_pad), log.(osim_norm_pad)

                # calculate the distances
                euc_log = euclidean(data_rf_norm_pad_log, osim_norm_pad_log)
            end
        end
        
        println("Starting optimizations.")

        # run the opt
        iSRange = [(0.0, 1.2),(0.0,1.2)];
        iopts = bbsetup(ibb_fn; Method = :adaptive_de_rand_1_bin_radiuslimited, SearchRange = iSRange, 
                    NumDimensions = 2, MaxTime = 86400.0, PopulationSize = 50, TraceMode=:silent); # 1 day: 86400.0
        ires = bboptimize(iopts);
        iopt_par = best_candidate(ires);
        ifitness = best_fitness(ires);

        # store the pars
        θ_vec_all[i] = iopt_par[1]; γ_vec_all[i] = iopt_par[2]; fitness_all[i] = ifitness;
        @printf "Calibration %d out of %d now done, giving parameters θ = %.3f, γ = %.3f with fitness %.3f." i length(import_data) θ_vec_all[i] γ_vec_all[i] ifitness
        println("")
    end
end


df_out = DataFrame(theta=θ_vec_all, gamma=γ_vec_all, fitness = fitness_all);
CSV.write("./data/cell-calibs-results-$(data_ind).csv", df_out)
pwd()
