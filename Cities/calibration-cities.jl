# import the relevant packages
using DataFrames, CSV, BlackBoxOptim, Distances, StatsBase, Printf, DelimitedFiles, LinearAlgebra, Interpolations
include("<path to diversity simulator>")
using .dsim

"""
Set-up: import and select data.
"""

# import the big city data
raw_data = DataFrame(CSV.File("./data/largest-cities.csv"));
# remove missing values
c_data = disallowmissing!(raw_data[completecases(raw_data), :]);
each_met_area = groupby(c_data, :AREA);
# the biggest 20 cities
sel_met_workers = []
for (i,e) in enumerate(each_met_area)
    push!(sel_met_workers,sort(e[!,:TOT_EMP],rev=true));
end
sel_met_workers = map(x->floor.(Int,x),sel_met_workers);

# find pop of each of the selected cities.
workers = [sum(x) for x in sel_met_workers];
import_data = sel_met_workers./20;# cities reduced by 20.
Ns = map(x->floor.(Int,x./20),workers);# numbers and diversities
Ds = length.(import_data);

# import the Waco city data.
read_w = DataFrame(CSV.File("/data/waco-data.csv"));
waco_vec = sort(read_w[!,:TOT_EMP],rev=true);
small_city = waco_vec./20;

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

# the vectors that will store the required values of θ, γ
θ_vec_all = Vector{Union{Missing, Float64}}(undef, length(import_data)); θ_vec_all .= 0.0;
γ_vec_all = Vector{Union{Missing, Float64}}(undef, length(import_data)); γ_vec_all .= 0.0;
fitness_all = Vector{Union{Missing, Float64}}(undef, length(import_data)); fitness_all .= 0.0;

# padding function
function pad_rdist(d1, d2)
    min_not_0 = minimum(filter!(x->x!=0,vcat(d1, d2))); # so that the padded values are not exactly zero
    n = max(length(d1), length(d2));
    d1_padded = vcat(d1, ones(n - length(d1))*min_not_0/10); # divide by 10, so as not to punish the distance too badly (could make things too stochastic).
    d2_padded = vcat(d2, ones(n - length(d2))*min_not_0/10);
    return d1_padded, d2_padded, abs(length(d1) - length(d2));
end;
            
println("Storage vectors initialized.")

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
        
        # Set up the data interpolation for comparsion of HD in simulation results.
        data_rf = import_data[i]
        
        function ibb_fn(S::Array{Float64,1})
            # assign values of model pars
            θ, γ = S[1:2];
            # simulate the organization and collect the norm rf distribution
            osim = dsim.create_cell(θ,γ,111.869,convert(Int64,Ns[i]), small_city, :pdf)
            df_norm = data_rf/sum(data_rf)
            osim_norm = osim/sum(osim)

            # need to pad the distributions 
            df_norm_pad, osim_norm_pad = pad_rdist(df_norm,osim_norm)
            # find the logarithm of elements for the distance calculation
            df_norm_pad_log, osim_norm_pad_log = log.(df_norm_pad), log.(osim_norm_pad)
            euc_log = euclidean(df_norm_pad_log, osim_norm_pad_log)
            return euc_log
        end
        
        println("Starting optimizations.")

        # run the opt
        iSRange = [(0.8,1.5),(0.5,1.2)];
        iopts = bbsetup(ibb_fn; Method = :adaptive_de_rand_1_bin_radiuslimited, SearchRange = iSRange, 
                    NumDimensions = 2, MaxTime = 86400.0, PopulationSize = 50, TraceMode=:silent); # 1 day: 86400.0
        ires = bboptimize(iopts);
        iopt_par = best_candidate(ires);
        ifitness = best_fitness(ires);

        # store the pars
        θ_vec_all[i] = iopt_par[1]; γ_vec_all[i] = iopt_par[2];fitness_all[i] = ifitness
        @printf "Calibration %d out of %d now done, giving parameters θ = %.3f, γ = %.3f with fitness %.3f." i length(import_data) θ_vec_all[i] γ_vec_all[i] ifitness
        println("")
    end
end

df_out = DataFrame(theta=θ_vec_all, gamma=γ_vec_all, p0 = p0_vec_all, fitness = fitness_all);
CSV.write("./calibration_results/cities-calibs-results.csv", df_out)
pwd()
