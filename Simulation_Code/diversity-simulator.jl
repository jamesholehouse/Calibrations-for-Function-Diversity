module dsim
    using StatsBase
    export create_org, create_city

    """
    pN: returns the probability of creating a new function.
    args:
    1. θ, the PA strength for creating a new function.
    2. SV, the state vector containing the number of people in each function.
    3. p0, the proportionality constant for new function creation.
    """
    function pN(θ,SV,p0)
        return p0*(sum(SV.^θ))^(-1)
    end


    """
    create_org: Returns an organization state vector containing the number of workers in each function (normalised). The initial state assumes a single worker in a single role.
    args:
    1. θ, the PA strength for creating a new function.
    2. γ, the PA strength to join an existing function.
    3. M, the size of the organization when the simulation stops.
    """
    function create_org(θ,γ,M,output = :pdf)
        # define the initial state vector, 1 worker in 1 role.
        SV = [1]
        # R counts the number of roles.
        R = 1
        # run M-1 iterations of worker addition to the company (org starts with 1).
        for i in 2:M
            # first calculate if new function added or not, p0 is 1
            p = pN(θ,SV,1)
            # rand_p = 1 then new function added, else =0 and add to existing
            rand_p = sample([0,1],Weights([1-p,p]))
            if rand_p == 1
                # add function with 1 worker
                append!(SV,1)
                # update number of functions
                R += 1
            else
                # add worker to existing function, first sample that function
                chosen_fn = sample(collect(1:1:R),Weights(SV.^γ))
                SV[chosen_fn] += 1
            end
        end
        sort!(SV,rev=true)

        if output == :pdf
            return SV./sum(SV)
        elseif output == :cdf
            rf = SV./sum(SV)
            return [sum(rf[1:j]) for j in 1:length(rf)]
        else
            error("choice of output should be either :cdf or :pdf.")
        end
    end


    """
    create_city: Returns an organization state vector containing the number of workers in each function (normalised) for a set of times (and also returns the corresponding times).
    args:
    1. θ, the PA strength for creating a new function.
    2. γ, the PA strength to join an existing function.
    3. p0, the proportionality constant for new function creation.
    4. M, the size of the organization when the simulation stops. IN TERMS OF BASE UNITS.
    5. SV0, the initial config of a city, if unspecified then takes IC of org model. MUST BE NORMALIZED FOR BASE UNIT.
    """
    function create_city(θ,γ,p0,M::Int64, SV0=[1], output = :pdf)
        # define the initial state vector, 1 worker in 1 role.
        SV = copy(SV0)
        # R counts the number of roles.
        R = length(SV0)
        # N0 is the initial number of workers
        N0 = sum(SV0)
        # stop at the maximal time M (see args)

        # calculate the first probability of new addition
        p = pN(θ,SV,p0)
        # error if p > 1
        if p > 1 || p < 0
            error("p0 should be chosen so that p(N) is ∈ [0,1]. Change p0 or the initial city structure.")
        end

        # run M-(N0+1) iterations of worker addition to the company.
        for i in N0+1:M
            # first calculate if new function added or not
            p = pN(θ,SV,p0)
            # rand_p = 1 then new function added, else = 0 and add to existing
            rand_p = sample([0,1],Weights([1-p,p]))
            if rand_p == 1
                # add function with 1 base unit
                append!(SV,1)
                # update number of functions
                R += 1
            else
                # add worker to existing function, first sample that function
                chosen_fn = sample(collect(1:1:R),Weights(SV.^γ))
                SV[chosen_fn] += 1
            end
        end
        sort!(SV,rev=true)
        
        if output == :pdf
#             rf = SV./sum(SV)
#             SV = Nothing # save RAM here and re-assign
#             SV0 = Nothing
#             GC.gc() # garbage collection
            return SV#, rf
        elseif output == :cdf
            rf = SV./sum(SV)
            SV = Nothing # save RAM here and re-assign
            SV0 = Nothing
            GC.gc() # garbage collection
            return [sum(rf[1:j]) for j in 1:length(rf)]
        else
            error("choice of output should be either :cdf or :pdf.")
        end
    end

end