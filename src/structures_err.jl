# ----- Error Type -----
struct RelativeError{N,T,D1,D2}
    true_soln::NTuple{N,D1}
    test_soln::NTuple{N,D2}
    t_range::T
    error::NTuple{N,Array{Float64,2}}
    # Constructors
    function RelativeError(true_solns::Tuple{Vararg{D1}}, test_solns::Tuple{Vararg{D2}}) where {D1<:SSAdata,D2<:MCMdata}
        true_totals = Tuple(dat.D for dat in true_solns)
        test_totals = Tuple(get_total(dat) for dat in test_solns)
        error = Tuple((test_totals[i] .- true_totals[i]) ./ true_totals[i] for i in 1:length(true_totals))
        tr = test_solns[1].t_range
        new{length(true_solns),typeof(tr),D1,D2}(true_solns, test_solns, tr, error)
    end
    function RelativeError(true_soln::D1, test_soln::D2) where {D1<:SSAdata,D2<:MCMdata} 
        true_total = true_soln.D
        test_total = get_total(test_soln)
        error = ((test_total .- true_total)./ true_total,)
        tr = true_soln.t_range
        new{1, typeof(tr), D1, D2}((true_soln,), (test_soln,), tr, error)
    end
    RelativeError(tests::Tuple{Vararg{D2}}, trues::Tuple{Vararg{D1}}) where {D1<:SSAdata,D2<:MCMdata} = RelativeError(trues, tests)
    RelativeError(tst::D2, tru::D1) where {D1<:SSAdata, D2<:MCMdata} = RelativeError(tru,tst)
end

Base.length(::RelativeError{N,T,D1,D2}) where {N,T,D1,D2} = N
rep_no(R::RelativeError) = R.true_soln[1].n