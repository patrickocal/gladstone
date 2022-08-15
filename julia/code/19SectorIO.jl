include("concordance.jl")

IOIGAs19=Array{Union{Nothing, String}}(nothing, length(IOIG));
for i in eachindex(IOIG);
    IOIGAs19[i] = IOIGTo19[IOIG[i]]
end

IOHalf = DataFrame(copy(transpose(IOSource[[4:117; 119; 121:126; 128], 3:116])),:auto);
insertcols!(IOHalf ,1, :Industry => IOIGAs19);
IOHalfSplitIndustry = groupby(IOHalf, :Industry);
IOHalf = combine(IOHalfSplitIndustry, valuecols(IOHalfSplitIndustry) .=> sum);

sort!(IOHalf)
IOHalf = select!(IOHalf, Not(:Industry));
IOHalf = transpose(Matrix(IOHalf[:,:]));

IO19 = DataFrame(copy([IOHalf[1:114, :] IOSource[4:117,117:126]]),:auto);
insertcols!(IO19 ,1, :Industry => IOIGAs19);
IO19SplitIndustry = groupby(IO19, :Industry);
IO19 = combine(IO19SplitIndustry, valuecols(IO19SplitIndustry) .=> sum);

FullIO19 = [IO19; [IOHalf[115:122,:] IOSource[[119; 121:126; 128], 117:126]]];


writedlm("data"*pathmark*"FullIO19Draft.csv", FullIO19, ',');