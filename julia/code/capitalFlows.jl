#import Pkg; Pkg.add("Tables", "XLSX", "ExcelReaders", "DataFrames",
#       "JuMP", "Ipopt", "NamedArrays")

include("concordance.jl")


using DataFrames, JuMP, Ipopt, DelimitedFiles;

#==============================================================================
Aggregate US capital flows table to Aus 19 sectors
==============================================================================#
flows97 = ExcelReaders.readxlsheet("data"*pathmark*"flow1997.xls", 
    "180x22Combined");
flows = DataFrame(flows97[4:182, 4:25], :auto);
dataTypeEx = flows[41,4]
for i in [1:1:179;]
    for j in [1:1:22;]
        if typeof(flows[i,j]) == typeof(dataTypeEx)
            flows[i,j] = Float64(0.0)
        end
    end
end

# Aggregating to 20 sector
# Combining all manufacturing into the same sector
flows.x5 = flows.x5+flows.x6+flows.x7;
flows = select!(flows, Not(:x6));
flows = select!(flows, Not(:x7));
# Combining Transportation and Warehouses into the one sector
flows.x10 = flows.x10 + flows.x11;
flows = select!(flows, Not(:x11));
# Combining "Professional and technical services" and "Management of companies 
# and enterprises" into the one sector
flows.x15 = flows.x15 + flows.x16;
flows = select!(flows, Not(:x16));

#print(flows)

# Renaming to Aus 19 Sectors codes
#= Because we begin with the 180x22 flows table, and combined a few above, we 
have 18 industries. This below is simply mapping these from the 18 combined 
industries into the Aus 19 industries categories (public sector is missing). 
This is done below by hand since it is so few sectors. It is effectively the 
concordance. The numbers "x1" etc correlate to the numberings of these 22 
industries of the BEA data, so the mapping below is meaningful =#
rename!(flows, :x1 => :A)
rename!(flows, :x2 => :B)
rename!(flows, :x3 => :D)
rename!(flows, :x4 => :E)
rename!(flows, :x5 => :C)
rename!(flows, :x8 => :F)
rename!(flows, :x9 => :G)
rename!(flows, :x10 => :I)
rename!(flows, :x12 => :J)
rename!(flows, :x13 => :K)
rename!(flows, :x14 => :L)
rename!(flows, :x15 => :M)
rename!(flows, :x17 => :N)
rename!(flows, :x18 => :P)
rename!(flows, :x19 => :Q)
rename!(flows, :x20 => :R)
rename!(flows, :x21 => :H)
rename!(flows, :x22 => :S)

rowCode = flows97[4:182,2];
containsSpace = findall( x -> occursin(" ", x), rowCode);
rowCode[containsSpace] = replace.(rowCode[containsSpace], " " => "");
rowCode20 = Array{Union{Nothing, String}}(nothing, length(rowCode));
for i in eachindex(rowCode);
    rowCode20[i] = Comm180To20[rowCode[i]];
end
insertcols!(flows ,1, :Industry => rowCode20);
splitIndustry = groupby(flows, :Industry);
flows = combine(splitIndustry, valuecols(splitIndustry) .=> sum);

rename!(flows, :A_sum => :A)
rename!(flows, :B_sum => :B)
rename!(flows, :C_sum => :C)
rename!(flows, :D_sum => :D)
rename!(flows, :E_sum => :E)
rename!(flows, :F_sum => :F)
rename!(flows, :G_sum => :G)
rename!(flows, :H_sum => :H)
rename!(flows, :I_sum => :I)
rename!(flows, :J_sum => :J)
rename!(flows, :K_sum => :K)
rename!(flows, :L_sum => :L)
rename!(flows, :M_sum => :M)
rename!(flows, :N_sum => :N)
rename!(flows, :P_sum => :P)
rename!(flows, :Q_sum => :Q)
rename!(flows, :R_sum => :R)
rename!(flows, :S_sum => :S)
#------------------------------------------------------------------------------
# what year are we importing?
#------------------------------------------------------------------------------
fyend = "2019"
#==============================================================================
Aggregate GFCF data
==============================================================================#

# Bring in GFCF data from excel
ausGFCFall = ExcelReaders.readxlsheet(
  "data"*pathmark*"5204064_GFCF_By_Industry_Asset.xls", "Data1"
);
ausGFCFall[1, 1] = "top corner";
# remove columns
ausGFCFall = ausGFCFall[:,
  Not(findall(x -> occursin("ALL INDUSTRIES ;", x), string.(ausGFCFall[1,:])))
];
# select current totals
ausGFCFcurrent = ausGFCFall[:,
  findall(x -> 
          occursin("Gross fixed capital formation: Current", x)
          | occursin("Dwellings: Current", x)
          | occursin("transfer costs: Current", x),
    string.(ausGFCFall[1,:]))
];
# Select dates column
ausGFCFdates = ausGFCFall[:,1];
ausGFCFcurrent = hcat(ausGFCFdates, ausGFCFcurrent);
ausGFCFrow= ausGFCFcurrent[
  findall(x -> occursin(fyend, x), string.(ausGFCFcurrent[:,1])), :];
# add ownership transfer costs to dwellings
ausGFCFrow[length(ausGFCFrow) - 1] = (ausGFCFrow[length(ausGFCFrow) - 1]
                                      + ausGFCFrow[length(ausGFCFrow)]);
# remove date column and ownership transfer costs
ausGFCFrow = ausGFCFrow[Not(1, length(ausGFCFrow))];

#==============================================================================
wrangle the IO table
==============================================================================#
# Import IO data
table8raw = DataFrame(CSV.File("data"*pathmark*fyend*pathmark*"table8.csv",
                              header=false));
table5raw = DataFrame(CSV.File("data"*pathmark*fyend*pathmark*"table5.csv",
                              header=false));
# check the two tables have the same year
commonwealthrow8 = findfirst(x -> occursin("Commonwealth", x), 
                            string.(table8raw[:, 2]));
commonwealthrow5 = findfirst(x -> occursin("Commonwealth", x), 
                            string.(table5raw[:, 2]));
# testing
(table8raw[commonwealthrow8, 2] != table5raw[commonwealthrow5, 2]
 && println("WARNING: io table release years don't match!"))
# find and set the number of sectors
(occursin("111 INDUSTRIES", string(table8raw[1:6, 1]))
  ? numioig = 111
  : numioig = 114)
# and the number of columns
numcol = numioig + 12;
#==============================================================================
a function for creating complete tables
==============================================================================#
function makeioig(table)
#table = table5raw
  # find the title row index
  titlerow = findfirst(x -> occursin("USE", x), string.(table[:,2]));
  # standardise the table
  table = table[:, 1:numcol];
  # make column titles
  coltitles = collect(values(table[titlerow, :]))
  coltitles[1] = "IOIG"
  coltitles[2] = "industry"
  morecoltitles = collect(values(table[titlerow + 2, numioig + 3 : numioig + 12]))
  coltitles[numioig + 3 : numioig + 12] = morecoltitles
  coltitles[ismissing.(coltitles)] .= "0"
  coltitles = filter.(x -> !isspace(x), coltitles)
  for i in [1:1:numioig;]
    coltitles[i + 2] = "IOIG"*values(coltitles[i + 2])
  end
  rename!(table, string.(coltitles), makeunique=true)
  # remove initial rows and tidy up
  table = table[Not(range(1, titlerow + 2)), :]
  table = dropmissing(table, :industry)
  table = table[Not(findall(x -> occursin("Commonwealth", x),
                            string.(table.industry))),
                :]
  (findall(occursin.("Total uses", string.(table[:, 2]))) == Int64[]
    ? (AProw = findfirst(occursin.("Australian Production",
                                  string.(table[:, 2])));
      table[AProw, 1] = "TU8AP5";
      GDProw = findfirst(occursin.("GDP",
                                  string.(table[:, 1])));
      table = table[Not(GDProw), :];)
    : (TUrow = findfirst(occursin.("Total uses", string.(table[:, 2])));
       table[TUrow, 1] = "TU8AP5")
  )
  GVArow = findfirst(occursin.("Value Added", string.(table[:, 2])));
  table[GVArow, 1] = "GVA";
  table = dropmissing(table, :IOIG)
  table.IOIG = string.(table.IOIG)
  numrow = nrow(table)
  table[numioig + 1 : numrow, 1] = "`".*table[numioig + 1 : numrow, 1]
  # convert tables to numbers
  table[:, 3:numcol] = filter.(x -> !isspace(x), table[:, 3:numcol])
  table[!, 3:numcol] = parse.(Float64, table[!, 3:numcol])
  sort!(table)
  return table
end
table8 = makeioig(table8raw)
table5 = makeioig(table5raw)

#==============================================================================
 take difference of the two tables
==============================================================================#
function makediff(t8, t5)
  (ncol(t8) < 100
     ? (firstcol = 2; numsec = ncol(t8) - 11)
     : (firstcol = 3; numsec = ncol(t8) - 12)
  )
  # make diff
  (tdiff = deepcopy(t8);
    for i in range(1, nrow(t8))
      for j in range(firstcol, ncol(t8))
        tdiff[i, j] = t8[i, j] - t5[i, j]
      end
    end
  )
  # identify negative values (expect intermediates to be positive)
  (negv = DataFrame(rowind = Int64[], colind = Int64[], val = Float64[]);
    for i in range(1, nrow(tdiff))
      for j in range(firstcol, ncol(tdiff))
        tdiff[i, j] < 0 && push!(negv, [i j tdiff[i, j]])
      end
    end;
   minimum(negv.val) < 0 &&
     println("See negvals for the negative values in
            the table of differences between table8 and table5");
  )
  # identify positive values in diff P6 (positive means re-export of imports)
  (rowP6 = findfirst(occursin.("P6", tdiff[:, 1]));
   posP6 = DataFrame(rowind = Int64[], colind = Int64[], val = Float64[]);
   for i in range(firstcol, numsec + 1)
     (tdiff[rowP6, i] > 0 && push!(posP6, [rowP6 i tdiff[rowP6, i]]))
   end;
   (maximum(posP6.val) > 0
    && println("See (anz)posdiffP6 for the positive values of competing imports
               in the table of differences between (anz)table8 and (anz)table5)"
              )
   );
  )
  (t8P6 = Vector(t8[findfirst(occursin.("P6", t8[:,1])),
                    firstcol: numsec + firstcol - 1]);
   tdiffTU8AP5 = Vector(tdiff[end, firstcol: numsec + firstcol - 1]);
   (t8P6 - tdiffTU8AP5 != zeros(length(t8P6))
    && println("WARNING: P6 of table8 is not equal to TU8AP5 of tablediff")  
   );
  )
 return tdiff, negv, posP6
end
(tablediff, negvals, posdiffP6) = makediff(table8, table5);
#==============================================================================
 transform table to 20 sectors
==============================================================================#
# map ioig to ANZSIC20
# the following function is impure in that it depends on dicts and functions ..
# in concordance.jl
function makeanztable(table, numsec=numioig)
  # create a many-to-one column: anzcode per ioig
  (# isolate ioig codes
   ioigcodes = string.(table.IOIG[1:numsec]);
   ioigcodesFloat = parse.(Float64, ioigcodes);
   ioigto20 = mapioig20(ioigcodesFloat);
   tmp = String[];
   for i in eachindex(ioigcodesFloat)
     push!(tmp, ioigto20[ioigcodesFloat[i]])
   end;
   # and the remaining terms in the table;
   for i in range(1, nrow(table) - numsec)
     push!(tmp, table[numsec + i, 1])
   end;
  anzperioigcol = tmp;
  )
  # collapse rows
  (anztable = table;
   insertcols!(anztable, 1, "ANZDIV" => anzperioigcol);
   anztable = (combine(groupby(anztable, :ANZDIV),
                        names(table, Between(:IOIG0101, :T6))
                          .=> sum .=> Between(:IOIG0101, :T6))
               );
   sort!(anztable);
  )
  # collapse cols
  (tmptable = transpose(Matrix(anztable[:, 2 : numioig + 1]));
   tmptable = DataFrame(tmptable, :auto);
   insertcols!(tmptable, 1, "ANZDIV" => anzperioigcol[1:numioig]);
   tmptable = combine(groupby(tmptable, :ANZDIV), Not(1) .=> sum);
   sort!(tmptable);
   tmptablet = transpose(Matrix(tmptable[:, Not(1)]));
   tmptablet = DataFrame(tmptablet, ANZSICDivByLetter);
   insertcols!(tmptablet, 1, "ANZDIV" => anztable[:,1]);
   anztable = anztable[:, Not(names(anztable, Between(:IOIG0101, :IOIG9502)))];
   anztable = innerjoin(tmptablet, anztable, on = :ANZDIV);
  )
  return anztable
end
anztable8 = makeanztable(table8)
anztable5 = makeanztable(table5)
(anztablediff, anznegvals, anzposdiffP6) = makediff(anztable8, anztable5)

#==============================================================================
ras tables
==============================================================================#
#make the io tables ras-ready
#table = anztable8
function makerr(table)
  gvarow = findfirst(occursin.("GVA", table.ANZDIV))
  t1row = findfirst(occursin.("T1", table.ANZDIV))
  table = table[Not(gvarow, t1row), :]
  table = table[:, Not([:T4, :T5])]
  table.Q3 += table.Q4 + table.Q5
  table = table[:, Not([:Q4, :Q5])]
  rename!(table, :Q3 => "Q345")
  return table
end 
anztable8rr = makerr(anztable8)
anztable5rr = makerr(anztable5)

# import local data
grxgrp = ExcelReaders.readxlsheet(
  "data"*pathmark*".xls", "Data1"
);
stop

# Isolate GFCF
Q3 = parse.(Float64, table8.Q3)
Q4 = parse.(Float64, table8.Q4)
Q5 = parse.(Float64, table8.Q5)

ioigGfcf = Q3 + Q4 + Q5;
ausprodrow = findfirst(x -> occursin("Australian Production", x),
                  string.(table8.industry));
T1row = findfirst(x -> occursin("T1", x),
                  string.(table8.IOIG));
ioigGfcftot = ioigGfcf[ausprodrow];
ioigGfcf = ioigGfcf[1 : numioig] / ioigGfcf[T1row] * ioigGfcftot;
ioigGfcf = DataFrame(:inv => ioigGfcf);
#
# Grouping GFCF receivable by 20 sectors (dwellings is all zero in capital)
insertcols!(ioigGfcf, 1, :anzcode => anzperioigcol);
splitIndustry = groupby(ioigGfcf, :anzcode);
anzdivgfcf = combine(splitIndustry, valuecols(splitIndustry) .=> sum);
sort!(anzdivgfcf);

# Balancing row and column sums to IO table 8 total
ausGFCFtot = sum(ausGFCFrow);
for i in eachindex(ausGFCFrow)
    ausGFCFrow[i] = ausGFCFrow[i] * ioigGfcftot / ausGFCFtot;
end;
# pull in proportionalised kapital flows to ras
y = DataFrame(CSV.File("data/propd-to-ras.csv", header=false))
y = Matrix(y)
# generate an initial table for the ras
#y = zeros(length(ausGFCFrow), length(anzdivgfcf.inv_sum))
#for i in eachindex(ausGFCFrow), j in eachindex(anzdivgfcf);
#  y[i, j] = anzdivgfcf[i] / ioigGfcftot * ausGFCFrow[j]
#end  

#==============================================================================
Make prior scaled from Aus Data
==============================================================================#
# Make vector of the proportion of each row sum element as a fraction of the 
# total
rowSumsProps = ones(20);
ausGFCFtot = sum(ausGFCFrow);
for i in 1:20
    rowSumsProps[i] = ausGFCFrow[i] / ausGFCFtot;
end
# Make prior
ausPropPrior = ones(20, 20);
for i in 1:20
    for j in 1:20
        ausPropPrior[i,j] = rowSumsProps[j] *  anzdivgfcf.inv_sum[i];
    end
end

#==============================================================================
Make prior from US Data
==============================================================================#
# Adding empty rows
push!(flows, ["A" zeros(1, ncol(flows) - 1)])
push!(flows, ["D" zeros(1, ncol(flows) - 1)])
push!(flows, ["H" zeros(1, ncol(flows) - 1)])
push!(flows, ["K" zeros(1, ncol(flows) - 1)])
push!(flows, ["N" zeros(1, ncol(flows) - 1)])
push!(flows, ["P" zeros(1, ncol(flows) - 1)])
push!(flows, ["Q" zeros(1, ncol(flows) - 1)])
push!(flows, ["R" zeros(1, ncol(flows) - 1)])
push!(flows, ["S" zeros(1, ncol(flows) - 1)])
# Sorting Rows by industry index
sort!(flows)

# Take column sum
flowsTemp = deepcopy(flows);
flowsColSum = sum(eachcol(select!(flowsTemp, Not(:Industry))));

# Adding dwellings column
flows[!, :T] = flowsColSum *rowSumsProps[20];

# Adding public admin column
flows[!, :O] = flowsColSum *rowSumsProps[15];

# Sorting columns
flows = (flows[!, [:Industry, :A, :B, :C, :D, :E, :F, :G, :H, :I, :J, :K, :L, 
    :M, :N, :O, :P, :Q, :R, :S, :T]])

# Scale to Aus Data
flowsTemp = deepcopy(flows);
flowsTemp = select!(flowsTemp, Not(:Industry));
flowsTempSum = sum(Matrix(flowsTemp));
print(flowsTemp)
for i in 1:ncol(flowsTemp)
    for j in 1:nrow(flowsTemp)
        flowsTemp[j,i] = flowsTemp[j,i] / flowsTempSum * ausGFCFtot;
    end
end
print(flowsTemp)
flowsTemp[!, :Industry] = flows.Industry
flows = deepcopy(flowsTemp)

# Re -Sorting columns
flows = (flows[!, [:Industry, :A, :B, :C, :D, :E, :F, :G, :H, :I, :J, :K, :L, 
    :M, :N, :O, :P, :Q, :R, :S, :T]])

#==============================================================================
RAS
==============================================================================#

# Begin Optimisation
modCap = Model(Ipopt.Optimizer);
# Should be equal dimensions
@variable(modCap, x[1:length(anzdivgfcf.inv_sum),
                    1:length(ausGFCFrow)] >= 0.0001);
# Max entropy objective (or min relative to uniform)
@NLobjective(modCap,
             Min,
             sum((x[i, j] - y[i, j]) ^ 2
                 for i in eachindex(anzdivgfcf.inv_sum),
                 j in eachindex(ausGFCFrow)
                )
            );

# Row-sums constraint - must be equal to the GFCF totals
for j in eachindex(ausGFCFrow)
    @constraint(modCap, sum(x[:,j]) <= ausGFCFrow[j] + 1);
    @constraint(modCap, sum(x[:,j]) >= ausGFCFrow[j] - 1);
end;
# Col-sums constraint - must be equal to the IO totals
for k in eachindex(anzdivgfcf.inv_sum)
    @constraint(modCap, sum(x[k,:]) <= anzdivgfcf.inv_sum[k] + 1);
    @constraint(modCap, sum(x[k,:]) >= anzdivgfcf.inv_sum[k] - 1);
end;
optimize!(modCap);

# Add titles etc. for ease of reading -J works but this is clumsy fix 
ySol2 = DataFrame(value.(x), ANZSICDivShort);
insertcols!(ySol2, 1, :Divisions => ANZSICDivShort);

# Export tables as CSV
CSV.write("data"*pathmark*"capitalFlowsRAS2.csv", ySol2);
CSV.write("data"*pathmark*"UScapitalFlows.csv", flows);
