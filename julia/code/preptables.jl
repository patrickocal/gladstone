#import Pkg; Pkg.add("Tables", "XLSX", "ExcelReaders", "DataFrames",
#       "JuMP", "Ipopt", "NamedArrays")

include("concordance.jl")

#filepath cross system compatability code
if Sys.KERNEL === :NT || Sys.KERNEL === :Windows
    sep = "\\";
else
    sep = "/";
end

# how many sectors?
numdiv = length(ANZdiv);

(using DataFrames, JuMP, Ipopt, DelimitedFiles, Random, Formatting,
       PrettyTables, XLSX, ExcelReaders);

# present working directory is the parent folder "gladstone" or "maiwar"
datadir = "julia"*sep*"data"*sep
outputdir = "julia"*sep*"output"*sep
#=============================================================================#
# function to get list of row indices for given sectors etc. Used only as an 
# input to the getsubframe function at this stage
function getrows(dfColumn, lsOfSectors)
    rowVec = [];
    for i in 1:length(lsOfSectors)
        append!(rowVec, findall(x -> x == lsOfSectors[i], dfColumn)[1])
    end
    return rowVec
end
# general function to get a sub-frame based on row and column names
# makes the code for extracting each parameter much more readable
function getsubframe(df, lsRows, lsCols, divn=divhdr);
    subdf = df[getrows(df[:, divn], lsRows), append!([divn], lsCols)]
    return subdf
end
#==============================================================================
Aggregate US capital flows table to Aus 19 sectors
==============================================================================#
flows97 = ExcelReaders.readxlsheet(datadir*"flow1997.xls", 
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

# Aggregating to "numdiv" sectors
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
rowCodediv = repeat(["div"], inner=length(rowCode), outer=1);
for i in eachindex(rowCode);
    rowCodediv[i] = Comm180Todiv[rowCode[i]];
end
insertcols!(flows ,1, :ANZcode => rowCodediv);
splitIndustry = groupby(flows, :ANZcode);
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
wrangle the IO table
==============================================================================#
# Import IO data
path = datadir*fyend*sep
table8raw = DataFrame(CSV.File(path*"table8.csv",
                              header=false));
table5raw = DataFrame(CSV.File(datadir*fyend*sep*"table5.csv",
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
      table[AProw, 1] = "T3";
      GDProw = findfirst(occursin.("GDP",
                                  string.(table[:, 1])));
      table = table[Not(GDProw), :];)
    : (TUrow = findfirst(occursin.("Total uses", string.(table[:, 2])));
       table[TUrow, 1] = "T3")
  )
  GVArow = findfirst(occursin.("Value Added", string.(table[:, 2])));
  table[GVArow, 1] = "T2";
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
 transform table to "numdiv" sectors
==============================================================================#
# map ioig to ANZSICnumdiv
# the following function is impure in that it depends on dicts and functions ..
# in concordance.jl
function makeaustable(table, numsec=numioig)
  # create a many-to-one column: anzcode per ioig
  (# isolate ioig codes
   ioigcodes = string.(table.IOIG[1:numsec]);
   ioigcodesFloat = parse.(Float64, ioigcodes);
   ioigtonumdiv = mapioigdiv(ioigcodesFloat);
   tmp = String[];
   for i in eachindex(ioigcodesFloat)
     push!(tmp, ioigtonumdiv[ioigcodesFloat[i]])
   end;
   # and the remaining terms in the table;
   for i in range(1, nrow(table) - numsec)
     push!(tmp, table[numsec + i, 1])
   end;
  anzperioigcol = tmp;
  )
  # collapse rows
  (austable = table;
   insertcols!(austable, 1, "ANZcode" => anzperioigcol);
   austable = (combine(groupby(austable, :ANZcode),
                        names(table, Between(:IOIG0101, :T6))
                          .=> sum .=> Between(:IOIG0101, :T6))
               );
   sort!(austable);
  )
  # collapse cols
  (tmptable = transpose(Matrix(austable[:, 2 : numioig + 1]));
   tmptable = DataFrame(tmptable, :auto);
   insertcols!(tmptable, 1, "ANZcode" => anzperioigcol[1:numioig]);
   tmptable = combine(groupby(tmptable, :ANZcode), Not(1) .=> sum);
   sort!(tmptable);
   tmptablet = transpose(Matrix(tmptable[:, Not(1)]));
   tmptablet = DataFrame(tmptablet, ANZcode);
   insertcols!(tmptablet, 1, "ANZcode" => austable[:,1]);
   austable = austable[:, Not(names(austable, Between(:IOIG0101, :IOIG9502)))];
   austable = innerjoin(tmptablet, austable, on = :ANZcode);
  )
  return austable
end
austable8 = makeaustable(table8);
austable5 = makeaustable(table5);
(austablediff, anznegvals, anzposdiffP6) = makediff(austable8, austable5);

#==============================================================================
regionalise tables
==============================================================================#
#make the io tables ras-ready
#table = austable8
function makerr(table)
  gvarow = findfirst(occursin.("T2", table.ANZcode));
  t1row = findfirst(occursin.("T1", table.ANZcode));
  p5row = findfirst(occursin.("P5", table.ANZcode));
  table = table[Not(gvarow, t1row), :];
  table.Q3 += table.Q6;
  table.Q6 .= 0;
  table.Q4 += table.Q5;
  table.Q5 .= 0;
  table = table[:, Not([:T4, :T5])];
  #table[:, Not(:ANZcode)] = round.(table[:, Not(:ANZcode)]);
  rename!(table, :Q4 => "Q4");
  return table
end 
austable8rr = makerr(austable8)
austable5rr = makerr(austable5)
austablediffrr = makerr(austablediff)
# Export tables as CSV
CSV.write(outputdir * "austable8rr.csv", austable8rr)
CSV.write(outputdir * "austable5rr.csv", austable5rr)
for i in range(1, numdiv)
  for j in range(1, numdiv + 2)
    (austablediffrr[i, j + 1] < 0) && (austablediffrr[i, j + 1] = 0);
  end
end;
CSV.write(outputdir * "austablediffrr.csv", austablediffrr);
# import regional data from remplan directory
# what region?
region = "GLD"
reference = "AUS"
source = "remplan"
#source = "john"
dir = datadir*source*sep
#=============================================================================#
(source == "remplan"
 ? (outpersecraw = XLSX.readdata(dir*"outpersec.xlsx", "Report", "A9:C28");
    labpersecraw = XLSX.readdata(dir*"labpersec.xlsx", "Report", "A9:C28");
    gvapersecraw = XLSX.readdata(dir*"gvapersec.xlsx", "Report", "A9:C28");
    salpersecraw = XLSX.readdata(dir*"salpersec.xlsx", "Report", "A9:C28");
    ysapersecraw = XLSX.readdata(dir*"ysapersec.xlsx", "Report", "A9:C28");
    exapersecraw = XLSX.readdata(dir*"exapersec.xlsx", "Report", "A9:C28");
    grxgrpraw = XLSX.readdata(dir*"grxgrp2.xlsx", "Report", "A9:C20");
   )
 # otherwise john's data
 : (outpersecraw = XLSX.readdata(dir*"GRP estimates Gladstone.xlsx",
                                 "Gladstone", "F2:H21");
    grxgrpraw = XLSX.readdata(dir*"GRP estimates Gladstone.xlsx",
                              "Gladstone", "L2:M13");
    for i in [" ", "\$", "M", ",", ]
      outpersecraw[:, 2:3] = strip.(replace.(outpersecraw[:, 2:3], i => ""));
      grxgrpraw[1:9, 2] = strip.(replace.(grxgrpraw[1:9, 2], i => ""))
      outpersecraw[:, 2:end] = parse.(Float64, outpersecraw[:, 2:end])
    end;
   )
)
#=============================================================================#
labpersecraw[:, 2:end] = Int.(round.(labpersecraw[:, 2:end]))

# prepare regional output column
function makebasiclq(raw, reg=region, ref=reference)
  df = DataFrame(raw, ["ANZdiv", reg, ref]);
  (typeof(raw[1, 2]) == Int64
   && (println("This should be labour data");
      );
  );
  (typeof(raw[1, 2]) == Float64
   && (println("\nThis should be output data");
       df[1:end, 2:end] = Matrix(df[1:end, 2:end]) / 1e+6;
      );
  );
#  df[end, 1] = "Ownership of Dwellings"
  tmp = DataFrame("ANZcode" => ANZcode, "ANZdiv" => ANZdiv)
  df = innerjoin(tmp, df, on = "ANZdiv")
  sort!(df)
  print(df)
  refgrowth = (df[end, ref] - austable8[end, "T4"]) / austable8[end, "T4"];
  ((austable8[numdiv, 1] == "T") & (df[5, 1] == "E") & (austable8[5, 1] == "E")
   ? (aggregshr = df[end, reg] / df[end, ref];
      Erefshr = austable8[end, "E"] / austable8[end, "T4"]; 
      Trefshr = austable8[end, "T"] / austable8[end, "T4"];
      origEval = Vector(df[5, Cols(reg, ref)]);
      df[5, Cols(reg, ref)] = (Erefshr / (Erefshr + Trefshr) * origEval);
      df[end, Cols(reg, ref)] = (Trefshr / (Erefshr + Trefshr) * origEval);
     )
   : println("check austable8 E and T are not in place")
  )
  return (df, refgrowth)
end;
labpersec = makebasiclq(labpersecraw)[1];
# save regional labour to julia output folder as csv in ampl-ready form
CSV.write(outputdir * "labsecreg.csv", labpersec[1:numdiv, [1, 3, 4]]);
outrowpersec = makebasiclq(outpersecraw)[1];
(outpersec, outgrowth) = makebasiclq(outpersecraw);
# more regional data
(source == "remplan"
 ? (grxgrp = DataFrame(grxgrpraw, ["Aggregates", region, reference]);
    rng = 2:3;
  (grxgrp[1:9, rng]
   = abs.(Matrix(grxgrp[1:9, rng])) / 1e+6);
   )
 : (grxgrp = DataFrame(grxgrpraw, ["Aggregates", region]);
    rng = 2;
    (grxgrp[1:9, rng]
     = abs.(Vector(grxgrp[1:9, rng])) / 1e+6);
   )
);
(p1aggshrgva = austable8[numdiv+1, end]
 / sum(austable8[[numdiv+1, numdiv+2, numdiv+4], end]));
(p2aggshrgva = austable8[numdiv+2, end]
 / sum(austable8[[numdiv+1, numdiv+2, numdiv+4], end]));
(p4aggshrgva = austable8[numdiv+4, end]
 / sum(austable8[[numdiv+1, numdiv+2, numdiv+4], end]));
(p4aggshrgvanetsal = austable8[numdiv+4, end]
 / sum(austable8[[numdiv+2, numdiv+4], end]));
# make P1 to P6
p = zeros(6, 2);
(source == "remplan"
 ? (p[1, :] = salpersecraw[end, rng] / 1e+6;
    p[2, :] = (gvapersecraw[end, rng]
               - salpersecraw[end, rng]) * (1 - p4aggshrgvanetsal) / 1e+6;
    p[3, :] = (outpersecraw[end, rng]
               * austable8[numdiv+3, end] / austable8[end-2, end]) / 1e+6;
    p[4, :] = (gvapersecraw[end, rng]
               - salpersecraw[end, rng]) * p4aggshrgvanetsal / 1e+6;
   )
 : (for i in range(1, 4)
      p[i, :] = (outpersecraw[end, rng]
                 * austable8[numdiv+i, end] / austable8[end-2, end]) / 1e+6;
    end
   )
);
p[5, :] = zeros(2);
p[6, :] = sum(Matrix(grxgrp[7:8, rng]), dims=1);
# now for Q1 to Q7
q = zeros(7, 2);
q[[1, 2, 3, 4, 7], :] = Matrix(grxgrp[[1, 2, 3, 4, 6], rng]);
q[5, :] = zeros(2);
aggoutlq = outpersecraw[end, 2] / outpersecraw[end, 3];
q[6, :] = 0 * austable8[end, "Q6"] * [aggoutlq , 1];

gap = sum(q, dims=1) - sum(p, dims=1);
gapgvadbn = p[[1, 2, 4], :] / sum(p[[1, 2, 4], :], dims=1) .* gap;
p[[1, 2, 4], :] = p[[1, 2, 4], :] + gapgvadbn;
outperdiv = outpersec;
CSV.write(outputdir * "outperdiv.csv", outperdiv);
outpersec = outpersec[:, Not("ANZdiv")];
for i in range(1, 6)
  push!(outpersec, [austable8[numdiv + i, 1] p[i, 1] p[i, 2]])
end;
labperdiv = labpersec;
outrowpersec = outrowpersec[:, Not("ANZdiv")];
for i in range(1, 7)
  push!(outrowpersec, ["`Q"*string(i) q[i, 1] q[i, 2]])
end;
# 
println("")
println(sum(outrowpersec[:, 2]) - sum(outpersec[:, 2]));
println(sum(outrowpersec[:, 3]) - sum(outpersec[:, 3]));

#==============================================================================
RAS
==============================================================================#
numrow = size(austable8rr, 1) - 1;
Q7col = columnindex(austable8rr, "Q7");
Acol = columnindex(austable8rr, "A");
numcol = size(austable8rr, 2) - 2; 
newrowtot = austable8rr[end, Between("A", "Q7")];
newcoltot = austable8rr[1:numrow, end];
newrowtot = outrowpersec[1:numcol, 2];
newcoltot = outpersec[1:numrow, 2];

# with current data, we need to RAS 5 first, so imports are excluded from flows
tn = "5";
# make sure we have a valid table and, if so, set it
println("Table ", tn, " coming up.");
table5origdf = austable5rr[:, Between("A", "T6")];
tabledf = deepcopy(table5origdf);
for i in range(1, numrow)
  tabledf[i, "T6"] = newcoltot[i];
  for j in range(1, numcol)
    tabledf[end, j] = newrowtot[j];
    (tabledf[i, j] != 0
     && (tabledf[i, j] = tabledf[i, j] * tabledf[end, j] / table5origdf[end, j])
    );
  end;
end;
(sum(tabledf[end, Between("A", "Q7")]) - sum(tabledf[1:numrow, "T6"]) == 0
 ? tabledf[end, end] = sum(tabledf[end, Between("A", "Q7")])
 : println("WARNING: row and column sums don't match!")
)
table = Matrix(tabledf[:, Between("A", "T6")])
# instantiate an optimisation problem
ras5 = Model(Ipopt.Optimizer);
# Should be equal dimensions
lbx = 1e-1
@variable(ras5, x[1:numrow+1, 1:numcol+1] >= lbx);
set_start_value.(x, table)
# Max entropy objective (or min relative to uniform)
eps = 050e-2;
shr = zeros(numrow+1, numcol+1);
for i in range(1, numrow+1)
  for j in range(1, numcol+1)
    shr[i, j] = (1 / count(!iszero, table[:, j]) * 11
                 * ((rand(MersenneTwister(i + j)) - 50e-2) * 10e-2 + 1) 
                ) ^ (- 1 / eps);
  end
end;
rho = round(1 + 1 / eps);
rhohat = 1 / rho;
@NLobjective(ras5,
             Max,
             - 1e-8 * (sum(sum(shr[i, j] * (x[i, j] - table[i, j]) ^ rho
                        for i in range(1, numrow)
                       ) ^ (rhohat * 12)
                    for j in range(1, numcol)
                   )
               )
            );
tol = 1e+0;
rowP6 = numdiv + 6;
colQ1 = columnindex(tabledf, :Q1);
colQ7 = columnindex(tabledf, :Q7);
#------------------------------------------------------------------------------
# the following are constraints for table 5
#------------------------------------------------------------------------------
# Col-sums constraint - must be equal to the IO totals
for i in range(1, numrow)
    #@constraint(ras5, sum(x[i,:]) <= newcoltot[i] + tol);
    #@constraint(ras5, sum(x[i,:]) >= newcoltot[i] - tol);
    @constraint(ras5, sum(x[i, :]) == 2 * x[i, end]);
    @constraint(ras5, x[i, end] == (table[i, end]
                                    + (x[rowP6, colQ7] * table[i, end]
                                       / sum(table[1:end-1, end]))
                                   )
               );
end;
# Row-sums constraint - must be equal to the GFCF totals
for j in range(1, numcol)
    #@constraint(ras5, sum(x[:, j]) <= newrowtot[j] + tol);
    #@constraint(ras5, sum(x[:, j]) >= newrowtot[j] - tol);
    @constraint(ras5, sum(x[:, j]) == 2 * x[end, j]);
    @constraint(ras5, x[end, j] == (table[end, j]
                                    + (x[rowP6, colQ7] * table[end, j]
                                       / sum(table[end, 1:end-1]))
                                   )
               );
end;
@constraint(ras5, x[rowP6, colQ7] == lbx);
for j in range(1, numdiv)
    q1tableshr = table5origdf[j, :Q1] / sum(table5origdf[1:numdiv, :Q1])
    @constraint(ras5, x[j, colQ1] ==  q1tableshr * sum(x[1:numdiv, colQ1]));
end;
colC = columnindex(tabledf, :C);
@constraint(ras5, x[4, colC] >= 600);
@constraint(ras5, x[3, colC] >= 500);
#@constraint(ras5, x[numdiv+4, colC] == -100);
for i in range(1, numrow)
  for j in range(1, numcol)
    (table[i, j] == 0 && @constraint(ras5, x[i, j] == lbx)
    );
    (table[i, j] > 0 && @constraint(ras5, x[i, j]
                                    >= tol * rand(MersenneTwister(i + j)))
    );
  end
end;

#==============================================================================#
# and finally ras table 5
#==============================================================================#
set_optimizer_attribute(ras5, "max_iter", 2000);
optimize!(ras5);
rawsol5 = value.(x);
for i in range(1, numrow)
  for j in range(1, numcol)
    (abs(rawsol5[i, j]) < 1e-7) && (rawsol5[i, j] = 0);
  end
end;
println("biggest movement is ", maximum(abs.(rawsol5 - table))
        , "with index ", argmax(abs.(rawsol5 - table)));
rowerr = round.(sum(rawsol5[1:end-1,:], dims=1)' - table[end, :]; digits=2);
colerr = round.(sum(rawsol5[:, 1:end-1], dims=2) - table[:, end]; digits=2);
println("and the errors (rowerr and colerr) are: ")
println(rowerr')
println(colerr')
sol5 = DataFrame(rawsol5, names(tabledf));
#insertcols!(sol, numcol + 1, :T6 => newcoltot);
insertcols!(sol5, 1, :ANZcode => austable8rr[:, 1]);
#push!(sol, [austable5rr[end, 1] newrowtot' sum(newrowtot)]);
println("and the new table", tn, " is:")
print(pretty_table(sol5, nosubheader=true, formatters=ft_round(1)));
# Export tables as CSV
CSV.write(outputdir*"newtable"*tn*".csv", sol5);

#------------------------------------------------------------------------------
# now for table 8
#------------------------------------------------------------------------------
tn = "8";
# make sure we have a valid table and, if so, set it
println("Table ", tn, " coming up.");
table8origdf = austable8rr[:, Between("A", "T6")];
tabledf = deepcopy(table8origdf);
for i in range(1, numrow)
  tabledf[i, "T6"] = newcoltot[i];
  for j in range(1, numcol)
    tabledf[end, j] = newrowtot[j];
    (tabledf[i, j] != 0
     && (tabledf[i, j] = tabledf[i, j] * tabledf[end, j] / table8origdf[end, j])
    );
  end;
end;
(sum(tabledf[end, Between("A", "Q7")]) - sum(tabledf[1:numrow, "T6"]) == 0
 ? tabledf[end, end] = sum(tabledf[end, Between("A", "Q7")])
 : println("WARNING: row and column sums don't match!")
)
table = Matrix(tabledf[:, Between("A", "T6")])
# instantiate an optimisation problem
ras8 = Model(Ipopt.Optimizer);
# Should be equal dimensions
@variable(ras8, x[1:numrow+1, 1:numcol+1] >= lbx);
set_start_value.(x, table)
# Max entropy objective (or min relative to uniform)
eps = 050e-2;
@NLobjective(ras8,
             Max,
             - 1e-0 * (sum(sum(shr[i, j] * (x[i, j] - table[i, j]) ^ rho
                        for i in range(1, numrow)
                       ) ^ (rhohat * 12)
                    for j in range(1, numcol)
                   )
               )
            );
tol = 1e+1;
#------------------------------------------------------------------------------
# the following are constraints for table 8
#------------------------------------------------------------------------------
# Col-sums constraint - must be equal to the IO totals
for i in range(1, numrow)
    #@constraint(ras8, sum(x[i,:]) <= newcoltot[i] + tol);
    #@constraint(ras8, sum(x[i,:]) >= newcoltot[i] - tol);
    @constraint(ras8, sum(x[i, :]) == 2 * x[i, end]);
end;
# Row-sums constraint - must be equal to the GFCF totals
for j in range(1, numcol)
    #@constraint(ras8, sum(x[:, j]) <= newrowtot[j] + tol);
    #@constraint(ras8, sum(x[:, j]) >= newrowtot[j] - tol);
    @constraint(ras8, sum(x[:, j]) == 2 * x[end, j]);
    # gva and t1 have to be in proportion
    gvarows = [numdiv + 1, numdiv + 2, numdiv + 4];
    @constraint(ras8, sum(x[1:numdiv, j]) <= sum(x[gvarows, j]) + tol)
    @constraint(ras8, sum(x[1:numdiv, j]) >= sum(x[gvarows, j]) - tol)
end;
colC = columnindex(tabledf, :C);
colQ1 = columnindex(tabledf, :Q1);
colQ3 = columnindex(tabledf, :Q3);
colQ4 = columnindex(tabledf, :Q4);
rowG = 7; # the retail row
for i in range(1, numdiv)
    q1tableshr = table8origdf[i, :Q1] / sum(table8origdf[1:numdiv, :Q1])
    @constraint(ras8, x[i, colQ1] ==  q1tableshr * sum(x[1:numdiv, colQ1]));
    @constraint(ras8, x[i, colQ3] >= rawsol5[i, colQ3])
    @constraint(ras8, x[i, colQ4] >= rawsol5[i, colQ4])
end;
# mining imports of bauxite by the refineries
@constraint(ras8, x[2, colC] >= 800);
# alumina purchases by the smelter
@constraint(ras8, x[3, colC] >= 500);
# electricity purchases by the smelter
@constraint(ras8, x[4, colC] >= 600);
# retail has to be more import dependent than australia
#@constraint(ras8, x[rowG, end] == rawsol5[rowG, end] * 200e-2);
#@constraint(ras8, x[numdiv+4, colC] == -100);
for i in range(1, numrow)
  for j in range(1, numcol)
    (table[i, j] == 0 && @constraint(ras8, x[i, j] == lbx)
    );
    (table[i, j] > 0 && @constraint(ras8, x[i, j]
                                    >= tol * rand(MersenneTwister(i + j)))
    );
  end
end;
for i in range(1, numrow+1)
  for j in range(1, numcol+1)
    #@constraint(ras8, x[rowG, j] == rawsol5[rowG, j]);
    (((i <= numdiv) & (j <= numdiv)) 
     ? (@constraint(ras8, x[i, j] >= rawsol5[i, j]))
     : ((-1 * tol < round(table8origdf[i, j] - table5origdf[i, j]) < 1 * tol
         && @constraint(ras8, x[i, j] == rawsol5[i, j])
        );
        (round(table8origdf[i, j] - table5origdf[i, j]) >= 1 * tol
         && @constraint(ras8, x[i, j] >=  rawsol5[i, j] + 1 * tol)
        );
        (round(table8origdf[i, j] - table5origdf[i, j]) <= -1 * tol
         && @constraint(ras8, x[i, j] <=  rawsol5[i, j])
        );
       );
    )
  end
end;
# row column balance for sectors
for i in range(1, numdiv)
    @constraint(ras8, x[i, end] == x[end, i])
end;
for i in range(numdiv+1, numrow)
  @constraint(ras8, x[i, end] == rawsol5[i, end]);
end;
for i in range(numdiv+1, numcol)
  @constraint(ras8, x[end, i] == rawsol5[end, i]);
  i != numcol && @constraint(ras8, x[end-1, i] == lbx);
end;
@constraint(ras8, x[end, end] == rawsol5[end, end] + rawsol5[26, end]);
@constraint(ras8, sum(x[1:end-1, end]) == x[end, end]);
@constraint(ras8, sum(x[end, 1:end-1]) == x[end, end]);

optimize!(ras8);

rawsol8 = value.(x);
for i in range(1, numrow)
  for j in range(1, numcol)
    (abs(rawsol8[i, j]) < 1e-7) && (rawsol8[i, j] = 0);
  end
end;
println("biggest movement is ", maximum(abs.(rawsol8 - table))
        , "with index ", argmax(abs.(rawsol8 - table)));
err = round.(rawsol8[1:numdiv, end] - rawsol8[end, 1:numdiv])
println("and the error (diff between rowsum and colsum) is: ")
println(err')
sol8 = DataFrame(rawsol8, names(tabledf));
#insertcols!(sol, numcol + 1, :T6 => newcoltot);
insertcols!(sol8, 1, :ANZcode => austable8rr[:, 1]);
#push!(sol, [austable5rr[end, 1] newrowtot' sum(newrowtot)]);
println("and the new table", tn, " is:")
println(pretty_table(sol8, nosubheader=true, formatters=ft_round(1)));
# Export table as CSV
CSV.write(outputdir*"newtable"*tn*".csv", sol8);
rawsoldiff = rawsol8 - rawsol5;
soldiff = deepcopy(sol8);
soldiff[:, 2:end] = rawsol8 - rawsol5;
println("and the new table of differences is: ")
println(pretty_table(soldiff, nosubheader=true, formatters=ft_round(1)));
# Export table as CSV
for i in range(1, numdiv)
  for j in range(1, numdiv)
    (soldiff[i, j + 1]) < 0 && (soldiff[i, j + 1] = 0);
  end
end;
CSV.write(outputdir * "soldiff.csv", soldiff);


gldgva8 =  sum(Matrix(getsubframe(sol8,
                                 ["`P1", "`P2", "`P4"],
                                 sectorcodes))[:, 2:end],
              dims=1)
gldt18 =  sum(Matrix(getsubframe(sol8,
                                sectorcodes,
                                sectorcodes))[:, 2:end],
             dims=1)
#==============================================================================
kapital ras prep
==============================================================================#
#table 8
gfcfioig = (austable8rr[1:numdiv, :Q3]
            + austable8rr[1:numdiv, :Q4] + austable8rr[1:numdiv, :Q5]);
gfcfioigtot = sum(gfcfioig)
#==============================================================================
Aggregate GFCF data
==============================================================================#

# Bring in GFCF data from excel
ausGFCFall = ExcelReaders.readxlsheet(
  datadir*"5204064_GFCF_By_Industry_Asset.xls", "Data1"
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
gfcfbiba= ausGFCFcurrent[
  findall(x -> occursin(fyend, x), string.(ausGFCFcurrent[:,1])), :];
# exclude ownership transfer costs to dwellings
#gfcfbiba[length(gfcfbiba) - 1] = (gfcfbiba[length(gfcfbiba) - 1]
#                                      + gfcfbiba[length(gfcfbiba)]);
# remove date column and ownership transfer costs
gfcfbiba = gfcfbiba[Not(1, length(gfcfbiba))];

# no ownership of dwellings
gfcfbiba[12] += gfcfbiba[numdiv + 1];
gfcfbiba = gfcfbiba[Not(numdiv + 1)];
gfcfbibatot = sum(gfcfbiba);

# Balancing row and column sums to IO table 8 total
println("Before balancing, the difference between the two gfcf totals is: ",
        gfcfioigtot - gfcfbibatot)
for i in eachindex(gfcfbiba)
    gfcfbiba[i] = gfcfbiba[i] / gfcfbibatot * gfcfioigtot;
end;
println("gfcfbiba has been normalised to match the gfcfioig total.");

# prior via excel:
# pull in proportionalised kapital flows to ras
#y = DataFrame(CSV.File(datadir*"propd-to-ras.csv", header=false))
#y = Matrix(y)
# generate an initial table for the ras
#y = zeros(length(gfcfbiba), length(anzdivgfcf.inv_sum))
#for i in eachindex(gfcfbiba), j in eachindex(anzdivgfcf);
#  y[i, j] = anzdivgfcf[i] / gfcfioigtot * gfcfbiba[j]
#end  

#=============================================================================#
# Prepare US data
#=============================================================================#
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
flowsColSum = sum(eachcol(flowsTemp[:, Not(:ANZcode)]));


# Adding dwellings column
#flows[!, :T] = flowsColSum *rowSumsProps[numdiv];

# Adding public admin column
colO = columnindex(austable8rr[:, Not(:ANZcode)], :O);
flows[!, :O] = flowsColSum * gfcfbiba[colO] / gfcfioigtot;
# Sorting columns
flows = permutedims(sort(permutedims(flows, 1), :ANZcode), 1);
# write out us flows
CSV.write(outputdir*"uskapflw.csv", flows);

usrowsum = sum(eachcol(permutedims(flows, 1)[:, Not(:ANZcode)]));
#=============================================================================# 
# Back to Aus Data
#=============================================================================# 
ausshr = 80e-2
# Make prior from Aus Data
# Make vector of the proportion of each row sum element as a fraction of the 
# total
ausprior = deepcopy(flows)
for j in range(1, numdiv)
  # In the initial prior, every column has the same cost structure
  ausprior[:, j + 1] = gfcfioig  / gfcfioigtot * ausshr;
  # rescale us flows
  flows[:, j + 1] .*= 1 / usrowsum[j] * (1 - ausshr);
end
# the new prior is the convex combination per ausshr
ausprior[:, Not(:ANZcode)] .+= flows[:, Not(:ANZcode)];
# now scale each column back up per gfcfbiba
for i in range(1, numdiv)
  ausprior[:, i + 1] .*= gfcfbiba[i];
end;


raskap8 = Model(Ipopt.Optimizer);
# Should be equal dimensions
numcol = numdiv;
numrow = numdiv;
@variable(raskap8, x[1:numrow, 1:numcol] >= lbx);
set_start_value.(x, flows[:, Not(:ANZcode)]);
# Max entropy objective (or min relative to uniform)
eps = 050e-2;
@NLobjective(raskap8,
             Min,
             (sum(sum((x[i, j] - ausprior[i, j + 1]) ^ 2
                      for i in range(1, numrow)
                     )
                  for j in range(1, numcol)
                 )
             )
            );
tol = 1e+0;
#------------------------------------------------------------------------------
# the following are constraints for table 8
#------------------------------------------------------------------------------
# Col-sums constraint - must be equal to the IO totals
for i in range(1, numrow)
    @constraint(raskap8, sum(x[i,:]) <= gfcfioig[i] + tol);
    @constraint(raskap8, sum(x[i,:]) >= gfcfioig[i] - tol);
    #@constraint(raskap8, sum(x[i, :]) == 2 * x[i, end]);
end;
# Row-sums constraint - must be equal to the GFCF totals
for j in range(1, numcol)
    @constraint(raskap8, sum(x[:, j]) <= gfcfbiba[j] + tol);
    @constraint(raskap8, sum(x[:, j]) >= gfcfbiba[j] - tol);
    #@constraint(raskap8, sum(x[:, j]) == 2 * x[end, j]);
end;

optimize!(raskap8)

rawsolkap8 = value.(x);
solkap8 = deepcopy(flows)
solkap8[:, Not(:ANZcode)] = rawsolkap8;
println(pretty_table(solkap8, nosubheader=true, formatters=ft_round(1)));
CSV.write(outputdir * "auskapflw8.csv", solkap8);

#------------------------------------------------------------------------------
# regional kapital
#------------------------------------------------------------------------------
#rasregkap8 = Model(Ipopt.Optimizer);
#@variable(rasregkap8, x[1:numrow, 1:numcol] >= lbx);
#set_start_value.(x[1:numrow], solkap8[1:numrow, Not(:ANZcode)]);
#@NLobjective(rasregkap8,
#             Min,
#             sum(sum((x[i, j] - solkap8[i, j + 1]) ^ 2
#                     for i in range(1, numrow)
#                    ) ^ 2
#                 for j in range(1, numcol)
#                )
#            );
#regcolsum = sol8[1:numrow, :Q3] + sol8[1:numrow, :Q4];
#regcolsumtot = sum(regcolsum);
#regrowsum = gfcfbiba * regcolsumtot / gfcfbibatot;
## Col-sums constraint - must be equal to the IO totals
#for i in range(1, numrow)
#  @constraint(rasregkap8, sum(x[i, :]) <= regcolsum[i] + tol);
#  @constraint(rasregkap8, sum(x[i, :]) >= regcolsum[i] - tol);
#end; 
## Row-sums constraint - must be equal to the GFCF totals
#for j in range(1, numcol)
#  @constraint(rasregkap8, sum(x[:, j]) <= regrowsum[j] + tol);
#  @constraint(rasregkap8, sum(x[:, j]) >= regrowsum[j] - tol);
#    #@constraint(raskap8, sum(x[:, j]) == 2 * x[end, j]);
#end;
#
#optimize!(rasregkap8)
#
#rawsolregkap8 = value.(x);
#solregkap8 = deepcopy(flows)
#solregkap8[:, Not(:ANZcode)] = rawsolregkap8[1:numrow, :];
#insertcols!(solregkap8, ncol(solregkap8) + 1, :T8 => regcolsum)
#push!(solregkap8, ["T7" regrowsum' regcolsumtot])
#println(pretty_table(solregkap8, nosubheader=true, formatters=ft_round(1)));
## Export tables as CSV
#CSV.write(outputdir*"regkapflw8.csv", solregkap8);
