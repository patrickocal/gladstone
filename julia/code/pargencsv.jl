# include("preptables-pretable8specific.jl")
# for the time being this should be run with pwd being the gladstone parent file

using DataFrames, CSV, DelimitedFiles;

# other code directory (relative to present file)
codedir = "./"
include(codedir * "concordance.jl")

# Directory for output (and input while pulling template from .csv)
outputdir = "julia/output/"
# Name of column that contains the row indices
divhdr = "ANZdiv"
tbl8 = "anztable8rr"
tbl5 = "anztable5rr"
tbl8 = "newtable8"
tbl5 = "newtable5"
# Pull in tables from .csv files, ultimately it will come from the RAS code
newtable8 = CSV.read(outputdir * tbl8 * ".csv", DataFrame)
newtable5 = CSV.read(outputdir * tbl5 * ".csv", DataFrame)
# Quickly mock up newtablediff from tables 5 and 8, ultimately it will come 
# from the RAS code #include -d at the top of this file
newtablediff = select(newtable8, Not(divhdr))
for i in 1:nrow(newtablediff)
  for j in 1:ncol(newtablediff)
    newtablediff[i,j] = (newtablediff[i,j] - select(newtable5, Not(divhdr))[i,j])
  end
end

newtablediff[!, divhdr] = newtable8[:, divhdr]

# list of names of sectors, ultimately it will come from the RAS code 
# include -d at the top of this file
sectorcodes = ANZcode 
# vector of region names
region = ["GLD"]
# ["a", "b", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", 
#    "n", "o", "P", "Q", "R", "S", "T"];
numdiv = length(sectorcodes);
#=============================================================================#
# general code to add Sector titles and convert to a df from a matrix 
# dataset if that's what the input is (not the case right now)

# here need to just get the data, no headings etc. in matrix form
# comment out if they are already in this form
# newtable8 = Matrix(newtable8)
# newtable5 = Matrix(newtable5)
# newtablediff = Matrix(newtablediff)

# function addTitles(table)
#     df = dataFrame(table, sectorcodes);
#     insertcols!(df, 1, :Sectors => SectorCodes);
#     return df;
# end

# addtitles(newtable8)
# addtitles(newtable5)
# addtitles(newtablediff)

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
# General Function to get a sub-frame based
function prepsubframe(df, lsRows, lsCols, reg=region, divn=divhdr);
  subdf = getsubframe(df, lsRows, lsCols)
  m=length(lsRows);
  n=length(lsCols);
  (n > 1
   ? (dfNew = DataFrame(index0 = repeat(reg, inner=(m * n), outer=1),
                        index1 = repeat(lsRows, inner=n, outer=1),
                        index2 = repeat(lsCols, inner=1, outer=m),
                        Value = ones(m * n),
                       );
     )
   : (dfNew = DataFrame(index0 = repeat(reg, inner=(m * n), outer=1),
                        index1 = repeat(lsRows, inner=n, outer=1),
                        Value = ones(m * n),
                       );
     );
  );
  for i in 1:m
      for j in 1:n
        dfNew[j + (i - 1) * n, :Value] = select(subdf, Not(divn))[i, j]
      end
  end;
  return dfNew
end
# Output a parameter as a CSV with the correct format for AMPL
function outputForm(paramName, subFrame);
    colNames = names(subFrame)
    colNames[length(colNames)] = paramName
    rename!(subFrame, colNames)
    CSV.write(outputdir*paramName*".csv", subFrame)
end
#=============================================================================#
# Calculate and export RAW parameters as .csv files

# RAW_CON_FLW - all sectors in the Q1 col of table 8
outputForm("RAW_CON_FLW", prepsubframe(newtable8, 
    sectorcodes, ["Q1"]));

# RAW_MED_FLW - all inter-sector flows table 8
outputForm("RAW_MED_FLW", prepsubframe(newtable8, 
    sectorcodes, sectorcodes));

# RAW_KAP_OUT - all sectors in the P2 row of table 8
# First example of a row - hence the permutedims
outputForm("RAW_KAP_OUT",
           prepsubframe(permutedims(newtable8, 1), sectorcodes, ["`P2"]));

# RAW_LAB_OUT all sectors in the P1 row of table 8
outputForm("RAW_LAB_OUT",
           prepsubframe(permutedims(newtable8, 1), sectorcodes, ["`P1"]));

# RAW_MED_OUT - all sectors in the T1 row of table 8 (maybe use sum instead)
T1vec = vec(sum(Matrix(newtable8[1:numdiv, 2:numdiv+1]), dims=1));
T1 = DataFrame(ANZdiv = sectorcodes, T1 = T1vec);
T1 = prepsubframe(T1, sectorcodes, ["T1"]);
outputForm("RAW_MED_OUT", T1);

# RAW_DOM_CCON - all sectors in the Q1 col of table 5
outputForm("RAW_DOM_CCON", prepsubframe(newtable5, sectorcodes, ["Q1"]));

# RAW_YSA_CCON - all sectors flows of table 8-5
outputForm("RAW_YSA_CCON", prepsubframe(newtablediff, sectorcodes, ["Q1"]));

# RAW_DOM_CINV 
#Q3 = prepsubframe(newtable5, sectorcodes, ["Q3"]);
#Q4 = prepsubframe(newtable5, sectorcodes, ["Q4"]);
#Q5 = prepsubframe(newtable5, sectorcodes, ["Q5"]);
#Q345 = deepcopy(Q3);
#Q345[:, end] += Q4[:, end] + Q5[:, end];
#outputForm("RAW_DOM_CINV", Q3);

# RAW_YSA_CINV - all sectors flows of kapflows table - not yet imported
# CSV.write(outputdir*"RAW_YSA_CINV.csv", prepsubframe(kapflows, 
#     sectorcodes, sectorcodes));

# RAW_DOM_CMED - all sectors flows of table 5
outputForm("RAW_DOM_CMED", prepsubframe(newtable5, sectorcodes, sectorcodes));

# RAW_YSA_CMED - all sectors flows of table 8-5
outputForm("RAW_YSA_CMED", prepsubframe(newtablediff, 
sectorcodes, sectorcodes));

# RAW_EXO_JOUT -  all sectors in the Q1 col of table 8
outputForm("RAW_EXO_JOUT", prepsubframe(newtable8, sectorcodes, ["Q7"]));

# RAW_DOM_JOUT - all sectors in the (T6 - Q7) col of table 8
Q7 = prepsubframe(newtable8, sectorcodes, ["Q7"]);
T6 = prepsubframe(newtable8, sectorcodes, ["T6"]);
T5 = deepcopy(T6);
T5[:, end] -= Q7[:, end];
outputForm("RAW_DOM_JOUT", T5);

println("done saving to julia/output/*.csv")
