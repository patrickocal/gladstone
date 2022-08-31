# include("preptables-pretable8specific.jl")
# for the time being this should be run with pwd being the gladstone parent file

using DataFrames, CSV, DelimitedFiles;

# other code directory (relative to present file)
codedir = "./"
include(codedir * "concordance.jl")

# Directory for output (and input while pulling template from .csv)
outputdir = "julia/output/"
region = ["GLD"]
refreg = ["AUS"]
# Name of column that contains the row indices
divhdr = "ANZcode"

table8name = "austable8rr"
table5name = "austable5rr"
tablediffname = "austablediffrr"

#table8name = "newtable8"
#table5name = "newtable5"
#tablediffname = "soldiff"

tablekapname = "auskapflw8"
tablelabname = "labsecreg"
tablepopname = "popreg"
tableoutname = "outperdiv"
#tbl8 = "tbl8"
#tbl5 = "tbl5"
# Pull in tables from .csv files, ultimately it will come from the RAS code
tbl8 = CSV.read(outputdir * table8name * ".csv", DataFrame)
tbl5 = CSV.read(outputdir * table5name * ".csv", DataFrame)
tblkap = CSV.read(outputdir * tablekapname * ".csv", DataFrame)
tbllab = CSV.read(outputdir * tablelabname * ".csv", DataFrame)
tbldiff = CSV.read(outputdir * tablediffname * ".csv", DataFrame)
tbloutperdiv = CSV.read(outputdir * tableoutname * ".csv", DataFrame)
#tblpop = CSV.read(outputdir * tablepopname * ".csv", DataFrame)

# Quickly mock up tbldiff from tables 5 and 8, ultimately it will come 
# from the RAS code #include -d at the top of this file
#tbldiff = select(tbl8, Not(divhdr))
#for i in 1:nrow(tbldiff)
#  for j in 1:ncol(tbldiff)
#    tbldiff[i,j] = (tbldiff[i,j] - select(tbl5, Not(divhdr))[i,j])
#  end
#end

#tbldiff[!, divhdr] = tbl8[:, divhdr]

# list of names of sectors, ultimately it will come from the RAS code 
# include -d at the top of this file
sectorcodes = ANZcode 
# vector of region names
# ["a", "b", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", 
#    "n", "o", "P", "Q", "R", "S", "T"];
numdiv = length(sectorcodes);
#=============================================================================#
# general code to add Sector titles and convert to a df from a matrix 
# dataset if that's what the input is (not the case right now)

# here need to just get the data, no headings etc. in matrix form
# comment out if they are already in this form
# tbl8 = Matrix(tbl8)
# tbl5 = Matrix(tbl5)
# tbldiff = Matrix(tbldiff)

# function addTitles(table)
#     df = dataFrame(table, sectorcodes);
#     insertcols!(df, 1, :Sectors => SectorCodes);
#     return df;
# end

# addtitles(tbl8)
# addtitles(tbl5)
# addtitles(tbldiff)

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
function outtocsv(paramName, subFrame);
    colNames = names(subFrame)
    colNames[length(colNames)] = paramName
    rename!(subFrame, colNames)
    CSV.write(outputdir*paramName*".csv", subFrame)
end
#=============================================================================#
# Calculate and export RAW parameters as .csv files
#=============================================================================#

#-----------------------------------------------------------------------------#
# three indices: regions sectors and sectors: flows between sectors
#-----------------------------------------------------------------------------#
# RAW_INV_FLW: from kapflws
outtocsv("RAW_INV_FLW", prepsubframe(tblkap, sectorcodes, sectorcodes));
# RAW_MED_FLW - all inter-sector flows table 8
outtocsv("RAW_MED_FLW", prepsubframe(tbl8, sectorcodes, sectorcodes));
# RAW_DOM_CMED - all sectors flows of table 5
outtocsv("RAW_DOM_CMED", prepsubframe(tbl5, sectorcodes, sectorcodes));
# RAW_YPO_CMED - all sectors flows of table 8-5
outtocsv("RAW_YPO_CMED", prepsubframe(tbldiff, sectorcodes, sectorcodes));

#-----------------------------------------------------------------------------#
# two indices: regions and sectors
#-----------------------------------------------------------------------------#
# RAW_CON_FLW - all sectors in the Q1 col of table 8
Q12 = getsubframe(tbl8, sectorcodes, ["Q1", "Q2"]);
Q12[:, "Q1"] += Q12[:, "Q2"]; 
outtocsv("RAW_CON_FLW", prepsubframe(Q12, sectorcodes, ["Q1"]));
# RAW_DOM_CCON - all sectors in the Q1 col of table 5
Q12 = getsubframe(tbl5, sectorcodes, ["Q1", "Q2"]);
Q12[:, "Q1"] += Q12[:, "Q2"]; 
outtocsv("RAW_DOM_CCON", prepsubframe(Q12, sectorcodes, ["Q1"]));
# RAW_YPO_CCON - all sectors flows of table 8-5
Q12 = getsubframe(tbldiff, sectorcodes, ["Q1", "Q2"]);
Q12[:, "Q1"] += Q12[:, "Q2"]; 
outtocsv("RAW_YPO_CCON", prepsubframe(Q12, sectorcodes, ["Q1"]));
#
# RAW_KAP_OUT - all sectors in the P2 row of table 8
# P2 is a row - hence the permutedims
(table8name == "austable8rr"
 ?  outtocsv("RAW_KAP_OUT",
           prepsubframe(permutedims(tbl8, 1), sectorcodes, ["`P2"]))
 :  outtocsv("RAW_KAP_OUT",
           prepsubframe(permutedims(tbl5, 1), sectorcodes, ["`P2"]))
);
# RAW_LAB_OUT all sectors in the P1 row of table 8
(table8name == "austable8rr"
 ?  outtocsv("RAW_LAB_OUT",
           prepsubframe(permutedims(tbl8, 1), sectorcodes, ["`P1"]))
 :  outtocsv("RAW_LAB_OUT",
           prepsubframe(permutedims(tbl5, 1), sectorcodes, ["`P1"]))
);
# RAW_MED_OUT - all sectors in the T1 row of table 8 (maybe use sum instead)
(table8name == "austable8rr"
 ?  T1 = getsubframe(permutedims(tbl8, 1), sectorcodes, sectorcodes)
 :  T1 = getsubframe(permutedims(tbl5, 1), sectorcodes, sectorcodes)
);
for i in names(tbl5, Between(:B, :S))
  T1[:, "A"] += T1[:, i]
end;
T1 = T1[:, Between(:ANZcode, :A)];
rename!(T1, ["ANZcode", "T1"]);
outtocsv("RAW_MED_OUT", prepsubframe(T1, sectorcodes, ["T1"]));



# RAW_XPO_JOUT -  all sectors in the Q1 col of table 8
outtocsv("RAW_XPO_JOUT", prepsubframe(tbl8, sectorcodes, ["Q7"]));
# RAW_DOM_JOUT - all sectors in the (T6 - Q7) col of table 8
T5 = getsubframe(tbl8, sectorcodes, ["Q7", "T6"]);
T5[:, "T6"] -= T5[:, "Q7"];
outtocsv("RAW_DOM_JOUT", prepsubframe(T5, sectorcodes, ["T6"]));

# reference REF_LAB
outtocsv("REF_LAB", prepsubframe(tbllab, sectorcodes, refreg));
# regional REG_LAB
outtocsv("REG_LAB", prepsubframe(tbllab, sectorcodes, region));
# regional AUS_POP
#outtocsv("AUS_POP", prepsubframe(tblpop, sectorcodes, refreg));
# regional REG_POP
#outtocsv("REG_POP", prepsubframe(tblpop, sectorcodes, region[1]));
# RAW_DOM_CINV 
#Q3 = prepsubframe(tbl5, sectorcodes, ["Q3"]);
#Q4 = prepsubframe(tbl5, sectorcodes, ["Q4"]);
#Q5 = prepsubframe(tbl5, sectorcodes, ["Q5"]);
#Q345 = deepcopy(Q3);
#Q345[:, end] += Q4[:, end] + Q5[:, end];
#outtocsv("RAW_DOM_CINV", Q3);

# RAW_YPO_CINV - all sectors flows of kapflows table - not yet imported
# CSV.write(outputdir*"RAW_YPO_CINV.csv", prepsubframe(kapflows, 
#     sectorcodes, sectorcodes));

rowshr = getsubframe(permutedims(tbl8, 1), sectorcodes, sectorcodes);
colshr = getsubframe(tbl8, sectorcodes, sectorcodes);
for i in 1:numdiv
  rowshr[:, i + 1] *= 1 / tbl8[i, end];
  colshr[:, i + 1] *= 1 / tbl8[end, i + 1];
end 
rowshr = permutedims(rowshr, 1);
outtocsv("SHR_MED_ROW", prepsubframe(rowshr, sectorcodes, sectorcodes));
outtocsv("SHR_MED_COL", prepsubframe(colshr, sectorcodes, sectorcodes));


# RAW_OUT_REG_SEC output per region and division
outtocsv("RAW_OUT_REG_SEC", prepsubframe(tbloutperdiv, sectorcodes, region)); 

println("done saving to julia/output/*.csv")
