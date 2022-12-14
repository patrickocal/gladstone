#concordance dict between ANZSIC divisions (19 sectors plus Ownership of
#Dwellings) and various other industry classifications
using XLSX, ExcelReaders, DataFrames, Tables, JuMP, Ipopt, NamedArrays, DelimitedFiles, CSV, Tables;

#filepath cross system compatability code
if Sys.KERNEL === :NT || Sys.KERNEL === :Windows
    sep = "\\";
else
    sep = "/";
end
datadir = "julia"*sep*"data"*sep
#20 Sector to 4 Sector
fromdivTo4 = Dict(
  "A"=> "Primary",
  "B" => "Primary",
  "C" => "Secondary",
  "D" => "Secondary",
  "E" => "Secondary",
  "F" => "Tertiary",
  "G" => "Tertiary",
  "H" => "Tertiary",
  "I" => "Tertiary",
  "J" => "Tertiary",
  "K" => "Tertiary",
  "L" => "Tertiary",
  "M" => "Tertiary",
  "N" => "Tertiary",
  "O" => "Tertiary",
  "P" => "Tertiary",
  "Q" => "Tertiary",
  "R" => "Tertiary",
  "S" => "Tertiary",
#  "T" => "Ownership"
 )

#IOIG to 20 Sector
IOSource = ExcelReaders.readxlsheet(datadir*"5209055001DO001_201819.xls", "Table 8");
IOIG = IOSource[4:117, 1];
ANZSICDiv = ["Agriculture, Forestry and Fishing",
             "Mining",
             "Manufacturing",
             "Electricity, Gas, Water and Waste Services",
             "Construction",
             "Wholesale Trade",
             "Retail Trade",
             "Accomodation and Food Services",
             "Transport, Postal and Warehousing",
             "Information Media and Telecommunications",
             "Financial and Insurance Services",
             "Rental, Hiring and Real Estate Services",
             "Professional, Scientific and Technical Services",
             "Administrative and Support Services",
             "Public Administration and Safety",
             "Education and Training",
             "Health Care and Social Assistance",
             "Arts and Recreation Services",
             "Other Services",
#             "Ownership of Dwellings",
            ];
ANZdiv = ["Agriculture, Forestry & Fishing",
                "Mining",
                "Manufacturing",
                "Electricity, Gas, Water & Waste Services",
                "Construction",
                "Wholesale Trade",
                "Retail Trade",
                "Accommodation & Food Services",
                "Transport, Postal & Warehousing",
                "Information Media & Telecommunications",
                "Financial & Insurance Services",
                "Rental, Hiring & Real Estate Services",
                "Professional, Scientific & Technical Services",
                "Administrative & Support Services",
                "Public Administration & Safety",
                "Education & Training",
                "Health Care & Social Assistance",
                "Arts & Recreation Services",
                "Other Services",
#                "Ownership of Dwellings",
            ];
ANZSICDivShort = ["AgrForestFish", "Mining", "Manufacturing", "Utilities",
                  "Construction", "Wholesale", "Retail", "AccomFoodServ",
                  "Transport&Ware", "Communications", "Finance&Insur",
                  "RealEstate", "BusinessServ", "Admin", "PublicAdminSafe",
                  "Education", "Health&Social", "Arts&Rec", "OtherServices",
#                  "Dwellings"
                 ];
ANZcode =["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O", "P","Q",
          "R","S",
          # "T"
         ];
fromdivToInd = Dict{String, Int64}();
for i in eachindex(ANZcode);
    fromdivToInd[ANZcode[i]] = Int(i);
end

function mapioigdiv(ioigcol)
  tmpdct = Dict{Float64, String}();
  for i in [1:1:length(ioigcol);]
    test = ioigcol[i] / 100
      if 1 <= test < 6
          tmpdct[ioigcol[i]]="A"
      elseif 6 <= test < 11
          tmpdct[ioigcol[i]]="B"
      elseif 11 <= test < 26
          tmpdct[ioigcol[i]]="C"
      elseif 26 <= test < 30
          tmpdct[ioigcol[i]]="D"
      elseif 30 <= test < 33
          tmpdct[ioigcol[i]]="E"
      elseif 33 <= test < 39
          tmpdct[ioigcol[i]]="F"
      elseif 39 <= test < 44
          tmpdct[ioigcol[i]]="G"
      elseif 44 <= test < 46
          tmpdct[ioigcol[i]]="H"
      elseif 46 <= test < 54
          tmpdct[ioigcol[i]]="I"
      elseif 54 <= test < 62
          tmpdct[ioigcol[i]]="J"
      elseif 62 <= test < 66
          tmpdct[ioigcol[i]]="K"
#        elseif (66 <= test < 67 || 67.02 <= test < 69) #ownership of dwellings
      elseif (66 <= test < 69) # no ownership of dwellings
          tmpdct[ioigcol[i]]="L"
      elseif 69 <= test < 72
          tmpdct[ioigcol[i]]="M"
      elseif 72 <= test < 75
          tmpdct[ioigcol[i]]="N"
      elseif 75 <= test < 80
          tmpdct[ioigcol[i]]="O"
      elseif 80 <= test < 84
          tmpdct[ioigcol[i]]="P"
      elseif 84 <= test < 89
          tmpdct[ioigcol[i]]="Q"
      elseif 89 <= test < 94
          tmpdct[ioigcol[i]]="R"
      elseif 94 <= test < 96
          tmpdct[ioigcol[i]]="S"
#      elseif 67 <= test < 67.02 # ownership of dwellings
#          tmpdct[ioigcol[i]]="T" # ownership of dwellings
      else
          print("ERROR: An input has fallen outside of the range of categories")
      end
  end
  return tmpdct
end
IOIGTodiv = mapioigdiv(IOIG)
#ISIC 4.0 To div
ANZSICISICSource = CSV.read(datadir*"ANZSIC06-ISIC3pt1.csv", DataFrame);
ANZSICdiv = ANZSICISICSource[6:1484, 1][findall(x -> typeof(x)<:String15, ANZSICISICSource[6:1484, 4])];
ISIC = ANZSICISICSource[6:1484, 4][findall(x -> typeof(x)<:String15, ANZSICISICSource[6:1484, 4])];
for i in eachindex(ISIC);
    ISIC[i]=strip(ISIC[i], ['p']);
end
ISICTodiv = Dict(ISIC .=> ANZSICdiv);

#NAIC2007 To div via ISIC 4.0
NAICSISICSource = ExcelReaders.readxlsheet(datadir*"2007_NAICS_to_ISIC_4.xls", "NAICS 07 to ISIC 4 technical");
NAICS = string.(Int.(NAICSISICSource[4:1768,1]));
ISICAsANZSIC = NAICSISICSource[4:1768,3];
ISICAsANZSIC = string.(ISICAsANZSIC);
containsX = findall( x -> occursin("X", x), ISICAsANZSIC);
ISICAsANZSIC[containsX] = replace.(ISICAsANZSIC[containsX], "X" => "1");
ISICAsANZSIC = parse.(Float64, ISICAsANZSIC);
NAICSANZSICdiv = string.(zeros(length(ISICAsANZSIC)));
for i in eachindex(ISICAsANZSIC);
    NAICSANZSICdiv[i] = ISICTodiv[lpad(Int(ISICAsANZSIC[i]),4,"0")];
end
NAICS07Todiv = Dict(NAICS .=> NAICSANZSICdiv);

#NAIC2002 To div via NAIC2007
NAICS02To07 = CSV.read(datadir*"2002_to_2007_NAICS.csv", DataFrame);
NAICS02To0702 = string.(NAICS02To07[3:1202, 1]);
NAICS02To0707 = string.(NAICS02To07[3:1202, 3]);
NAICS07Asdiv = string.(zeros(length(NAICS02To0707)));
for i in eachindex(NAICS02To0707);
    NAICS07Asdiv[i] = NAICS07Todiv[NAICS02To0707[i]];
end
NAICS02Todiv = Dict(NAICS02To0702 .=> NAICS07Asdiv);

#NAIC1997 To div via NAIC2002
NAICS97To02 = CSV.read(datadir*"1997_NAICS_to_2002_NAICS.csv", DataFrame);
NAICS97To0297 = string.(NAICS97To02[1:1355, 1]);
NAICS97To0202 = string.(NAICS97To02[1:1355, 3]);
NAICS02Asdiv = string.(zeros(length(NAICS97To0202)));
for i in eachindex(NAICS97To0202);
    NAICS02Asdiv[i] = NAICS02Todiv[NAICS97To0202[i]];
end
NAICS97Todiv = Dict(NAICS97To0297 .=> NAICS02Asdiv);
NAICS97To0297Trunc = first.(string.(NAICS97To02[1:1355, 1]),4);
NAICS97TodivTrunc = Dict(NAICS97To0297Trunc .=> NAICS02Asdiv);

#Comm180 To div via NAIC 1997
NAICS97ToComm180 = CSV.read(datadir*"NAICS_to_Comm180.csv", DataFrame);
NAICS97ToComm18097 = first.([NAICS97ToComm180[1:90,4];NAICS97ToComm180[1:89,9]],4);
containsStar = findall( x -> occursin("*", x), NAICS97ToComm18097);
NAICS97ToComm18097[containsStar] = replace.(NAICS97ToComm18097[containsStar], "*" => "");
tooShort = findall( x -> occursin(",", x), NAICS97ToComm18097);
NAICS97ToComm18097[tooShort] = first.(NAICS97ToComm18097[tooShort],2);
NAICS97ToComm180180 = [NAICS97ToComm180[1:90,2];NAICS97ToComm180[1:89,7]];
containsStar = findall( x -> occursin("*", x), NAICS97ToComm180180);
NAICS97ToComm180180[containsStar] = replace.(NAICS97ToComm180180[containsStar], "*" => "");
containsSpace = findall( x -> occursin(" ", x), NAICS97ToComm180180);
NAICS97ToComm180180[containsSpace] = replace.(NAICS97ToComm180180[containsSpace], " " => "");
NAICS97Asdiv = string.(zeros(length(NAICS97ToComm18097)));
for i in eachindex(NAICS97ToComm18097);
    NAICS97ToComm18097[i] = rpad(parse(Int64, NAICS97ToComm18097[i]),4,"1");
end
Invalid4Dig = findall( x -> occursin("2311", x), NAICS97ToComm18097);
NAICS97ToComm18097[Invalid4Dig].=["2331"];
for i in eachindex(NAICS97ToComm18097);
    NAICS97Asdiv[i] = NAICS97TodivTrunc[NAICS97ToComm18097[i]];
end
Comm180Todiv=Dict(NAICS97ToComm180180 .=> NAICS97Asdiv);

#Final concordance
finalConcordance = [NAICS97ToComm180180 NAICS97Asdiv];
writedlm(datadir*"Comm180TodivConcordance.csv", finalConcordance, ',');
