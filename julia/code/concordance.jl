#concordance dict between ANZSIC divisions (19 sectors plus Ownership of
#Dwellings) and various other industry classifications
using XLSX, ExcelReaders, DataFrames, Tables, JuMP, Ipopt, NamedArrays, DelimitedFiles, CSV, Tables;

#filepath cross system compatability code
if Sys.KERNEL === :NT || Sys.KERNEL === :Windows
    pathmark = "\\";
else
    pathmark = "/";
end

#20 Sector to 4 Sector
from20To4 = Dict(
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
  "T" => "Ownership"
 )

#IOIG to 20 Sector
IOSource = ExcelReaders.readxlsheet("data"*pathmark*"5209055001DO001_201819.xls", "Table 8");
IOIG = IOSource[4:117, 1];
ANZSICDiv = ["Agriculture,
             forestry and fishing",
             "Mining",
             "Manufacturing",
             "Electricity, gas, water and waste services",
             "Construction",
             "Wholesale trade",
             "Retail trade",
             "Accomodation and food services",
             "Transport, postal and warehousing",
             "Information media and telecommunications",
             "Financial and insurance services",
             "Rental, hiring and real estate services",
             "Professional, scientific and technical services",
             "Administrative and support services",
             "Public administration and safety",
             "Education and training",
             "Health care and social assistance",
             "Arts and recreation services",
             "Other services",
             "Ownership of Dwellings",
            ];
ANZSICDivShort = ["AgrForestFish", "Mining", "Manufacturing", "Utilities",
                  "Construction", "Wholesale", "Retail", "AccomFoodServ",
                  "Transport&Ware", "Communications", "Finance&Insur",
                  "RealEstate", "BusinessServ", "Admin", "PublicAdminSafe",
                  "Education", "Health&Social", "Arts&Rec", "OtherServices",
                  "Dwellings"];
ANZSICDivByLetter =["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O",
                    "P","Q","R","S", "T"];
from20ToInd = Dict{String, Int64}();
for i in eachindex(ANZSICDivByLetter);
    from20ToInd[ANZSICDivByLetter[i]] = Int(i);
end

function mapioig20(ioigcol)
  tmpdct = Dict{Float64, String}();
  for i in [1:1:length(ioigcol);]
    test = ioigcol[i] / 100
      if 1 <= test < 6
          tmpdct[ioigcol[i]]="A"
      elseif 6 <= test < 10
          tmpdct[ioigcol[i]]="B"
      elseif 10 <= test < 26
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
        elseif (66 <= test < 67 || 67.02 <= test < 69) #ownership of dwellings
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
      elseif 67 <= test < 67.02
          tmpdct[ioigcol[i]]="T"
      else
          print("ERROR: An input has fallen outside of the range of categories")
      end
  end
  return tmpdct
end
IOIGTo20 = mapioig20(IOIG)
#ISIC 4.0 To 20 Sectors
ANZSICISICSource = CSV.read("data"*pathmark*"ANZSIC06-ISIC3pt1.csv", DataFrame);
ANZSIC20 = ANZSICISICSource[6:1484, 1][findall(x -> typeof(x)<:String15, ANZSICISICSource[6:1484, 4])];
ISIC = ANZSICISICSource[6:1484, 4][findall(x -> typeof(x)<:String15, ANZSICISICSource[6:1484, 4])];
for i in eachindex(ISIC);
    ISIC[i]=strip(ISIC[i], ['p']);
end
ISICTo20 = Dict(ISIC .=> ANZSIC20);

#NAIC2007 To 20 Sectors via ISIC 4.0
NAICSISICSource = ExcelReaders.readxlsheet("data"*pathmark*"2007_NAICS_to_ISIC_4.xls", "NAICS 07 to ISIC 4 technical");
NAICS = string.(Int.(NAICSISICSource[4:1768,1]));
ISICAsANZSIC = NAICSISICSource[4:1768,3];
ISICAsANZSIC = string.(ISICAsANZSIC);
containsX = findall( x -> occursin("X", x), ISICAsANZSIC);
ISICAsANZSIC[containsX] = replace.(ISICAsANZSIC[containsX], "X" => "1");
ISICAsANZSIC = parse.(Float64, ISICAsANZSIC);
NAICSANZSIC20 = string.(zeros(length(ISICAsANZSIC)));
for i in eachindex(ISICAsANZSIC);
    NAICSANZSIC20[i] = ISICTo20[lpad(Int(ISICAsANZSIC[i]),4,"0")];
end
NAICS07To20 = Dict(NAICS .=> NAICSANZSIC20);

#NAIC2002 To 20 Sectors via NAIC2007
NAICS02To07 = CSV.read("data"*pathmark*"2002_to_2007_NAICS.csv", DataFrame);
NAICS02To0702 = string.(NAICS02To07[3:1202, 1]);
NAICS02To0707 = string.(NAICS02To07[3:1202, 3]);
NAICS07As20 = string.(zeros(length(NAICS02To0707)));
for i in eachindex(NAICS02To0707);
    NAICS07As20[i] = NAICS07To20[NAICS02To0707[i]];
end
NAICS02To20 = Dict(NAICS02To0702 .=> NAICS07As20);

#NAIC1997 To 20 Sectors via NAIC2002
NAICS97To02 = CSV.read("data"*pathmark*"1997_NAICS_to_2002_NAICS.csv", DataFrame);
NAICS97To0297 = string.(NAICS97To02[1:1355, 1]);
NAICS97To0202 = string.(NAICS97To02[1:1355, 3]);
NAICS02As20 = string.(zeros(length(NAICS97To0202)));
for i in eachindex(NAICS97To0202);
    NAICS02As20[i] = NAICS02To20[NAICS97To0202[i]];
end
NAICS97To20 = Dict(NAICS97To0297 .=> NAICS02As20);
NAICS97To0297Trunc = first.(string.(NAICS97To02[1:1355, 1]),4);
NAICS97To20Trunc = Dict(NAICS97To0297Trunc .=> NAICS02As20);

#Comm180 To 20 Sectors via NAIC 1997
NAICS97ToComm180 = CSV.read("data"*pathmark*"NAICS_to_Comm180.csv", DataFrame);
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
NAICS97As20 = string.(zeros(length(NAICS97ToComm18097)));
for i in eachindex(NAICS97ToComm18097);
    NAICS97ToComm18097[i] = rpad(parse(Int64, NAICS97ToComm18097[i]),4,"1");
end
Invalid4Dig = findall( x -> occursin("2311", x), NAICS97ToComm18097);
NAICS97ToComm18097[Invalid4Dig].=["2331"];
for i in eachindex(NAICS97ToComm18097);
    NAICS97As20[i] = NAICS97To20Trunc[NAICS97ToComm18097[i]];
end
Comm180To20=Dict(NAICS97ToComm180180 .=> NAICS97As20);

#Final concordance
finalConcordance = [NAICS97ToComm180180 NAICS97As20];
writedlm("data"*pathmark*"Comm180To20Concordance.csv", finalConcordance, ',');
