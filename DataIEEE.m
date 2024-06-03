function [busdata,linedata,GenRestric]=DataIEEE(stmbus)

if stmbus==14
    [busdata,linedata,GenRestric] = IEEE14BusS; % Bus System
elseif stmbus==30
    [busdata,linedata,GenRestric] = IEEE30Bus; % Bus System
elseif stmbus==35
    [busdata,linedata,GenRestric] = IEEE35Bus; % Bus System SULSEL
end


